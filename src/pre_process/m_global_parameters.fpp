!>
!! @file m_global_parameters.f90
!! @brief Contains module m_global_parameters

!> @brief This module contains all of the parameters characterizing the
!!              computational domain, simulation algorithm, initial condition
!!              and the stiffened equation of state.
module m_global_parameters

    ! Dependencies =============================================================
#ifdef MFC_MPI
    use mpi                     ! Message passing interface (MPI) module
#endif

    use m_derived_types         ! Definitions of the derived types

    ! ==========================================================================

    implicit none

    ! Logistics ================================================================
    integer :: num_procs            !< Number of processors
    integer, parameter :: num_stcls_min = 5    !< Mininum # of stencils
    integer, parameter :: path_len = 400  !< Maximum path length
    integer, parameter :: name_len = 50   !< Maximum name length
    real(kind(0d0)), parameter :: dflt_real = -1d6 !< Default real value
    integer, parameter :: dflt_int = -100 !< Default integer value
    character(LEN=path_len) :: case_dir             !< Case folder location
    logical :: old_grid             !< Use existing grid data
    logical :: old_ic               !< Use existing IC data
    integer :: t_step_old           !< Existing IC/grid folder
    ! ==========================================================================

    ! Computational Domain Parameters ==========================================

    integer :: proc_rank !< Rank of the local processor

    integer :: m
    integer :: n

    integer :: m_glb, n_glb

    integer :: num_dims !< Number of spatial dimensions

    real(kind(0d0)), allocatable, dimension(:) :: x_cc, y_cc !<
    !! Locations of cell-centers (cc) in x-, y- and z-directions, respectively

    real(kind(0d0)), allocatable, dimension(:) :: x_cb, y_cb!<
    !! Locations of cell-boundaries (cb) in x-, y- and z-directions, respectively

    real(kind(0d0)) :: dx, dy !<
    !! Minimum cell-widths in the x-, y- and z-coordinate directions

    type(bounds_info) :: x_domain, y_domain !<
    !! Locations of the domain bounds in the x-, y- and z-coordinate directions

    logical :: stretch_x, stretch_y !<
    !! Grid stretching flags for the x-, y- and z-coordinate directions

    ! Parameters of the grid stretching function for the x-, y- and z-coordinate
    ! directions. The "a" parameters are a measure of the rate at which the grid
    ! is stretched while the remaining parameters are indicative of the location
    ! on the grid at which the stretching begins.
    real(kind(0d0)) :: a_x, a_y
    integer :: loops_x, loops_y
    real(kind(0d0)) :: x_a, y_a
    real(kind(0d0)) :: x_b, y_b

    ! ==========================================================================

    ! Simulation Algorithm Parameters ==========================================
    integer :: num_fluids      !< Number of different fluids present in the flow
    integer :: sys_size        !< Number of unknowns in the system of equations
    integer :: weno_order      !< Order of accuracy for the WENO reconstruction

    ! Annotations of the structure, i.e. the organization, of the state vectors
    type(int_bounds_info) :: cont_idx                   !< Indexes of first & last continuity eqns.
    type(int_bounds_info) :: mom_idx                    !< Indexes of first & last momentum eqns.
    integer :: E_idx                      !< Index of total energy equation
    type(int_bounds_info) :: adv_idx                    !< Indexes of first & last advection eqns.

    type(int_bounds_info) :: bc_x, bc_y !<
    !! Boundary conditions in the x-, y- and z-coordinate directions

    logical :: parallel_io !< Format of the data files
    integer :: precision !< Precision of output files


    integer, allocatable, dimension(:) :: proc_coords !<
    !! Processor coordinates in MPI_CART_COMM

    integer, allocatable, dimension(:) :: start_idx !<
    !! Starting cell-center index of local processor in global grid

#ifdef MFC_MPI

    type(mpi_io_var), public :: MPI_IO_DATA

    character(LEN=name_len) :: mpiiofs
    integer :: mpi_info_int !<
    !! MPI info for parallel IO with Lustre file systems

#endif

    integer, private :: ierr
    ! ==========================================================================

    ! Initial Condition Parameters =============================================
    integer :: num_patches     !< Number of patches composing initial condition

    type(ic_patch_parameters), dimension(num_patches_max) :: patch_icpp !<
    !! Database of the initial condition patch parameters (icpp) for each of the
    !! patches employed in the configuration of the initial condition. Note that
    !! the maximum allowable number of patches, num_patches_max, may be changed
    !! in the module m_derived_types.f90.
    ! ==========================================================================

    ! Fluids Physical Parameters ===============================================
    type(physical_parameters), dimension(num_fluids_max) :: fluid_pp !<
    !! Database of the physical parameters of each of the fluids that is present
    !! in the flow. These include the stiffened gas equation of state parameters,
    !! the Reynolds numbers and the Weber numbers.

    ! ==========================================================================


    !> @name Index variables used for m_variables_conversion
    !> @{
    integer :: momxb, momxe
    integer :: advxb, advxe
    integer :: contxb, contxe
    integer :: intxb, intxe
    !> @}

    integer, allocatable, dimension(:, :) :: logic_grid

    ! Mathematical and Physical Constants ======================================
    real(kind(0d0)), parameter :: pi = 3.141592653589793d0
    ! ==========================================================================

contains

    !>  Assigns default values to user inputs prior to reading
        !!              them in. This allows for an easier consistency check of
        !!              these parameters once they are read from the input file.
    subroutine s_assign_default_values_to_user_inputs() ! ------------------

        integer :: i !< Generic loop operator

        ! Logistics
        case_dir = ' '
        old_grid = .false.
        old_ic = .false.
        t_step_old = dflt_int

        ! Computational domain parameters
        m = dflt_int; n = 0; 

        x_domain%beg = dflt_real
        x_domain%end = dflt_real
        y_domain%beg = dflt_real
        y_domain%end = dflt_real

        stretch_x = .false.
        stretch_y = .false.

        a_x = dflt_real
        a_y = dflt_real
        loops_x = 1
        loops_y = 1
        x_a = dflt_real
        x_b = dflt_real
        y_a = dflt_real
        y_b = dflt_real

        ! Simulation algorithm parameters
        num_fluids = dflt_int
        weno_order = dflt_int

        bc_x%beg = dflt_int
        bc_x%end = dflt_int
        bc_y%beg = dflt_int
        bc_y%end = dflt_int

        parallel_io = .false.
        precision = 2

        ! Initial condition parameters
        num_patches = dflt_int

        do i = 1, num_patches_max
            patch_icpp(i)%geometry = dflt_int
            patch_icpp(i)%x_centroid = dflt_real
            patch_icpp(i)%y_centroid = dflt_real
            patch_icpp(i)%length_x = dflt_real
            patch_icpp(i)%length_y = dflt_real
            patch_icpp(i)%radius = dflt_real
            patch_icpp(i)%epsilon = dflt_real
            patch_icpp(i)%beta = dflt_real
            patch_icpp(i)%normal = dflt_real
            patch_icpp(i)%radii = dflt_real
            patch_icpp(i)%alter_patch = .false.
            patch_icpp(i)%alter_patch(0) = .true.
            patch_icpp(i)%smoothen = .false.
            patch_icpp(i)%smooth_patch_id = i
            patch_icpp(i)%smooth_coeff = dflt_real
            patch_icpp(i)%alpha_rho = dflt_real
            patch_icpp(i)%rho = dflt_real
            patch_icpp(i)%vel = dflt_real
            patch_icpp(i)%pres = dflt_real
            patch_icpp(i)%alpha = dflt_real
            patch_icpp(i)%gamma = dflt_real
            patch_icpp(i)%pi_inf = dflt_real
        end do


        ! Fluids physical parameters
        do i = 1, num_fluids_max
            fluid_pp(i)%gamma = dflt_real
            fluid_pp(i)%pi_inf = dflt_real
        end do

    end subroutine s_assign_default_values_to_user_inputs ! ----------------

    !> Computation of parameters, allocation procedures, and/or
        !! any other tasks needed to properly setup the module
    subroutine s_initialize_global_parameters_module() ! ----------------------

        integer :: i

        ! Determining the layout of the state vectors and overall size of
        ! the system of equations, given the dimensionality and choice of
        ! the equations of motion

        ! Annotating structure of the state and flux vectors belonging
        ! to the system of equations defined by the selected number of
        ! spatial dimensions and the volume fraction model
        cont_idx%beg = 1
        cont_idx%end = num_fluids
        mom_idx%beg = cont_idx%end + 1
        mom_idx%end = cont_idx%end + num_dims
        E_idx = mom_idx%end + 1
        adv_idx%beg = E_idx + 1
        adv_idx%end = E_idx + num_fluids

        sys_size = adv_idx%end
        ! ==================================================================
        momxb = mom_idx%beg
        momxe = mom_idx%end
        advxb = adv_idx%beg
        advxe = adv_idx%end
        contxb = cont_idx%beg
        contxe = cont_idx%end

        ! ==================================================================

#ifdef MFC_MPI

        allocate (MPI_IO_DATA%view(1:sys_size))
        allocate (MPI_IO_DATA%var(1:sys_size))

        do i = 1, sys_size
            allocate (MPI_IO_DATA%var(i)%sf(0:m, 0:n))
            MPI_IO_DATA%var(i)%sf => null()
        end do

#endif

        ! Allocating grid variables for the x-direction
        allocate (x_cc(0:m), x_cb(-1:m))
        ! Allocating grid variables for the y- and z-directions
        if (n > 0) then
            allocate (y_cc(0:n), y_cb(-1:n))
        end if

        allocate (logic_grid(0:m, 0:n))

    end subroutine s_initialize_global_parameters_module ! --------------------


    subroutine s_initialize_parallel_io() ! --------------------------------

        num_dims = 1 + min(1, n)

        allocate (proc_coords(1:num_dims))

        if (parallel_io .neqv. .true.) return

#ifdef MFC_MPI

        ! Option for Lustre file system (Darter/Comet/Stampede)
        write (mpiiofs, '(A)') '/lustre_'
        mpiiofs = trim(mpiiofs)
        call MPI_INFO_CREATE(mpi_info_int, ierr)
        call MPI_INFO_SET(mpi_info_int, 'romio_ds_write', 'disable', ierr)

        ! Option for UNIX file system (Hooke/Thomson)
        ! WRITE(mpiiofs, '(A)') '/ufs_'
        ! mpiiofs = TRIM(mpiiofs)
        ! mpi_info_int = MPI_INFO_NULL

        allocate (start_idx(1:num_dims))

#endif

    end subroutine s_initialize_parallel_io ! ------------------------------

    subroutine s_finalize_global_parameters_module() ! ------------------------

        integer :: i

        ! Deallocating grid variables for the x-direction
        deallocate (x_cc, x_cb)
        ! Deallocating grid variables for the y- and z-directions
        if (n > 0) then
            deallocate (y_cc, y_cb)
        end if

        deallocate (proc_coords)

#ifdef MFC_MPI

        if (parallel_io) then
            deallocate (start_idx)
            do i = 1, sys_size
                MPI_IO_DATA%var(i)%sf => null()
            end do

            deallocate (MPI_IO_DATA%var)
            deallocate (MPI_IO_DATA%view)
        end if

#endif

    end subroutine s_finalize_global_parameters_module ! ----------------------


end module m_global_parameters
