!>
!! @file m_global_parameters.f90
!! @brief Contains module m_global_parameters

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief The module contains all of the parameters describing the program
!!              logistics, the computational domain and the simulation algorithm.
!!              Additionally, for the volume fraction model, physical parameters
!!              of each of the fluids present in the flow are located here. They
!!              include stiffened gas equation of state parameters, the Reynolds
!!              numbers and the Weber numbers.
module m_global_parameters

    ! Dependencies =============================================================
#ifdef MFC_MPI
    use mpi                    !< Message passing interface (MPI) module
#endif

    use m_derived_types        !< Definitions of the derived types

#ifdef _OPENACC
    use openacc
#endif

    ! ==========================================================================

    implicit none

    ! Logistics ================================================================
    integer :: num_procs             !< Number of processors
    integer, parameter :: num_stcls_min = 5     !< Mininum # of stencils
    integer, parameter :: path_len = 400   !< Maximum path length
    integer, parameter :: name_len = 50    !< Maximum name length
    character, parameter :: dflt_char = ' '   !< Default string value
    real(kind(0d0)), parameter :: dflt_real = -1d6  !< Default real value
    integer, parameter :: dflt_int = -100  !< Default integer value
    real(kind(0d0)), parameter :: sgm_eps = 1d-16 !< Segmentation tolerance
    character(LEN=path_len) :: case_dir              !< Case folder location
    logical :: run_time_info         !< Run-time output flag
    integer :: t_step_old            !< Existing IC/grid folder
    ! ==========================================================================

    ! Computational Domain Parameters ==========================================
    integer :: proc_rank !< Rank of the local processor

    !> @name Number of cells in the x-, y- and z-directions, respectively
    !> @{
    integer :: m, n
    !> @}

    !> @name Global number of cells in each direction
    !> @{
    integer :: m_glb, n_glb
    !> @}

    !> @name Cell-boundary (CB) locations in the x-, y- and z-directions, respectively
    !> @{
    real(kind(0d0)), target, allocatable, dimension(:) :: x_cb, y_cb
    !> @}

    !> @name Cell-center (CC) locations in the x-, y- and z-directions, respectively
    !> @{
    real(kind(0d0)), target, allocatable, dimension(:) :: x_cc, y_cc
    !> @}

    !> @name Cell-width distributions in the x-, y- and z-directions, respectively
    !> @{
    real(kind(0d0)), target, allocatable, dimension(:) :: dx, dy
    !> @}

    real(kind(0d0)) :: dt !< Size of the time-step

    !$acc declare create(x_cb, y_cb, x_cc, y_cc, dx, dy, dt, m, n)

    !> @name Starting time-step iteration, stopping time-step iteration and the number
    !! of time-step iterations between successive solution backups, respectively
    !> @{
    integer :: t_step_start, t_step_stop, t_step_save
    !> @}

    ! ==========================================================================

    ! Simulation Algorithm Parameters ==========================================
    #:if MFC_CASE_OPTIMIZATION
        integer, parameter :: num_dims = ${num_dims}$       !< Number of spatial dimensions
    #:else
        integer :: num_dims       !< Number of spatial dimensions
    #:endif
    integer :: num_fluids     !< Number of fluids in the flow
    integer :: time_stepper   !< Time-stepper algorithm

    #:if MFC_CASE_OPTIMIZATION
        integer, parameter :: weno_polyn = ${weno_polyn}$ !< Degree of the WENO polynomials (polyn)
        integer, parameter :: weno_order = ${weno_order}$ !< Order of the WENO reconstruction
    #:else
        integer :: weno_polyn     !< Degree of the WENO polynomials (polyn)
        integer :: weno_order     !< Order of the WENO reconstruction
    #:endif

    real(kind(0d0)) :: weno_eps       !< Binding for the WENO nonlinear weights
    logical :: weno_Re_flux   !< WENO reconstruct velocity gradients for viscous stress tensor

    integer :: cpu_start, cpu_end, cpu_rate

    #:if not MFC_CASE_OPTIMIZATION
        !$acc declare create(num_dims, weno_polyn, weno_order)
    #:endif

!$acc declare create(num_fluids, weno_eps)

    !> @name Boundary conditions (BC) in the x-, y- and z-directions, respectively
    !> @{
    type(int_bounds_info) :: bc_x, bc_y
    type(int_bounds_info) :: bc_x_glb, bc_y_glb
    !> @}

    logical :: parallel_io !< Format of the data files
    integer :: precision !< Precision of output files

    integer, allocatable, dimension(:) :: proc_coords !<
    !! Processor coordinates in MPI_CART_COMM

    integer, allocatable, dimension(:) :: start_idx !<
    !! Starting cell-center index of local processor in global grid

    type(mpi_io_var), public :: MPI_IO_DATA

    !> @name MPI info for parallel IO with Lustre file systems
    !> @{
    character(LEN=name_len) :: mpiiofs
    integer :: mpi_info_int
    !> @}

    integer, private :: ierr

    !> @name Annotations of the structure of the state and flux vectors in terms of the
    !! size and the configuration of the system of equations to which they belong
    !> @{
    integer :: sys_size                  !< Number of unknowns in system of eqns.
    type(int_bounds_info) :: cont_idx                  !< Indexes of first & last continuity eqns.
    type(int_bounds_info) :: mom_idx                   !< Indexes of first & last momentum eqns.
    integer :: E_idx                     !< Index of energy equation
    type(int_bounds_info) :: adv_idx                   !< Indexes of first & last advection eqns.
    integer :: gamma_idx                 !< Index of specific heat ratio func. eqn.
    integer :: pi_inf_idx                !< Index of liquid stiffness func. eqn.
    !> @}

    !> @name The number of fluids, along with their identifying indexes, respectively,
    !! for which viscous effects, e.g. the shear and/or the volume Reynolds (Re)
    !! numbers, will be non-negligible.
    !> @{
    integer, dimension(2) :: Re_size
    integer, allocatable, dimension(:, :) :: Re_idx
    !> @}
!$acc declare create(Re_size, Re_idx)

    !> @name The coordinate direction indexes and flags (flg), respectively, for which
    !! the configurations will be determined with respect to a working direction
    !! and that will be used to isolate the contributions, in that direction, in
    !! the dimensionally split system of equations.
    !> @{
    integer, dimension(2) :: dir_idx
    real(kind(0d0)), dimension(2) :: dir_flg
    !> @}
!$acc declare create(dir_idx, dir_flg)

    integer :: buff_size !<
    !! The number of cells that are necessary to be able to store enough boundary
    !! conditions data to march the solution in the physical computational domain
    !! to the next time-step.

    integer :: startx, starty

!$acc declare create(sys_size, buff_size, startx, starty, E_idx, gamma_idx, pi_inf_idx)

    ! END: Simulation Algorithm Parameters =====================================

    ! Fluids Physical Parameters ===============================================

    type(physical_parameters), dimension(num_fluids_max) :: fluid_pp !<
    !! Database of the physical parameters of each of the fluids that is present
    !! in the flow. These include the stiffened gas equation of state parameters,
    !! the Reynolds numbers and the Weber numbers.

    ! ==========================================================================

    integer :: fd_order !<
    !! The order of the finite-difference (fd) approximations of the first-order
    !! derivatives that need to be evaluated when the CoM or flow probe data
    !! files are to be written at each time step

    integer :: fd_number !<
    !! The finite-difference number is given by MAX(1, fd_order/2). Essentially,
    !! it is a measure of the half-size of the finite-difference stencil for the
    !! selected order of accuracy.

    logical :: probe_wrt
    integer :: num_probes
    type(probe_parameters), dimension(num_probes_max) :: probe

     integer :: momxb, momxe
     integer :: advxb, advxe
     integer :: contxb, contxe
     !$acc declare create(momxb, momxe, advxb, advxe, contxb, contxe)

    real(kind(0d0)), allocatable, dimension(:) :: gammas, pi_infs
    !$acc declare create(gammas, pi_infs)


    real(kind(0d0)) :: mytime       !< Current simulation time
    real(kind(0d0)) :: finaltime    !< Final simulation time

    logical :: weno_flat, riemann_flat, cu_mpi

    ! ======================================================================

    ! Mathematical and Physical Constants ======================================
    real(kind(0d0)), parameter :: pi = 3.141592653589793d0 !< Pi

    ! ==========================================================================

contains

    !> Assigns default values to the user inputs before reading
        !!  them in. This enables for an easier consistency check of
        !!  these parameters once they are read from the input file.
    subroutine s_assign_default_values_to_user_inputs() ! ------------------

        integer :: i !< Generic loop iterator

        ! Logistics
        case_dir = dflt_char
        run_time_info = .false.
        t_step_old = dflt_int

        ! Computational domain parameters
        m = dflt_int
        n = 0

        dt = dflt_real

        t_step_start = dflt_int
        t_step_stop = dflt_int
        t_step_save = dflt_int

        ! Simulation algorithm parameters
        num_fluids = dflt_int
        time_stepper = dflt_int
        weno_eps = dflt_real
        weno_Re_flux = .false.
        parallel_io = .false.
        precision = 2
        cu_mpi = .false.

        bc_x%beg = dflt_int; bc_x%end = dflt_int
        bc_y%beg = dflt_int; bc_y%end = dflt_int

        ! Fluids physical parameters
        do i = 1, num_fluids_max
            fluid_pp(i)%gamma = dflt_real
            fluid_pp(i)%pi_inf = dflt_real
            fluid_pp(i)%Re(:) = dflt_real
        end do

        #:if not MFC_CASE_OPTIMIZATION
            weno_order = dflt_int
        #:endif

        fd_order = dflt_int
        probe_wrt = .false.
        num_probes = dflt_int

        do i = 1, num_probes_max
            probe(i)%x = dflt_real
            probe(i)%y = dflt_real
        end do

    end subroutine s_assign_default_values_to_user_inputs ! ----------------

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_global_parameters_module() ! -------------------

        integer :: i, k

        #:if not MFC_CASE_OPTIMIZATION
            ! Determining the degree of the WENO polynomials
            weno_polyn = (weno_order - 1)/2
!$acc update device(weno_polyn)
        #:endif

        ! Initializing the number of fluids for which viscous effects will
        ! be non-negligible, the number of distinctive material interfaces
        ! for which surface tension will be important and also, the number
        ! of fluids for which the physical and geometric curvatures of the
        ! interfaces will be computed
        Re_size = 0

        cont_idx%beg = 1
        cont_idx%end = num_fluids
        mom_idx%beg = cont_idx%end + 1
        mom_idx%end = cont_idx%end + num_dims
        E_idx = mom_idx%end + 1
        adv_idx%beg = E_idx + 1
        adv_idx%end = E_idx + num_fluids

        sys_size = adv_idx%end

        ! Determining the number of fluids for which the shear and the
        ! volume Reynolds numbers, e.g. viscous effects, are important
        do i = 1, num_fluids
            if (fluid_pp(i)%Re(1) > 0) Re_size(1) = Re_size(1) + 1
            if (fluid_pp(i)%Re(2) > 0) Re_size(2) = Re_size(2) + 1
        end do

        ! Bookkeeping the indexes of any viscous fluids and any pairs of
        ! fluids whose interface will support effects of surface tension
        if (any(Re_size > 0)) then

            @:ALLOCATE(Re_idx(1:2, 1:maxval(Re_size)))

            k = 0
            do i = 1, num_fluids
                if (fluid_pp(i)%Re(1) > 0) then
                    k = k + 1; Re_idx(1, k) = i
                end if
            end do

            k = 0
            do i = 1, num_fluids
                if (fluid_pp(i)%Re(2) > 0) then
                    k = k + 1; Re_idx(2, k) = i
                end if
            end do

        end if

        ! END: Volume Fraction Model =======================================

        allocate (MPI_IO_DATA%view(1:sys_size))
        allocate (MPI_IO_DATA%var(1:sys_size))

        do i = 1, sys_size
            allocate (MPI_IO_DATA%var(i)%sf(0:m, 0:n))
            MPI_IO_DATA%var(i)%sf => null()
        end do

!$acc update device(Re_size)
        ! Determining the number of cells that are needed in order to store
        ! sufficient boundary conditions data as to iterate the solution in
        ! the physical computational domain from one time-step iteration to
        ! the next one
        if (any(Re_size > 0)) then
            buff_size = 2*weno_polyn + 2
        else
            buff_size = weno_polyn + 2
        end if

        ! Configuring Coordinate Direction Indexes =========================
        if (probe_wrt) then
            fd_number = max(1, fd_order/2)
            buff_size = buff_size + fd_number
        end if

        startx = -buff_size
        starty = 0
        if (n > 0) then
            starty = -buff_size
        end if

!$acc update device(startx, starty)

        momxb = mom_idx%beg
        momxe = mom_idx%end
        advxb = adv_idx%beg
        advxe = adv_idx%end
        contxb = cont_idx%beg
        contxe = cont_idx%end

!$acc update device(momxb, momxe, advxb, advxe, contxb, contxe, sys_size, buff_size, E_idx)

        ! Allocating grid variables for the x-, y- and z-directions
        @:ALLOCATE(x_cb(-1 - buff_size:m + buff_size))
        @:ALLOCATE(x_cc(-buff_size:m + buff_size))
        @:ALLOCATE(dx(-buff_size:m + buff_size))

        if (n == 0) return;
        
        @:ALLOCATE(y_cb(-1 - buff_size:n + buff_size))
        @:ALLOCATE(y_cc(-buff_size:n + buff_size))
        @:ALLOCATE(dy(-buff_size:n + buff_size))

    end subroutine s_initialize_global_parameters_module ! -----------------


    !> Initializes parallel infrastructure
    subroutine s_initialize_parallel_io() ! --------------------------------

        #:if not MFC_CASE_OPTIMIZATION
            num_dims = 1 + min(1, n) 
        #:endif

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

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_global_parameters_module() ! ---------------------

        ! Deallocating the variables bookkeeping the indexes of any viscous
        ! fluids and any pairs of fluids whose interfaces supported effects
        ! of surface tension
        if (any(Re_size > 0)) then
            @:DEALLOCATE(Re_idx)
        end if

        ! Deallocating grid variables for the x-, y- and z-directions
        @:DEALLOCATE(x_cb, x_cc, dx)
        
        if (n == 0) return;
        @:DEALLOCATE(y_cb, y_cc, dy)

    end subroutine s_finalize_global_parameters_module ! -------------------

end module m_global_parameters
