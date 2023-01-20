!>
!! @file m_start_up.f90
!! @brief Contains module m_start_up

#:include 'case.fpp'

!> @brief This module contains subroutines that read, and check consistency
!!              of, the user provided inputs and data.
module m_start_up

    ! Dependencies =============================================================
    use m_derived_types          !< Definitions of the derived types

    use m_global_parameters      !< Global parameters for the code

    use m_mpi_proxy              !< Message passing interface (MPI) module proxy

    use m_data_output            !< Procedures to write the grid data and the
                                 !! conservative variables to files

#ifdef MFC_MPI
    use mpi                      !< Message passing interface (MPI) module
#endif

    use m_compile_specific

    use m_check_patches
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_start_up_module, &
         s_read_input_file, &
         s_check_input_file, &
         s_read_grid_data_files, &
         s_read_ic_data_files, &
         s_read_serial_grid_data_files, &
         s_read_serial_ic_data_files, &
         s_read_parallel_grid_data_files, &
         s_read_parallel_ic_data_files, &
         s_check_grid_data_files, &
         s_finalize_start_up_module

    abstract interface ! ===================================================

        subroutine s_read_abstract_grid_data_files(dflt_int)! ----------

            integer, intent(IN) :: dflt_int

        end subroutine s_read_abstract_grid_data_files ! ---------------

        subroutine s_read_abstract_ic_data_files(q_cons_vf) ! -----------

            import :: scalar_field, sys_size

            ! Conservative variables
            type(scalar_field), &
                dimension(sys_size), &
                intent(INOUT) :: q_cons_vf

        end subroutine s_read_abstract_ic_data_files ! -----------------

    end interface ! ========================================================

    character(LEN=path_len + name_len) :: proc_rank_dir !<
    !! Location of the folder associated with the rank of the local processor

    character(LEN=path_len + 2*name_len), private :: t_step_dir !<
    !! Possible location of time-step folder containing preexisting grid and/or
    !! conservative variables data to be used as starting point for pre-process

    procedure(s_read_abstract_grid_data_files), pointer :: s_read_grid_data_files => null()
    procedure(s_read_abstract_ic_data_files), pointer :: s_read_ic_data_files => null()

contains

    !>  Reads the configuration file pre_process.inp, in order to
        !!      populate the parameters in module m_global_parameters.f90
        !!      with the user provided inputs
    subroutine s_read_input_file() ! ---------------------------------------

        character(LEN=name_len) :: file_loc  !<
            !! Generic string used to store the address of a particular file

        logical :: file_check !<
            !! Generic logical used for the purpose of asserting whether a file
            !! is or is not present in the designated location

        ! Namelist for all of the parameters to be inputed by the user
        namelist /user_inputs/ case_dir, old_grid, old_ic, &
            t_step_old, m, n, x_domain, y_domain, &
            stretch_x, stretch_y, a_x, a_y, &
            x_a, y_a, x_b, y_b, &
            model_eqns, num_fluids, &
            adv_alphan, mpp_lim, &
            weno_order, bc_x, bc_y, num_patches, &
            patch_icpp, fluid_pp, &
            precision, parallel_io, &
            fluid_rho, &
            cyl_coord, loops_x, loops_y

        ! Inquiring the status of the pre_process.inp file
        file_loc = 'pre_process.inp'
        inquire (FILE=trim(file_loc), EXIST=file_check)

        ! Checking whether the input file is there. If it is, the input file
        ! is read. If not, the program is terminated.
        if (file_check) then
            open (1, FILE=trim(file_loc), FORM='formatted', &
                  STATUS='old', ACTION='read')
            read (1, NML=user_inputs)
            close (1)
            m_glb = m
            n_glb = n
        else
            print '(A)', 'File pre_process.inp is missing. Exiting ...'
            call s_mpi_abort()
        end if

    end subroutine s_read_input_file ! -------------------------------------

    !>  Checking that the user inputs make sense, i.e. that the
    !!      individual choices are compatible with the code's options
    !!      and that the combination of these choices results into a
    !!      valid configuration for the pre-process
    subroutine s_check_input_file() ! --------------------------------------

        character(LEN=len_trim(case_dir)) :: file_loc !<
            !! Generic string used to store the address of a particular file

        logical :: dir_check !<
            !! Logical variable used to test the existence of folders

        integer :: i !<
            !! Generic loop iterator

        ! Checking the existence of the case folder
        case_dir = adjustl(case_dir)

        file_loc = trim(case_dir)//'/.'

        call my_inquire(file_loc, dir_check)

        ! Constraint on the location of the case directory
        if (dir_check .neqv. .true.) then
            print '(A)', 'Unsupported choice for the value of '// &
                'case_dir.'
            print '(A)', 'WARNING: Ensure that compiler flags/choices in Makefiles match your compiler! '
            print '(A)', 'WARNING: Ensure that preprocessor flags are enabled! '
            call s_mpi_abort()

            ! Constraints on the use of a preexisting grid and initial condition
        elseif ((old_grid .neqv. .true.) .and. old_ic) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for old_grid and old_ic. Exiting ...'
            call s_mpi_abort()

        elseif ((old_grid .or. old_ic) .and. t_step_old == dflt_int) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for old_grid and old_ic and t_step_old. Exiting ...'
            call s_mpi_abort()

            ! Constraints on dimensionality and the number of cells for the grid
        elseif (m <= 0) then
            print '(A)', 'Unsupported choice for the value of m. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (n < 0) then
            print '(A)', 'Unsupported choice for the value of n. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif ((m + 1)*(n + 1) &
                < &
                2**(min(1, m) + min(1, n))*num_procs) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for num_procs, m, n. '// &
                'Exiting ...'
            call s_mpi_abort()

            ! Constraints on domain boundaries locations in the x-direction
        elseif ((old_grid .and. x_domain%beg /= dflt_real) &
                .or. &
                ((old_grid .neqv. .true.) .and. &
                 x_domain%beg == dflt_real)) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for old_grid and x_domain%beg. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif ((old_grid .and. x_domain%end /= dflt_real) &
                .or. &
                ((old_grid .neqv. .true.) .and. &
                 x_domain%end == dflt_real)) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for old_grid and x_domain%end. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif ((old_grid .neqv. .true.) &
                .and. &
                x_domain%beg >= x_domain%end) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for old_grid, x_domain%beg and '// &
                'x_domain%end. Exiting ...'
            call s_mpi_abort()
        end if

        if (cyl_coord .neqv. .true.) then ! Cartesian coordinates

            ! Constraints on domain boundaries locations in the y-direction
            if ((n == 0 .and. y_domain%beg /= dflt_real) &
                .or. &
                (n > 0 &
                 .and. &
                 ((old_grid .and. y_domain%beg /= dflt_real) &
                  .or. &
                  ((old_grid .neqv. .true.) .and. &
                   y_domain%beg == dflt_real)))) then
                print '(A)', 'Unsupported choice of the combination of '// &
                    'values for old_grid, n and y_domain%beg. '// &
                    'Exiting ...'
                call s_mpi_abort()
            elseif ((n == 0 .and. y_domain%end /= dflt_real) &
                    .or. &
                    (n > 0 &
                     .and. &
                     ((old_grid .and. y_domain%end /= dflt_real) &
                      .or. &
                      ((old_grid .neqv. .true.) .and. &
                       y_domain%end == dflt_real)))) then
                print '(A)', 'Unsupported choice of the combination of '// &
                    'values for old_grid, n and y_domain%end. '// &
                    'Exiting ...'
                call s_mpi_abort()
            elseif (n > 0 &
                    .and. &
                    (old_grid .neqv. .true.) &
                    .and. &
                    y_domain%beg >= y_domain%end) then
                print '(A)', 'Unsupported choice of the combination of '// &
                    'values for old_grid, n, y_domain%beg and '// &
                    'y_domain%end. Exiting ...'
                call s_mpi_abort()
            end if
        end if

        ! Constraints on the grid stretching in the x-direction
        if (old_grid .and. stretch_x) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for old_grid and stretch_x. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (stretch_x .and. a_x == dflt_real) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for stretch_x and a_x. Exiting ...'
            call s_mpi_abort()
        elseif (stretch_x .and. x_a == dflt_real) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for stretch_x and x_a. Exiting ...'
            call s_mpi_abort()
        elseif (stretch_x .and. x_b == dflt_real) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for stretch_x and x_b. Exiting ...'
            call s_mpi_abort()
        elseif (stretch_x .and. x_a >= x_b) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for stretch_x, x_a and x_b. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (stretch_x &
                .and. &
                (a_x + log(cosh(a_x*(x_domain%beg - x_a))) &
                 + log(cosh(a_x*(x_domain%beg - x_b))) &
                 - 2d0*log(cosh(0.5d0*a_x*(x_b - x_a))))/a_x <= 0d0) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for x_domain%beg, stretch_x, a_x, '// &
                'x_a, and x_b. Exiting ...'
            call s_mpi_abort()
        elseif (stretch_x &
                .and. &
                (a_x + log(cosh(a_x*(x_domain%end - x_a))) &
                 + log(cosh(a_x*(x_domain%end - x_b))) &
                 - 2d0*log(cosh(0.5d0*a_x*(x_b - x_a))))/a_x <= 0d0) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for x_domain%end, stretch_x, a_x, '// &
                'x_a, and x_b. Exiting ...'
            call s_mpi_abort()
        elseif (loops_x < 1) then
            print '(A)', 'Unsupported choice for the value of loops_x. '// &
                'Exiting ...'
            call s_mpi_abort()

            ! Constraints on the grid stretching in the y-direction
        elseif (old_grid .and. stretch_y) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for old_grid and stretch_y. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (n == 0 .and. stretch_y) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for n and stretch_y. Exiting ...'
            call s_mpi_abort()
        elseif (stretch_y .and. a_y == dflt_real) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for stretch_y and a_y. Exiting ...'
            call s_mpi_abort()
        elseif (stretch_y .and. y_a == dflt_real) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for stretch_y and y_a. Exiting ...'
            call s_mpi_abort()
        elseif (stretch_y .and. y_b == dflt_real) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for stretch_y and y_b. Exiting ...'
            call s_mpi_abort()
        elseif (stretch_y .and. y_a >= y_b) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for stretch_y, y_a and y_b. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (stretch_y &
                .and. &
                (a_y + log(cosh(a_y*(y_domain%beg - y_a))) &
                 + log(cosh(a_y*(y_domain%beg - y_b))) &
                 - 2d0*log(cosh(0.5d0*a_y*(y_b - y_a))))/a_y <= 0d0) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for y_domain%beg, stretch_y, a_y, '// &
                'y_a, and y_b. Exiting ...'
            call s_mpi_abort()
        elseif (stretch_y &
                .and. &
                (a_y + log(cosh(a_y*(y_domain%end - y_a))) &
                 + log(cosh(a_y*(y_domain%end - y_b))) &
                 - 2d0*log(cosh(0.5d0*a_y*(y_b - y_a))))/a_y <= 0d0) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for y_domain%end, stretch_y, a_y, '// &
                'y_a, and y_b. Exiting ...'
            call s_mpi_abort()
        elseif (loops_y < 1) then
            print '(A)', 'Unsupported choice for the value of loops_y. '// &
                'Exiting ...'
            call s_mpi_abort()

            ! Constraints on the grid stretching in the z-direction
        elseif (all(model_eqns /= (/ 2 /))) then
            print '(A)', 'Unsupported value of model_eqns. Exiting ...'
            call s_mpi_abort()
        elseif (num_fluids /= dflt_int &
                .and. &
                (num_fluids < 1 .or. num_fluids > num_fluids)) then
            print '(A)', 'Unsupported value of num_fluids. Exiting ...'
            call s_mpi_abort()
        elseif (model_eqns == 1 .and. adv_alphan) then
            print '(A)', 'Unsupported combination of values of '// &
                'model_eqns and adv_alphan. '// &
                'Exiting ...'
            call s_mpi_abort()
            ! Constraints on the order of the WENO scheme
        elseif (weno_order /= 1 .and. weno_order /= 3 &
                .and. &
                weno_order /= 5) then
            print '(A)', 'Unsupported choice for the value of '// &
                'weno_order. Exiting ...'
            call s_mpi_abort()
        elseif (m + 1 < weno_order) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for m and weno_order. Exiting ...'
            call s_mpi_abort()
        elseif (n > 0 .and. n + 1 < weno_order) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for n and weno_order. Exiting ...'
            call s_mpi_abort()
        elseif ((m + 1)*(n + 1) &
                < &
                weno_order**(min(1, m) + min(1, n))*num_procs) &
            then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for num_procs, m, n, and '// &
                'weno_order. Exiting ...'
            call s_mpi_abort()

            ! Constraints on the boundary conditions in the x-direction
        elseif (bc_x%beg < -12 .or. bc_x%beg > -1) then
            print '(A)', 'Unsupported choice for the value of '// &
                'bc_x%beg. Exiting ...'
            call s_mpi_abort()
        elseif (bc_x%end < -12 .or. bc_x%end > -1) then
            print '(A)', 'Unsupported choice for the value of '// &
                'bc_x%end. Exiting ...'
            call s_mpi_abort()
        elseif ((bc_x%beg == -1 .and. bc_x%end /= -1) &
                .or. &
                (bc_x%end == -1 .and. bc_x%beg /= -1)) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for bc_x%beg and bc_x%end. '// &
                'Exiting ...'
            call s_mpi_abort()
        end if


        ! Constraints on the boundary conditions in the y-direction
        if (bc_y%beg /= dflt_int &
            .and. &
            (bc_y%beg < -12 .or. bc_y%beg > -1)) then
            print '(A)', 'Unsupported choice for the value of '// &
                'bc_y%beg. Exiting ...'
            call s_mpi_abort()
        elseif (bc_y%end /= dflt_int &
                .and. &
                (bc_y%end < -12 .or. bc_y%end > -1)) then
            print '(A)', 'Unsupported choice for the value of '// &
                'bc_y%end. Exiting ...'
            call s_mpi_abort()
        elseif ((n == 0 .and. bc_y%beg /= dflt_int) &
                .or. &
                (n > 0 .and. bc_y%beg == dflt_int)) then
            print '(A)', 'Unsupported choice for the value of n and '// &
                'bc_y%beg. Exiting ...'
            call s_mpi_abort()
        elseif ((n == 0 .and. bc_y%end /= dflt_int) &
                .or. &
                (n > 0 .and. bc_y%end == dflt_int)) then
            print '(A)', 'Unsupported choice for the value of n and '// &
                'bc_y%end. Exiting ...'
            call s_mpi_abort()
        elseif (n > 0 &
                .and. &
                ((bc_y%beg == -1 .and. bc_y%end /= -1) &
                 .or. &
                 (bc_y%end == -1 .and. bc_y%beg /= -1))) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for n, bc_y%beg and bc_y%end. '// &
                'Exiting ...'
            call s_mpi_abort()
        end if

        ! Constraints on number of patches making up the initial condition
        if (num_patches < 0 .or. num_patches > num_patches .or. &
            (num_patches == 0 .and. t_step_old == dflt_int)) then
            print '(A)', 'Unsupported choice for the value of '// &
                'num_patches. Exiting ...'
            call s_mpi_abort()
        end if


        ! check all the patch properties
        call s_check_patches()

        ! Constraints on the stiffened equation of state fluids parameters
        do i = 1, num_fluids

            if (fluid_pp(i)%gamma /= dflt_real &
                .and. &
                fluid_pp(i)%gamma <= 0d0) then
                print '(A,I0,A)', 'Unsupported value of '// &
                    'fluid_pp(', i, ')%'// &
                    'gamma. Exiting ...'
                call s_mpi_abort()
            elseif ((i <= num_fluids .and. fluid_pp(i)%gamma <= 0d0) &
                    .or. &
                    (i > num_fluids .and. fluid_pp(i)%gamma /= dflt_real)) &
                then
                print '(A,I0,A)', 'Unsupported combination '// &
                    'of values of num_fluids '// &
                    'and fluid_pp(', i, ')%'// &
                    'gamma. Exiting ...'
                call s_mpi_abort()
            elseif (fluid_pp(i)%pi_inf /= dflt_real &
                    .and. &
                    fluid_pp(i)%pi_inf < 0d0) then
                print '(A,I0,A)', 'Unsupported value of '// &
                    'fluid_pp(', i, ')%'// &
                    'pi_inf. Exiting ...'
                call s_mpi_abort()
            elseif (model_eqns == 1 &
                    .and. &
                    fluid_pp(i)%pi_inf /= dflt_real) then
                print '(A,I0,A)', 'Unsupported combination '// &
                    'of values of model_eqns '// &
                    'and fluid_pp(', i, ')%'// &
                    'pi_inf. Exiting ...'
                call s_mpi_abort()
            elseif ((i <= num_fluids .and. fluid_pp(i)%pi_inf < 0d0) &
                    .or. &
                    (i > num_fluids .and. fluid_pp(i)%pi_inf /= dflt_real)) &
                then
                print '(A,I0,A)', 'Unsupported combination '// &
                    'of values of num_fluids '// &
                    'and fluid_pp(', i, ')%'// &
                    'pi_inf. Exiting ...'
                call s_mpi_abort()
            end if

        end do

    end subroutine s_check_input_file ! ------------------------------------


    !> The goal of this subroutine is to read in any preexisting
        !!      grid data as well as based on the imported grid, complete
        !!      the necessary global computational domain parameters.
        !! @param dflt_int Default null integer
    subroutine s_read_serial_grid_data_files(dflt_int) ! ---

        integer, intent(IN) :: dflt_int

        ! Generic string used to store the address of a particular file
        character(LEN=len_trim(case_dir) + 3*name_len) :: file_loc

        ! Logical variable used to test the existence of folders
        logical :: dir_check

        ! Generic logical used for the purpose of asserting whether a file
        ! is or is not present in the designated location
        logical :: file_check

        ! Setting address of the local processor rank and time-step directory
        write (proc_rank_dir, '(A,I0)') '/p', proc_rank
        proc_rank_dir = trim(case_dir)//trim(proc_rank_dir)

        write (t_step_dir, '(A,I0)') '/', t_step_old
        t_step_dir = trim(proc_rank_dir)//trim(t_step_dir)

        ! Inquiring as to the existence of the time-step directory
        file_loc = trim(t_step_dir)//'/.'
        call my_inquire(file_loc, dir_check)

        ! If the time-step directory is missing, the pre-process exits
        if (dir_check .neqv. .true.) then
            print '(A)', 'Time-step folder '//trim(t_step_dir)// &
                ' is missing. Exiting ...'
            call s_mpi_abort()
        end if

        ! Reading the Grid Data File for the x-direction ===================

        ! Checking whether x_cb.dat exists
        file_loc = trim(t_step_dir)//'/x_cb.dat'
        inquire (FILE=trim(file_loc), EXIST=file_check)

        ! If it exists, x_cb.dat is read
        if (file_check) then
            open (1, FILE=trim(file_loc), FORM='unformatted', &
                  STATUS='old', ACTION='read')
            read (1) x_cb(-1:m)
            close (1)
        else
            print '(A)', 'File x_cb.dat is missing in '// &
                trim(t_step_dir)//'. Exiting ...'
            call s_mpi_abort()
        end if

        ! Computing cell-center locations
        x_cc(0:m) = (x_cb(0:m) + x_cb(-1:(m - 1)))/2d0

        ! Computing minimum cell-width
        dx = minval(x_cb(0:m) - x_cb(-1:m - 1))
        if (num_procs > 1) call s_mpi_reduce_min(dx)

        ! Setting locations of domain bounds
        x_domain%beg = x_cb(-1)
        x_domain%end = x_cb(m)

        ! ==================================================================

        ! Reading the Grid Data File for the y-direction ===================

        if (n > 0) then

            ! Checking whether y_cb.dat exists
            file_loc = trim(t_step_dir)//'/y_cb.dat'
            inquire (FILE=trim(file_loc), EXIST=file_check)

            ! If it exists, y_cb.dat is read
            if (file_check) then
                open (1, FILE=trim(file_loc), FORM='unformatted', &
                      STATUS='old', ACTION='read')
                read (1) y_cb(-1:n)
                close (1)
            else
                print '(A)', 'File y_cb.dat is missing in '// &
                    trim(t_step_dir)//'. Exiting ...'
                call s_mpi_abort()
            end if

            ! Computing cell-center locations
            y_cc(0:n) = (y_cb(0:n) + y_cb(-1:(n - 1)))/2d0

            ! Computing minimum cell-width
            dy = minval(y_cb(0:n) - y_cb(-1:n - 1))
            if (num_procs > 1) call s_mpi_reduce_min(dy)

            ! Setting locations of domain bounds
            y_domain%beg = y_cb(-1)
            y_domain%end = y_cb(n)

            ! ==================================================================

        end if

        ! ==================================================================

        ! If only the preexisting grid data files are read in and there will
        ! not be any preexisting initial condition data files imported, then
        ! the directory associated with the rank of the local processor may
        ! be cleaned to make room for the new pre-process data. In addition,
        ! the time-step directory that will contain the new grid and initial
        ! condition data are also generated.
        if (old_ic .neqv. .true.) then
            call s_delete_directory(trim(proc_rank_dir)//'/*')
            call s_create_directory(trim(proc_rank_dir)//'/0')
        end if

    end subroutine s_read_serial_grid_data_files ! --------------------------------

    !> Cell-boundary data are checked for consistency by looking
        !!      at the (non-)uniform cell-width distributions for all the
        !!      active coordinate directions and making sure that all of
        !!      the cell-widths are positively valued
    subroutine s_check_grid_data_files() ! -----------------

        ! Cell-boundary Data Consistency Check in x-direction ==============

        if (any(x_cb(0:m) - x_cb(-1:m - 1) <= 0d0)) then
            print '(A)', 'x_cb.dat in '//trim(t_step_dir)// &
                ' contains non-positive cell-spacings. Exiting ...'
            call s_mpi_abort()
        end if

        ! ==================================================================

        ! Cell-boundary Data Consistency Check in y-direction ==============

        if (n > 0) then

            if (any(y_cb(0:n) - y_cb(-1:n - 1) <= 0d0)) then
                print '(A)', 'y_cb.dat in '//trim(t_step_dir)// &
                    ' contains non-positive cell-spacings. '// &
                    'Exiting ...'
                call s_mpi_abort()
            end if

        end if

        ! ==================================================================

    end subroutine s_check_grid_data_files ! -------------------------------

    !> The goal of this subroutine is to read in any preexisting
        !!      initial condition data files so that they may be used by
        !!      the pre-process as a starting point in the creation of an
        !!      all new initial condition.
        !! @param q_cons_vf Conservative variables
    subroutine s_read_serial_ic_data_files(q_cons_vf) ! ---------------------------

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: q_cons_vf

        character(LEN=len_trim(case_dir) + 3*name_len) :: file_loc !<
        ! Generic string used to store the address of a particular file

        character(LEN= &
                  int(floor(log10(real(sys_size, kind(0d0))))) + 1) :: file_num !<
            !! Used to store the variable position, in character form, of the
            !! currently manipulated conservative variable file

        logical :: file_check !<
            !! Generic logical used for the purpose of asserting whether a file
            !! is or is not present in the designated location

        integer :: i !< Generic loop iterator

        ! Reading the Conservative Variables Data Files ====================
        do i = 1, sys_size

            ! Checking whether data file associated with variable position
            ! of the currently manipulated conservative variable exists
            write (file_num, '(I0)') i
            file_loc = trim(t_step_dir)//'/q_cons_vf'// &
                       trim(file_num)//'.dat'
            inquire (FILE=trim(file_loc), EXIST=file_check)

            ! If it exists, the data file is read
            if (file_check) then
                open (1, FILE=trim(file_loc), FORM='unformatted', &
                      STATUS='old', ACTION='read')
                read (1) q_cons_vf(i)%sf
                close (1)
            else
                print '(A)', 'File q_cons_vf'//trim(file_num)// &
                    '.dat is missing in '//trim(t_step_dir)// &
                    '. Exiting ...'
                call s_mpi_abort()
            end if

        end do

        ! ==================================================================

        ! Since the preexisting grid and initial condition data files have
        ! been read in, the directory associated with the rank of the local
        ! process may be cleaned out to make room for new pre-process data.
        ! In addition, the time-step folder that will contain the new grid
        ! and initial condition data are also generated.
        call s_create_directory(trim(proc_rank_dir)//'/*')
        call s_create_directory(trim(proc_rank_dir)//'/0')

    end subroutine s_read_serial_ic_data_files ! ----------------------------------

    !> Cell-boundary data are checked for consistency by looking
        !!      at the (non-)uniform cell-width distributions for all the
        !!      active coordinate directions and making sure that all of
        !!      the cell-widths are positively valued
        !! @param dflt_int Default null integer
    subroutine s_read_parallel_grid_data_files(dflt_int)

        integer, intent(IN) :: dflt_int

#ifdef MFC_MPI

        real(kind(0d0)), allocatable, dimension(:) :: x_cb_glb, y_cb_glb

        integer :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE) :: status

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist

        allocate (x_cb_glb(-1:m_glb))
        allocate (y_cb_glb(-1:n_glb))

        ! Read in cell boundary locations in x-direction
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'x_cb.dat'
        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (file_exist) then
            data_size = m_glb + 2
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
            call MPI_FILE_READ_ALL(ifile, x_cb_glb, data_size, MPI_DOUBLE_PRECISION, status, ierr)
            call MPI_FILE_CLOSE(ifile, ierr)
        else
            print '(A)', 'File ', trim(file_loc), ' is missing. Exiting... '
            call s_mpi_abort()
        end if

        ! Assigning local cell boundary locations
        x_cb(-1:m) = x_cb_glb((start_idx(1) - 1):(start_idx(1) + m))
        ! Computing cell center locations
        x_cc(0:m) = (x_cb(0:m) + x_cb(-1:(m - 1)))/2d0
        ! Computing minimum cell width
        dx = minval(x_cb(0:m) - x_cb(-1:(m - 1)))
        if (num_procs > 1) call s_mpi_reduce_min(dx)
        ! Setting locations of domain bounds
        x_domain%beg = x_cb(-1)
        x_domain%end = x_cb(m)

        if (n > 0) then
            ! Read in cell boundary locations in y-direction
            file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'y_cb.dat'
            inquire (FILE=trim(file_loc), EXIST=file_exist)

            if (file_exist) then
                data_size = n_glb + 2
                call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                call MPI_FILE_READ_ALL(ifile, y_cb_glb, data_size, MPI_DOUBLE_PRECISION, status, ierr)
                call MPI_FILE_CLOSE(ifile, ierr)
            else
                print '(A)', 'File ', trim(file_loc), ' is missing. Exiting... '
                call s_mpi_abort()
            end if

            ! Assigning local cell boundary locations
            y_cb(-1:n) = y_cb_glb((start_idx(2) - 1):(start_idx(2) + n))
            ! Computing cell center locations
            y_cc(0:n) = (y_cb(0:n) + y_cb(-1:(n - 1)))/2d0
            ! Computing minimum cell width
            dy = minval(y_cb(0:n) - y_cb(-1:(n - 1)))
            if (num_procs > 1) call s_mpi_reduce_min(dy)
            ! Setting locations of domain bounds
            y_domain%beg = y_cb(-1)
            y_domain%end = y_cb(n)

        end if

        deallocate (x_cb_glb, y_cb_glb)

#endif

    end subroutine s_read_parallel_grid_data_files ! -----------------------

    !> The goal of this subroutine is to read in any preexisting
        !!      initial condition data files so that they may be used by
        !!      the pre-process as a starting point in the creation of an
        !!      all new initial condition.
        !! @param q_cons_vf Conservative variables
    subroutine s_read_parallel_ic_data_files(q_cons_vf) ! ------------------

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: q_cons_vf

#ifdef MFC_MPI

        integer :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(KIND=MPI_OFFSET_KIND) :: disp
        integer(KIND=MPI_OFFSET_KIND) :: m_MOK, n_MOK
        integer(KIND=MPI_OFFSET_KIND) :: WP_MOK, var_MOK, str_MOK
        integer(KIND=MPI_OFFSET_KIND) :: NVARS_MOK
        integer(KIND=MPI_OFFSET_KIND) :: MOK

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist

        integer :: i

        ! Open the file to read
        write (file_loc, '(I0,A)') t_step_old, '.dat'
        file_loc = trim(restart_dir)//trim(mpiiofs)//trim(file_loc)
        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (file_exist) then
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

            ! Initialize MPI data I/O
            call s_initialize_mpi_data(q_cons_vf)

            ! Size of local arrays
            data_size = (m + 1)*(n + 1)

            ! Resize some integers so MPI can read even the biggest files
            m_MOK = int(m_glb + 1, MPI_OFFSET_KIND)
            n_MOK = int(n_glb + 1, MPI_OFFSET_KIND)
            WP_MOK = int(8d0, MPI_OFFSET_KIND)
            MOK = int(1d0, MPI_OFFSET_KIND)
            str_MOK = int(name_len, MPI_OFFSET_KIND)
            NVARS_MOK = int(sys_size, MPI_OFFSET_KIND)

            ! Read the data for each variable
            do i = 1, sys_size
                var_MOK = int(i, MPI_OFFSET_KIND)

                ! Initial displacement to skip at beginning of file
                disp = m_MOK*max(MOK, n_MOK)*WP_MOK*(var_MOK - 1)

                call MPI_FILE_SET_VIEW(ifile, disp, MPI_DOUBLE_PRECISION, MPI_IO_DATA%view(i), &
                                       'native', mpi_info_int, ierr)
                call MPI_FILE_READ(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                   MPI_DOUBLE_PRECISION, status, ierr)
            end do

            call s_mpi_barrier()

            call MPI_FILE_CLOSE(ifile, ierr)

        else
            print '(A)', 'File ', trim(file_loc), ' is missing. Exiting... '
            call s_mpi_abort()
        end if
        call s_mpi_barrier()
        if (proc_rank == 0) call s_create_directory(trim(file_loc))

#endif

    end subroutine s_read_parallel_ic_data_files ! -------------------------

    subroutine s_initialize_start_up_module() !-----------------------------

        if (parallel_io .neqv. .true.) then
            s_read_grid_data_files => s_read_serial_grid_data_files
            s_read_ic_data_files => s_read_serial_ic_data_files
        else
            s_read_grid_data_files => s_read_parallel_grid_data_files
            s_read_ic_data_files => s_read_parallel_ic_data_files
        end if

    end subroutine s_initialize_start_up_module ! --------------------------

    subroutine s_finalize_start_up_module() ! ------------------------------

        s_read_grid_data_files => null()
        s_read_ic_data_files => null()

    end subroutine s_finalize_start_up_module ! ----------------------------

end module m_start_up
