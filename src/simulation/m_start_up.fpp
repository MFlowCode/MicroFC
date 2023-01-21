!>
!! @file m_start_up.f90
!! @brief Contains module m_start_up

#:include 'case.fpp'

!> @brief The purpose of the module is primarily to read in the files that
!!              contain the inputs, the initial condition data and the grid data
!!              that are provided by the user. The module is additionally tasked
!!              with verifying the consistency of the user inputs and completing
!!              the grid variablesThe purpose of the module is primarily to read
!!              in the files that
!!              contain the inputs, the initial condition data and the grid data
!!              that are provided by the user. The module is additionally tasked
!!              with verifying the consistency of the user inputs and completing
!!              the grid variables.
module m_start_up

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_compile_specific
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_start_up_module, &
 s_read_input_file, &
 s_check_input_file, &
 s_read_data_files, &
 s_read_serial_data_files, &
 s_read_parallel_data_files, &
 s_populate_grid_variables_buffers, &
 s_finalize_start_up_module

    abstract interface ! ===================================================

        !! @param q_cons_vf  Conservative variables
        subroutine s_read_abstract_data_files(q_cons_vf) ! -----------

            import :: scalar_field, sys_size

            type(scalar_field), &
                dimension(sys_size), &
                intent(INOUT) :: q_cons_vf

        end subroutine s_read_abstract_data_files ! -----------------

    end interface ! ========================================================

    type(scalar_field), allocatable, dimension(:) :: grad_x_vf, grad_y_vf, norm_vf

    procedure(s_read_abstract_data_files), pointer :: s_read_data_files => null()

contains

    !>  The purpose of this procedure is to first verify that an
        !!      input file has been made available by the user. Provided
        !!      that this is so, the input file is then read in.
    subroutine s_read_input_file() ! ---------------------------------------

        ! Relative path to the input file provided by the user
        character(LEN=name_len) :: file_path = './simulation.inp'

        logical :: file_exist !<
            !! Logical used to check the existence of the input file

        ! Namelist of the global parameters which may be specified by user
        namelist /user_inputs/ case_dir, run_time_info, m, n, dt, &
            t_step_start, t_step_stop, t_step_save, &
            num_fluids, adv_alphan, &
            time_stepper, &
            weno_eps, weno_flat, riemann_flat, cu_mpi, &
            wave_speeds, avg_state, &
            bc_x, bc_y, &
            fluid_pp, probe_wrt, &
            fd_order, probe, num_probes, t_step_old, &
            weno_Re_flux, &
#:if not MFC_CASE_OPTIMIZATION
            weno_order, &
#:endif
            precision, parallel_io 


        ! Checking that an input file has been provided by the user. If it
        ! has, then the input file is read in, otherwise, simulation exits.
        inquire (FILE=trim(file_path), EXIST=file_exist)

        if (file_exist) then
            open (1, FILE=trim(file_path), &
                  FORM='formatted', &
                  ACTION='read', &
                  STATUS='old')
            read (1, NML=user_inputs); close (1)

            ! Store BC information into global BC
            bc_x_glb%beg = bc_x%beg
            bc_x_glb%end = bc_x%end
            bc_y_glb%beg = bc_y%beg
            bc_y_glb%end = bc_y%end

            ! Store m,n,p into global m,n,p
            m_glb = m
            n_glb = n

        else
            print '(A)', trim(file_path)//' is missing. Exiting ...'
            call s_mpi_abort()
        end if

    end subroutine s_read_input_file ! -------------------------------------

    !> The goal of this procedure is to verify that each of the
    !!      user provided inputs is valid and that their combination
    !!      consitutes a meaningful configuration for the simulation.
    subroutine s_check_input_file() ! --------------------------------------

        ! Relative path to the current directory file in the case directory
        character(LEN=path_len) :: file_path

        ! Logical used to check the existence of the current directory file
        logical :: file_exist

        ! Generic loop iterators
        integer :: i, j


        ! Logistics ========================================================
        file_path = trim(case_dir)//'/.'

        call my_inquire(file_path, file_exist)

        if (file_exist .neqv. .true.) then
            print '(A)', trim(file_path)//' is missing. Exiting ...'
            call s_mpi_abort()
        end if
        ! ==================================================================

#if !(defined(_OPENACC) && defined(__PGI))
        if (cu_mpi) then
            print '(A)', 'Unsupported value of cu_mpi. Exiting ...'
            call s_mpi_abort()
        end if
#endif

        ! Computational Domain Parameters ==================================
        if (m <= 0) then
            print '(A)', 'Unsupported value of m. Exiting ...'
            call s_mpi_abort()
        elseif (n < 0) then
            print '(A)', 'Unsupported value of n. Exiting ...'
            call s_mpi_abort()
        elseif (dt <= 0) then
            print '(A)', 'Unsupported value of dt. Exiting ...'
            call s_mpi_abort()
        elseif (t_step_start < 0) then
            print '(A)', 'Unsupported value of t_step_start. Exiting ...'
            call s_mpi_abort()
        elseif (t_step_stop <= t_step_start) then
            print '(A)', 'Unsupported combination of values of '// &
                't_step_start and t_step_stop. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (t_step_save > t_step_stop - t_step_start) then
            print '(A)', 'Unsupported combination of values of '// &
                't_step_start, t_step_stop and '// &
                't_step_save. Exiting ...'
            call s_mpi_abort()
        end if
        ! ==================================================================

        ! Simulation Algorithm Parameters ==================================
        if (num_fluids /= dflt_int &
                .and. &
                (num_fluids < 1 .or. num_fluids > num_fluids)) then
            print '(A)', 'Unsupported value of num_fluids. Exiting ...'
            call s_mpi_abort()
        elseif (time_stepper < 1 .or. time_stepper > 5) then
            if (time_stepper /= 23) then
                print '(A)', 'Unsupported value of time_stepper. Exiting ...'
                call s_mpi_abort()
            end if
        elseif (all(weno_order /= (/1, 3, 5/))) then
            print '(A)', 'Unsupported value of weno_order. Exiting ...'
            call s_mpi_abort()
        elseif (m + 1 < num_stcls_min*weno_order) then
            print '(A)', 'Unsupported combination of values of '// &
                'm and weno_order. Exiting ...'
            call s_mpi_abort()
        elseif (n + 1 < min(1, n)*num_stcls_min*weno_order) then
            print '(A)', 'Unsupported combination of values of '// &
                'n and weno_order. Exiting ...'
            call s_mpi_abort()
        elseif (weno_eps <= 0d0 .or. weno_eps > 1d-6) then
            print '(A)', 'Unsupported value of weno_eps. Exiting ...'
            call s_mpi_abort()
        elseif (all(wave_speeds /= (/dflt_int, 1, 2/))) then
            print '(A)', 'Unsupported value of wave_speeds. Exiting ...'
            call s_mpi_abort()
        elseif (all(avg_state /= (/dflt_int, 1, 2/))) then
            print '(A)', 'Unsupported value of avg_state. Exiting ...'
            call s_mpi_abort()
        elseif (bc_x%beg < -12 .or. bc_x%beg > -1) then
            print '(A)', 'Unsupported value of bc_x%beg. Exiting ...'
            call s_mpi_abort()
        elseif (bc_x%end < -12 .or. bc_x%end > -1) then
            print '(A)', 'Unsupported value of bc_x%end. Exiting ...'
            call s_mpi_abort()
        elseif ((bc_x%beg == -1 .and. bc_x%end /= -1) &
                .or. &
                (bc_x%end == -1 .and. bc_x%beg /= -1)) then
            print '(A)', 'Unsupported combination of values of '// &
                'bc_x%beg and bc_x%end. Exiting ...'
            call s_mpi_abort()
        elseif (bc_y%end /= dflt_int &
                .and. &
                (bc_y%end < -12 .or. bc_y%end > -1)) then
            print '(A)', 'Unsupported value of bc_y%end. Exiting ...'
            call s_mpi_abort()
        elseif ((n == 0 .and. bc_y%beg /= dflt_int) &
                .or. &
                (n > 0 .and. bc_y%beg == dflt_int)) then
            print '(A)', 'Unsupported combination of values of '// &
                'n and bc_y%beg. Exiting ...'
            call s_mpi_abort()
        elseif ((n == 0 .and. bc_y%end /= dflt_int) &
                .or. &
                (n > 0 .and. bc_y%end == dflt_int)) then
            print '(A)', 'Unsupported combination of values of '// &
                'n and bc_y%end. Exiting ...'
            call s_mpi_abort()
        elseif ((bc_y%beg == -1 .and. bc_y%end /= -1) &
                .or. &
                (bc_y%end == -1 .and. bc_y%beg /= -1)) then
            print '(A)', 'Unsupported combination of values of '// &
                'bc_y%beg and bc_y%end. Exiting ...'
            call s_mpi_abort()
        end if
        ! END: Simulation Algorithm Parameters =============================

        ! Finite Difference Parameters =====================================
        if (fd_order /= dflt_int &
            .and. &
            fd_order /= 1 .and. fd_order /= 2 .and. fd_order /= 4) then
            print '(A)', 'Unsupported choice for the value of '// &
                'fd_order. Exiting ...'
            call s_mpi_abort()
        elseif (probe_wrt .and. fd_order == dflt_int) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for probe_wrt, and fd_order. '// &
                'Exiting ...'
            call s_mpi_abort()
        end if
        ! END: Finite Difference Parameters ================================

        ! Fluids Physical Parameters =======================================
        do i = 1, num_fluids

            if (fluid_pp(i)%gamma /= dflt_real &
                .and. &
                fluid_pp(i)%gamma <= 0d0) then
                print '(A,I0,A)', 'Unsupported value of '// &
                    'fluid_pp(', i, ')%'// &
                    'gamma. Exiting ...'
                call s_mpi_abort()
            elseif (fluid_pp(i)%pi_inf /= dflt_real &
                    .and. &
                    fluid_pp(i)%pi_inf < 0d0) then
                print '(A,I0,A)', 'Unsupported value of '// &
                    'fluid_pp(', i, ')%'// &
                    'pi_inf. Exiting ...'
                call s_mpi_abort()
            end if

            do j = 1, 2
                if (fluid_pp(i)%Re(j) /= dflt_real &
                    .and. &
                    fluid_pp(i)%Re(j) <= 0d0) then
                    print '(A,I0,A,I0,A)', 'Unsupported value of '// &
                        'fluid_pp(', i, ')%'// &
                        'Re(', j, '). Exiting ...'
                    call s_mpi_abort()
                end if

                if (i > num_fluids &
                    .and. &
                    fluid_pp(i)%Re(j) /= dflt_real) then
                    print '(A,I0,A,I0,A)', 'Unsupported combination '// &
                        'of values of num_fluids '// &
                        'and fluid_pp(', i, ')%'// &
                        'Re(', j, '). Exiting ...'
                    call s_mpi_abort()
                end if

            end do

        end do
        ! END: Fluids Physical Parameters ==================================

    end subroutine s_check_input_file ! ------------------------------------

        !!              initial condition and grid data files. The cell-average
        !!              conservative variables constitute the former, while the
        !!              cell-boundary locations in x-, y- and z-directions make
        !!              up the latter. This procedure also calculates the cell-
        !!              width distributions from the cell-boundary locations.
        !! @param q_cons_vf Cell-averaged conservative variables
    subroutine s_read_serial_data_files(q_cons_vf) ! ------------------------------

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf

        character(LEN=path_len + 2*name_len) :: t_step_dir !<
            !! Relative path to the starting time-step directory

        character(LEN=path_len + 3*name_len) :: file_path !<
            !! Relative path to the grid and conservative variables data files

        logical :: file_exist !<
        ! Logical used to check the existence of the data files

        integer :: i !< Generic loop iterator

        ! Confirming that the directory from which the initial condition and
        ! the grid data files are to be read in exists and exiting otherwise
        write (t_step_dir, '(A,I0,A,I0)') &
            trim(case_dir)//'/p_all/p', proc_rank, '/', t_step_start

        file_path = trim(t_step_dir)//'/.'
        call my_inquire(file_path, file_exist)

        if (file_exist .neqv. .true.) then
            print '(A)', trim(file_path)//' is missing. Exiting ...'
            call s_mpi_abort()
        end if

        ! Cell-boundary Locations in x-direction ===========================
        file_path = trim(t_step_dir)//'/x_cb.dat'

        inquire (FILE=trim(file_path), EXIST=file_exist)

        if (file_exist) then
            open (2, FILE=trim(file_path), &
                  FORM='unformatted', &
                  ACTION='read', &
                  STATUS='old')
            read (2) x_cb(-1:m); close (2)
        else
            print '(A)', trim(file_path)//' is missing. Exiting ...'
            call s_mpi_abort()
        end if

        dx(0:m) = x_cb(0:m) - x_cb(-1:m - 1)
        x_cc(0:m) = x_cb(-1:m - 1) + dx(0:m)/2d0
        ! ==================================================================

        ! Cell-boundary Locations in y-direction ===========================
        if (n > 0) then

            file_path = trim(t_step_dir)//'/y_cb.dat'

            inquire (FILE=trim(file_path), EXIST=file_exist)

            if (file_exist) then
                open (2, FILE=trim(file_path), &
                      FORM='unformatted', &
                      ACTION='read', &
                      STATUS='old')
                read (2) y_cb(-1:n); close (2)
            else
                print '(A)', trim(file_path)//' is missing. Exiting ...'
                call s_mpi_abort()
            end if

            dy(0:n) = y_cb(0:n) - y_cb(-1:n - 1)
            y_cc(0:n) = y_cb(-1:n - 1) + dy(0:n)/2d0

        end if
        ! ==================================================================

        ! ==================================================================

        ! Cell-average Conservative Variables ==============================
        do i = 1, adv_idx%end
            write (file_path, '(A,I0,A)') &
                trim(t_step_dir)//'/q_cons_vf', i, '.dat'
            inquire (FILE=trim(file_path), EXIST=file_exist)
            if (file_exist) then
                open (2, FILE=trim(file_path), &
                      FORM='unformatted', &
                      ACTION='read', &
                      STATUS='old')
                read (2) q_cons_vf(i)%sf(0:m, 0:n); close (2)
            else
                print '(A)', trim(file_path)//' is missing. Exiting ...'
                call s_mpi_abort()
            end if
        end do
        ! ==================================================================

    end subroutine s_read_serial_data_files ! -------------------------------------

        !! @param q_cons_vf Conservative variables
    subroutine s_read_parallel_data_files(q_cons_vf) ! ---------------------------

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: q_cons_vf

#ifdef MFC_MPI

        real(kind(0d0)), allocatable, dimension(:) :: x_cb_glb, y_cb_glb

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

        allocate (x_cb_glb(-1:m_glb))
        allocate (y_cb_glb(-1:n_glb))

        ! Read in cell boundary locations in x-direction
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'x_cb.dat'
        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (file_exist) then
            data_size = m_glb + 2
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
            call MPI_FILE_READ(ifile, x_cb_glb, data_size, MPI_DOUBLE_PRECISION, status, ierr)
            call MPI_FILE_CLOSE(ifile, ierr)
        else
            print '(A)', 'File ', trim(file_loc), ' is missing. Exiting...'
            call s_mpi_abort()
        end if

        ! Assigning local cell boundary locations
        x_cb(-1:m) = x_cb_glb((start_idx(1) - 1):(start_idx(1) + m))
        ! Computing the cell width distribution
        dx(0:m) = x_cb(0:m) - x_cb(-1:m - 1)
        ! Computing the cell center locations
        x_cc(0:m) = x_cb(-1:m - 1) + dx(0:m)/2d0

        if (n > 0) then
            ! Read in cell boundary locations in y-direction
            file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'y_cb.dat'
            inquire (FILE=trim(file_loc), EXIST=file_exist)

            if (file_exist) then
                data_size = n_glb + 2
                call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                call MPI_FILE_READ(ifile, y_cb_glb, data_size, MPI_DOUBLE_PRECISION, status, ierr)
                call MPI_FILE_CLOSE(ifile, ierr)
            else
                print '(A)', 'File ', trim(file_loc), ' is missing. Exiting...'
                call s_mpi_abort()
            end if

            ! Assigning local cell boundary locations
            y_cb(-1:n) = y_cb_glb((start_idx(2) - 1):(start_idx(2) + n))
            ! Computing the cell width distribution
            dy(0:n) = y_cb(0:n) - y_cb(-1:n - 1)
            ! Computing the cell center locations
            y_cc(0:n) = y_cb(-1:n - 1) + dy(0:n)/2d0

        end if

        ! Open the file to read conservative variables
        write (file_loc, '(I0,A)') t_step_start, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (file_exist) then
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

            ! Initialize MPI data I/O
            call s_initialize_mpi_data(q_cons_vf)

            ! Size of local arrays
            data_size = (m + 1)*(n + 1)

            ! Resize some integers so MPI can read even the biggest file
            m_MOK = int(m_glb + 1, MPI_OFFSET_KIND)
            n_MOK = int(n_glb + 1, MPI_OFFSET_KIND)
            WP_MOK = int(8d0, MPI_OFFSET_KIND)
            MOK = int(1d0, MPI_OFFSET_KIND)
            str_MOK = int(name_len, MPI_OFFSET_KIND)
            NVARS_MOK = int(sys_size, MPI_OFFSET_KIND)

            ! Read the data for each variable
            do i = 1, adv_idx%end
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
            print '(A)', 'File ', trim(file_loc), ' is missing. Exiting...'
            call s_mpi_abort()
        end if

        deallocate (x_cb_glb, y_cb_glb)

#endif

    end subroutine s_read_parallel_data_files ! -------------------------------

    !> The purpose of this subroutine is to populate the buffers
        !!          of the grid variables, which are constituted of the cell-
        !!          boundary locations and cell-width distributions, based on
        !!          the boundary conditions.
    subroutine s_populate_grid_variables_buffers() ! -----------------------

        integer :: i !< Generic loop iterator

        ! Population of Buffers in x-direction =============================

        ! Populating cell-width distribution buffer, at the beginning of the
        ! coordinate direction, based on the selected boundary condition. In
        ! order, these are the ghost-cell extrapolation, symmetry, periodic,
        ! and processor boundary conditions.
        if (bc_x%beg <= -3) then
            do i = 1, buff_size
                dx(-i) = dx(0)
            end do
        elseif (bc_x%beg == -2) then
            do i = 1, buff_size
                dx(-i) = dx(i - 1)
            end do
        elseif (bc_x%beg == -1) then
            do i = 1, buff_size
                dx(-i) = dx(m - (i - 1))
            end do
        else
            call s_mpi_sendrecv_grid_variables_buffers(1, -1)
        end if

        ! Computing the cell-boundary locations buffer, at the beginning of
        ! the coordinate direction, from the cell-width distribution buffer
        do i = 1, buff_size
            x_cb(-1 - i) = x_cb(-i) - dx(-i)
        end do
        ! Computing the cell-center locations buffer, at the beginning of
        ! the coordinate direction, from the cell-width distribution buffer
        do i = 1, buff_size
            x_cc(-i) = x_cc(1 - i) - (dx(1 - i) + dx(-i))/2d0
        end do

        ! Populating the cell-width distribution buffer, at the end of the
        ! coordinate direction, based on desired boundary condition. These
        ! include, in order, ghost-cell extrapolation, symmetry, periodic,
        ! and processor boundary conditions.
        if (bc_x%end <= -3) then
            do i = 1, buff_size
                dx(m + i) = dx(m)
            end do
        elseif (bc_x%end == -2) then
            do i = 1, buff_size
                dx(m + i) = dx(m - (i - 1))
            end do
        elseif (bc_x%end == -1) then
            do i = 1, buff_size
                dx(m + i) = dx(i - 1)
            end do
        else
            call s_mpi_sendrecv_grid_variables_buffers(1, 1)
        end if

        ! Populating the cell-boundary locations buffer, at the end of the
        ! coordinate direction, from buffer of the cell-width distribution
        do i = 1, buff_size
            x_cb(m + i) = x_cb(m + (i - 1)) + dx(m + i)
        end do
        ! Populating the cell-center locations buffer, at the end of the
        ! coordinate direction, from buffer of the cell-width distribution
        do i = 1, buff_size
            x_cc(m + i) = x_cc(m + (i - 1)) + (dx(m + (i - 1)) + dx(m + i))/2d0
        end do

        ! END: Population of Buffers in x-direction ========================

        ! Population of Buffers in y-direction =============================

        ! Populating cell-width distribution buffer, at the beginning of the
        ! coordinate direction, based on the selected boundary condition. In
        ! order, these are the ghost-cell extrapolation, symmetry, periodic,
        ! and processor boundary conditions.
        if (n == 0) then
            return
        elseif (bc_y%beg <= -3 .and. bc_y%beg /= -13) then
            do i = 1, buff_size
                dy(-i) = dy(0)
            end do
        elseif (bc_y%beg == -2 .or. bc_y%beg == -13) then
            do i = 1, buff_size
                dy(-i) = dy(i - 1)
            end do
        elseif (bc_y%beg == -1) then
            do i = 1, buff_size
                dy(-i) = dy(n - (i - 1))
            end do
        else
            call s_mpi_sendrecv_grid_variables_buffers(2, -1)
        end if

        ! Computing the cell-boundary locations buffer, at the beginning of
        ! the coordinate direction, from the cell-width distribution buffer
        do i = 1, buff_size
            y_cb(-1 - i) = y_cb(-i) - dy(-i)
        end do
        ! Computing the cell-center locations buffer, at the beginning of
        ! the coordinate direction, from the cell-width distribution buffer
        do i = 1, buff_size
            y_cc(-i) = y_cc(1 - i) - (dy(1 - i) + dy(-i))/2d0
        end do

        ! Populating the cell-width distribution buffer, at the end of the
        ! coordinate direction, based on desired boundary condition. These
        ! include, in order, ghost-cell extrapolation, symmetry, periodic,
        ! and processor boundary conditions.
        if (bc_y%end <= -3) then
            do i = 1, buff_size
                dy(n + i) = dy(n)
            end do
        elseif (bc_y%end == -2) then
            do i = 1, buff_size
                dy(n + i) = dy(n - (i - 1))
            end do
        elseif (bc_y%end == -1) then
            do i = 1, buff_size
                dy(n + i) = dy(i - 1)
            end do
        else
            call s_mpi_sendrecv_grid_variables_buffers(2, 1)
        end if

        ! Populating the cell-boundary locations buffer, at the end of the
        ! coordinate direction, from buffer of the cell-width distribution
        do i = 1, buff_size
            y_cb(n + i) = y_cb(n + (i - 1)) + dy(n + i)
        end do
        ! Populating the cell-center locations buffer, at the end of the
        ! coordinate direction, from buffer of the cell-width distribution
        do i = 1, buff_size
            y_cc(n + i) = y_cc(n + (i - 1)) + (dy(n + (i - 1)) + dy(n + i))/2d0
        end do

        ! END: Population of Buffers in y-direction ========================

    end subroutine s_populate_grid_variables_buffers ! ---------------------

    subroutine s_initialize_start_up_module() !-----------------------------

        type(int_bounds_info) :: ix, iy, iz

        integer :: i !< Generic loop iterator

        if (parallel_io .neqv. .true.) then
            s_read_data_files => s_read_serial_data_files
        else
            s_read_data_files => s_read_parallel_data_files
        end if

    end subroutine s_initialize_start_up_module ! --------------------------

    subroutine s_finalize_start_up_module() ! ------------------------------

        integer :: i !< Generic loop interator

        s_read_data_files => null()

    end subroutine s_finalize_start_up_module ! ----------------------------

end module m_start_up
