!>
!! @file m_data_output.f90
!! @brief Contains module m_data_output

#:include 'macros.fpp'

!> @brief The primary purpose of this module is to output the grid and the
!!              conservative variables data at the chosen time-step interval. In
!!              addition, this module is also in charge of outputting a run-time
!!              information file which summarizes the time-dependent behavior !of
!!              the stability criteria. The latter include the inviscid Courant–
!!              Friedrichs–Lewy (ICFL), viscous CFL (VCFL), capillary CFL (CCFL)
!!              and cell Reynolds (Rc) numbers.
module m_data_output

    !  Dependencies ============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_compile_specific
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_data_output_module, &
 s_open_run_time_information_file, &
 s_open_probe_files, &
 s_write_run_time_information, &
 s_write_data_files, &
 s_write_serial_data_files, &
 s_write_parallel_data_files, &
 s_write_probe_files, &
 s_close_run_time_information_file, &
 s_close_probe_files, &
 s_finalize_data_output_module

    abstract interface ! ===================================================

        !> Write data files
        !! @param q_cons_vf Conservative variables
        !! @param t_step Current time step
        subroutine s_write_abstract_data_files(q_cons_vf, t_step)

            import :: scalar_field, sys_size

            type(scalar_field), &
                dimension(sys_size), &
                intent(IN) :: q_cons_vf

            integer, intent(IN) :: t_step

        end subroutine s_write_abstract_data_files ! -------------------
    end interface ! ========================================================

    real(kind(0d0)), allocatable, dimension(:, :) :: icfl_sf  !< ICFL stability criterion
    real(kind(0d0)), allocatable, dimension(:, :) :: vcfl_sf  !< VCFL stability criterion
    real(kind(0d0)), allocatable, dimension(:, :) :: ccfl_sf  !< CCFL stability criterion
    real(kind(0d0)), allocatable, dimension(:, :) :: Rc_sf  !< Rc stability criterion

!$acc declare create(icfl_sf, vcfl_sf, ccfl_sf, Rc_sf)

    real(kind(0d0)) :: icfl_max_loc, icfl_max_glb !< ICFL stability extrema on local and global grids
    real(kind(0d0)) :: vcfl_max_loc, vcfl_max_glb !< VCFL stability extrema on local and global grids
    real(kind(0d0)) :: ccfl_max_loc, ccfl_max_glb !< CCFL stability extrema on local and global grids
    real(kind(0d0)) :: Rc_min_loc, Rc_min_glb !< Rc   stability extrema on local and global grids

!$acc declare create(icfl_max_loc, icfl_max_glb, vcfl_max_loc, vcfl_max_glb, ccfl_max_loc, ccfl_max_glb, Rc_min_loc, Rc_min_glb)

    !> @name ICFL, VCFL, CCFL and Rc stability criteria extrema over all the time-steps
    !> @{
    real(kind(0d0)) :: icfl_max !< ICFL criterion maximum
    real(kind(0d0)) :: vcfl_max !< VCFL criterion maximum
    real(kind(0d0)) :: ccfl_max !< CCFL criterion maximum
    real(kind(0d0)) :: Rc_min !< Rc criterion maximum
    !> @}


    procedure(s_write_abstract_data_files), pointer :: s_write_data_files => null()

contains

    !>  The purpose of this subroutine is to open a new or pre-
        !!          existing run-time information file and append to it the
        !!      basic header information relevant to current simulation.
        !!      In general, this requires generating a table header for
        !!      those stability criteria which will be written at every
        !!      time-step.
    subroutine s_open_run_time_information_file() ! ------------------------

        character(LEN=name_len) :: file_name = 'run_time.inf' !<
            !! Name of the run-time information file

        character(LEN=path_len + name_len) :: file_path !<
            !! Relative path to a file in the case directory

        character(LEN=8) :: file_date !<
            !! Creation date of the run-time information file

        logical :: file_exist !<
            !! Logical used to check existence of run-time information file

        ! Opening the run-time information file
        file_path = trim(case_dir)//'/'//trim(file_name)

        inquire (FILE=trim(file_path), EXIST=file_exist)

        open (1, FILE=trim(file_path), &
              FORM='formatted', &
              POSITION='append', &
              STATUS='unknown')

        ! Generating file header for a new run-time information file
        if (file_exist .neqv. .true.) then

            write (1, '(A)') 'Description: Stability information at '// &
                'each time-step of the simulation. This'
            write (1, '(13X,A)') 'data is composed of the inviscid '// &
                'Courant–Friedrichs–Lewy (ICFL)'
            write (1, '(13X,A)') 'number, the viscous CFL (VCFL) number, '// &
                'the capillary CFL (CCFL)'
            write (1, '(13X,A)') 'number and the cell Reynolds (Rc) '// &
                'number. Please note that only'
            write (1, '(13X,A)') 'those stability conditions pertinent '// &
                'to the physics included in'
            write (1, '(13X,A)') 'the current computation are displayed.'

            call date_and_time(DATE=file_date)

            write (1, '(A)') 'Date: '//file_date(5:6)//'/'// &
                file_date(7:8)//'/'// &
                file_date(3:4)

        end if

        write (1, '(A)') ''; write (1, '(A)') ''

        ! Generating table header for the stability criteria to be outputted
        if (any(Re_size > 0)) then
            write (1, '(A)') '==== Time-steps ====== Time ======= ICFL '// &
                'Max ==== VCFL Max ====== Rc Min ======='
        else
            write (1, '(A)') '=========== Time-steps ============== Time '// &
                '============== ICFL Max ============='
        end if

    end subroutine s_open_run_time_information_file ! ----------------------

    !>  This opens a formatted data file where the root processor
        !!      can write out flow probe information
    subroutine s_open_probe_files() ! --------------------------------------

        character(LEN=path_len + 3*name_len) :: file_path !<
            !! Relative path to the probe data file in the case directory

        integer :: i !< Generic loop iterator

        do i = 1, num_probes
            ! Generating the relative path to the data file
            write (file_path, '(A,I0,A)') '/D/probe', i, '_prim.dat'
            file_path = trim(case_dir)//trim(file_path)

            ! Creating the formatted data file and setting up its
            ! structure
            open (i + 30, FILE=trim(file_path), &
                  FORM='formatted', &
                  STATUS='unknown')
        end do

    end subroutine s_open_probe_files ! ------------------------------------

    !>  The goal of the procedure is to output to the run-time
        !!      information file the stability criteria extrema in the
        !!      entire computational domain and at the given time-step.
        !!      Moreover, the subroutine is also in charge of tracking
        !!      these stability criteria extrema over all time-steps.
        !!  @param q_prim_vf Cell-average primitive variables
        !!  @param t_step Current time step
    subroutine s_write_run_time_information(q_prim_vf, t_step) ! -----------

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        integer, intent(IN) :: t_step

        real(kind(0d0)), dimension(num_fluids) :: alpha_rho  !< Cell-avg. partial density
        real(kind(0d0)) :: rho        !< Cell-avg. density
        real(kind(0d0)), dimension(num_dims) :: vel        !< Cell-avg. velocity
        real(kind(0d0)) :: pres       !< Cell-avg. pressure
        real(kind(0d0)), dimension(num_fluids) :: alpha      !< Cell-avg. volume fraction
        real(kind(0d0)) :: gamma      !< Cell-avg. sp. heat ratio
        real(kind(0d0)) :: pi_inf     !< Cell-avg. liquid stiffness function
        real(kind(0d0)) :: c          !< Cell-avg. sound speed
        real(kind(0d0)), dimension(2) :: Re         !< Cell-avg. Reynolds numbers

        ! ICFL, VCFL, CCFL and Rc stability criteria extrema for the current
        ! time-step and located on both the local (loc) and the global (glb)
        ! computational domains

        real(kind(0d0)) :: blkmod1, blkmod2 !<
            !! Fluid bulk modulus for Woods mixture sound speed

        integer :: i, j, k !< Generic loop iterators

        integer :: Nfq

        ! Computing Stability Criteria at Current Time-step ================
!$acc parallel loop collapse(3) gang vector default(present) private(alpha_rho, vel, alpha, Re)
            do k = 0, n
                do j = 0, m

                    do i = 1, num_fluids
                        alpha_rho(i) = q_prim_vf(i)%sf(j, k)
                        alpha(i) = q_prim_vf(E_idx + i)%sf(j, k)
                    end do

                    call s_convert_species_to_mixture_variables_acc(rho, gamma, pi_inf, alpha, alpha_rho, Re, j, k)

                    do i = 1, num_dims
                        vel(i) = q_prim_vf(contxe + i)%sf(j, k)
                    end do

                    pres = q_prim_vf(E_idx)%sf(j, k)

                    ! Compute mixture sound speed
                    c = (((gamma + 1d0)*pres + pi_inf)/(gamma*rho))

                    c = sqrt(c)

                    if (n > 0) then
                        !2D
                        icfl_sf(j, k) = dt/min(dx(j)/(abs(vel(1)) + c), &
                                                  dy(k)/(abs(vel(2)) + c))

                        if (any(Re_size > 0)) then

                            vcfl_sf(j, k) = maxval(dt/Re)/min(dx(j), dy(k))**2d0

                            Rc_sf(j, k) = min(dx(j)*(abs(vel(1)) + c), &
                                                 dy(k)*(abs(vel(2)) + c)) &
                                             /maxval(1d0/Re)

                        end if

                    else
                        !1D
                        icfl_sf(j, k) = (dt/dx(j))*(abs(vel(1)) + c)

                        if (any(Re_size > 0)) then

                            vcfl_sf(j, k) = maxval(dt/Re)/dx(j)**2d0

                            Rc_sf(j, k) = dx(j)*(abs(vel(1)) + c)/maxval(1d0/Re)

                        end if

                    end if

                end do
            end do
        ! END: Computing Stability Criteria at Current Time-step ===========

        ! Determining local stability criteria extrema at current time-step

        !$acc kernels
        icfl_max_loc = maxval(icfl_sf)
        !$acc end kernels

        if (any(Re_size > 0)) then
            !$acc kernels
            vcfl_max_loc = maxval(vcfl_sf)
            Rc_min_loc = minval(Rc_sf)
            !$acc end kernels
        end if

        !$acc update host(icfl_max_loc, vcfl_max_loc, Rc_min_loc)

        ! Determining global stability criteria extrema at current time-step
        if (num_procs > 1) then
            call s_mpi_reduce_stability_criteria_extrema(icfl_max_loc, &
                                                         vcfl_max_loc, &
                                                         ccfl_max_loc, &
                                                         Rc_min_loc, &
                                                         icfl_max_glb, &
                                                         vcfl_max_glb, &
                                                         ccfl_max_glb, &
                                                         Rc_min_glb)
        else
            icfl_max_glb = icfl_max_loc
            if (any(Re_size > 0)) vcfl_max_glb = vcfl_max_loc
            if (any(Re_size > 0)) Rc_min_glb = Rc_min_loc
        end if

        ! Determining the stability criteria extrema over all the time-steps
        if (icfl_max_glb > icfl_max) icfl_max = icfl_max_glb

        if (any(Re_size > 0)) then
            if (vcfl_max_glb > vcfl_max) vcfl_max = vcfl_max_glb
            if (Rc_min_glb < Rc_min) Rc_min = Rc_min_glb
        end if

        ! Outputting global stability criteria extrema at current time-step
        if (proc_rank == 0) then
            if (any(Re_size > 0)) then
                write (1, '(6X,I8,6X,F10.6,6X,F9.6,6X,F9.6,6X,F10.6)') &
                    t_step, t_step*dt, icfl_max_glb, &
                    vcfl_max_glb, &
                    Rc_min_glb
            else
                write (1, '(13X,I8,14X,F10.6,13X,F9.6)') &
                    t_step, t_step*dt, icfl_max_glb
            end if

            if (icfl_max_glb /= icfl_max_glb) then
                print '(A)', 'ICFL is NaN. Exiting ...'
                ! print*, (dt/dx(:)),ABS(vel(1)),c

                call s_mpi_abort()
            elseif (icfl_max_glb > 1d0) then
                print '(A)', 'ICFL is greater than 1.0. Exiting ...'
                print *, 'icfl', icfl_max_glb
                call s_mpi_abort()
            end if
        end if

        call s_mpi_barrier()

    end subroutine s_write_run_time_information ! --------------------------

    !>  The goal of this subroutine is to output the grid and
        !!      conservative variables data files for given time-step.
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param t_step Current time-step
    subroutine s_write_serial_data_files(q_cons_vf, t_step) ! ---------------------

        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        integer, intent(IN) :: t_step

        character(LEN=path_len + 2*name_len) :: t_step_dir !<
            !! Relative path to the current time-step directory

        character(LEN=path_len + 3*name_len) :: file_path !<
            !! Relative path to the grid and conservative variables data files

        logical :: file_exist !<
            !! Logical used to check existence of current time-step directory

        character(LEN=15) :: FMT

        integer :: i, j, k, ii !< Generic loop iterators

        real(kind(0d0)) :: gamma, lit_gamma, pi_inf     !< Temporary EOS params
        real(kind(0d0)) :: rho                          !< Temporary density
        real(kind(0d0)), dimension(2) :: Re !< Temporary Reynolds number

        ! Creating or overwriting the time-step root directory
        write (t_step_dir, '(A,I0,A,I0)') trim(case_dir)//'/p_all'

        ! Creating or overwriting the current time-step directory
        write (t_step_dir, '(A,I0,A,I0)') trim(case_dir)//'/p_all/p', &
            proc_rank, '/', t_step

        file_path = trim(t_step_dir)//'/.'
        call my_inquire(file_path, file_exist)
        if (file_exist) call s_delete_directory(trim(t_step_dir))
        call s_create_directory(trim(t_step_dir))

        ! Writing the grid data file in the x-direction
        file_path = trim(t_step_dir)//'/x_cb.dat'

        open (2, FILE=trim(file_path), &
              FORM='unformatted', &
              STATUS='new')
        write (2) x_cb(-1:m); close (2)

        ! Writing the grid data files in the y- and z-directions
        if (n > 0) then

            file_path = trim(t_step_dir)//'/y_cb.dat'

            open (2, FILE=trim(file_path), &
                  FORM='unformatted', &
                  STATUS='new')
            write (2) y_cb(-1:n); close (2)


        end if

        ! Writing the conservative variables data files
        do i = 1, sys_size
            write (file_path, '(A,I0,A)') trim(t_step_dir)//'/q_cons_vf', &
                i, '.dat'

            open (2, FILE=trim(file_path), &
                  FORM='unformatted', &
                  STATUS='new')

            write (2) q_cons_vf(i)%sf(0:m, 0:n); close (2)
        end do

        gamma = fluid_pp(1)%gamma
        lit_gamma = 1d0/fluid_pp(1)%gamma + 1d0
        pi_inf = fluid_pp(1)%pi_inf

        if (precision == 1) then
            FMT = "(2F30.3)"
        else
            FMT = "(2F40.14)"
        end if

        ! writting an output directory
        write (t_step_dir, '(A,I0,A,I0)') trim(case_dir)//'/D'
        file_path = trim(t_step_dir)//'/.'

        inquire (FILE=trim(file_path), EXIST=file_exist)

        if (.not. file_exist) call s_create_directory(trim(t_step_dir))

        !1D
        if (n == 0) then

            do i = 1, sys_size
                write (file_path, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/prim.', i, '.', proc_rank, '.', t_step, '.dat'

                open (2, FILE=trim(file_path))
                do j = 0, m
                    call s_convert_to_mixture_variables(q_cons_vf, j, 0, rho, gamma, pi_inf, Re )
                    lit_gamma = 1d0/gamma + 1d0

                    if (((i >= cont_idx%beg) .and. (i <= cont_idx%end)) &
                        .or. &
                        ((i >= adv_idx%beg) .and. (i <= adv_idx%end)) &
                        ) then
                        write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0)
                    else if (i == mom_idx%beg) then !u
                        write (2, FMT) x_cb(j), q_cons_vf(mom_idx%beg)%sf(j, 0)/rho
                    else if (i == E_idx) then !p
                        !Stiffened gas pressure from energy
                        write (2, FMT) x_cb(j), &
                            ( &
                            q_cons_vf(E_idx)%sf(j, 0) - &
                            0.5d0*(q_cons_vf(mom_idx%beg)%sf(j, 0)**2.d0)/rho - &
                            pi_inf &
                            )/gamma
                    end if
                end do
                close (2)
            end do

            do i = 1, sys_size
                write (file_path, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/cons.', i, '.', proc_rank, '.', t_step, '.dat'

                open (2, FILE=trim(file_path))
                do j = 0, m
                    write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0)
                end do
                close (2)
            end do
        end if

        if (precision == 1) then
            FMT = "(3F30.7)"
        else
            FMT = "(3F40.14)"
        end if

        ! 2D
        if (n > 0) then
            do i = 1, sys_size
                write (file_path, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/cons.', i, '.', proc_rank, '.', t_step, '.dat'
                open (2, FILE=trim(file_path))
                do j = 0, m
                    do k = 0, n
                        write (2, FMT) x_cb(j), y_cb(k), q_cons_vf(i)%sf(j, k)
                    end do
                    write (2, *)
                end do
                close (2)
            end do
        end if

        if (precision == 1) then
            FMT = "(4F30.7)"
        else
            FMT = "(4F40.14)"
        end if


    end subroutine s_write_serial_data_files ! ------------------------------------

    !>  The goal of this subroutine is to output the grid and
        !!      conservative variables data files for given time-step.
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param t_step Current time-step
    subroutine s_write_parallel_data_files(q_cons_vf, t_step) ! --

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_cons_vf

        integer, intent(IN) :: t_step

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

        integer :: i !< Generic loop iterator

        ! Initialize MPI data I/O
        call s_initialize_mpi_data(q_cons_vf)

        ! Open the file to write all flow variables
        write (file_loc, '(I0,A)') t_step, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
        inquire (FILE=trim(file_loc), EXIST=file_exist)
        if (file_exist .and. proc_rank == 0) then
            call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
        end if
        call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                           mpi_info_int, ifile, ierr)

        ! Size of local arrays
        data_size = (m + 1)*(n + 1)

        ! Resize some integers so MPI can write even the biggest files
        m_MOK = int(m_glb + 1, MPI_OFFSET_KIND)
        n_MOK = int(n_glb + 1, MPI_OFFSET_KIND)
        WP_MOK = int(8d0, MPI_OFFSET_KIND)
        MOK = int(1d0, MPI_OFFSET_KIND)
        str_MOK = int(name_len, MPI_OFFSET_KIND)
        NVARS_MOK = int(sys_size, MPI_OFFSET_KIND)

        do i = 1, sys_size !TODO: check if correct (sys_size
            var_MOK = int(i, MPI_OFFSET_KIND)

            ! Initial displacement to skip at beginning of file
            disp = m_MOK*max(MOK, n_MOK)*WP_MOK*(var_MOK - 1)

            call MPI_FILE_SET_VIEW(ifile, disp, MPI_DOUBLE_PRECISION, MPI_IO_DATA%view(i), &
                                   'native', mpi_info_int, ierr)
            call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                    MPI_DOUBLE_PRECISION, status, ierr)
        end do

        call MPI_FILE_CLOSE(ifile, ierr)

#endif

    end subroutine s_write_parallel_data_files ! ---------------------------


    !>  This writes a formatted data file for the flow probe information
        !!  @param t_step Current time-step
        !!  @param q_cons_vf Conservative variables
    subroutine s_write_probe_files(t_step, q_cons_vf) ! -----------

        integer, intent(IN) :: t_step
        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf

        real(kind(0d0)), dimension(-1:m) :: distx
        real(kind(0d0)), dimension(-1:n) :: disty

        ! The cell-averaged partial densities, density, velocity, pressure,
        ! volume fractions, specific heat ratio function, liquid stiffness
        ! function, and sound speed.
        real(kind(0d0)) :: rho
        real(kind(0d0)), dimension(num_dims) :: vel
        real(kind(0d0)) :: pres
        real(kind(0d0)), dimension(num_fluids) :: alpha
        real(kind(0d0)) :: gamma
        real(kind(0d0)) :: pi_inf
        real(kind(0d0)) :: c
        real(kind(0d0)), dimension(2) :: Re

        integer :: i, j, k, s !< Generic loop iterator

        real(kind(0d0)) :: nondim_time !< Non-dimensional time

        real(kind(0d0)) :: tmp !<
            !! Temporary variable to store quantity for mpi_allreduce

        real(kind(0d0)) :: blkmod1, blkmod2 !<
            !! Fluid bulk modulus for Woods mixture sound speed

        ! Non-dimensional time calculation
        if (time_stepper == 23) then
            nondim_time = mytime
        else
            if (t_step_old /= dflt_int) then
                nondim_time = real(t_step + t_step_old, kind(0d0))*dt
            else
                nondim_time = real(t_step, kind(0d0))*dt !*1.d-5/10.0761131451d0
            end if
        end if

        do i = 1, num_probes
            ! Zeroing out flow variables for all processors
            rho = 0d0
            do s = 1, num_dims
                vel(s) = 0d0
            end do
            pres = 0d0
            gamma = 0d0
            pi_inf = 0d0
            c = 0d0

            ! Find probe location in terms of indices on a
            ! specific processor
            if (n == 0) then ! 1D simulation
                if ((probe(i)%x >= x_cb(-1)) .and. (probe(i)%x <= x_cb(m))) then
                    do s = -1, m
                        distx(s) = x_cb(s) - probe(i)%x
                        if (distx(s) < 0d0) distx(s) = 1000d0
                    end do
                    j = minloc(distx, 1)
                    if (j == 1) j = 2 ! Pick first point if probe is at edge
                    k = 0

                    ! Computing/Sharing necessary state variables
                    call s_convert_to_mixture_variables(q_cons_vf, j - 2, k, &
                                                        rho, gamma, pi_inf, &
                                                        Re)
                    do s = 1, num_dims
                        vel(s) = q_cons_vf(cont_idx%end + s)%sf(j - 2, k)/rho
                    end do

                    !Stiffened gas pressure from energy
                    pres = ( &
                           q_cons_vf(E_idx)%sf(j - 2, k) - &
                           0.5d0*(q_cons_vf(2)%sf(j - 2, k)**2.d0)/q_cons_vf(1)%sf(j - 2, k) - &
                           pi_inf &
                           )/gamma

                    ! Compute mixture sound speed
                    c = (((gamma + 1d0)*pres + pi_inf)/(gamma*rho))

                    c = sqrt(c)

                end if
            else ! 2D simulation
                if ((probe(i)%x >= x_cb(-1)) .and. (probe(i)%x <= x_cb(m))) then
                    if ((probe(i)%y >= y_cb(-1)) .and. (probe(i)%y <= y_cb(n))) then
                        do s = -1, m
                            distx(s) = x_cb(s) - probe(i)%x
                            if (distx(s) < 0d0) distx(s) = 1000d0
                        end do
                        do s = -1, n
                            disty(s) = y_cb(s) - probe(i)%y
                            if (disty(s) < 0d0) disty(s) = 1000d0
                        end do
                        j = minloc(distx, 1)
                        k = minloc(disty, 1)
                        if (j == 1) j = 2 ! Pick first point if probe is at edge
                        if (k == 1) k = 2 ! Pick first point if probe is at edge

                        ! Computing/Sharing necessary state variables
                        call s_convert_to_mixture_variables(q_cons_vf, j - 2, k - 2, & 
                                                            rho, gamma, pi_inf, &
                                                            Re)
                        do s = 1, num_dims
                            vel(s) = q_cons_vf(cont_idx%end + s)%sf(j - 2, k - 2)/rho
                        end do


                        !Stiffened gas pressure from energy
                        pres = ( &
                               q_cons_vf(E_idx)%sf(j - 2, k - 2) - &
                               0.5d0*((q_cons_vf(2)%sf(j - 2, k - 2)**2.d0 + &
                                       q_cons_vf(3)%sf(j - 2, k - 2)**2.d0)/q_cons_vf(1)%sf(j - 2, k - 2)) - &
                               pi_inf &
                               )/gamma


                        ! Compute mixture sound speed
                        c = (((gamma + 1d0)*pres + pi_inf)/(gamma*rho))
                        c = sqrt(c)

                    end if
                end if
            end if

            if (num_procs > 1) then
                #:for VAR in ['rho','pres','gamma','pi_inf','c']
                    tmp = ${VAR}$
                    call s_mpi_allreduce_sum(tmp, ${VAR}$)
                #:endfor

                do s = 1, num_dims
                    tmp = vel(s)
                    call s_mpi_allreduce_sum(tmp, vel(s))
                end do
            end if

       end do

    end subroutine s_write_probe_files ! -----------------------------------

    !>  The goal of this subroutine is to write to the run-time
        !!      information file basic footer information applicable to
        !!      the current computation and to close the file when done.
        !!      The footer contains the stability criteria extrema over
        !!      all of the time-steps and the simulation run-time.
    subroutine s_close_run_time_information_file() ! -----------------------

        real(kind(0d0)) :: run_time !< Run-time of the simulation

        ! Writing the footer of and closing the run-time information file
        write (1, '(A)') '----------------------------------------'// &
            '----------------------------------------'
        write (1, '(A)') ''

        write (1, '(A,F9.6)') 'ICFL Max: ', icfl_max
        if (any(Re_size > 0)) write (1, '(A,F9.6)') 'VCFL Max: ', vcfl_max
        if (any(Re_size > 0)) write (1, '(A,F10.6)') 'Rc Min: ', Rc_min

        call cpu_time(run_time)

        write (1, '(A)') ''
        write (1, '(A,I0,A)') 'Run-time: ', int(anint(run_time)), 's'
        write (1, '(A)') '========================================'// &
            '========================================'
        close (1)

    end subroutine s_close_run_time_information_file ! ---------------------

    !> Closes probe files
    subroutine s_close_probe_files() ! -------------------------------------

        integer :: i !< Generic loop iterator

        do i = 1, num_probes
            close (i + 30)
        end do

    end subroutine s_close_probe_files ! -----------------------------------

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_data_output_module() ! -------------------------

        type(int_bounds_info) :: ix, iy, iz

        integer :: i !< Generic loop iterator


        ! Allocating/initializing ICFL, VCFL, CCFL and Rc stability criteria
        @:ALLOCATE(icfl_sf(0:m, 0:n))
        icfl_max = 0d0
        
        if (any(Re_size > 0)) then
            @:ALLOCATE(vcfl_sf(0:m, 0:n))
            @:ALLOCATE(Rc_sf  (0:m, 0:n))
            
            vcfl_max = 0d0
            Rc_min   = 1d3
        end if

        ! Associating the procedural pointer to the appropriate subroutine
        ! that will be utilized in the conversion to the mixture variables

        s_convert_to_mixture_variables => &
            s_convert_species_to_mixture_variables


        if (parallel_io .neqv. .true.) then
            s_write_data_files => s_write_serial_data_files
        else
            s_write_data_files => s_write_parallel_data_files
        end if

    end subroutine s_initialize_data_output_module ! -----------------------

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_data_output_module() ! ---------------------------

        integer :: i !< Generic loop iterator

        ! Deallocating the ICFL, VCFL, CCFL, and Rc stability criteria
        @:DEALLOCATE(icfl_sf)
        if (any(Re_size > 0)) then
            @:DEALLOCATE(vcfl_sf, Rc_sf)
        end if


        ! Disassociating the pointer to the procedure that was utilized to
        ! to convert mixture or species variables to the mixture variables
        s_convert_to_mixture_variables => null()
        s_write_data_files => null()

    end subroutine s_finalize_data_output_module ! -------------------------

end module m_data_output
