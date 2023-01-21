!>
!! @file m_start_up.f90
!! @brief  Contains module m_start_up

!> @brief This module contains the subroutines that read in and check the
!!              consistency of the user provided inputs.
module m_start_up

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_compile_specific
    ! ==========================================================================

    implicit none

contains

    !>  Reads the configuration file post_process.inp, in order
        !!      to populate parameters in module m_global_parameters.f90
        !!      with the user provided inputs
    subroutine s_read_input_file() ! ---------------------------------------

        character(LEN=name_len) :: file_loc !<
            !! Generic string used to store the address of a particular file

        logical :: file_check !<
            !! Generic logical used for the purpose of asserting whether a file
            !! is or is not present in the designated location

        ! Namelist for all of the parameters to be inputed by the user
        namelist /user_inputs/ case_dir, m, n, t_step_start, &
            t_step_stop, t_step_save, &
            num_fluids, weno_order, bc_x, &
            bc_y, fluid_pp, format, precision, &
            cons_vars_wrt, prim_vars_wrt, &
            omega_wrt, schlieren_wrt, schlieren_alpha, fd_order, &
            parallel_io

        ! Inquiring the status of the post_process.inp file
        file_loc = 'post_process.inp'
        inquire (FILE=trim(file_loc), EXIST=file_check)

        ! Checking whether the input file is there. If it is, the input file
        ! is read. If not, the program is terminated.
        if (file_check) then
            open (1, FILE=trim(file_loc), FORM='formatted', &
                  STATUS='old', ACTION='read')
            read (1, NML=user_inputs)
            close (1)
            ! Store m,n,p into global m,n,p
            m_glb = m
            n_glb = n
        else
            print '(A)', 'File post_process.inp is missing. Exiting ...'
            call s_mpi_abort()
        end if

    end subroutine s_read_input_file ! -------------------------------------

    !>  Checking that the user inputs make sense, i.e. that the
        !!      individual choices are compatible with the code's options
        !!      and that the combination of these choices results into a
        !!      valid configuration for the post-process
    subroutine s_check_input_file() ! --------------------------------------

        character(LEN=len_trim(case_dir)) :: file_loc !<
            !! Generic string used to store the address of a particular file

        logical :: dir_check !<
            !! Logical variable used to test the existence of folders

        integer :: i  !< Generic loop iterator


        ! Checking the existence of the case folder
        case_dir = adjustl(case_dir)

        file_loc = trim(case_dir)//'/.'

        call my_inquire(file_loc, dir_check)

        ! Constraint on the location of the case directory
        if (dir_check .neqv. .true.) then
            print '(A)', 'Unsupported choice for the value of '// &
                'case_dir. Exiting ...'
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

            ! Constraints on the time-stepping parameters
        elseif (t_step_start < 0) then
            print '(A)', 'Unsupported choice for the value of '// &
                't_step_start. Exiting ...'
            call s_mpi_abort()
        elseif (t_step_stop < t_step_start) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for t_step_start and t_step_stop. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (t_step_save > t_step_stop - t_step_start) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for t_step_start, t_step_stop and '// &
                't_step_save. Exiting ...'
            call s_mpi_abort()
        elseif (num_fluids /= dflt_int &
                .and. &
                (num_fluids < 1 .or. num_fluids > num_fluids)) then
            print '(A)', 'Unsupported value of num_fluids. Exiting ...'
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
                'values for num_procs, m, n and '// &
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

        ! Constraints on the stiffened equation of state fluids parameters
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

        end do

        ! Constraints on the format of the formatted database file(s)
        if (format /= 1 .and. format /= 2) then
            print '(A)', 'Unsupported choice for the value of format. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif ((precision /= 2) .and. (parallel_io .neqv. .false.)) then
            print '(A)', 'Unsupported combination of precision and parallel IO. '// &
                'Please use precision == 2 when enabling parallel_io.  Exiting ...'
            call s_mpi_abort()

            ! Constraints on the precision of the formatted database file(s)
        elseif (precision /= 1 .and. precision /= 2) then
            print '(A)', 'Unsupported choice for the value of '// &
                'precision. Exiting ...'
            call s_mpi_abort()
        end if

        ! Constraints on the post-processing of the vorticity
        if (n == 0 .and. omega_wrt(1)) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for n and omega_wrt(1). Exiting ...'
            call s_mpi_abort()
        elseif (n == 0 .and. omega_wrt(2)) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for n and omega_wrt(2). Exiting ...'
            call s_mpi_abort()
        elseif (n == 0 .and. schlieren_wrt) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for n and schlieren_wrt. Exiting ...'
            call s_mpi_abort()

            ! Constraints on post-processing combination of flow variables
        elseif ((any((/ cons_vars_wrt, &
                       prim_vars_wrt, &
                       schlieren_wrt/)) .neqv. .true.) &
                .and. &
                (any(omega_wrt) .neqv. .true.)) then
            print '(A)', 'None of the flow variables have been '// &
                'selected for post-process. Exiting ...'
            call s_mpi_abort()
        end if

        ! Constraints on the coefficients of numerical Schlieren function
        do i = 1, num_fluids
            if (schlieren_alpha(i) /= dflt_real &
                .and. &
                schlieren_alpha(i) <= 0d0) then
                print '(A,I0,A)', 'Unsupported choice for the value of '// &
                    'schlieren_alpha(', i, '). Exiting ...'
                call s_mpi_abort()
            elseif (((i > num_fluids .or. (schlieren_wrt .neqv. .true.)) &
                     .and. &
                     schlieren_alpha(i) /= dflt_real) &
                    .or. &
                    ((i <= num_fluids .and. schlieren_wrt) &
                     .and. &
                     schlieren_alpha(i) <= 0d0)) then
                print '(A,I0,A)', 'Unsupported choice of the '// &
                    'combination of values for '// &
                    'num_fluids, schlieren_wrt and '// &
                    'schlieren_alpha(', i, '). Exiting ...'
                call s_mpi_abort()
            end if
        end do

        ! Constraints on the order of the finite difference scheme
        if (fd_order /= dflt_int &
            .and. &
            fd_order /= 1 .and. fd_order /= 2 .and. fd_order /= 4) then
            print '(A)', 'Unsupported choice for the value of '// &
                'fd_order. Exiting ...'
            call s_mpi_abort()
        elseif ((any(omega_wrt) .or. schlieren_wrt) &
                .and. &
                fd_order == dflt_int) then
            print '(A)', 'BB Unsupported choice of the combination of '// &
                'values for omega_wrt, schlieren_wrt and '// &
                'fd_order. Exiting ...'
            call s_mpi_abort()
        end if

    end subroutine s_check_input_file ! ------------------------------------

end module m_start_up
