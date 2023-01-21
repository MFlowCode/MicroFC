!>
!! @file p_main.f90
!! @brief Contains program p_main

!> @brief The post-process restructures raw unformatted data, outputted by
!!              the simulation, into a formatted database, Silo-HDF5 or Binary,
!!              chosen by the user. The user may also specify which variables to
!!              include in the database. The choices range from any one of the
!!              primitive and conservative variables, as well as quantities that
!!              can be derived from those such as the unadvected volume fraction,
!!              specific heat ratio, liquid stiffness, speed of sound, vorticity
!!              and the numerical Schlieren function.
program p_main

    ! Dependencies =============================================================
    use m_derived_types         !< Definitions of the derived types

    use m_global_parameters     !< Global parameters for the code

    use m_mpi_proxy             !< Message passing interface (MPI) module proxy

    use m_variables_conversion  !< Subroutines to change the state variables from
                                !! one form to another

    use m_start_up              !< Subroutines that read in and check consistency
                                !! of the user provided inputs and data

    use m_data_input            !< Procedures reading raw simulation data to fill
                                !! the conservative, primitive and grid variables

    use m_data_output           !< Procedures that write the grid and chosen flow
                                !! variable(s) to the formatted database file(s)

    use m_derived_variables     !< Procedures used to compute quantites derived
                                !! from the conservative and primitive variables
    ! ==========================================================================

    implicit none

    integer :: t_step !< Iterator for the main time-stepping loop

    character(LEN=name_len) :: varname !<
    !! Generic storage for the name(s) of the flow variable(s) that will be added
    !! to the formatted database file(s)

    !> @name Generic loop iterator
    !> @{
    integer :: i, j, k
    !> @}

    real(kind(0d0)) :: total_volume !<
    !! Variable for the total volume of the second volume fraction
    !! to later on track the evolution of the radius of a bubble over time

    ! Initialization of the MPI environment
    call s_mpi_initialize()

    ! Processor with rank 0 assigns default user input values prior to reading
    ! those in from the input file. Next, the user inputs are read in and their
    ! consistency is checked. The detection of any inconsistencies automatically
    ! leads to the termination of the post-process.
    if (proc_rank == 0) then
        call s_assign_default_values_to_user_inputs()
        call s_read_input_file()
        call s_check_input_file()
        
        print '(" Post-processing a "I0"x"I0"x"I0" case on "I0" rank(s)")', m, n, num_procs
    end if

    ! Broadcasting the user inputs to all of the processors and performing the
    ! parallel computational domain decomposition. Neither procedure has to be
    ! carried out if the simulation is in fact not truly executed in parallel.
    call s_mpi_bcast_user_inputs()
    call s_initialize_parallel_io()
    call s_mpi_decompose_computational_domain()

    ! Computation of parameters, allocation procedures, and/or any other tasks
    ! needed to properly setup the modules
    call s_initialize_global_parameters_module()
    if (num_procs > 1) call s_initialize_mpi_proxy_module()
    call s_initialize_variables_conversion_module()
    call s_initialize_data_input_module()
    call s_initialize_derived_variables_module()
    call s_initialize_data_output_module()

    ! Associate pointers for serial or paralle I/O
    if (parallel_io .neqv. .true.) then
        s_read_data_files => s_read_serial_data_files
    else
        s_read_data_files => s_read_parallel_data_files
    end if

    ! Setting the time-step iterator to the first time step to be post-processed
    t_step = t_step_start

    ! Time-Marching Loop =======================================================
    do
        if (proc_rank == 0) then
            print '(" ["I3"%]  Save "I8" of "I0" @ t_step = "I0"")',            &
                  int(100*(real(t_step + 1)/(t_step_stop - t_step_start + 1))), &
                  (t_step      - t_step_start)/t_step_save + 1,                 &
                  (t_step_stop - t_step_start)/t_step_save + 1,                 &
                  t_step
        end if

        ! Populating the grid and conservative variables
        call s_read_data_files(t_step)
        ! Populating the buffer regions of the grid variables
        if (buff_size > 0) then
            call s_populate_grid_variables_buffer_regions()
        end if

        ! Populating the buffer regions of the conservative variables
        if (buff_size > 0) then
            call s_populate_conservative_variables_buffer_regions()
        end if

        ! Converting the conservative variables to the primitive ones
        call s_convert_conservative_to_primitive_variables(q_cons_vf, q_prim_vf)
        ! Opening a new formatted database file
        call s_open_formatted_database_file(t_step)

        ! Adding the grid to the formatted database file
        call s_write_grid_to_formatted_database_file(t_step)

        ! Computing centered finite-difference coefficients in x-direction
        if (omega_wrt(2) .or. schlieren_wrt) then
            call s_compute_finite_difference_coefficients(m, offset_x, x_cc, &
                                                          fd_coeff_x)
        end if

        ! Computing centered finite-difference coefficients in y-direction
        if (omega_wrt(1) .or. (n > 0 .and. schlieren_wrt)) then
            call s_compute_finite_difference_coefficients(n, offset_y, y_cc, &
                                                          fd_coeff_y)
        end if

        ! Adding the partial densities to the formatted database file ----------
        do i = 1, num_fluids
            if (cons_vars_wrt .or. prim_vars_wrt) then

                q_sf = q_cons_vf(i)%sf(-offset_x%beg:m + offset_x%end, &
                                       -offset_y%beg:n + offset_y%end )

                write (varname, '(A,I0)') 'alpha_rho', i
                call s_write_variable_to_formatted_database_file(varname, t_step)

                varname(:) = ' '
            end if
        end do
        ! ----------------------------------------------------------------------

        ! Adding the density to the formatted database file --------------------
        if (cons_vars_wrt .or. prim_vars_wrt) then

            q_sf = rho_sf(-offset_x%beg:m + offset_x%end, &
                          -offset_y%beg:n + offset_y%end)

            write (varname, '(A)') 'rho'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '

        end if
        ! ----------------------------------------------------------------------

        ! Adding the momentum to the formatted database file -------------------
        do i = 1, E_idx - mom_idx%beg
            if (cons_vars_wrt) then

                q_sf = q_cons_vf(i + cont_idx%end)%sf( &
                       -offset_x%beg:m + offset_x%end, &
                       -offset_y%beg:n + offset_y%end )

                write (varname, '(A,I0)') 'mom', i
                call s_write_variable_to_formatted_database_file(varname, t_step)

                varname(:) = ' '

            end if
        end do
        ! ----------------------------------------------------------------------

        ! Adding the velocity to the formatted database file -------------------
        do i = 1, E_idx - mom_idx%beg
            if (prim_vars_wrt) then

                q_sf = q_prim_vf(i + cont_idx%end)%sf( &
                       -offset_x%beg:m + offset_x%end, &
                       -offset_y%beg:n + offset_y%end )

                write (varname, '(A,I0)') 'vel', i
                call s_write_variable_to_formatted_database_file(varname, t_step)

                varname(:) = ' '

            end if
        end do
        ! ----------------------------------------------------------------------

        ! Adding the energy to the formatted database file ---------------------
        if (cons_vars_wrt) then

            q_sf = q_cons_vf(E_idx)%sf(-offset_x%beg:m + offset_x%end, &
                                       -offset_y%beg:n + offset_y%end )

            write (varname, '(A)') 'E'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '

        end if
        ! ----------------------------------------------------------------------

        ! Adding the pressure to the formatted database file -------------------
        if (prim_vars_wrt) then

            q_sf = q_prim_vf(E_idx)%sf(-offset_x%beg:m + offset_x%end, &
                                       -offset_y%beg:n + offset_y%end )

            write (varname, '(A)') 'pres'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '

        end if
        ! ----------------------------------------------------------------------

        ! Adding the volume fraction(s) to the formatted database file ---------
        do i = 1, num_fluids 
            if (cons_vars_wrt .or. prim_vars_wrt) then

                q_sf = q_cons_vf(i + E_idx)%sf( &
                       -offset_x%beg:m + offset_x%end, &
                       -offset_y%beg:n + offset_y%end )

                write (varname, '(A,I0)') 'alpha', i
                call s_write_variable_to_formatted_database_file(varname, t_step)

                varname(:) = ' '

            end if
        end do


        ! ----------------------------------------------------------------------

        ! Adding the vorticity to the formatted database file ------------------
        if (n > 0) then
            do i = 1, E_idx - cont_idx%end
                if (omega_wrt(i)) then

                    call s_derive_vorticity_component(i, q_prim_vf, q_sf)

                    write (varname, '(A,I0)') 'omega', i
                    call s_write_variable_to_formatted_database_file(varname, t_step)

                    varname(:) = ' '
                end if
            end do
        end if
        ! ----------------------------------------------------------------------

        ! Adding numerical Schlieren function to formatted database file -------
        if (schlieren_wrt) then

            call s_derive_numerical_schlieren_function(q_cons_vf, rho_sf, q_sf)

            write (varname, '(A)') 'schlieren'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '

        end if
        ! ----------------------------------------------------------------------

        ! Closing the formatted database file
        call s_close_formatted_database_file()

        ! Modifies the time-step iterator so that it may reach the final time-
        ! step to be post-processed, in the case that this one is not originally
        ! attainable through constant incrementation from the first time-step.
        ! This modification is performed upon reaching the final time-step. In
        ! case that it is not needed, the post-processor is done and may exit.
        if ((t_step_stop - t_step) < t_step_save .and. t_step_stop /= t_step) then
            t_step = t_step_stop - t_step_save
        elseif (t_step == t_step_stop) then
            exit
        end if

        ! Incrementing time-step iterator to next time-step to be post-processed
        t_step = t_step + t_step_save

    end do
    ! END: Time-Marching Loop ==================================================

    close (11)

    ! Disassociate pointers for serial and parallel I/O
    s_read_data_files => null()

    ! Deallocation procedures for the modules
    call s_finalize_data_output_module()
    call s_finalize_derived_variables_module()
    call s_finalize_data_input_module()
    call s_finalize_variables_conversion_module()
    if (num_procs > 1) call s_finalize_mpi_proxy_module()
    call s_finalize_global_parameters_module()

    ! Finalizing the MPI environment
    call s_mpi_finalize()

end program p_main
