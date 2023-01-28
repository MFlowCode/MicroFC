!>
!! @file m_mpi_proxy.f90
!! @brief Contains module m_mpi_proxy

!> @brief This module serves as a proxy to the parameters and subroutines
!!              available in the MPI implementation's MPI module. Specifically,
!!              the role of the proxy is to harness basic MPI commands into more
!!              complex procedures as to achieve the required pre-processing
!!              communication goals.
module m_mpi_proxy

    ! Dependencies =============================================================
#ifdef MFC_MPI
    use mpi                     !< Message passing interface (MPI) module
#endif

    use m_derived_types         !< Definitions of the derived types

    use m_global_parameters     !< Global parameters for the code

    use m_mpi_common
    ! ==========================================================================

    implicit none

    integer, private :: ierr !<
    !! Generic flags used to identify and report MPI errors

contains


    !> Since only processor with rank 0 is in charge of reading
        !!       and checking the consistency of the user provided inputs,
         !!       these are not available to the remaining processors. This
        !!       subroutine is then in charge of broadcasting the required
        !!       information.
    subroutine s_mpi_bcast_user_inputs() ! ---------------------------------

#ifdef MFC_MPI

        ! Generic loop iterator
        integer :: i

        ! Logistics
        call MPI_BCAST(case_dir, len(case_dir), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

        #:for VAR in ['t_step_old', 'm', 'n', 'm_glb', 'n_glb',  &
            & 'loops_x', 'loops_y',  'num_fluids',     &
            & 'weno_order', 'precision', & 
            & 'num_patches' ]
            call MPI_BCAST(${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in [ 'old_grid','old_ic','stretch_x','stretch_y', 'parallel_io'  ]
            call MPI_BCAST(${VAR}$, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        #:endfor


        #:for VAR in [ 'x_domain%beg', 'x_domain%end', 'y_domain%beg',         &
            & 'y_domain%end', 'a_x', 'a_y', &
            & 'x_a', 'x_b', 'y_a', 'y_b', 'bc_x%beg',     &
            & 'bc_x%end', 'bc_y%beg', 'bc_y%end' ]
            call MPI_BCAST(${VAR}$, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        do i = 1, num_patches_max
            call MPI_BCAST(patch_icpp(i)%geometry, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_icpp(i)%smooth_patch_id, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

            call MPI_BCAST(patch_icpp(i)%smoothen, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_icpp(i)%alter_patch(0), num_patches_max, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

            #:for VAR in [ 'x_centroid', 'y_centroid', &
                & 'length_x', 'length_y',  'radius', 'epsilon',     &
                & 'beta', 'smooth_coeff', 'rho',        &
                & 'pres', 'gamma', 'pi_inf',  ]
                call MPI_BCAST(patch_icpp(i)%${VAR}$, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            call MPI_BCAST(patch_icpp(i)%normal(1), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_icpp(i)%radii(1), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_icpp(i)%vel(1), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_icpp(i)%alpha_rho(1), num_fluids_max, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_icpp(i)%alpha(1), num_fluids_max - 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        end do

        ! Fluids physical parameters
        do i = 1, num_fluids_max
            #:for VAR in [ 'gamma','pi_inf' ]
                call MPI_BCAST(fluid_pp(i)%${VAR}$, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do
#endif

    end subroutine s_mpi_bcast_user_inputs ! -------------------------------



    !> Description: This subroutine takes care of efficiently distributing
        !!              the computational domain among the available processors
        !!             as well as recomputing some of the global parameters so
        !!              that they reflect the configuration of sub-domain that is
        !!              overseen by the local processor.
    subroutine s_mpi_decompose_computational_domain() ! --------------------

#ifdef MFC_MPI

        ! # of processors in the x-, y- and z-coordinate directions
        integer :: num_procs_x, num_procs_y

        ! Temporary # of processors in x-, y- and z-coordinate directions
        ! used during the processor factorization optimization procedure
        real(kind(0d0)) :: tmp_num_procs_x, tmp_num_procs_y

        ! Processor factorization (fct) minimization parameter
        real(kind(0d0)) :: fct_min

        ! Cartesian processor topology communicator
        integer :: MPI_COMM_CART

        ! Number of remaining cells for a particular coordinate direction
        ! after the bulk has evenly been distributed among the available
        ! processors for that coordinate direction
        integer :: rem_cells

        ! Generic loop iterators
        integer :: i

        if (num_procs == 1 .and. parallel_io) then
            do i = 1, num_dims
                start_idx(i) = 0
            end do
            return
        end if

        if (n > 0) then


            ! Initial values of the processor factorization optimization
            num_procs_x = 1
            num_procs_y = num_procs
            ierr = -1

            ! Computing minimization variable for these initial values
            tmp_num_procs_x = num_procs_x
            tmp_num_procs_y = num_procs_y
            fct_min = 10d0*abs((m + 1)/tmp_num_procs_x &
                               - (n + 1)/tmp_num_procs_y)

            ! Searching for optimal computational domain distribution
            do i = 1, num_procs

                if (mod(num_procs, i) == 0 &
                    .and. &
                    (m + 1)/i >= num_stcls_min*weno_order) then

                    tmp_num_procs_x = i
                    tmp_num_procs_y = num_procs/i

                    if (fct_min >= abs((m + 1)/tmp_num_procs_x &
                                       - (n + 1)/tmp_num_procs_y) &
                        .and. &
                        (n + 1)/tmp_num_procs_y &
                        >= &
                        num_stcls_min*weno_order) then

                        num_procs_x = i
                        num_procs_y = num_procs/i
                        fct_min = abs((m + 1)/tmp_num_procs_x &
                                      - (n + 1)/tmp_num_procs_y)
                        ierr = 0

                    end if

                end if

            end do

            ! Checking whether the decomposition of the computational
            ! domain was successful
            if (proc_rank == 0 .and. ierr == -1) then
                print '(A)', 'Unable to decompose computational '// &
                    'domain for selected number of '// &
                    'processors. Exiting ...'
                call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
            end if

            ! Creating a new communicator using Cartesian topology
            call MPI_CART_CREATE(MPI_COMM_WORLD, 2, (/num_procs_x, &
                                                      num_procs_y/), (/.true., &
                                                                       .true./), .false., MPI_COMM_CART, &
                                 ierr)

            ! Finding corresponding Cartesian coordinates of the local
            ! processor rank in newly declared cartesian communicator
            call MPI_CART_COORDS(MPI_COMM_CART, proc_rank, 2, &
                                 proc_coords, ierr)


            ! END: Generating 2D Cartesian Processor Topology ==================

            ! Sub-domain Global Parameters in y-direction ======================

            ! Number of remaining cells after majority has been distributed
            rem_cells = mod(n + 1, num_procs_y)

            ! Preliminary uniform cell-width spacing
            if (old_grid .neqv. .true.) then
                dy = (y_domain%end - y_domain%beg)/real(n + 1, kind(0d0))
            end if

            ! Optimal number of cells per processor
            n = (n + 1)/num_procs_y - 1

            ! Distributing any remaining cells
            do i = 1, rem_cells
                if (proc_coords(2) == i - 1) then
                    n = n + 1
                    exit
                end if
            end do

            ! Beginning and end sub-domain boundary locations
            if (parallel_io .neqv. .true.) then
                if (old_grid .neqv. .true.) then
                    if (proc_coords(2) < rem_cells) then
                        y_domain%beg = y_domain%beg + dy*real((n + 1)* &
                                                              proc_coords(2))
                        y_domain%end = y_domain%end - dy*real((n + 1)* &
                                                              (num_procs_y - proc_coords(2) - 1) &
                                                              - (num_procs_y - rem_cells))
                    else
                        y_domain%beg = y_domain%beg + dy*real((n + 1)* &
                                                              proc_coords(2) + rem_cells)
                        y_domain%end = y_domain%end - dy*real((n + 1)* &
                                                              (num_procs_y - proc_coords(2) - 1))
                    end if
                end if
            else
                if (proc_coords(2) < rem_cells) then
                    start_idx(2) = (n + 1)*proc_coords(2)
                else
                    start_idx(2) = (n + 1)*proc_coords(2) + rem_cells
                end if
            end if

            ! ==================================================================

            ! Generating 1D Cartesian Processor Topology =======================

        else

            ! Number of processors in the coordinate direction is equal to
            ! the total number of processors available
            num_procs_x = num_procs

            ! Creating a new communicator using Cartesian topology
            call MPI_CART_CREATE(MPI_COMM_WORLD, 1, (/num_procs_x/), &
                                 (/.true./), .false., MPI_COMM_CART, &
                                 ierr)

            ! Finding the corresponding Cartesian coordinates of the local
            ! processor rank in the newly declared cartesian communicator
            call MPI_CART_COORDS(MPI_COMM_CART, proc_rank, 1, &
                                 proc_coords, ierr)

        end if

        ! ==================================================================

        ! Sub-domain Global Parameters in x-direction ======================

        ! Number of remaining cells after majority has been distributed
        rem_cells = mod(m + 1, num_procs_x)

        ! Preliminary uniform cell-width spacing
        if (old_grid .neqv. .true.) then
            dx = (x_domain%end - x_domain%beg)/real(m + 1, kind(0d0))
        end if

        ! Optimal number of cells per processor
        m = (m + 1)/num_procs_x - 1

        ! Distributing any remaining cells
        do i = 1, rem_cells
            if (proc_coords(1) == i - 1) then
                m = m + 1
                exit
            end if
        end do

        ! Beginning and end sub-domain boundary locations
        if (parallel_io .neqv. .true.) then
            if (old_grid .neqv. .true.) then
                if (proc_coords(1) < rem_cells) then
                    x_domain%beg = x_domain%beg + dx*real((m + 1)* &
                                                          proc_coords(1))
                    x_domain%end = x_domain%end - dx*real((m + 1)* &
                                                          (num_procs_x - proc_coords(1) - 1) &
                                                          - (num_procs_x - rem_cells))
                else
                    x_domain%beg = x_domain%beg + dx*real((m + 1)* &
                                                          proc_coords(1) + rem_cells)
                    x_domain%end = x_domain%end - dx*real((m + 1)* &
                                                          (num_procs_x - proc_coords(1) - 1))
                end if
            end if
        else
            if (proc_coords(1) < rem_cells) then
                start_idx(1) = (m + 1)*proc_coords(1)
            else
                start_idx(1) = (m + 1)*proc_coords(1) + rem_cells
            end if
        end if

        ! ==================================================================


#endif

    end subroutine s_mpi_decompose_computational_domain ! ------------------


end module m_mpi_proxy
