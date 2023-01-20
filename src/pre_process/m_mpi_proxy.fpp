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
    ! ==========================================================================

    implicit none

    integer, private :: err_code, ierr !<
    !! Generic flags used to identify and report MPI errors

contains

    !>  The subroutine intializes the MPI environment and queries
        !!      both the number of processors that will be available for
        !!      the job as well as the local processor rank.
    subroutine s_mpi_initialize() ! ----------------------------

#ifndef MFC_MPI

        ! Serial run only has 1 processor
        num_procs = 1
        ! Local processor rank is 0
        proc_rank = 0

#else

        ! Establishing the MPI environment
        call MPI_INIT(ierr)

        ! Checking whether the MPI environment has been properly intialized
        if (ierr /= MPI_SUCCESS) then
            print '(A)', 'Unable to initialize MPI environment. Exiting ...'
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end if

        ! Querying number of processors available for the job
        call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)

        ! Identifying the rank of the local processor
        call MPI_COMM_RANK(MPI_COMM_WORLD, proc_rank, ierr)

#endif

    end subroutine s_mpi_initialize ! --------------------------

    !> The subroutine terminates the MPI execution environment.
    subroutine s_mpi_abort() ! ---------------------------------------------

#ifndef MFC_MPI

        stop 1

#else

        ! Terminating the MPI environment
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)

#endif

    end subroutine s_mpi_abort ! -------------------------------------------

        !! @param q_cons_vf Conservative variables
    subroutine s_initialize_mpi_data(q_cons_vf) ! --------------------------

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_cons_vf

#ifdef MFC_MPI

        integer, dimension(num_dims) :: sizes_glb, sizes_loc
        integer :: ierr

        ! Generic loop iterator
        integer :: i

        do i = 1, sys_size
            MPI_IO_DATA%var(i)%sf => q_cons_vf(i)%sf
        end do

        ! Define global(g) and local(l) sizes for flow variables
        sizes_glb(1) = m_glb + 1; sizes_loc(1) = m + 1
        if (n > 0) then
            sizes_glb(2) = n_glb + 1; sizes_loc(2) = n + 1
        end if

        ! Define the view for each variable
        do i = 1, sys_size
            call MPI_TYPE_CREATE_SUBARRAY(num_dims, sizes_glb, sizes_loc, start_idx, &
                                          MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, MPI_IO_DATA%view(i), ierr)
            call MPI_TYPE_COMMIT(MPI_IO_DATA%view(i), ierr)
        end do

#endif

    end subroutine s_initialize_mpi_data ! ---------------------------------

    !> Halts all processes until all have reached barrier.
    subroutine s_mpi_barrier() ! -------------------------------------------

#ifdef MFC_MPI

        ! Calling MPI_BARRIER
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

#endif

    end subroutine s_mpi_barrier ! -----------------------------------------

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

        #:for VAR in [ 'old_grid','old_ic','stretch_x','stretch_y',&
            & 'adv_alphan','mpp_lim', 'parallel_io'  ]
            call MPI_BCAST(${VAR}$, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        #:endfor
        call MPI_BCAST(fluid_rho(1), num_fluids_max, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)


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

    subroutine mpi_bcast_time_step_values(proc_time, time_avg)

        real(kind(0d0)), dimension(0:num_procs - 1), intent(INOUT) :: proc_time
        real(kind(0d0)), intent(INOUT) :: time_avg

#ifdef MFC_MPI

        integer :: j

        call MPI_GATHER(time_avg, 1, MPI_DOUBLE_PRECISION, proc_time(0), 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

#endif

    end subroutine mpi_bcast_time_step_values

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
        integer :: i, j

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

    !>  The following subroutine takes the inputted variable and
        !!      determines its minimum value on the entire computational
        !!      domain. The result is stored back into inputted variable.
        !!  @param var_loc holds the local value to be reduced among
        !!      all the processors in communicator. On output, the variable holds
        !!      the minimum value, reduced amongst all of the local values.
    subroutine s_mpi_reduce_min(var_loc) ! ---------------------------------

        real(kind(0d0)), intent(INOUT) :: var_loc

#ifdef MFC_MPI

        ! Temporary storage variable that holds the reduced minimum value
        real(kind(0d0)) :: var_glb

        ! Performing reduction procedure and eventually storing its result
        ! into the variable that was initially inputted into the subroutine
        call MPI_REDUCE(var_loc, var_glb, 1, MPI_DOUBLE_PRECISION, &
                        MPI_MIN, 0, MPI_COMM_WORLD, ierr)

        call MPI_BCAST(var_glb, 1, MPI_DOUBLE_PRECISION, &
                       0, MPI_COMM_WORLD, ierr)

        var_loc = var_glb

#endif

    end subroutine s_mpi_reduce_min ! --------------------------------------

    !> Finalization of all MPI related processes
    subroutine s_mpi_finalize() ! ------------------------------

#ifdef MFC_MPI

        ! Terminating the MPI environment
        call MPI_FINALIZE(ierr)

#endif

    end subroutine s_mpi_finalize ! ----------------------------

end module m_mpi_proxy
