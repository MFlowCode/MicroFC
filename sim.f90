    subroutine s_mpi_decompose_computational_domain() ! --------------------

#ifdef MFC_MPI

        integer :: num_procs_x, num_procs_y

        real(kind(0d0)) :: tmp_num_procs_x, tmp_num_procs_y

        real(kind(0d0)) :: fct_min

        integer :: MPI_COMM_CART

        integer :: rem_cells

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

            call MPI_CART_CREATE(MPI_COMM_WORLD, 2, (/num_procs_x, &
                                                      num_procs_y/), (/.true., &
                                                                       .true./), .false., MPI_COMM_CART, &
                                 ierr)

            ! Finding corresponding Cartesian coordinates of the local
            ! processor rank in newly declared cartesian communicator
            call MPI_CART_COORDS(MPI_COMM_CART, proc_rank, 2, &
                                 proc_coords, ierr)



            ! Sub-domain Global Parameters in y-direction ======================

            ! Number of remaining cells after majority has been distributed
            rem_cells = mod(n + 1, num_procs_y)

            ! Optimal number of cells per processor
            n = (n + 1)/num_procs_y - 1

            ! Distributing any remaining cells
            do i = 1, rem_cells
                if (proc_coords(2) == i - 1) then
                    n = n + 1
                    exit
                end if
            end do

            ! Boundary condition at the beginning
            if (proc_coords(2) > 0 .or. bc_y%beg == -1) then
                proc_coords(2) = proc_coords(2) - 1
                call MPI_CART_RANK(MPI_COMM_CART, proc_coords, bc_y%beg, &
                                   ierr)
                proc_coords(2) = proc_coords(2) + 1
            end if

            ! Ghost zone at the beginning
            if (proc_coords(2) > 0 .and. format == 1) then
                offset_y%beg = 2
            else
                offset_y%beg = 0
            end if

            ! Boundary condition at the end
            if (proc_coords(2) < num_procs_y - 1 .or. bc_y%end == -1) then
                proc_coords(2) = proc_coords(2) + 1
                call MPI_CART_RANK(MPI_COMM_CART, proc_coords, bc_y%end, &
                                   ierr)
                proc_coords(2) = proc_coords(2) - 1
            end if

            ! Ghost zone at the end
            if (proc_coords(2) < num_procs_y - 1 .and. format == 1) then
                offset_y%end = 2
            else
                offset_y%end = 0
            end if

            if (parallel_io) then
                if (proc_coords(2) < rem_cells) then
                    start_idx(2) = (n + 1)*proc_coords(2)
                else
                    start_idx(2) = (n + 1)*proc_coords(2) + rem_cells
                end if
            end if
            ! ==================================================================


        else

            ! Number of processors in the coordinate direction is equal to
            ! the total number of processors available
            num_procs_x = num_procs

            ! Number of cells in undecomposed computational domain needed
            ! for sub-domain reassembly during formatted data output
            m_root = m

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

        ! Optimal number of cells per processor
        m = (m + 1)/num_procs_x - 1

        ! Distributing any remaining cells
        do i = 1, rem_cells
            if (proc_coords(1) == i - 1) then
                m = m + 1
                exit
            end if
        end do

        ! Boundary condition at the beginning
        if (proc_coords(1) > 0 .or. bc_x%beg == -1) then
            proc_coords(1) = proc_coords(1) - 1
            call MPI_CART_RANK(MPI_COMM_CART, proc_coords, bc_x%beg, ierr)
            proc_coords(1) = proc_coords(1) + 1
        end if

        ! Ghost zone at the beginning
        if (proc_coords(1) > 0 .and. format == 1 .and. n > 0) then
            offset_x%beg = 2
        else
            offset_x%beg = 0
        end if

        ! Boundary condition at the end
        if (proc_coords(1) < num_procs_x - 1 .or. bc_x%end == -1) then
            proc_coords(1) = proc_coords(1) + 1
            call MPI_CART_RANK(MPI_COMM_CART, proc_coords, bc_x%end, ierr)
            proc_coords(1) = proc_coords(1) - 1
        end if

        ! Ghost zone at the end
        if (proc_coords(1) < num_procs_x - 1 .and. format == 1 .and. n > 0) then
            offset_x%end = 2
        else
            offset_x%end = 0
        end if

        if (parallel_io) then
            if (proc_coords(1) < rem_cells) then
                start_idx(1) = (m + 1)*proc_coords(1)
            else
                start_idx(1) = (m + 1)*proc_coords(1) + rem_cells
            end if
        end if
        ! ==================================================================

#endif

    end subroutine s_mpi_decompose_computational_domain ! ------------------
