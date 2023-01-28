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

            num_procs_x = 1
            num_procs_y = num_procs
            ierr = -1

            ! Benchmarking the quality of this initial guess
            tmp_num_procs_x = num_procs_x
            tmp_num_procs_y = num_procs_y
            fct_min = 10d0*abs((m + 1)/tmp_num_procs_x &
                               - (n + 1)/tmp_num_procs_y)

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

            ! Verifying that a valid decomposition of the computational
            ! domain has been established. If not, the simulation exits.
            if (proc_rank == 0 .and. ierr == -1) then
                print '(A)', 'Unsupported combination of values '// &
                    'of num_procs, m, n and '// &
                    'weno_order. Exiting ...'
                    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
            end if

            call MPI_CART_CREATE(MPI_COMM_WORLD, 2, (/num_procs_x, &
                                                      num_procs_y/), (/.true., &
                                                                       .true./), .false., MPI_COMM_CART, &
                                 ierr)

            ! Finding the Cartesian coordinates of the local process
            call MPI_CART_COORDS(MPI_COMM_CART, proc_rank, 2, &
                                 proc_coords, ierr)


            ! Global Parameters for y-direction ================================

            ! Number of remaining cells
            rem_cells = mod(n + 1, num_procs_y)

            ! Optimal number of cells per processor
            n = (n + 1)/num_procs_y - 1

            ! Distributing the remaining cells
            do i = 1, rem_cells
                if (proc_coords(2) == i - 1) then
                    n = n + 1; exit
                end if
            end do

            ! Boundary condition at the beginning
            if (proc_coords(2) > 0 .or. bc_y%beg == -1) then
                proc_coords(2) = proc_coords(2) - 1
                call MPI_CART_RANK(MPI_COMM_CART, proc_coords, &
                                   bc_y%beg, ierr)
                proc_coords(2) = proc_coords(2) + 1
            end if

            ! Boundary condition at the end
            if (proc_coords(2) < num_procs_y - 1 .or. bc_y%end == -1) then
                proc_coords(2) = proc_coords(2) + 1
                call MPI_CART_RANK(MPI_COMM_CART, proc_coords, &
                                   bc_y%end, ierr)
                proc_coords(2) = proc_coords(2) - 1
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

            num_procs_x = num_procs

            call MPI_CART_CREATE(MPI_COMM_WORLD, 1, (/num_procs_x/), &
                                 (/.true./), .false., MPI_COMM_CART, &
                                 ierr)

            ! Finding the Cartesian coordinates of the local process
            call MPI_CART_COORDS(MPI_COMM_CART, proc_rank, 1, &
                                 proc_coords, ierr)

        end if
        ! ==================================================================

        ! Global Parameters for x-direction ================================

        ! Number of remaining cells
        rem_cells = mod(m + 1, num_procs_x)

        ! Optimal number of cells per processor
        m = (m + 1)/num_procs_x - 1

        ! Distributing the remaining cells
        do i = 1, rem_cells
            if (proc_coords(1) == i - 1) then
                m = m + 1; exit
            end if
        end do

        ! Boundary condition at the beginning
        if (proc_coords(1) > 0 .or. bc_x%beg == -1) then
            proc_coords(1) = proc_coords(1) - 1
            call MPI_CART_RANK(MPI_COMM_CART, proc_coords, bc_x%beg, ierr)
            proc_coords(1) = proc_coords(1) + 1
        end if

        ! Boundary condition at the end
        if (proc_coords(1) < num_procs_x - 1 .or. bc_x%end == -1) then
            proc_coords(1) = proc_coords(1) + 1
            call MPI_CART_RANK(MPI_COMM_CART, proc_coords, bc_x%end, ierr)
            proc_coords(1) = proc_coords(1) - 1
        end if

        if (parallel_io) then
            if (proc_coords(1) < rem_cells) then
                start_idx(1) = (m + 1)*proc_coords(1)
            else
                start_idx(1) = (m + 1)*proc_coords(1) + rem_cells
            end if
        end if
        ! ==================================================================

        if (proc_rank == 0) then
            print *, m, n
        end if

#endif

    end subroutine s_mpi_decompose_computational_domain ! ------------------

