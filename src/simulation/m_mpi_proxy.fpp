!>
!! @file m_mpi_proxy.f90
!! @brief Contains module m_mpi_proxy

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief The module serves as a proxy to the parameters and subroutines
!!          available in the MPI implementation's MPI module. Specifically,
!!          the purpose of the proxy is to harness basic MPI commands into
!!          more complicated procedures as to accomplish the communication
!!          goals for the simulation.
module m_mpi_proxy

    ! Dependencies =============================================================
#ifdef MFC_MPI
    use mpi                    !< Message passing interface (MPI) module
#endif

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters
    ! ==========================================================================

    implicit none

    real(kind(0d0)), private, allocatable, dimension(:) :: q_cons_buff_send !<
    !! This variable is utilized to pack and send the buffer of the cell-average
    !! conservative variables, for a single computational domain boundary at the
    !! time, to the relevant neighboring processor.

    real(kind(0d0)), private, allocatable, dimension(:) :: q_cons_buff_recv !<
    !! q_cons_buff_recv is utilized to receive and unpack the buffer of the cell-
    !! average conservative variables, for a single computational domain boundary
    !! at the time, from the relevant neighboring processor.

    !> @name Generic flags used to identify and report MPI errors
    !> @{
    integer, private :: err_code, ierr
    !> @}

!$acc declare create(q_cons_buff_send, q_cons_buff_recv)

    !real :: s_time, e_time
    !real :: compress_time, mpi_time, decompress_time
    !integer :: nCalls_time = 0

contains

    !> The subroutine intializes the MPI execution environment
        !!      and queries both the number of processors which will be
        !!      available for the job and the local processor rank.
    subroutine s_mpi_initialize() ! ----------------------------------------

#ifndef MFC_MPI

        ! Serial run only has 1 processor
        num_procs = 1
        ! Local processor rank is 0
        proc_rank = 0

#else

        ! Initializing the MPI environment
        call MPI_INIT(ierr)

        ! Checking whether the MPI environment has been properly intialized
        if (ierr /= MPI_SUCCESS) then
            print '(A)', 'Unable to initialize MPI environment. Exiting ...'
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end if

        ! Querying the number of processors available for the job
        call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)

        ! Querying the rank of the local processor
        call MPI_COMM_RANK(MPI_COMM_WORLD, proc_rank, ierr)

#endif

    end subroutine s_mpi_initialize ! --------------------------------------

    !> The subroutine terminates the MPI execution environment.
    subroutine s_mpi_abort() ! ---------------------------------------------

#ifndef MFC_MPI

        stop 1

#else

        ! Terminating the MPI environment
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)

#endif

    end subroutine s_mpi_abort ! -------------------------------------------

    !> The subroutine that initializes MPI data structures
        !!  @param q_cons_vf Conservative variables
    subroutine s_initialize_mpi_data(q_cons_vf) ! --------------------------

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_cons_vf

        integer, dimension(num_dims) :: sizes_glb, sizes_loc
        integer :: ierr

        integer :: i !< Generic loop iterator

#ifdef MFC_MPI

        do i = 1, sys_size
            MPI_IO_DATA%var(i)%sf => q_cons_vf(i)%sf(0:m, 0:n)
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

    !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_mpi_proxy_module() ! ---------------------------

#ifdef MFC_MPI

        ! Allocating q_cons_buff_send and q_cons_buff_recv. Please note that
        ! for the sake of simplicity, both variables are provided sufficient
        ! storage to hold the largest buffer in the computational domain.

        if (n > 0) then
                @:ALLOCATE(q_cons_buff_send(0:-1 + buff_size*sys_size* &
                                         & (max(m, n) + 2*buff_size + 1)))
        else
            @:ALLOCATE(q_cons_buff_send(0:-1 + buff_size*sys_size))
        end if

        @:ALLOCATE(q_cons_buff_recv(0:ubound(q_cons_buff_send, 1)))

#endif

    end subroutine s_initialize_mpi_proxy_module ! -------------------------

    !>  Since only the processor with rank 0 reads and verifies
        !!      the consistency of user inputs, these are initially not
        !!      available to the other processors. Then, the purpose of
        !!      this subroutine is to distribute the user inputs to the
        !!      remaining processors in the communicator.
    subroutine s_mpi_bcast_user_inputs() ! ---------------------------------

#ifdef MFC_MPI

        integer :: i, j !< Generic loop iterator

        call MPI_BCAST(case_dir, len(case_dir), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

        #:for VAR in ['t_step_old', 'm', 'n', 'm_glb', 'n_glb',   &
            & 't_step_start','t_step_stop','t_step_save',         &
            & 'num_fluids','time_stepper',  & 
            & 'wave_speeds', 'avg_state', 'precision', 'bc_x%beg', 'bc_x%end', & 
            & 'bc_y%beg', 'bc_y%end', 'bc_x_glb%beg',  & 
            & 'bc_x_glb%end', 'bc_y_glb%beg', 'bc_y_glb%end',  &
            &  'fd_order', 'num_probes' ]
            call MPI_BCAST(${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in [ 'run_time_info', 'cu_mpi', 'weno_Re_flux', 'parallel_io', 'probe_wrt' ]
            call MPI_BCAST(${VAR}$, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in [ 'dt','weno_eps' ]
            call MPI_BCAST(${VAR}$, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:if not MFC_CASE_OPTIMIZATION
            call MPI_BCAST(weno_order, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        #:endif

        do i = 1, num_fluids_max
            #:for VAR in [ 'gamma','pi_inf' ]
                call MPI_BCAST(fluid_pp(i)%${VAR}$, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            call MPI_BCAST(fluid_pp(i)%Re(1), 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        end do

        do j = 1, num_probes_max
            #:for VAR in [ 'x','y' ]
                call MPI_BCAST(probe(j)%${VAR}$, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
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

    !>  The purpose of this procedure is to optimally decompose
        !!      the computational domain among the available processors.
        !!      This is performed by attempting to award each processor,
        !!      in each of the coordinate directions, approximately the
        !!      same number of cells, and then recomputing the affected
        !!      global parameters.
    subroutine s_mpi_decompose_computational_domain() ! --------------------

#ifdef MFC_MPI

        integer :: num_procs_x, num_procs_y !<
            !! Optimal number of processors in the x-, y- and z-directions

        real(kind(0d0)) :: tmp_num_procs_x, tmp_num_procs_y !<
            !! Non-optimal number of processors in the x-, y- and z-directions

        real(kind(0d0)) :: fct_min !<
            !! Processor factorization (fct) minimization parameter

        integer :: MPI_COMM_CART !<
            !! Cartesian processor topology communicator

        integer :: rem_cells !<
            !! Remaining number of cells, in a particular coordinate direction,
            !! after the majority is divided up among the available processors

        integer :: i, j !< Generic loop iterators

        if (num_procs == 1 .and. parallel_io) then
            do i = 1, num_dims
                start_idx(i) = 0
            end do
            return
        end if

        if (n > 0) then

            ! Initial estimate of optimal processor topology
            num_procs_x = 1
            num_procs_y = num_procs
            ierr = -1

            ! Benchmarking the quality of this initial guess
            tmp_num_procs_x = num_procs_x
            tmp_num_procs_y = num_procs_y
            fct_min = 10d0*abs((m + 1)/tmp_num_procs_x &
                               - (n + 1)/tmp_num_procs_y)

            ! Optimization of the initial processor topology
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
                call s_mpi_abort()
            end if

            ! Creating new communicator using the Cartesian topology
            call MPI_CART_CREATE(MPI_COMM_WORLD, 2, (/num_procs_x, &
                                                      num_procs_y/), (/.true., &
                                                                       .true./), .false., MPI_COMM_CART, &
                                 ierr)

            ! Finding the Cartesian coordinates of the local process
            call MPI_CART_COORDS(MPI_COMM_CART, proc_rank, 2, &
                                 proc_coords, ierr)

            ! END: 2D Cartesian Processor Topology =============================

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

            ! 1D Cartesian Processor Topology ==================================
        else

            ! Optimal processor topology
            num_procs_x = num_procs

            ! Creating new communicator using the Cartesian topology
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

    !>  The goal of this procedure is to populate the buffers of
        !!      the grid variables by communicating with the neighboring
        !!      processors. Note that only the buffers of the cell-width
        !!      distributions are handled in such a way. This is because
        !!      the buffers of cell-boundary locations may be calculated
        !!      directly from those of the cell-width distributions.
        !!  @param mpi_dir MPI communication coordinate direction
        !!  @param pbc_loc Processor boundary condition (PBC) location
    subroutine s_mpi_sendrecv_grid_variables_buffers(mpi_dir, pbc_loc) ! ---

        integer, intent(IN) :: mpi_dir
        integer, intent(IN) :: pbc_loc

#ifdef MFC_MPI

        ! MPI Communication in x-direction =================================
        if (mpi_dir == 1) then

            if (pbc_loc == -1) then      ! PBC at the beginning

                if (bc_x%end >= 0) then      ! PBC at the beginning and end

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        dx(m - buff_size + 1), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_x%end, 0, &
                        dx(-buff_size), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_x%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                else                        ! PBC at the beginning only

                    ! Send/receive buffer to/from bc_x%beg/bc_x%beg
                    call MPI_SENDRECV( &
                        dx(0), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_x%beg, 1, &
                        dx(-buff_size), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_x%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                end if

            else                        ! PBC at the end

                if (bc_x%beg >= 0) then      ! PBC at the end and beginning

                    ! Send/receive buffer to/from bc_x%beg/bc_x%end
                    call MPI_SENDRECV( &
                        dx(0), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_x%beg, 1, &
                        dx(m + 1), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_x%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                else                        ! PBC at the end only

                    ! Send/receive buffer to/from bc_x%end/bc_x%end
                    call MPI_SENDRECV( &
                        dx(m - buff_size + 1), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_x%end, 0, &
                        dx(m + 1), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_x%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                end if

            end if
            ! END: MPI Communication in x-direction ============================

            ! MPI Communication in y-direction =================================
        elseif (mpi_dir == 2) then

            if (pbc_loc == -1) then      ! PBC at the beginning

                if (bc_y%end >= 0) then      ! PBC at the beginning and end

                    ! Send/receive buffer to/from bc_y%end/bc_y%beg
                    call MPI_SENDRECV( &
                        dy(n - buff_size + 1), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_y%end, 0, &
                        dy(-buff_size), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_y%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                else                        ! PBC at the beginning only

                    ! Send/receive buffer to/from bc_y%beg/bc_y%beg
                    call MPI_SENDRECV( &
                        dy(0), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_y%beg, 1, &
                        dy(-buff_size), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_y%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                end if

            else                        ! PBC at the end

                if (bc_y%beg >= 0) then      ! PBC at the end and beginning

                    ! Send/receive buffer to/from bc_y%beg/bc_y%end
                    call MPI_SENDRECV( &
                        dy(0), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_y%beg, 1, &
                        dy(n + 1), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_y%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                else                        ! PBC at the end only

                    ! Send/receive buffer to/from bc_y%end/bc_y%end
                    call MPI_SENDRECV( &
                        dy(n - buff_size + 1), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_y%end, 0, &
                        dy(n + 1), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_y%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                end if

            end if
            ! END: MPI Communication in y-direction ============================

        end if
        ! END: MPI Communication in z-direction ============================

#endif

    end subroutine s_mpi_sendrecv_grid_variables_buffers ! -----------------

    !>  The goal of this subroutine is to determine the global
        !!      extrema of the stability criteria in the computational
        !!      domain. This is performed by sifting through the local
        !!      extrema of each stability criterion. Note that each of
        !!      the local extrema is from a single process, within its
        !!      assigned section of the computational domain. Finally,
        !!      note that the global extrema values are only bookkeept
        !!      on the rank 0 processor.
        !!  @param icfl_max_loc Local maximum ICFL stability criterion
        !!  @param vcfl_max_loc Local maximum VCFL stability criterion
        !!  @param ccfl_max_loc Local maximum CCFL stability criterion
        !!  @param Rc_min_loc Local minimum Rc stability criterion
        !!  @param icfl_max_glb Global maximum ICFL stability criterion
        !!  @param vcfl_max_glb Global maximum VCFL stability criterion
        !!  @param ccfl_max_glb Global maximum CCFL stability criterion
        !!  @param Rc_min_glb Global minimum Rc stability criterion
    subroutine s_mpi_reduce_stability_criteria_extrema(icfl_max_loc, & ! --
                                                       vcfl_max_loc, &
                                                       ccfl_max_loc, &
                                                       Rc_min_loc, &
                                                       icfl_max_glb, &
                                                       vcfl_max_glb, &
                                                       ccfl_max_glb, &
                                                       Rc_min_glb)

        real(kind(0d0)), intent(IN) :: icfl_max_loc
        real(kind(0d0)), intent(IN) :: vcfl_max_loc
        real(kind(0d0)), intent(IN) :: ccfl_max_loc
        real(kind(0d0)), intent(IN) :: Rc_min_loc

        real(kind(0d0)), intent(OUT) :: icfl_max_glb
        real(kind(0d0)), intent(OUT) :: vcfl_max_glb
        real(kind(0d0)), intent(OUT) :: ccfl_max_glb
        real(kind(0d0)), intent(OUT) :: Rc_min_glb

#ifdef MFC_MPI

        ! Reducing local extrema of ICFL, VCFL, CCFL and Rc numbers to their
        ! global extrema and bookkeeping the results on the rank 0 processor
        call MPI_REDUCE(icfl_max_loc, icfl_max_glb, 1, &
                        MPI_DOUBLE_PRECISION, MPI_MAX, 0, &
                        MPI_COMM_WORLD, ierr)

        if (any(Re_size > 0)) then
            call MPI_REDUCE(vcfl_max_loc, vcfl_max_glb, 1, &
                            MPI_DOUBLE_PRECISION, MPI_MAX, 0, &
                            MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(Rc_min_loc, Rc_min_glb, 1, &
                            MPI_DOUBLE_PRECISION, MPI_MIN, 0, &
                            MPI_COMM_WORLD, ierr)
        end if

#endif

    end subroutine s_mpi_reduce_stability_criteria_extrema ! ---------------

    !>  The following subroutine takes the input local variable
        !!      from all processors and reduces to the sum of all
        !!      values. The reduced variable is recorded back onto the
        !!      original local variable on each processor.
        !!  @param var_loc Some variable containing the local value which should be
        !!  reduced amongst all the processors in the communicator.
        !!  @param var_glb The globally reduced value
    subroutine s_mpi_allreduce_sum(var_loc, var_glb) ! ---------------------

        real(kind(0d0)), intent(IN) :: var_loc
        real(kind(0d0)), intent(OUT) :: var_glb

#ifdef MFC_MPI

        ! Performing the reduction procedure
        call MPI_ALLREDUCE(var_loc, var_glb, 1, MPI_DOUBLE_PRECISION, &
                           MPI_SUM, MPI_COMM_WORLD, ierr)

#endif

    end subroutine s_mpi_allreduce_sum ! -----------------------------------

    !>  The following subroutine takes the input local variable
        !!      from all processors and reduces to the minimum of all
        !!      values. The reduced variable is recorded back onto the
        !!      original local variable on each processor.
        !!  @param var_loc Some variable containing the local value which should be
        !!  reduced amongst all the processors in the communicator.
        !!  @param var_glb The globally reduced value
    subroutine s_mpi_allreduce_min(var_loc, var_glb) ! ---------------------

        real(kind(0d0)), intent(IN) :: var_loc
        real(kind(0d0)), intent(OUT) :: var_glb

#ifdef MFC_MPI

        ! Performing the reduction procedure
        call MPI_ALLREDUCE(var_loc, var_glb, 1, MPI_DOUBLE_PRECISION, &
                           MPI_MIN, MPI_COMM_WORLD, ierr)

#endif

    end subroutine s_mpi_allreduce_min ! -----------------------------------

    !>  The following subroutine takes the input local variable
        !!      from all processors and reduces to the maximum of all
        !!      values. The reduced variable is recorded back onto the
        !!      original local variable on each processor.
        !!  @param var_loc Some variable containing the local value which should be
        !!  reduced amongst all the processors in the communicator.
        !!  @param var_glb The globally reduced value
    subroutine s_mpi_allreduce_max(var_loc, var_glb) ! ---------------------

        real(kind(0d0)), intent(IN) :: var_loc
        real(kind(0d0)), intent(OUT) :: var_glb

#ifdef MFC_MPI

        ! Performing the reduction procedure
        call MPI_ALLREDUCE(var_loc, var_glb, 1, MPI_DOUBLE_PRECISION, &
                           MPI_MAX, MPI_COMM_WORLD, ierr)

#endif

    end subroutine s_mpi_allreduce_max ! -----------------------------------

    !>  The goal of this procedure is to populate the buffers of
        !!      the cell-average conservative variables by communicating
        !!      with the neighboring processors.
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param mpi_dir MPI communication coordinate direction
        !!  @param pbc_loc Processor boundary condition (PBC) location
    subroutine s_mpi_sendrecv_conservative_variables_buffers(q_cons_vf, &
                                                             mpi_dir, &
                                                             pbc_loc)

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        integer, intent(IN) :: mpi_dir
        integer, intent(IN) :: pbc_loc

        integer :: i, j, k, r !< Generic loop iterators

#ifdef MFC_MPI

        !nCalls_time = nCalls_time + 1

        ! MPI Communication in x-direction =================================
        if (mpi_dir == 1) then

            if (pbc_loc == -1) then      ! PBC at the beginning

                if (bc_x%end >= 0) then      ! PBC at the beginning and end

                    ! Packing buffer to be sent to bc_x%end
!$acc parallel loop collapse(3) gang vector default(present) private(r)
                        do k = 0, n
                            do j = m - buff_size + 1, m
                                do i = 1, sys_size
                                    r = (i - 1) + sys_size* &
                                        ((j - m - 1) + buff_size*(k + 1))
                                    q_cons_buff_send(r) = q_cons_vf(i)%sf(j, k)
                                end do
                            end do
                        end do

                    !call MPI_Barrier(MPI_COMM_WORLD, ierr)

#if defined(_OPENACC) && defined(__PGI)
                    if (cu_mpi) then
!$acc host_data use_device( q_cons_buff_recv, q_cons_buff_send )

                        ! Send/receive buffer to/from bc_x%end/bc_x%beg
                        call MPI_SENDRECV( &
                            q_cons_buff_send(0), &
                            buff_size*sys_size*(n + 1), &
                            MPI_DOUBLE_PRECISION, bc_x%end, 0, &
                            q_cons_buff_recv(0), &
                            buff_size*sys_size*(n + 1), &
                            MPI_DOUBLE_PRECISION, bc_x%beg, 0, &
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

!$acc end host_data
!$acc wait
                    else
#endif // #if defined(_OPENACC) && defined(__PGI)
                        !call cpu_time(s_time)
!$acc update host(q_cons_buff_send)
                        !call cpu_time(e_time)
                        !compress_time = compress_time + (e_time - s_time)

                        !call cpu_time(s_time)

! Send/receive buffer to/from bc_x%end/bc_x%beg
                        call MPI_SENDRECV( &
                            q_cons_buff_send(0), &
                            buff_size*sys_size*(n + 1), &
                            MPI_DOUBLE_PRECISION, bc_x%end, 0, &
                            q_cons_buff_recv(0), &
                            buff_size*sys_size*(n + 1), &
                            MPI_DOUBLE_PRECISION, bc_x%beg, 0, &
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                        !call cpu_time(e_time)
                        !mpi_time = mpi_time + (e_time - s_time)

#if defined(_OPENACC) && defined(__PGI)
                    end if
#endif

                else                        ! PBC at the beginning only

                    ! Packing buffer to be sent to bc_x%beg
!$acc parallel loop collapse(3) gang vector default(present) private(r)
                        do k = 0, n
                            do j = 0, buff_size - 1
                                do i = 1, sys_size
                                    r = (i - 1) + sys_size* &
                                        (j + buff_size*k)
                                    q_cons_buff_send(r) = q_cons_vf(i)%sf(j, k)
                                end do
                            end do
                        end do

                    !call MPI_Barrier(MPI_COMM_WORLD, ierr)

#if defined(_OPENACC) && defined(__PGI)
                    if (cu_mpi) then
!$acc host_data use_device( q_cons_buff_recv, q_cons_buff_send )

                        ! Send/receive buffer to/from bc_x%end/bc_x%beg
                        call MPI_SENDRECV( &
                            q_cons_buff_send(0), &
                            buff_size*sys_size*(n + 1), &
                            MPI_DOUBLE_PRECISION, bc_x%beg, 1, &
                            q_cons_buff_recv(0), &
                            buff_size*sys_size*(n + 1), &
                            MPI_DOUBLE_PRECISION, bc_x%beg, 0, &
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

!$acc end host_data
!$acc wait
                    else
#endif
                        !call cpu_time(s_time)
!$acc update host(q_cons_buff_send)
                        !call cpu_time(e_time)
                        !compress_time = compress_time + (e_time - s_time)

                        !call cpu_time(s_time)
                        ! Send/receive buffer to/from bc_x%end/bc_x%beg
                        call MPI_SENDRECV( &
                            q_cons_buff_send(0), &
                            buff_size*sys_size*(n + 1), &
                            MPI_DOUBLE_PRECISION, bc_x%beg, 1, &
                            q_cons_buff_recv(0), &
                            buff_size*sys_size*(n + 1), &
                            MPI_DOUBLE_PRECISION, bc_x%beg, 0, &
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                        !call cpu_time(e_time)
                        !mpi_time = mpi_time + (e_time - s_time)

#if defined(_OPENACC) && defined(__PGI)
                    end if
#endif

                end if

#if defined(_OPENACC) && defined(__PGI)
                if (cu_mpi .eqv. .false.) then
                    !call cpu_time(s_time)
!$acc update device(q_cons_buff_recv)
                    !call cpu_time(e_time)
                    !decompress_time = decompress_time + (e_time - s_time)
                end if
#endif

                ! Unpacking buffer received from bc_x%beg
!$acc parallel loop collapse(3) gang vector default(present) private(r)
                    do k = 0, n
                        do j = -buff_size, -1
                            do i = 1, sys_size
                                r = (i - 1) + sys_size * &
                                    (j + buff_size*((k + 1)))
                                q_cons_vf(i)%sf(j, k) = q_cons_buff_recv(r)
                            end do
                        end do
                    end do

            else                        ! PBC at the end

                if (bc_x%beg >= 0) then      ! PBC at the end and beginning

!$acc parallel loop collapse(3) gang vector default(present) private(r)
                    ! Packing buffer to be sent to bc_x%beg
                        do k = 0, n
                            do j = 0, buff_size - 1
                                do i = 1, sys_size
                                    r = (i - 1) + sys_size* &
                                        (j + buff_size*(k))
                                    q_cons_buff_send(r) = q_cons_vf(i)%sf(j, k)
                                end do
                            end do
                        end do

                    !call MPI_Barrier(MPI_COMM_WORLD, ierr)

#if defined(_OPENACC) && defined(__PGI)
                    if (cu_mpi) then
!$acc host_data use_device( q_cons_buff_recv, q_cons_buff_send )

                        ! Send/receive buffer to/from bc_x%end/bc_x%beg
                        call MPI_SENDRECV( &
                            q_cons_buff_send(0), &
                            buff_size*sys_size*(n + 1), &
                            MPI_DOUBLE_PRECISION, bc_x%beg, 1, &
                            q_cons_buff_recv(0), &
                            buff_size*sys_size*(n + 1), &
                            MPI_DOUBLE_PRECISION, bc_x%end, 1, &
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

!$acc end host_data
!$acc wait
                    else
#endif
                        !call cpu_time(s_time)
!$acc update host(q_cons_buff_send)
                        !call cpu_time(e_time)
                        !compress_time = compress_time + (e_time - s_time)

                        !call cpu_time(s_time)
                        ! Send/receive buffer to/from bc_x%end/bc_x%beg
                        call MPI_SENDRECV( &
                            q_cons_buff_send(0), &
                            buff_size*sys_size*(n + 1), &
                            MPI_DOUBLE_PRECISION, bc_x%beg, 1, &
                            q_cons_buff_recv(0), &
                            buff_size*sys_size*(n + 1), &
                            MPI_DOUBLE_PRECISION, bc_x%end, 1, &
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                        !call cpu_time(e_time)
                        !mpi_time = mpi_time + (e_time - s_time)

#if defined(_OPENACC) && defined(__PGI)
                    end if
#endif

                else                        ! PBC at the end only

                    ! Packing buffer to be sent to bc_x%end
!$acc parallel loop collapse(3) gang vector default(present) private(r)
                        do k = 0, n
                            do j = m - buff_size + 1, m
                                do i = 1, sys_size
                                    r = (i - 1) + sys_size* &
                                        ((j - m - 1) + buff_size*((k + 1)))
                                    q_cons_buff_send(r) = q_cons_vf(i)%sf(j, k)
                                end do
                            end do
                        end do

                    !call MPI_Barrier(MPI_COMM_WORLD, ierr)

#if defined(_OPENACC) && defined(__PGI)
                    if (cu_mpi) then
!$acc host_data use_device( q_cons_buff_recv, q_cons_buff_send )

                        ! Send/receive buffer to/from bc_x%end/bc_x%beg
                        call MPI_SENDRECV( &
                            q_cons_buff_send(0), &
                            buff_size*sys_size*(n + 1), &
                            MPI_DOUBLE_PRECISION, bc_x%end, 0, &
                            q_cons_buff_recv(0), &
                            buff_size*sys_size*(n + 1), &
                            MPI_DOUBLE_PRECISION, bc_x%end, 1, &
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

!$acc end host_data
!$acc wait
                    else
#endif
                        !call cpu_time(s_time)
!$acc update host(q_cons_buff_send)
                        !call cpu_time(e_time)
                        !compress_time = compress_time + (e_time - s_time)

                        !call cpu_time(s_time)
                        ! Send/receive buffer to/from bc_x%end/bc_x%beg
                        call MPI_SENDRECV( &
                            q_cons_buff_send(0), &
                            buff_size*sys_size*(n + 1), &
                            MPI_DOUBLE_PRECISION, bc_x%end, 0, &
                            q_cons_buff_recv(0), &
                            buff_size*sys_size*(n + 1), &
                            MPI_DOUBLE_PRECISION, bc_x%end, 1, &
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                        !call cpu_time(e_time)
                        !mpi_time = mpi_time + (e_time - s_time)

#if defined(_OPENACC) && defined(__PGI)
                    end if
#endif

                end if

                if (cu_mpi .eqv. .false.) then
                    !call cpu_time(s_time)
!$acc update device(q_cons_buff_recv)
                    !call cpu_time(e_time)
                    !decompress_time = decompress_time + (e_time - s_time)
                end if

                ! Unpacking buffer received from bc_x%end
!$acc parallel loop collapse(3) gang vector default(present) private(r)
                    do k = 0, n
                        do j = m + 1, m + buff_size
                            do i = 1, sys_size
                                r = (i - 1) + sys_size* &
                                    ((j - m - 1) + buff_size*(k))
                                q_cons_vf(i)%sf(j, k) = q_cons_buff_recv(r)
                            end do
                        end do
                    end do

            end if
            ! END: MPI Communication in x-direction ============================

            ! MPI Communication in y-direction =================================
        elseif (mpi_dir == 2) then

            if (pbc_loc == -1) then      ! PBC at the beginning

                if (bc_y%end >= 0) then      ! PBC at the beginning and end

                    ! Packing buffer to be sent to bc_y%end
!$acc parallel loop collapse(3) gang vector default(present) private(r)
                    do i = 1, sys_size
                            do k = n - buff_size + 1, n
                                do j = -buff_size, m + buff_size
                                    r = (i - 1) + sys_size* &
                                        ((j + buff_size) + (m + 2*buff_size + 1)* &
                                         ((k - n + buff_size - 1)))
                                    q_cons_buff_send(r) = q_cons_vf(i)%sf(j, k)
                                end do
                            end do
                    end do

                    !call MPI_Barrier(MPI_COMM_WORLD, ierr)

#if defined(_OPENACC) && defined(__PGI)
                    if (cu_mpi) then
!$acc host_data use_device( q_cons_buff_recv, q_cons_buff_send )

                        ! Send/receive buffer to/from bc_x%end/bc_x%beg
                        call MPI_SENDRECV( &
                            q_cons_buff_send(0), &
                            buff_size*sys_size*(m + 2*buff_size + 1), &
                            MPI_DOUBLE_PRECISION, bc_y%end, 0, &
                            q_cons_buff_recv(0), &
                            buff_size*sys_size*(m + 2*buff_size + 1), &
                            MPI_DOUBLE_PRECISION, bc_y%beg, 0, &
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

!$acc end host_data
!$acc wait
                    else
#endif
!$acc update host(q_cons_buff_send)
                        call MPI_SENDRECV( &
                            q_cons_buff_send(0), &
                            buff_size*sys_size*(m + 2*buff_size + 1), &
                            MPI_DOUBLE_PRECISION, bc_y%end, 0, &
                            q_cons_buff_recv(0), &
                            buff_size*sys_size*(m + 2*buff_size + 1), &
                            MPI_DOUBLE_PRECISION, bc_y%beg, 0, &
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

#if defined(_OPENACC) && defined(__PGI)
                    end if
#endif

                else                        ! PBC at the beginning only

                    ! Packing buffer to be sent to bc_y%beg
!$acc parallel loop collapse(3) gang vector default(present) private(r)
                    do i = 1, sys_size
                            do k = 0, buff_size - 1
                                do j = -buff_size, m + buff_size
                                    r = (i - 1) + sys_size* &
                                        ((j + buff_size) + (m + 2*buff_size + 1)* &
                                         (k))
                                    q_cons_buff_send(r) = q_cons_vf(i)%sf(j, k)
                                end do
                            end do
                    end do

                    !call MPI_Barrier(MPI_COMM_WORLD, ierr)

#if defined(_OPENACC) && defined(__PGI)
                    if (cu_mpi) then
!$acc host_data use_device( q_cons_buff_recv, q_cons_buff_send )

                        ! Send/receive buffer to/from bc_x%end/bc_x%beg
                        call MPI_SENDRECV( &
                            q_cons_buff_send(0), &
                            buff_size*sys_size*(m + 2*buff_size + 1), &
                            MPI_DOUBLE_PRECISION, bc_y%beg, 1, &
                            q_cons_buff_recv(0), &
                            buff_size*sys_size*(m + 2*buff_size + 1), &
                            MPI_DOUBLE_PRECISION, bc_y%beg, 0, &
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

!$acc end host_data
!$acc wait
                    else
#endif
                        !call cpu_time(s_time)
!$acc update host(q_cons_buff_send)
                        !call cpu_time(e_time)
                        !compress_time = compress_time + (e_time - s_time)

                        !call cpu_time(s_time)
                        ! Send/receive buffer to/from bc_x%end/bc_x%beg
                        call MPI_SENDRECV( &
                            q_cons_buff_send(0), &
                            buff_size*sys_size*(m + 2*buff_size + 1), &
                            MPI_DOUBLE_PRECISION, bc_y%beg, 1, &
                            q_cons_buff_recv(0), &
                            buff_size*sys_size*(m + 2*buff_size + 1), &
                            MPI_DOUBLE_PRECISION, bc_y%beg, 0, &
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                        !call cpu_time(e_time)
                        !mpi_time = mpi_time + (e_time - s_time)

#if defined(_OPENACC) && defined(__PGI)
                    end if
#endif

                end if

#if defined(_OPENACC) && defined(__PGI)
                if (cu_mpi .eqv. .false.) then
                    !call cpu_time(s_time)
!$acc update device(q_cons_buff_recv)
                    !call cpu_time(e_time)
                    !decompress_time = decompress_time + (e_time - s_time)
                end if
#endif

                ! Unpacking buffer received from bc_y%beg
!$acc parallel loop collapse(3) gang vector default(present) private(r)
                do i = 1, sys_size
                        do k = -buff_size, -1
                            do j = -buff_size, m + buff_size
                                r = (i - 1) + sys_size* &
                                    ((j + buff_size) + (m + 2*buff_size + 1)* &
                                     ((k + buff_size)))
                                q_cons_vf(i)%sf(j, k) = q_cons_buff_recv(r)
                            end do
                        end do
                end do

            else                        ! PBC at the end

                if (bc_y%beg >= 0) then      ! PBC at the end and beginning

                    ! Packing buffer to be sent to bc_y%beg
!$acc parallel loop collapse(3) gang vector default(present) private(r)
                    do i = 1, sys_size
                            do k = 0, buff_size - 1
                                do j = -buff_size, m + buff_size
                                    r = (i - 1) + sys_size* &
                                        ((j + buff_size) + (m + 2*buff_size + 1)* &
                                         (k))
                                    q_cons_buff_send(r) = q_cons_vf(i)%sf(j, k)
                                end do
                            end do
                    end do

                    !call MPI_Barrier(MPI_COMM_WORLD, ierr)

#if defined(_OPENACC) && defined(__PGI)
                    if (cu_mpi) then
!$acc host_data use_device( q_cons_buff_recv, q_cons_buff_send )

                        ! Send/receive buffer to/from bc_x%end/bc_x%beg
                        call MPI_SENDRECV( &
                            q_cons_buff_send(0), &
                            buff_size*sys_size*(m + 2*buff_size + 1), &
                            MPI_DOUBLE_PRECISION, bc_y%beg, 1, &
                            q_cons_buff_recv(0), &
                            buff_size*sys_size*(m + 2*buff_size + 1), &
                            MPI_DOUBLE_PRECISION, bc_y%end, 1, &
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

!$acc end host_data
!$acc wait
                    else
#endif
                        !call cpu_time(s_time)
!$acc update host(q_cons_buff_send)
                        !call cpu_time(e_time)
                        !compress_time = compress_time + (e_time - s_time)

                        !call cpu_time(s_time)
                        ! Send/receive buffer to/from bc_x%end/bc_x%beg
                        call MPI_SENDRECV( &
                            q_cons_buff_send(0), &
                            buff_size*sys_size*(m + 2*buff_size + 1), &
                            MPI_DOUBLE_PRECISION, bc_y%beg, 1, &
                            q_cons_buff_recv(0), &
                            buff_size*sys_size*(m + 2*buff_size + 1), &
                            MPI_DOUBLE_PRECISION, bc_y%end, 1, &
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                        !call cpu_time(e_time)
                        !mpi_time = mpi_time + (e_time - s_time)

#if defined(_OPENACC) && defined(__PGI)
                    end if
#endif

                else                        ! PBC at the end only

                    ! Packing buffer to be sent to bc_y%end
!$acc parallel loop collapse(3) gang vector default(present) private(r)
                    do i = 1, sys_size
                            do k = n - buff_size + 1, n
                                do j = -buff_size, m + buff_size
                                    r = (i - 1) + sys_size* &
                                        ((j + buff_size) + (m + 2*buff_size + 1)* &
                                         ((k - n + buff_size - 1)))
                                    q_cons_buff_send(r) = q_cons_vf(i)%sf(j, k)
                                end do
                            end do
                    end do

                    !call MPI_Barrier(MPI_COMM_WORLD, ierr)

#if defined(_OPENACC) && defined(__PGI)
                    if (cu_mpi) then
!$acc host_data use_device( q_cons_buff_recv, q_cons_buff_send )

                        ! Send/receive buffer to/from bc_x%end/bc_x%beg
                        call MPI_SENDRECV( &
                            q_cons_buff_send(0), &
                            buff_size*sys_size*(m + 2*buff_size + 1), &
                            MPI_DOUBLE_PRECISION, bc_y%end, 0, &
                            q_cons_buff_recv(0), &
                            buff_size*sys_size*(m + 2*buff_size + 1), &
                            MPI_DOUBLE_PRECISION, bc_y%end, 1, &
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

!$acc end host_data
!$acc wait
                    else
#endif
                        !call cpu_time(s_time)
!$acc update host(q_cons_buff_send)
                        !call cpu_time(e_time)
                        !compress_time = compress_time + (e_time - s_time)

                        !call cpu_time(s_time)
                        ! Send/receive buffer to/from bc_x%end/bc_x%beg
                        call MPI_SENDRECV( &
                            q_cons_buff_send(0), &
                            buff_size*sys_size*(m + 2*buff_size + 1), &
                            MPI_DOUBLE_PRECISION, bc_y%end, 0, &
                            q_cons_buff_recv(0), &
                            buff_size*sys_size*(m + 2*buff_size + 1), &
                            MPI_DOUBLE_PRECISION, bc_y%end, 1, &
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                        !call cpu_time(e_time)
                        !mpi_time = mpi_time + (e_time - s_time)

#if defined(_OPENACC) && defined(__PGI)
                    end if
#endif

                end if

#if defined(_OPENACC) && defined(__PGI)
                if (cu_mpi .eqv. .false.) then
                    !call cpu_time(s_time)
!$acc update device(q_cons_buff_recv)
                    !call cpu_time(e_time)
                    !decompress_time = decompress_time + (e_time - s_time)
                end if
#endif

                ! Unpacking buffer received form bc_y%end
!$acc parallel loop collapse(3) gang vector default(present) private(r)
                do i = 1, sys_size
                    do k = n + 1, n + buff_size
                        do j = -buff_size, m + buff_size
                            r = (i - 1) + sys_size* &
                                ((j + buff_size) + (m + 2*buff_size + 1)* &
                                 ((k - n - 1)))
                            q_cons_vf(i)%sf(j, k) = q_cons_buff_recv(r)
                        end do
                    end do
                end do

            end if
            ! END: MPI Communication in y-direction ============================

        end if
        ! END: MPI Communication in z-direction ============================

#endif

    end subroutine s_mpi_sendrecv_conservative_variables_buffers ! ---------

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_mpi_proxy_module() ! -----------------------------

#ifdef MFC_MPI

        ! Deallocating q_cons_buff_send and q_cons_buff_recv
        @:DEALLOCATE(q_cons_buff_send, q_cons_buff_recv)

#endif

    end subroutine s_finalize_mpi_proxy_module ! ---------------------------

    !> The subroutine finalizes the MPI execution environment.
    subroutine s_mpi_finalize() ! ------------------------------------------

#ifdef MFC_MPI

        ! Finalizing the MPI environment
        call MPI_FINALIZE(ierr)

#endif

    end subroutine s_mpi_finalize ! ----------------------------------------

end module m_mpi_proxy
