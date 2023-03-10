!>
!! @file m_rhs.f90
!! @brief Contains module m_rhs

#:include 'macros.fpp'

!> @brief The module contains the subroutines used to calculate the right-
!!              hand-side (RHS) in the quasi-conservative, shock- and interface-
!!              capturing finite-volume framework for the multicomponent Navier-
!!              Stokes equations supplemented by appropriate advection equations
!!              used to capture the material interfaces. The system of equations
!!              is closed by the stiffened gas equation of state, as well as any
!!              required mixture relationships. Capillarity effects are included
!!              and are modeled by the means of a volume force acting across the
!!              diffuse material interface region. The implementation details of
!!              surface tension may be found in Perigaud and Saurel (2005). Note
!!              that both viscous and surface tension effects are only available
!!              in the volume fraction model.
module m_rhs

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_weno                 !< Weighted and essentially non-oscillatory (WENO)
                               !! schemes for spatial reconstruction of variables
    use m_riemann_solvers      !< Exact and approximate Riemann problem solvers

    use m_viscous
    
    use m_nvtx
    
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_rhs_module, &
 s_compute_rhs, &
 s_finalize_rhs_module


    type(vector_field) :: q_cons_qp 
    type(vector_field) :: q_prim_qp

    type(vector_field), allocatable, dimension(:) :: qL_cons_n
    type(vector_field), allocatable, dimension(:) :: qR_cons_n

    type(vector_field), allocatable, dimension(:) :: qL_prim_n
    type(vector_field), allocatable, dimension(:) :: qR_prim_n

    type(vector_field) :: dq_prim_dx_qp
    type(vector_field) :: dq_prim_dy_qp
    type(vector_field) :: gm_vel_qp

    !> @name The left and right WENO-reconstructed cell-boundary values of the cell-
    !! average first-order spatial derivatives of the primitive variables. The
    !! cell-average of the first-order spatial derivatives may be found in the
    !! variables dq_prim_ds_qp, where s = x, y or z.
    !> @{
    type(vector_field), allocatable, dimension(:) :: dqL_prim_dx_n
    type(vector_field), allocatable, dimension(:) :: dqL_prim_dy_n
    type(vector_field), allocatable, dimension(:) :: dqR_prim_dx_n
    type(vector_field), allocatable, dimension(:) :: dqR_prim_dy_n
    !> @}

    !> @name The cell-boundary values of the fluxes (src - source). 
    !! These are computed by applying the chosen Riemann problem solver
    !! on the left and right cell-boundary values of the primitive variables,
    !! qK_prim_n, the first-order spatial derivatives, dqK_prim_ds_n, as
    !! well as the curvature of volume fractions, kappaK_n.
    !> @{
    type(vector_field), allocatable, dimension(:) :: flux_n
    type(vector_field), allocatable, dimension(:) :: flux_src_n
    !> @}

    type(vector_field), allocatable, dimension(:) :: qL_prim, qR_prim

    type(int_bounds_info) :: iv !< Vector field indical bounds

    !> @name Indical bounds in the x-, y- and z-directions
    !> @{
    type(int_bounds_info) :: ix, iy
    !> @}

    type(int_bounds_info) :: is1, is2

    type(int_bounds_info) :: ixt, iyt

    type(vector_field), allocatable, dimension(:) :: myflux_vf, myflux_src_vf

    real(kind(0d0)), allocatable, dimension(:, :, :) :: qL_rsx_vf, qL_rsy_vf, qR_rsx_vf, qR_rsy_vf
    real(kind(0d0)), allocatable, dimension(:, :, :) :: dqL_rsx_vf, dqL_rsy_vf, dqR_rsx_vf, dqR_rsy_vf

   

    real(kind(0d0)), allocatable, dimension(:) :: gamma_min, pres_inf
!$acc declare create(gamma_min, pres_inf)

    real(kind(0d0)), allocatable, dimension(:, :) :: Res
!$acc declare create(Res)

!$acc declare create(q_cons_qp,q_prim_qp,qL_cons_n,qR_cons_n,qL_prim_n,qR_prim_n,  &
!$acc   dq_prim_dx_qp,dq_prim_dy_qp,gm_vel_qp,dqL_prim_dx_n,dqL_prim_dy_n, &
!$acc   dqR_prim_dx_n,dqR_prim_dy_n,       &
!$acc   flux_n,flux_src_n,       &
!$acc   qL_prim, qR_prim, iv,ix, iy,is1,is2, &
!$acc   myflux_vf, myflux_src_vf, &
!$acc   qL_rsx_vf, qL_rsy_vf, qR_rsx_vf, qR_rsy_vf, &
!$acc   dqL_rsx_vf, dqL_rsy_vf, dqR_rsx_vf, dqR_rsy_vf, &
!$acc   ixt, iyt)


contains

    !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_rhs_module() ! ---------------------------------

        integer :: i, j, l !< Generic loop iterators

        ! Configuring Coordinate Direction Indexes =========================
        ix%beg = -buff_size; iy%beg = 0

        if (n > 0) iy%beg = -buff_size

        ix%end = m - ix%beg; iy%end = n - iy%beg
        ! ==================================================================

        !$acc update device(ix, iy)
        
        ixt = ix; iyt = iy

        @:ALLOCATE(q_cons_qp%vf(1:sys_size))
        @:ALLOCATE(q_prim_qp%vf(1:sys_size))

        do l = 1, sys_size
            @:ALLOCATE(q_cons_qp%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end))
        end do

        do l = mom_idx%beg, E_idx
            @:ALLOCATE(q_prim_qp%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end))
        end do

        do l = adv_idx%end + 1, sys_size
            @:ALLOCATE(q_prim_qp%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end))
        end do

        do l = 1, cont_idx%end
            q_prim_qp%vf(l)%sf => &
                q_cons_qp%vf(l)%sf
!$acc enter data attach(q_prim_qp%vf(l)%sf)
        end do

        do l = adv_idx%beg, adv_idx%end
            q_prim_qp%vf(l)%sf => &
                q_cons_qp%vf(l)%sf
!$acc enter data attach(q_prim_qp%vf(l)%sf)
        end do

        ! ==================================================================


        ! Allocation/Association of qK_cons_n and qK_prim_n ==========
        @:ALLOCATE(qL_cons_n(1:num_dims))
        @:ALLOCATE(qR_cons_n(1:num_dims))
        @:ALLOCATE(qL_prim_n(1:num_dims))
        @:ALLOCATE(qR_prim_n(1:num_dims))

        @:ALLOCATE(qL_prim(1:num_dims))
        @:ALLOCATE(qR_prim(1:num_dims))

        do i = 1, num_dims
            @:ALLOCATE(qL_prim(i)%vf(1:sys_size))
            @:ALLOCATE(qR_prim(i)%vf(1:sys_size))
        end do

        if (weno_Re_flux) then

            do i = 1, num_dims
                do l = mom_idx%beg, mom_idx%end
                    @:ALLOCATE(qL_prim(i)%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end))
                    @:ALLOCATE(qR_prim(i)%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end))
                end do
            end do
        end if
        @:ALLOCATE(myflux_vf(1:num_dims))
        @:ALLOCATE(myflux_src_vf(1:num_dims))


        do i = 1, num_dims
            @:ALLOCATE(qL_cons_n(i)%vf(1:sys_size))
            @:ALLOCATE(qR_cons_n(i)%vf(1:sys_size))
            @:ALLOCATE(qL_prim_n(i)%vf(1:sys_size))
            @:ALLOCATE(qR_prim_n(i)%vf(1:sys_size))
            @:ALLOCATE(myflux_vf(i)%vf(1:sys_size))
            @:ALLOCATE(myflux_src_vf(i)%vf(1:sys_size))
        end do
        ! END: Allocation/Association of qK_cons_n and qK_prim_n =====

        @:ALLOCATE(qL_rsx_vf(ix%beg:ix%end, &
                                 iy%beg:iy%end, 1:sys_size))
        @:ALLOCATE(qR_rsx_vf(ix%beg:ix%end, &
                                 iy%beg:iy%end, 1:sys_size))

        if (n > 0) then

            @:ALLOCATE(qL_rsy_vf(iy%beg:iy%end, &
                                     ix%beg:ix%end, 1:sys_size))
            @:ALLOCATE(qR_rsy_vf(iy%beg:iy%end, &
                                     ix%beg:ix%end, 1:sys_size))
        else
            @:ALLOCATE(qL_rsy_vf(ix%beg:ix%end, &
                                     iy%beg:iy%end, 1:sys_size))
            @:ALLOCATE(qR_rsy_vf(ix%beg:ix%end, &
                                     iy%beg:iy%end, 1:sys_size))
        end if


        ! Allocation of dq_prim_ds_qp ======================================

        if (any(Re_size > 0)) then

            @:ALLOCATE(dq_prim_dx_qp%vf(1:sys_size))
            @:ALLOCATE(dq_prim_dy_qp%vf(1:sys_size))
            @:ALLOCATE(gm_vel_qp%vf(1:sys_size))
            
            if (any(Re_size > 0)) then

                do l = mom_idx%beg, mom_idx%end
                    @:ALLOCATE(dq_prim_dx_qp%vf(l)%sf( &
                              & ix%beg:ix%end, &
                              & iy%beg:iy%end ))
                    @:ALLOCATE(gm_vel_qp%vf(l)%sf( &
                              & ix%beg:ix%end, &
                              & iy%beg:iy%end ))
                end do

                if (n > 0) then

                    do l = mom_idx%beg, mom_idx%end
                        @:ALLOCATE(dq_prim_dy_qp%vf(l)%sf( &
                                 & ix%beg:ix%end, &
                                 & iy%beg:iy%end ))
                    end do

                end if

            end if

        end if
        ! END: Allocation of dq_prim_ds_qp =================================

        ! Allocation/Association of dqK_prim_ds_n =======================
        @:ALLOCATE(dqL_prim_dx_n(1:num_dims))
        @:ALLOCATE(dqL_prim_dy_n(1:num_dims))
        @:ALLOCATE(dqR_prim_dx_n(1:num_dims))
        @:ALLOCATE(dqR_prim_dy_n(1:num_dims))

        if (any(Re_size > 0)) then
            do i = 1, num_dims
                @:ALLOCATE(dqL_prim_dx_n(i)%vf(1:sys_size))
                @:ALLOCATE(dqL_prim_dy_n(i)%vf(1:sys_size))
                @:ALLOCATE(dqR_prim_dx_n(i)%vf(1:sys_size))
                @:ALLOCATE(dqR_prim_dy_n(i)%vf(1:sys_size))
                
                if (any(Re_size > 0)) then

                    do l = mom_idx%beg, mom_idx%end
                        @:ALLOCATE(dqL_prim_dx_n(i)%vf(l)%sf( &
                                 & ix%beg:ix%end, &
                                 & iy%beg:iy%end ))
                        @:ALLOCATE(dqR_prim_dx_n(i)%vf(l)%sf( &
                                 & ix%beg:ix%end, &
                                 & iy%beg:iy%end ))
                    end do

                    if (n > 0) then
                        do l = mom_idx%beg, mom_idx%end
                            @:ALLOCATE(dqL_prim_dy_n(i)%vf(l)%sf( &
                                     & ix%beg:ix%end, &
                                     & iy%beg:iy%end ))
                            @:ALLOCATE(dqR_prim_dy_n(i)%vf(l)%sf( &
                                     & ix%beg:ix%end, &
                                     & iy%beg:iy%end ))
                        end do
                    end if

                end if

            end do
        end if
        ! END: Allocation/Association of d K_prim_ds_n ==================

        if (any(Re_size > 0)) then
            if (weno_Re_flux) then
                @:ALLOCATE(dqL_rsx_vf(ix%beg:ix%end, &
                                          iy%beg:iy%end, mom_idx%beg:mom_idx%end))
                @:ALLOCATE(dqR_rsx_vf(ix%beg:ix%end, &
                                          iy%beg:iy%end, mom_idx%beg:mom_idx%end))

                if (n > 0) then

                    @:ALLOCATE(dqL_rsy_vf(iy%beg:iy%end, &
                                              ix%beg:ix%end, mom_idx%beg:mom_idx%end))
                    @:ALLOCATE(dqR_rsy_vf(iy%beg:iy%end, &
                                              ix%beg:ix%end, mom_idx%beg:mom_idx%end))
                else
                    @:ALLOCATE(dqL_rsy_vf(ix%beg:ix%end, &
                                              iy%beg:iy%end, mom_idx%beg:mom_idx%end))
                    @:ALLOCATE(dqR_rsy_vf(ix%beg:ix%end, &
                                              iy%beg:iy%end, mom_idx%beg:mom_idx%end))

                end if

            end if
        end if

        ! Allocation/Association of flux_n, flux_src_n===
        @:ALLOCATE(flux_n(1:num_dims))
        @:ALLOCATE(flux_src_n(1:num_dims))

        do i = 1, num_dims

            @:ALLOCATE(flux_n(i)%vf(1:sys_size))
            @:ALLOCATE(flux_src_n(i)%vf(1:sys_size))

            if (i == 1) then

                do l = 1, sys_size
                    @:ALLOCATE(flux_n(i)%vf(l)%sf( &
                             & ix%beg:ix%end, &
                             & iy%beg:iy%end ))
                end do

                if (any(Re_size > 0)) then
                    do l = mom_idx%beg, E_idx
                        @:ALLOCATE(flux_src_n(i)%vf(l)%sf( &
                                 & ix%beg:ix%end, &
                                 & iy%beg:iy%end ))
                    end do
                end if

                @:ALLOCATE(flux_src_n(i)%vf(adv_idx%beg)%sf( &
                         & ix%beg:ix%end, &
                         & iy%beg:iy%end ))

                do l = adv_idx%beg + 1, adv_idx%end
                    flux_src_n(i)%vf(l)%sf => &
                        flux_src_n(i)%vf(adv_idx%beg)%sf
                    !$acc enter data attach(flux_src_n(i)%vf(l)%sf(ix%beg:ix%end,iy%beg:iy%end))
                end do

            else
                do l = 1, sys_size
                    flux_n(i)%vf(l)%sf => &
                        flux_n(1)%vf(l)%sf
                    flux_src_n(i)%vf(l)%sf => &
                        flux_src_n(1)%vf(l)%sf

                    !$acc enter data attach(flux_n(i)%vf(l)%sf,flux_src_n(i)%vf(l)%sf)
                end do

            end if
        end do

        ! END: Allocation/Association of flux_n, flux_src_n ===
        @:ALLOCATE(gamma_min(1:num_fluids), pres_inf(1:num_fluids))

        do i = 1, num_fluids
            gamma_min(i) = 1d0/fluid_pp(i)%gamma + 1d0
            pres_inf(i) = fluid_pp(i)%pi_inf/(1d0 + fluid_pp(i)%gamma)
        end do
!$acc update device(gamma_min, pres_inf)

        if (any(Re_size > 0)) then
            @:ALLOCATE(Res(1:2, 1:maxval(Re_size)))
        end if

        if (any(Re_size > 0)) then
            do i = 1, 2
                do j = 1, Re_size(i)
                    Res(i, j) = fluid_pp(Re_idx(i, j))%Re(i)
                end do
            end do
!$acc update device(Res, Re_idx, Re_size)
        end if


        ! Associating the procedural pointer to the appropriate subroutine
        ! that will be utilized in the conversion to the mixture variables
        s_convert_to_mixture_variables => s_convert_species_to_mixture_variables

    end subroutine s_initialize_rhs_module ! -------------------------------

    ! [SHB]: This is a 'pruned' version of s_compute_rhs
    !   (see compute_rhs_full below for full version)
    !   it exercises all the key things, but gets rid of some of the extraneous
    !   calls that might hold back progress
    subroutine s_compute_rhs(q_cons_vf, q_prim_vf, rhs_vf, t_step) ! -------

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(INOUT) :: rhs_vf
        integer, intent(IN) :: t_step

        integer :: i, j, k, l, q, id !< Generic loop iterators

        ! Configuring Coordinate Direction Indexes =========================
        ix%beg = -buff_size; iy%beg = 0

        if (n > 0) iy%beg = -buff_size

        ix%end = m - ix%beg; iy%end = n - iy%beg
        ! ==================================================================

        !$acc update device(ix, iy)

        ! Association/Population of Working Variables ======================
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, sys_size
                do k = iy%beg, iy%end
                    do j = ix%beg, ix%end
                        q_cons_qp%vf(i)%sf(j, k) = q_cons_vf(i)%sf(j, k)
                    end do
                end do
        end do

        call nvtxStartRange("RHS-MPI")
        call s_populate_conservative_variables_buffers()
        call nvtxEndRange

        
        ! ==================================================================

        ! Converting Conservative to Primitive Variables ==================


        call nvtxStartRange("RHS-CONVERT")
        call s_convert_conservative_to_primitive_variables( &
            q_cons_qp%vf, &
            q_prim_qp%vf )
        call nvtxEndRange


        
        if (t_step == t_step_stop) return
        ! ==================================================================

        call nvtxStartRange("Viscous")
        if (any(Re_size > 0)) call s_get_viscous(qL_rsx_vf, qL_rsy_vf, &
                                            dqL_prim_dx_n, dqL_prim_dy_n, &
                                            qL_prim, &
                                            qR_rsx_vf, qR_rsy_vf, &
                                            dqR_prim_dx_n, dqR_prim_dy_n, &
                                            qR_prim, &
                                            q_prim_qp, &
                                            dq_prim_dx_qp, dq_prim_dy_qp, &
                                            ix, iy)
        call nvtxEndRange()
        
        ! Dimensional Splitting Loop =======================================

        do id = 1, num_dims

            ! Configuring Coordinate Direction Indexes ======================
            ix%beg = -buff_size; iy%beg = 0

            if (n > 0) iy%beg = -buff_size

            ix%end = m - ix%beg; iy%end = n - iy%beg
            ! ===============================================================
            ! Reconstructing Primitive/Conservative Variables ===============
            
            if (all(Re_size == 0)) then
                    iv%beg = 1; iv%end = sys_size
                !call nvtxStartRange("RHS-WENO")
                call nvtxStartRange("RHS-WENO")
                call s_reconstruct_cell_boundary_values( &
                    q_prim_qp%vf(1:sys_size), &
                    qL_rsx_vf, qL_rsy_vf, &
                    qR_rsx_vf, qR_rsy_vf, &
                    id)
                call nvtxEndRange
            else
                call nvtxStartRange("RHS-WENO")
                iv%beg = 1; iv%end = contxe
                call s_reconstruct_cell_boundary_values( &
                    q_prim_qp%vf(iv%beg:iv%end), &
                    qL_rsx_vf, qL_rsy_vf, &
                    qR_rsx_vf, qR_rsy_vf, &
                    id)

                iv%beg = E_idx; iv%end = E_idx
                call s_reconstruct_cell_boundary_values( &
                    q_prim_qp%vf(iv%beg:iv%end), &
                    qL_rsx_vf, qL_rsy_vf, &
                    qR_rsx_vf, qR_rsy_vf, &
                    id)

                iv%beg = advxb; iv%end = advxe
                call s_reconstruct_cell_boundary_values( &
                    q_prim_qp%vf(iv%beg:iv%end), &
                    qL_rsx_vf, qL_rsy_vf, &
                    qR_rsx_vf, qR_rsy_vf, &
                    id)

                iv%beg = mom_idx%beg; iv%end = mom_idx%end
                if (weno_Re_flux) then
                    call s_reconstruct_cell_boundary_values_visc_deriv( &
                        dq_prim_dx_qp%vf(iv%beg:iv%end), &
                        dqL_rsx_vf, dqL_rsy_vf, &
                        dqR_rsx_vf, dqR_rsy_vf, &
                        id, dqL_prim_dx_n(id)%vf(iv%beg:iv%end), dqR_prim_dx_n(id)%vf(iv%beg:iv%end), &
                        ix, iy)
                    if (n > 0) then
                        call s_reconstruct_cell_boundary_values_visc_deriv( &
                            dq_prim_dy_qp%vf(iv%beg:iv%end), &
                            dqL_rsx_vf, dqL_rsy_vf, &
                            dqR_rsx_vf, dqR_rsy_vf, &
                            id, dqL_prim_dy_n(id)%vf(iv%beg:iv%end), dqR_prim_dy_n(id)%vf(iv%beg:iv%end), &
                            ix, iy)
                    end if
                end if
                call nvtxEndRange
            end if

            ! Configuring Coordinate Direction Indexes ======================
            if (id == 1) then
                ix%beg = -1; iy%beg = 0
            elseif (id == 2) then
                ix%beg = 0; iy%beg = -1
            end if
            ix%end = m; iy%end = n
            ! ===============================================================
            call nvtxStartRange("RHS-Riemann")

            ! Computing Riemann Solver Flux and Source Flux =================

            call s_hllc_riemann_solver( &
                                  qR_rsx_vf, &
                                  qR_rsy_vf, &
                                  dqR_prim_dx_n(id)%vf, &
                                  dqR_prim_dy_n(id)%vf, &
                                  qL_rsx_vf, &
                                  qL_rsy_vf, &
                                  dqL_prim_dx_n(id)%vf, &
                                  dqL_prim_dy_n(id)%vf, &
                                  flux_n(id)%vf, &
                                  flux_src_n(id)%vf, &
                                  id, ix, iy)

            call nvtxEndRange

            ! ===============================================================



            call nvtxStartRange("RHS_Flux_Add")
            if (id == 1) then

                !$acc parallel loop collapse(3) gang vector default(present)
                do j = 1, sys_size
                        do l = 0, n
                            do k = 0, m
                                rhs_vf(j)%sf(k, l) = 1d0/dx(k)* &
                                                        (flux_n(1)%vf(j)%sf(k - 1, l) &
                                                         - flux_n(1)%vf(j)%sf(k, l))
                            end do
                        end do
                end do

                !$acc parallel loop collapse(3) gang vector default(present)
                do j = advxb, advxe
                        do l = 0, n
                            do k = 0, m
                                rhs_vf(j)%sf(k, l) = &
                                    rhs_vf(j)%sf(k, l) + 1d0/dx(k)* &
                                    q_cons_qp%vf(j)%sf(k, l)* &
                                    (flux_src_n(1)%vf(j)%sf(k, l) &
                                     - flux_src_n(1)%vf(j)%sf(k - 1, l))
                            end do
                        end do
                end do

                if (any(Re_size > 0)) then
!$acc parallel loop collapse(2) gang vector default(present)
                        do k = 0, n
                            do j = 0, m
                                !$acc loop seq
                                do i = momxb, E_idx
                                    rhs_vf(i)%sf(j, k) = &
                                        rhs_vf(i)%sf(j, k) + 1d0/dx(j)* &
                                        (flux_src_n(1)%vf(i)%sf(j - 1, k) &
                                         - flux_src_n(1)%vf(i)%sf(j, k))
                                end do
                            end do
                    end do
                end if

            elseif (id == 2) then
                ! RHS Contribution in y-direction ===============================
                ! Applying the Riemann fluxes

                !$acc parallel loop collapse(3) gang vector default(present)
                do j = 1, sys_size
                        do k = 0, n
                            do q = 0, m
                                rhs_vf(j)%sf(q, k) = &
                                    rhs_vf(j)%sf(q, k) + 1d0/dy(k)* &
                                    (flux_n(2)%vf(j)%sf(q, k - 1) &
                                     - flux_n(2)%vf(j)%sf(q, k))
                            end do
                        end do
                end do
                ! Applying source terms to the RHS of the advection equations


                !$acc parallel loop collapse(3) gang vector default(present)
                do j = advxb, advxe
                    do k = 0, n
                        do q = 0, m
                            rhs_vf(j)%sf(q, k) = &
                                rhs_vf(j)%sf(q, k) + 1d0/dy(k)* &
                                q_cons_qp%vf(j)%sf(q, k)* &
                                (flux_src_n(2)%vf(j)%sf(q, k) &
                                 - flux_src_n(2)%vf(j)%sf(q, k - 1))
                        end do
                    end do
                end do



                if (any(Re_size > 0)) then
                    !$acc parallel loop collapse(2) gang vector default(present)
                    do k = 0, n
                        do j = 0, m
                            !$acc loop seq
                            do i = momxb, E_idx
                                rhs_vf(i)%sf(j, k) = &
                                    rhs_vf(i)%sf(j, k) + 1d0/dy(k)* &
                                    (flux_src_n(2)%vf(i)%sf(j, k - 1) &
                                     - flux_src_n(2)%vf(i)%sf(j, k))
                            end do
                        end do
                    end do
                end if
            end if  ! id loop
            call nvtxEndRange

        end do
        ! END: Dimensional Splitting Loop =================================

        if (run_time_info .or. probe_wrt) then

            ix%beg = -buff_size; iy%beg = 0
            if (n > 0) iy%beg = -buff_size; 
            ix%end = m - ix%beg; iy%end = n - iy%beg
            !$acc update device(ix, iy)

            !$acc parallel loop collapse(3) gang vector default(present)
            do i = 1, sys_size
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end
                            q_prim_vf(i)%sf(j, k) = q_prim_qp%vf(i)%sf(j, k)
                        end do
                    end do
            end do

        end if

        ! ==================================================================

    end subroutine s_compute_rhs ! -----------------------------------------



    !>  The purpose of this procedure is to populate the buffers
        !!      of the conservative variables, depending on the selected
        !!      boundary conditions.
    subroutine s_populate_conservative_variables_buffers() ! ---------------

        integer :: i, j, l, k !< Generic loop iterators

        ! Population of Buffers in x-direction =============================

        if (bc_x%beg <= -3) then         ! Ghost-cell extrap. BC at beginning

!$acc parallel loop collapse(3) gang vector default(present)
            do i = 1, sys_size
                    do k = 0, n
                        do j = 1, buff_size
                            q_cons_qp%vf(i)%sf(-j, k) = &
                                q_cons_qp%vf(i)%sf(0, k)
                        end do
                    end do
            end do

        elseif (bc_x%beg == -2) then     ! Symmetry BC at beginning

!$acc parallel loop collapse(2) gang vector default(present)
                do k = 0, n
                    do j = 1, buff_size
!$acc loop seq
                        do i = 1, contxe
                            q_cons_qp%vf(i)%sf(-j, k) = &
                                q_cons_qp%vf(i)%sf(j - 1, k)
                        end do

                        q_cons_qp%vf(momxb)%sf(-j, k) = &
                            -q_cons_qp%vf(momxb)%sf(j - 1, k)
!$acc loop seq
                        do i = momxb + 1, sys_size
                            q_cons_qp%vf(i)%sf(-j, k) = &
                                q_cons_qp%vf(i)%sf(j - 1, k)
                        end do
                    end do
            end do

        elseif (bc_x%beg == -1) then     ! Periodic BC at beginning

!$acc parallel loop collapse(3) gang vector default(present)
            do i = 1, sys_size
                    do k = 0, n
                        do j = 1, buff_size
                            q_cons_qp%vf(i)%sf(-j, k) = &
                                q_cons_qp%vf(i)%sf(m - (j - 1), k)
                        end do
                    end do
            end do

        else                            ! Processor BC at beginning

            call s_mpi_sendrecv_conservative_variables_buffers( &
                q_cons_qp%vf, 1, -1)

        end if

        if (bc_x%end <= -3) then         ! Ghost-cell extrap. BC at end

!$acc parallel loop collapse(3) gang vector default(present)
            do i = 1, sys_size
                    do k = 0, n
                        do j = 1, buff_size
                            q_cons_qp%vf(i)%sf(m + j, k) = &
                                q_cons_qp%vf(i)%sf(m, k)
                        end do
                    end do
            end do

        elseif (bc_x%end == -2) then     ! Symmetry BC at end

!$acc parallel loop collapse(2) default(present)
                do k = 0, n
                    do j = 1, buff_size

!$acc loop seq
                        do i = 1, contxe
                            q_cons_qp%vf(i)%sf(m + j, k) = &
                                q_cons_qp%vf(i)%sf(m - (j - 1), k)
                        end do

                        q_cons_qp%vf(momxb)%sf(m + j, k) = &
                            -q_cons_qp%vf(momxb)%sf(m - (j - 1), k)

!$acc loop seq
                        do i = momxb + 1, sys_size
                            q_cons_qp%vf(i)%sf(m + j, k) = &
                                q_cons_qp%vf(i)%sf(m - (j - 1), k)
                        end do

                    end do
                end do

        elseif (bc_x%end == -1) then     ! Periodic BC at end

!$acc parallel loop collapse(3) gang vector default(present)
            do i = 1, sys_size
                    do k = 0, n
                        do j = 1, buff_size
                            q_cons_qp%vf(i)%sf(m + j, k) = &
                                q_cons_qp%vf(i)%sf(j - 1, k)
                        end do
                    end do
            end do

        else                            ! Processor BC at end

            call s_mpi_sendrecv_conservative_variables_buffers( &
                q_cons_qp%vf, 1, 1)

        end if

        ! END: Population of Buffers in x-direction ========================

        ! Population of Buffers in y-direction =============================

        if (n == 0) then

            return

        elseif (bc_y%beg <= -3 .and. bc_y%beg /= -13) then     ! Ghost-cell extrap. BC at beginning

!$acc parallel loop collapse(3) gang vector default(present)
            do i = 1, sys_size
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            q_cons_qp%vf(i)%sf(l, -j) = &
                                q_cons_qp%vf(i)%sf(l, 0)
                        end do
                    end do
            end do

        elseif (bc_y%beg == -2) then     ! Symmetry BC at beginning
!$acc parallel loop collapse(2) gang vector default(present)
                do j = 1, buff_size
                    do l = -buff_size, m + buff_size
!$acc loop seq
                        do i = 1, momxb
                            q_cons_qp%vf(i)%sf(l, -j) = &
                                q_cons_qp%vf(i)%sf(l, j - 1)
                        end do

                        q_cons_qp%vf(momxb + 1)%sf(l, -j) = &
                            -q_cons_qp%vf(momxb + 1)%sf(l, j - 1)
!$acc loop seq
                        do i = momxb + 2, sys_size
                            q_cons_qp%vf(i)%sf(l, -j) = &
                                q_cons_qp%vf(i)%sf(l, j - 1)
                        end do
                    end do
                end do

        elseif (bc_y%beg == -1) then     ! Periodic BC at beginning
!$acc parallel loop collapse(3) gang vector default(present)
            do i = 1, sys_size
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            q_cons_qp%vf(i)%sf(l, -j) = &
                                q_cons_qp%vf(i)%sf(l, n - (j - 1))
                        end do
                    end do
            end do

        else                            ! Processor BC at beginning

            call s_mpi_sendrecv_conservative_variables_buffers( &
                q_cons_qp%vf, 2, -1)

        end if

        if (bc_y%end <= -3) then         ! Ghost-cell extrap. BC at end
!$acc parallel loop collapse(3) gang vector default(present)
            do i = 1, sys_size
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            q_cons_qp%vf(i)%sf(l, n + j) = &
                                q_cons_qp%vf(i)%sf(l, n)
                        end do
                    end do
            end do

        elseif (bc_y%end == -2) then     ! Symmetry BC at end

!$acc parallel loop collapse(2) gang vector default(present)
                do j = 1, buff_size
                    do l = -buff_size, m + buff_size
!$acc loop seq
                        do i = 1, momxb
                            q_cons_qp%vf(i)%sf(l, n + j) = &
                                q_cons_qp%vf(i)%sf(l, n - (j - 1))
                        end do

                        q_cons_qp%vf(momxb + 1)%sf(l, n + j) = &
                            -q_cons_qp%vf(momxb + 1)%sf(l, n - (j - 1))
!$acc loop seq
                        do i = momxb + 2, sys_size
                            q_cons_qp%vf(i)%sf(l, n + j) = &
                                q_cons_qp%vf(i)%sf(l, n - (j - 1))
                        end do
                    end do
                end do

        elseif (bc_y%end == -1) then     ! Periodic BC at end
!$acc parallel loop collapse(3) gang vector default(present)
            do i = 1, sys_size
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            q_cons_qp%vf(i)%sf(l, n + j) = &
                                q_cons_qp%vf(i)%sf(l, j - 1)
                        end do
                    end do
            end do

        else                            ! Processor BC at end

            call s_mpi_sendrecv_conservative_variables_buffers( &
                q_cons_qp%vf, 2, 1)

        end if
        ! END: Population of Buffers in y-direction ========================

    end subroutine s_populate_conservative_variables_buffers ! -------------

    !>  The purpose of this subroutine is to WENO-reconstruct the
        !!      left and the right cell-boundary values, including values
        !!      at the Gaussian quadrature points, from the cell-averaged
        !!      variables.
        !!  @param v_vf Cell-average variables
        !!  @param vL_qp Left WENO-reconstructed, cell-boundary values including
        !!          the values at the quadrature points, of the cell-average variables
        !!  @param vR_qp Right WENO-reconstructed, cell-boundary values including
        !!          the values at the quadrature points, of the cell-average variables
        !!  @param norm_dir Splitting coordinate direction
    subroutine s_reconstruct_cell_boundary_values(v_vf, vL_x, vL_y, vR_x, vR_y, & ! -
                                                      norm_dir)

        type(scalar_field), dimension(iv%beg:iv%end), intent(IN) :: v_vf
        real(kind(0d0)), dimension(startx:, starty:, 1:), intent(INOUT) :: vL_x, vL_y, vR_x, vR_y
        integer, intent(IN) :: norm_dir
        integer :: weno_dir !< Coordinate direction of the WENO reconstruction

        if (norm_dir == 1) then
            is1 = ix; is2 = iy
            weno_dir = 1; is1%beg = is1%beg + weno_polyn
            is1%end = is1%end - weno_polyn

        elseif (norm_dir == 2) then
            is1 = iy; is2 = ix
            weno_dir = 2; is1%beg = is1%beg + weno_polyn
            is1%end = is1%end - weno_polyn
        end if

        if (n > 0) then
            call s_weno(v_vf(iv%beg:iv%end), &
                vL_x(:, :, iv%beg:iv%end), vL_y(:, :, iv%beg:iv%end), vR_x(:, :, iv%beg:iv%end), vR_y(:, :, iv%beg:iv%end),  &
                weno_dir, &
                is1, is2)
        else
            call s_weno(v_vf(iv%beg:iv%end), &
                vL_x(:, :, iv%beg:iv%end), vL_y(:, :, :), vR_x(:, :, iv%beg:iv%end), vR_y(:, :, :), &
                weno_dir, &
                is1, is2)
        end if

        ! ==================================================================
    end subroutine s_reconstruct_cell_boundary_values ! --------------------

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_rhs_module() ! -----------------------------------

        integer :: i, j, l !< Generic loop iterators

        do j = cont_idx%beg, cont_idx%end
!$acc exit data detach(q_prim_qp%vf(j)%sf)
            nullify (q_prim_qp%vf(j)%sf)
        end do

        do j = adv_idx%beg, adv_idx%end
!$acc exit data detach(q_prim_qp%vf(j)%sf)
            nullify (q_prim_qp%vf(j)%sf)
        end do

        do j = mom_idx%beg, E_idx
            @:DEALLOCATE(q_cons_qp%vf(j)%sf)
            @:DEALLOCATE(q_prim_qp%vf(j)%sf)
        end do

        @:DEALLOCATE(q_cons_qp%vf, q_prim_qp%vf)

        @:DEALLOCATE(qL_rsx_vf, qR_rsx_vf)

        if (n > 0) then
            @:DEALLOCATE(qL_rsy_vf, qR_rsy_vf)
        end if

        if (weno_Re_flux) then
            @:DEALLOCATE(dqL_rsx_vf, dqR_rsx_vf)

            if (n > 0) then
                @:DEALLOCATE(dqL_rsy_vf, dqR_rsy_vf)
            end if

        end if

        do i = num_dims, 1, -1
            @:DEALLOCATE(qL_cons_n(i)%vf, qL_prim_n(i)%vf)
            @:DEALLOCATE(qR_cons_n(i)%vf, qR_prim_n(i)%vf)
        end do


        @:DEALLOCATE(qL_cons_n, qR_cons_n, qL_prim_n, qR_prim_n)

        if (any(Re_size > 0)) then
            do l = mom_idx%beg, mom_idx%end
                @:DEALLOCATE(dq_prim_dx_qp%vf(l)%sf)
                @:DEALLOCATE(gm_vel_qp%vf(l)%sf)
            end do

            if (n > 0) then

                do l = mom_idx%beg, mom_idx%end
                    @:DEALLOCATE(dq_prim_dy_qp%vf(l)%sf)
                end do

            end if

            @:DEALLOCATE(dq_prim_dx_qp%vf)
            @:DEALLOCATE(dq_prim_dy_qp%vf)
            @:DEALLOCATE(gm_vel_qp%vf)
        end if

        if (any(Re_size > 0)) then
            do i = num_dims, 1, -1
                if (any(Re_size > 0)) then

                    do l = mom_idx%beg, mom_idx%end
                        @:DEALLOCATE(dqL_prim_dx_n(i)%vf(l)%sf)
                        @:DEALLOCATE(dqR_prim_dx_n(i)%vf(l)%sf)
                    end do

                    if (n > 0) then
                        do l = mom_idx%beg, mom_idx%end
                            @:DEALLOCATE(dqL_prim_dy_n(i)%vf(l)%sf)
                            @:DEALLOCATE(dqR_prim_dy_n(i)%vf(l)%sf)
                        end do
                    end if

                end if

                @:DEALLOCATE(dqL_prim_dx_n(i)%vf)
                @:DEALLOCATE(dqL_prim_dy_n(i)%vf)
                @:DEALLOCATE(dqR_prim_dx_n(i)%vf)
                @:DEALLOCATE(dqR_prim_dy_n(i)%vf)
            end do
        end if

        @:DEALLOCATE(dqL_prim_dx_n, dqL_prim_dy_n)
        @:DEALLOCATE(dqR_prim_dx_n, dqR_prim_dy_n)


        do i = num_dims, 1, -1
            if (i /= 1) then
                do l = 1, sys_size
                    nullify (flux_n(i)%vf(l)%sf)
                    nullify (flux_src_n(i)%vf(l)%sf)
                end do
            else
                do l = 1, sys_size
                    @:DEALLOCATE(flux_n(i)%vf(l)%sf)
                end do

                if (any(Re_size > 0)) then
                    do l = mom_idx%beg, E_idx
                        @:DEALLOCATE(flux_src_n(i)%vf(l)%sf)
                    end do
                end if

                do l = adv_idx%beg + 1, adv_idx%end
                    nullify (flux_src_n(i)%vf(l)%sf)
                end do

                @:DEALLOCATE(flux_src_n(i)%vf(adv_idx%beg)%sf)
            end if

            @:DEALLOCATE(flux_n(i)%vf, flux_src_n(i)%vf)
        end do

        @:DEALLOCATE(flux_n, flux_src_n)

        s_convert_to_mixture_variables => null()

    end subroutine s_finalize_rhs_module ! ---------------------------------

end module m_rhs
