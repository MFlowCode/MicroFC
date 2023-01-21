!>
!! @file m_riemann_solvers.f90
!! @brief Contains module m_riemann_solvers

!> @brief This module features a database of approximate and exact Riemann
!!              problem solvers for the Navier-Stokes system of equations, which
!!              is supplemented by appropriate advection equations that are used
!!              to capture the material interfaces. The closure of the system is
!!              achieved by the stiffened gas equation of state and any required
!!              mixture relations. Surface tension effects are accounted for and
!!              are modeled by means of a volume force acting across the diffuse
!!              material interface region. The implementation details of viscous
!!              and capillary effects, into the Riemann solvers, may be found in
!!              Perigaud and Saurel (2005). Note that both effects are available
!!              only in the volume fraction model. At this time, the approximate
!!              and exact Riemann solvers that are listed below are available:
!!                  1) Harten-Lax-van Leer (HLL)
!!                  2) Harten-Lax-van Leer-Contact (HLLC)
!!                  3) Exact
#:include 'inline_riemann.fpp'

module m_riemann_solvers

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_riemann_solvers_module, &
                        s_riemann_solver, &
                        s_hllc_riemann_solver, &
                        s_finalize_riemann_solvers_module

    abstract interface ! =======================================================

        !> Abstract interface to the subroutines that are utilized to compute the
        !! Riemann problem solution. For additional information please reference:
        !!                        1) s_hll_riemann_solver
        !!                        2) s_hllc_riemann_solver
        !!                        3) s_exact_riemann_solver
        !!  @param qL_prim_vf The  left WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables
        !!  @param qR_prim_vf The right WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables
        !!  @param dqL_prim_dx_vf The  left WENO-reconstructed cell-boundary values of the
        !!      first-order x-dir spatial derivatives
        !!  @param dqL_prim_dy_vf The  left WENO-reconstructed cell-boundary values of the
        !!      first-order y-dir spatial derivatives
        !!  @param dqR_prim_dx_vf The right WENO-reconstructed cell-boundary values of the
        !!      first-order x-dir spatial derivatives
        !!  @param dqR_prim_dy_vf The right WENO-reconstructed cell-boundary values of the
        !!      first-order y-dir spatial derivatives
        !!  @param gm_alphaL_vf  Left averaged gradient magnitude
        !!  @param gm_alphaR_vf Right averaged gradient magnitude
        !!  @param flux_vf Intra-cell fluxes
        !!  @param flux_src_vf Intra-cell fluxes sources
        !!  @param norm_dir Dir. splitting direction
        !!  @param ix Index bounds in the x-dir
        !!  @param iy Index bounds in the y-dir
        !!  @param q_prim_vf Cell-averaged primitive variables
        subroutine s_abstract_riemann_solver(qL_prim_rsx_vf, qL_prim_rsy_vf, dqL_prim_dx_vf, &
                                             dqL_prim_dy_vf, &
                                             qL_prim_vf, &
                                             qR_prim_rsx_vf, qR_prim_rsy_vf, dqR_prim_dx_vf, &
                                             dqR_prim_dy_vf, &
                                             qR_prim_vf, &
                                             q_prim_vf, &
                                             flux_vf, flux_src_vf, &
                                             norm_dir, ix, iy)

            import :: scalar_field, int_bounds_info, sys_size, startx, starty

            real(kind(0d0)), dimension(startx:, starty:, 1:), intent(INOUT) :: qL_prim_rsx_vf, qL_prim_rsy_vf, qR_prim_rsx_vf, qR_prim_rsy_vf
            type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf

            type(scalar_field), allocatable, dimension(:), intent(INOUT) :: qL_prim_vf, qR_prim_vf

            type(scalar_field), &
                allocatable, dimension(:), &
                intent(INOUT) :: dqL_prim_dx_vf, dqR_prim_dx_vf, &
                                 dqL_prim_dy_vf, dqR_prim_dy_vf

            type(scalar_field), &
                dimension(sys_size), &
                intent(INOUT) :: flux_vf, flux_src_vf

            integer, intent(IN) :: norm_dir

            type(int_bounds_info), intent(IN) :: ix, iy

        end subroutine s_abstract_riemann_solver

        !> The abstract interface to the subroutines that are utilized to compute
        !! the wave speeds of the Riemann problem either directly or by the means
        !! of pressure-velocity estimates. For more information please refer to:
        !!      1) s_compute_direct_wave_speeds
        !!      2) s_compute_pressure_velocity_wave_speeds
        !!  @param i First coordinate location index
        !!  @param j Second coordinate location index
        !!  @param k Third coordinate location index
        subroutine s_compute_abstract_wave_speeds(i, j, k)

            integer, intent(IN) :: i, j, k

        end subroutine s_compute_abstract_wave_speeds

        !> The abstract interface to the subroutines that are utilized to compute
        !! the viscous source fluxes for either Cartesian or cylindrical geometries.
        !! For more information please refer to:
        !!      1) s_compute_cartesian_viscous_source_flux
        !!      2) s_compute_cylindrical_viscous_source_flux
        subroutine s_compute_abstract_viscous_source_flux(velL_vf, & ! -------------
                                                          dvelL_dx_vf, &
                                                          dvelL_dy_vf, &
                                                          velR_vf, &
                                                          dvelR_dx_vf, &
                                                          dvelR_dy_vf, &
                                                          flux_src_vf, &
                                                          norm_dir, &
                                                          ix, iy)

            import :: scalar_field, int_bounds_info, num_dims, sys_size

            type(scalar_field), &
                dimension(num_dims), &
                intent(IN) :: velL_vf, velR_vf, &
                              dvelL_dx_vf, dvelR_dx_vf, &
                              dvelL_dy_vf, dvelR_dy_vf

            type(scalar_field), &
                dimension(sys_size), &
                intent(INOUT) :: flux_src_vf

            integer, intent(IN) :: norm_dir

            type(int_bounds_info), intent(IN) :: ix, iy

        end subroutine s_compute_abstract_viscous_source_flux

    end interface ! ============================================================



    !> The cell-boundary values of the fluxes (src - source) that are computed
    !! through the chosen Riemann problem solver, and the direct evaluation of
    !! source terms, by using the left and right states given in qK_prim_rs_vf,
    !! dqK_prim_ds_vf and kappaK_rs_vf, where ds = dx, dy.
    !> @{
    real(kind(0d0)), allocatable, dimension(:, :, :) :: flux_rsx_vf, flux_src_rsx_vf
    real(kind(0d0)), allocatable, dimension(:, :, :) :: flux_rsy_vf, flux_src_rsy_vf

    !> @}

    !! The cell-boundary values of the geometrical source flux that are computed
    !! through the chosen Riemann problem solver by using the left and right
    !! states given in qK_prim_rs_vf. Currently 2D axisymmetric for inviscid only.

    ! The cell-boundary values of the velocity. vel_src_rs_vf is determined as
    ! part of Riemann problem solution and is used to evaluate the source flux.
    real(kind(0d0)), allocatable, dimension(:, :, :) :: vel_src_rsx_vf
    real(kind(0d0)), allocatable, dimension(:, :, :) :: vel_src_rsy_vf

    real(kind(0d0)), allocatable, dimension(:, :, :) :: mom_sp_rsx_vf
    real(kind(0d0)), allocatable, dimension(:, :, :) :: mom_sp_rsy_vf

    !> @name Left and right, WENO-reconstructed, cell-boundary values of cell-average
    !! partial densities, density, velocity, pressure, internal energy, energy, enthalpy, volume
    !! fractions, mass fractions, the specific heat ratio and liquid stiffness functions, speed
    !! of sound, shear and volume Reynolds numbers and the Weber numbers. These
    !! variables are left and right states of the Riemann problem obtained from
    !! qK_prim_rs_vf and kappaK_rs_vf.
    !> @{
    real(kind(0d0)), allocatable, dimension(:) :: alpha_rho_L, alpha_rho_R
    real(kind(0d0)) :: rho_L, rho_R
    real(kind(0d0)), allocatable, dimension(:) :: vel_L, vel_R
    real(kind(0d0)) :: pres_L, pres_R
    real(kind(0d0)) :: E_L, E_R
    real(kind(0d0)) :: H_L, H_R
    real(kind(0d0)), allocatable, dimension(:) :: alpha_L, alpha_R
    real(kind(0d0)) :: Y_L, Y_R
    real(kind(0d0)) :: gamma_L, gamma_R
    real(kind(0d0)) :: pi_inf_L, pi_inf_R
    real(kind(0d0)) :: c_L, c_R
    real(kind(0d0)), dimension(2) :: Re_L, Re_R
    real(kind(0d0)), allocatable, dimension(:) :: tau_e_L, tau_e_R
    real(kind(0d0)), allocatable, dimension(:) :: G_L, G_R

!$acc declare create(alpha_rho_L, alpha_rho_R,rho_L, rho_R,vel_L, vel_R,pres_L, pres_R, &
!$acc    E_L, E_R, H_L, H_R, alpha_L, alpha_R, Y_L, Y_R, gamma_L, gamma_R,pi_inf_L, pi_inf_R, &
!$acc    c_L, c_R,Re_L, Re_R,tau_e_L, tau_e_R, G_L, G_R)

    !> @}

    !> @name Roe or arithmetic average density, velocity, enthalpy, volume fractions,
    !! specific heat ratio function, speed of sound, shear and volume Reynolds
    !! numbers, Weber numbers and curvatures, at the cell-boundaries, computed
    !! from the left and the right states of the Riemann problem
    !> @{
    real(kind(0d0)) :: rho_avg
    real(kind(0d0)) :: H_avg
    type(scalar_field), allocatable, dimension(:) :: alpha_avg_rs_vf
    real(kind(0d0)) :: gamma_avg
    real(kind(0d0)) :: c_avg
    real(kind(0d0)), allocatable, dimension(:, :, :) :: Re_avg_rsx_vf
    real(kind(0d0)), allocatable, dimension(:, :, :) :: Re_avg_rsy_vf
!$acc declare create(rho_avg, H_avg, alpha_avg_rs_vf, gamma_avg, c_avg,  Re_avg_rsx_vf, Re_avg_rsy_vf)
    !> @}

    !> @name Left, right and star (S) region wave speeds
    !> @{
    real(kind(0d0)) :: s_L, s_R, s_S
    !> @}

    !> Minus (M) and plus (P) wave speeds
    !> @{
    real(kind(0d0)) :: s_M, s_P
    !> @}

    !> Minus and plus wave speeds functions
    !> @{
    real(kind(0d0)) :: xi_M, xi_P
    !> @}
    real(kind(0d0)) :: xi_L, xi_R

!$acc declare create(s_L, s_R, s_S, s_M, s_P, xi_M, xi_P, xi_L, xi_R)

    procedure(s_abstract_riemann_solver), &
        pointer :: s_riemann_solver => null() !<
    !! Pointer to the procedure that is utilized to calculate either the HLL,
    !! HLLC or exact intercell fluxes, based on the choice of Riemann solver

    procedure(s_compute_abstract_wave_speeds), &
        pointer :: s_compute_wave_speeds => null() !<
    !! Pointer to the subroutine that is utilized to compute the wave speeds of
    !! the Riemann problem either directly or by the means of pressure-velocity
    !! estimates, based on the selected method of estimation of the wave speeds

    procedure(s_compute_abstract_viscous_source_flux), &
        pointer :: s_compute_viscous_source_flux => null() !<
    !! Pointer to the subroutine that is utilized to compute the viscous source
    !! flux for either Cartesian or cylindrical geometries.

    !> @name Indical bounds in the s1-, s2- and s3-directions
    !> @{
    type(int_bounds_info) :: is1, is2
    type(int_bounds_info) :: isx, isy
    !> @}
!$acc declare create( &
!$acc    is1, is2, isx, isy)

!$acc declare create(&
!$acc    flux_rsx_vf, flux_src_rsx_vf, flux_rsy_vf, flux_src_rsy_vf, vel_src_rsx_vf, vel_src_rsy_vf )


    real(kind(0d0)), allocatable, dimension(:, :) :: Res
!$acc declare create(Res)

contains


    subroutine s_hllc_riemann_solver(qL_prim_rsx_vf, qL_prim_rsy_vf, dqL_prim_dx_vf, & ! ------
                                     dqL_prim_dy_vf, &
                                     qL_prim_vf, &
                                     qR_prim_rsx_vf, qR_prim_rsy_vf, dqR_prim_dx_vf, &
                                     dqR_prim_dy_vf, &
                                     qR_prim_vf, &
                                     q_prim_vf, &
                                     flux_vf, flux_src_vf, &
                                     norm_dir, ix, iy)

        real(kind(0d0)), dimension(startx:, starty:, 1:), intent(INOUT) :: qL_prim_rsx_vf, qL_prim_rsy_vf, qR_prim_rsx_vf, qR_prim_rsy_vf
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf

        type(scalar_field), allocatable, dimension(:), intent(INOUT) :: qL_prim_vf, qR_prim_vf

        type(scalar_field), &
            allocatable, dimension(:), &
            intent(INOUT) :: dqL_prim_dx_vf, dqR_prim_dx_vf, &
                             dqL_prim_dy_vf, dqR_prim_dy_vf

        ! Intercell fluxes
        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_vf, flux_src_vf

        integer, intent(IN) :: norm_dir
        type(int_bounds_info), intent(IN) :: ix, iy

        real(kind(0d0)), dimension(num_fluids) :: alpha_rho_L, alpha_rho_R
        real(kind(0d0)) :: rho_L, rho_R
        real(kind(0d0)), dimension(num_dims) :: vel_L, vel_R
        real(kind(0d0)) :: pres_L, pres_R
        real(kind(0d0)) :: E_L, E_R
        real(kind(0d0)) :: H_L, H_R
        real(kind(0d0)), dimension(num_fluids) :: alpha_L, alpha_R
        real(kind(0d0)) :: Y_L, Y_R
        real(kind(0d0)) :: gamma_L, gamma_R
        real(kind(0d0)) :: pi_inf_L, pi_inf_R
        real(kind(0d0)) :: c_L, c_R
        real(kind(0d0)), dimension(2) :: Re_L, Re_R

        real(kind(0d0)) :: rho_avg
        real(kind(0d0)) :: H_avg
        real(kind(0d0)) :: gamma_avg
        real(kind(0d0)) :: c_avg

        real(kind(0d0)) :: s_L, s_R, s_M, s_P, s_S
        real(kind(0d0)) :: xi_L, xi_R !< Left and right wave speeds functions
        real(kind(0d0)) :: xi_M, xi_P

        real(kind(0d0)) :: PbwR3Lbar, Pbwr3Rbar
        real(kind(0d0)) :: R3Lbar, R3Rbar
        real(kind(0d0)) :: R3V2Lbar, R3V2Rbar

        real(kind(0d0)) :: vel_L_rms, vel_R_rms, vel_avg_rms
        real(kind(0d0)) :: blkmod1, blkmod2
        real(kind(0d0)) :: pres_SL, pres_SR, Ms_L, Ms_R

        integer :: i, j, k, l, q !< Generic loop iterators
        integer :: idx1, idxi

        ! Populating the buffers of the left and right Riemann problem
        ! states variables, based on the choice of boundary conditions
        call s_populate_riemann_states_variables_buffers( &
            qL_prim_rsx_vf, qL_prim_rsy_vf, dqL_prim_dx_vf, &
            dqL_prim_dy_vf, &
            qL_prim_vf, &
            qR_prim_rsx_vf, qR_prim_rsy_vf, dqR_prim_dx_vf, &
            dqR_prim_dy_vf, &
            qR_prim_vf, &
            norm_dir, ix, iy)

        ! Reshaping inputted data based on dimensional splitting direction
        call s_initialize_riemann_solver( &
            q_prim_vf, &
            flux_vf, flux_src_vf, &
            norm_dir, ix, iy)
        #:for NORM_DIR, XYZ in [(1, 'x'), (2, 'y')]

            if (norm_dir == ${NORM_DIR}$) then


                !$acc parallel loop collapse(2) gang vector default(present) private(vel_L, vel_R, Re_L, Re_R, &
                !$acc rho_avg, h_avg, gamma_avg, s_L, s_R, s_S, vel_avg_rms)
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end
                            idx1 = 1; if (dir_idx(1) == 2) idx1 = 2; if (dir_idx(1) == 3) idx1 = 3

                            vel_L_rms = 0d0; vel_R_rms = 0d0
                            !$acc loop seq
                            do i = 1, num_dims
                                vel_L(i) = qL_prim_rs${XYZ}$_vf(j, k,contxe + i)
                                vel_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k,contxe + i)
                                vel_L_rms = vel_L_rms + vel_L(i)**2d0
                                vel_R_rms = vel_R_rms + vel_R(i)**2d0
                            end do

                            pres_L = qL_prim_rs${XYZ}$_vf(j, k,E_idx)
                            pres_R = qR_prim_rs${XYZ}$_vf(j + 1, k,E_idx)

                            rho_L = 0d0
                            gamma_L = 0d0
                            pi_inf_L = 0d0

                            rho_R = 0d0
                            gamma_R = 0d0
                            pi_inf_R = 0d0

                            !$acc loop seq
                            do i = 1, num_fluids
                                rho_L = rho_L + qL_prim_rs${XYZ}$_vf(j, k,i)
                                gamma_L = gamma_L + qL_prim_rs${XYZ}$_vf(j, k,E_idx + i)*gammas(i)
                                pi_inf_L = pi_inf_L + qL_prim_rs${XYZ}$_vf(j, k,E_idx + i)*pi_infs(i)

                                rho_R = rho_R + qR_prim_rs${XYZ}$_vf(j + 1, k,i)
                                gamma_R = gamma_R + qR_prim_rs${XYZ}$_vf(j + 1, k,E_idx + i)*gammas(i)
                                pi_inf_R = pi_inf_R + qR_prim_rs${XYZ}$_vf(j + 1, k,E_idx + i)*pi_infs(i)
                            end do

                            if (any(Re_size > 0)) then
                                !$acc loop seq
                                do i = 1, 2
                                    Re_L(i) = dflt_real

                                    if (Re_size(i) > 0) Re_L(i) = 0d0

                                    !$acc loop seq
                                    do q = 1, Re_size(i)
                                        Re_L(i) = qL_prim_rs${XYZ}$_vf(j, k,E_idx + Re_idx(i, q))/Res(i, q) &
                                                  + Re_L(i)
                                    end do

                                    Re_L(i) = 1d0/max(Re_L(i), sgm_eps)

                                end do

                                !$acc loop seq
                                do i = 1, 2
                                    Re_R(i) = dflt_real

                                    if (Re_size(i) > 0) Re_R(i) = 0d0

                                    !$acc loop seq
                                    do q = 1, Re_size(i)
                                        Re_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k,E_idx + Re_idx(i, q))/Res(i, q) &
                                                  + Re_R(i)
                                    end do

                                    Re_R(i) = 1d0/max(Re_R(i), sgm_eps)
                                end do
                            end if

                            E_L = gamma_L*pres_L + pi_inf_L + 5d-1*rho_L*vel_L_rms

                            E_R = gamma_R*pres_R + pi_inf_R + 5d-1*rho_R*vel_R_rms

                            H_L = (E_L + pres_L)/rho_L
                            H_R = (E_R + pres_R)/rho_R

                            @:compute_average_state()


                            c_avg = sqrt((H_avg - 5d-1*vel_avg_rms)/gamma_avg)

                            c_L = ((H_L - 5d-1*vel_L_rms)/gamma_L)

                            c_R = ((H_R - 5d-1*vel_R_rms)/gamma_R)

                            c_L = sqrt(c_L)
                            c_R = sqrt(c_R)

                            if (any(Re_size > 0)) then
                                !$acc loop seq
                                do i = 1, 2
                                    Re_avg_rs${XYZ}$_vf(j, k,i) = 2d0/(1d0/Re_L(i) + 1d0/Re_R(i))
                                end do
                            end if

                            if (wave_speeds == 1) then
                                s_L = min(vel_L(idx1) - c_L, vel_R(idx1) - c_R)
                                s_R = max(vel_R(idx1) + c_R, vel_L(idx1) + c_L)

                                s_S = (pres_R - pres_L + rho_L*vel_L(idx1)* &
                                       (s_L - vel_L(idx1)) - &
                                       rho_R*vel_R(idx1)* &
                                       (s_R - vel_R(idx1))) &
                                      /(rho_L*(s_L - vel_L(idx1)) - &
                                        rho_R*(s_R - vel_R(idx1)))

                            elseif (wave_speeds == 2) then
                                pres_SL = 5d-1*(pres_L + pres_R + rho_avg*c_avg* &
                                                (vel_L(idx1) - &
                                                 vel_R(idx1)))

                                pres_SR = pres_SL

                                Ms_L = max(1d0, sqrt(1d0 + ((5d-1 + gamma_L)/(1d0 + gamma_L))* &
                                                     (pres_SL/pres_L - 1d0)*pres_L/ &
                                                     ((pres_L + pi_inf_L/(1d0 + gamma_L)))))
                                Ms_R = max(1d0, sqrt(1d0 + ((5d-1 + gamma_R)/(1d0 + gamma_R))* &
                                                     (pres_SR/pres_R - 1d0)*pres_R/ &
                                                     ((pres_R + pi_inf_R/(1d0 + gamma_R)))))

                                s_L = vel_L(idx1) - c_L*Ms_L
                                s_R = vel_R(idx1) + c_R*Ms_R

                                s_S = 5d-1*((vel_L(idx1) + vel_R(idx1)) + &
                                            (pres_L - pres_R)/ &
                                            (rho_avg*c_avg))
                            end if

                            ! follows Einfeldt et al.
                            ! s_M/P = min/max(0.,s_L/R)
                            s_M = min(0d0, s_L); s_P = max(0d0, s_R)

                            ! goes with q_star_L/R = xi_L/R * (variable)
                            ! xi_L/R = ( ( s_L/R - u_L/R )/(s_L/R - s_star) )
                            xi_L = (s_L - vel_L(idx1))/(s_L - s_S)
                            xi_R = (s_R - vel_R(idx1))/(s_R - s_S)

                            ! goes with numerical velocity in x/y/z directions
                            ! xi_P/M = 0.5 +/m sgn(0.5,s_star)
                            xi_M = (5d-1 + sign(5d-1, s_S))
                            xi_P = (5d-1 - sign(5d-1, s_S))

                            !$acc loop seq
                            do i = 1, contxe
                                flux_rs${XYZ}$_vf(j, k,i) = &
                                    xi_M*qL_prim_rs${XYZ}$_vf(j, k,i) &
                                    *(vel_L(idx1) + s_M*(xi_L - 1d0)) &
                                    + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k,i) &
                                    *(vel_R(idx1) + s_P*(xi_R - 1d0))
                            end do

                            ! Momentum flux.
                            ! f = \rho u u + p I, q = \rho u, q_star = \xi * \rho*(s_star, v, w)

                            !$acc loop seq
                            do i = 1, num_dims
                                idxi = dir_idx(i)
                                flux_rs${XYZ}$_vf(j, k,contxe + idxi) = &
                                    xi_M*(rho_L*(vel_L(idx1)* &
                                                 vel_L(idxi) + &
                                                 s_M*(xi_L*(dir_flg(idxi)*s_S + &
                                                            (1d0 - dir_flg(idxi))* &
                                                            vel_L(idxi)) - vel_L(idxi))) + &
                                          dir_flg(idxi)*(pres_L)) &
                                    + xi_P*(rho_R*(vel_R(idx1)* &
                                                   vel_R(idxi) + &
                                                   s_P*(xi_R*(dir_flg(idxi)*s_S + &
                                                              (1d0 - dir_flg(idxi))* &
                                                              vel_R(idxi)) - vel_R(idxi))) + &
                                            dir_flg(idxi)*(pres_R))
                                ! if (j==0) print*, 'flux_rs_vf', flux_rs_vf(cont_idx%end+dir_idx(i))%sf(j,k,l)
                            end do

                            ! Energy flux.
                            ! f = u*(E+p), q = E, q_star = \xi*E+(s-u)(\rho s_star + p/(s-u))

                            flux_rs${XYZ}$_vf(j, k,E_idx) = &
                                xi_M*(vel_L(idx1)*(E_L + pres_L) + &
                                      s_M*(xi_L*(E_L + (s_S - vel_L(idx1))* &
                                                 (rho_L*s_S + pres_L/ &
                                                  (s_L - vel_L(idx1)))) - E_L)) &
                                + xi_P*(vel_R(idx1)*(E_R + pres_R) + &
                                        s_P*(xi_R*(E_R + (s_S - vel_R(idx1))* &
                                                   (rho_R*s_S + pres_R/ &
                                                    (s_R - vel_R(idx1)))) - E_R))

                            ! Volume fraction flux

                            !$acc loop seq
                            do i = advxb, advxe
                                flux_rs${XYZ}$_vf(j, k,i) = &
                                    xi_M*qL_prim_rs${XYZ}$_vf(j, k,i) &
                                    *(vel_L(idx1) + s_M*(xi_L - 1d0)) &
                                    + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k,i) &
                                    *(vel_R(idx1) + s_P*(xi_R - 1d0))
                            end do

                            ! Source for volume fraction advection equation
                            !$acc loop seq
                            do i = 1, num_dims
                                idxi = dir_idx(i)
                                vel_src_rs${XYZ}$_vf(j, k,idxi) = &
                                    xi_M*(vel_L(idxi) + &
                                          dir_flg(idxi)* &
                                          s_M*(xi_L - 1d0)) &
                                    + xi_P*(vel_R(idxi) + &
                                            dir_flg(idxi)* &
                                            s_P*(xi_R - 1d0))

                            end do
                            flux_src_rs${XYZ}$_vf(j, k,advxb) = vel_src_rs${XYZ}$_vf(j, k,idx1)
                        end do
                    end do
            end if
        #:endfor

        if (any(Re_size > 0)) then
            if (weno_Re_flux) then
                call s_compute_viscous_source_flux( &
                    qL_prim_vf(momxb:momxe), &
                    dqL_prim_dx_vf(momxb:momxe), &
                    dqL_prim_dy_vf(momxb:momxe), &
                    qR_prim_vf(momxb:momxe), &
                    dqR_prim_dx_vf(momxb:momxe), &
                    dqR_prim_dy_vf(momxb:momxe), &
                    flux_src_vf, norm_dir, ix, iy)
            else
                call s_compute_viscous_source_flux( &
                    q_prim_vf(momxb:momxe), &
                    dqL_prim_dx_vf(momxb:momxe), &
                    dqL_prim_dy_vf(momxb:momxe), &
                    q_prim_vf(momxb:momxe), &
                    dqR_prim_dx_vf(momxb:momxe), &
                    dqR_prim_dy_vf(momxb:momxe), &
                    flux_src_vf, norm_dir, ix, iy)
            end if
        end if

        call s_finalize_riemann_solver(flux_vf, flux_src_vf, &
                                       norm_dir, ix, iy)

    end subroutine s_hllc_riemann_solver


    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_riemann_solvers_module() ! ---------------------

        ! Allocating the variables that will be utilized to formulate the
        ! left, right, and average states of the Riemann problem, as well
        ! the Riemann problem solution
        integer :: i, j

        if (any(Re_size > 0)) then
            allocate (Res(1:2, 1:maxval(Re_size)))
        end if

        if (any(Re_size > 0)) then
            do i = 1, 2
                do j = 1, Re_size(i)
                    Res(i, j) = fluid_pp(Re_idx(i, j))%Re(i)
                end do
            end do
!$acc update device(Res, Re_idx, Re_size)
        end if

        ! Associating procedural pointer to the subroutine that will be
        ! utilized to calculate the solution of a given Riemann problem
        s_riemann_solver => s_hllc_riemann_solver

        ! Associating the procedural pointers to the procedures that will be
        ! utilized to compute the average state and estimate the wave speeds

        ! Associating procedural pointer to the subroutine that will be
        ! utilized to compute the viscous source flux
        s_compute_viscous_source_flux => s_compute_cartesian_viscous_source_flux

        ! Associating the procedural pointer to the appropriate subroutine
        ! that will be utilized in the conversion to the mixture variables
        s_convert_to_mixture_variables => s_convert_species_to_mixture_variables

        is1%beg = -1; is2%beg = 0
        is1%end = m; is2%end = n

        !allocate(qL_prim_rsx_vf(is1%beg:is1%end, is2%beg:is2%end,  1:sys_size))
        !allocate(qR_prim_rsx_vf(is1%beg + 1:is1%end + 1, is2%beg:is2%end,  1:sys_size))
        allocate (flux_rsx_vf(is1%beg:is1%end, &
                                   is2%beg:is2%end, &
                                    1:sys_size))
        allocate (flux_src_rsx_vf(is1%beg:is1%end, &
                                       is2%beg:is2%end, &
                                        advxb:sys_size))
        allocate (vel_src_rsx_vf(is1%beg:is1%end, &
                                      is2%beg:is2%end, &
                                       1:num_dims))

        if (any(Re_size > 0)) then
            allocate (Re_avg_rsx_vf(is1%beg:is1%end, &
                                         is2%beg:is2%end, &
                                          1:2))
        end if

        if (n == 0) return

        is1%beg = -1; is2%beg = 0; 
        is1%end = n; is2%end = m; 

        !allocate(qL_prim_rsy_vf(is1%beg:is1%end, is2%beg:is2%end,  1:sys_size))
        !allocate(qR_prim_rsy_vf(is1%beg + 1:is1%end + 1, is2%beg:is2%end,  1:sys_size))
        allocate (flux_rsy_vf(is1%beg:is1%end, &
                                   is2%beg:is2%end, &
                                    1:sys_size))
        allocate (flux_src_rsy_vf(is1%beg:is1%end, &
                                       is2%beg:is2%end, &
                                        advxb:sys_size))
        allocate (vel_src_rsy_vf(is1%beg:is1%end, &
                                      is2%beg:is2%end, &
                                       1:num_dims))


        if (any(Re_size > 0)) then
            allocate (Re_avg_rsy_vf(is1%beg:is1%end, &
                                         is2%beg:is2%end, &
                                          1:2))
        end if

    end subroutine s_initialize_riemann_solvers_module ! -------------------

    !>  The purpose of this subroutine is to populate the buffers
        !!      of the left and right Riemann states variables, depending
        !!      on the boundary conditions.
        !!  @param qL_prim_vf The  left WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables
        !!  @param qR_prim_vf The right WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables
        !!  @param dqL_prim_dx_vf The  left WENO-reconstructed cell-boundary values of the
        !!      first-order x-dir spatial derivatives
        !!  @param dqL_prim_dy_vf The  left WENO-reconstructed cell-boundary values of the
        !!      first-order y-dir spatial derivatives
        !!  @param dqR_prim_dx_vf The right WENO-reconstructed cell-boundary values of the
        !!      first-order x-dir spatial derivatives
        !!  @param dqR_prim_dy_vf The right WENO-reconstructed cell-boundary values of the
        !!      first-order y-dir spatial derivatives
        !!  @param gm_alphaL_vf  Left averaged gradient magnitude
        !!  @param gm_alphaR_vf Right averaged gradient magnitude
        !!  @param norm_dir Dir. splitting direction
        !!  @param ix Index bounds in the x-dir
        !!  @param iy Index bounds in the y-dir
    subroutine s_populate_riemann_states_variables_buffers( & ! ------------
        qL_prim_rsx_vf, qL_prim_rsy_vf, dqL_prim_dx_vf, &
        dqL_prim_dy_vf, &
        qL_prim_vf, &
        qR_prim_rsx_vf, qR_prim_rsy_vf, dqR_prim_dx_vf, &
        dqR_prim_dy_vf, &
        qR_prim_vf, &
        norm_dir, ix, iy)

        real(kind(0d0)), dimension(startx:, starty:, 1:), intent(INOUT) :: qL_prim_rsx_vf, qL_prim_rsy_vf, qR_prim_rsx_vf, qR_prim_rsy_vf

        type(scalar_field), &
            allocatable, dimension(:), &
            intent(INOUT) :: dqL_prim_dx_vf, dqR_prim_dx_vf, &
                             dqL_prim_dy_vf, dqR_prim_dy_vf, &
                             qL_prim_vf, qR_prim_vf

        integer, intent(IN) :: norm_dir

        type(int_bounds_info), intent(IN) :: ix, iy

        integer :: i, j, k, l !< Generic loop iterator

        if (norm_dir == 1) then
            is1 = ix; is2 = iy
            dir_idx = (/1, 2/); dir_flg = (/1d0, 0d0/)
        elseif (norm_dir == 2) then
            is1 = iy; is2 = ix
            dir_idx = (/2, 1/); dir_flg = (/0d0, 1d0/)
        end if


        isx = ix; isy = iy

        !$acc update device(is1, is2, dir_idx, dir_flg, isx, isy)

        ! Population of Buffers in x-direction =============================
        if (norm_dir == 1) then

            if (bc_x%beg == -4) then    ! Riemann state extrap. BC at beginning
                !$acc parallel loop collapse(2) gang vector default(present)
                do i = 1, sys_size
                        do k = is2%beg, is2%end
                            qL_prim_rsx_vf(-1, k,i) = &
                                qR_prim_rsx_vf(0, k,i)
                        end do
                end do

                if (any(Re_size > 0)) then
                    !$acc parallel loop collapse(2) gang vector default(present)
                    do i = momxb, momxe
                            do k = isy%beg, isy%end

                                dqL_prim_dx_vf(i)%sf(-1, k) = &
                                    dqR_prim_dx_vf(i)%sf(0, k)
                            end do
                    end do

                    if (n > 0) then
                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = momxb, momxe
                                do k = isy%beg, isy%end

                                    dqL_prim_dy_vf(i)%sf(-1, k) = &
                                        dqR_prim_dy_vf(i)%sf(0, k)
                                end do
                        end do

                    end if

                end if

            end if

            if (bc_x%end == -4) then    ! Riemann state extrap. BC at end

                !$acc parallel loop collapse(2) gang vector default(present)
                do i = 1, sys_size
                        do k = is2%beg, is2%end
                            qR_prim_rsx_vf(m + 1, k,i) = &
                                qL_prim_rsx_vf(m, k,i)
                        end do
                end do

                if (any(Re_size > 0)) then

                    !$acc parallel loop collapse(2) gang vector default(present)
                    do i = momxb, momxe
                            do k = isy%beg, isy%end

                                dqR_prim_dx_vf(i)%sf(m + 1, k) = &
                                    dqL_prim_dx_vf(i)%sf(m, k)
                            end do
                    end do

                    if (n > 0) then
                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = momxb, momxe
                                do k = isy%beg, isy%end

                                    dqR_prim_dy_vf(i)%sf(m + 1, k) = &
                                        dqL_prim_dy_vf(i)%sf(m, k)
                                end do
                        end do

                    end if

                end if

            end if
            ! END: Population of Buffers in x-direction ========================

            ! Population of Buffers in y-direction =============================
        elseif (norm_dir == 2) then

            if (bc_y%beg == -4) then    ! Riemann state extrap. BC at beginning
                !$acc parallel loop collapse(2) gang vector default(present)
                do i = 1, sys_size
                        do k = is2%beg, is2%end
                            qL_prim_rsy_vf(-1, k,i) = &
                                qR_prim_rsy_vf(0, k,i)
                        end do
                end do

                if (any(Re_size > 0)) then

                    !$acc parallel loop collapse(2) gang vector default(present)
                    do i = momxb, momxe
                            do j = isx%beg, isx%end
                                dqL_prim_dx_vf(i)%sf(j, -1) = &
                                    dqR_prim_dx_vf(i)%sf(j, 0)
                            end do
                    end do

                    !$acc parallel loop collapse(2) gang vector default(present)
                    do i = momxb, momxe
                            do j = isx%beg, isx%end
                                dqL_prim_dy_vf(i)%sf(j, -1) = &
                                    dqR_prim_dy_vf(i)%sf(j, 0)
                            end do
                    end do


                end if

            end if

            if (bc_y%end == -4) then    ! Riemann state extrap. BC at end

                !$acc parallel loop collapse(2) gang vector default(present)
                do i = 1, sys_size
                        do k = is2%beg, is2%end
                            qR_prim_rsy_vf(n + 1, k,i) = &
                                qL_prim_rsy_vf(n, k,i)
                        end do
                end do

                if (any(Re_size > 0)) then
                    !$acc parallel loop collapse(2) gang vector default(present)
                    do i = momxb, momxe
                            do j = isx%beg, isx%end
                                dqR_prim_dx_vf(i)%sf(j, n + 1) = &
                                    dqL_prim_dx_vf(i)%sf(j, n)
                            end do
                    end do

                    !$acc parallel loop collapse(2) gang vector default(present)
                    do i = momxb, momxe
                            do j = isx%beg, isx%end
                                dqR_prim_dy_vf(i)%sf(j, n + 1) = &
                                    dqL_prim_dy_vf(i)%sf(j, n)
                            end do
                    end do
                end if
            end if
            ! END: Population of Buffers in y-direction ========================
        end if

    end subroutine s_populate_riemann_states_variables_buffers ! -----------

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures needed to configure the chosen Riemann
        !!      solver algorithm.
        !!  @param qL_prim_vf The  left WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables
        !!  @param qR_prim_vf The right WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables
        !!  @param flux_vf Intra-cell fluxes
        !!  @param flux_src_vf Intra-cell fluxes sources
        !!  @param norm_dir Dir. splitting direction
        !!  @param ix Index bounds in the x-dir
        !!  @param iy Index bounds in the y-dir
        !!  @param q_prim_vf Cell-averaged primitive variables
    subroutine s_initialize_riemann_solver( &
        q_prim_vf, &
        flux_vf, flux_src_vf, &
        norm_dir, ix, iy)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_vf, flux_src_vf

        integer, intent(IN) :: norm_dir

        type(int_bounds_info), intent(IN) :: ix, iy

        integer :: i, j, k, l ! Generic loop iterators

        ! Reshaping Inputted Data in x-direction ===========================

        if (norm_dir == 1) then

            if (any(Re_size > 0)) then
!$acc parallel loop collapse(3) gang vector default(present)
                do i = momxb, E_idx
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                flux_src_vf(i)%sf(j, k) = 0d0
                            end do
                        end do
                end do
            end if

        elseif (norm_dir == 2) then
            ! Reshaping Inputted Data in y-direction ===========================

            if (any(Re_size > 0)) then
!$acc parallel loop collapse(3) gang vector default(present)
                do i = momxb, E_idx
                        do j = is1%beg, is1%end
                            do k = is2%beg, is2%end
                                flux_src_vf(i)%sf(k, j) = 0d0
                            end do
                        end do
                end do
            end if
        end if

    end subroutine s_initialize_riemann_solver ! ---------------------------

    !>  The goal of this subroutine is to evaluate and account
        !!      for the contribution of viscous stresses in the source
        !!      flux for the momentum and energy.
        !!  @param velL_vf  Left, WENO reconstructed, cell-boundary values of the velocity
        !!  @param velR_vf Right, WENO reconstructed, cell-boundary values of the velocity
        !!  @param dvelL_dx_vf  Left, WENO reconstructed cell-avg. x-dir derivative of the velocity
        !!  @param dvelL_dy_vf  Left, WENO reconstructed cell-avg. y-dir derivative of the velocity
        !!  @param dvelR_dx_vf Right, WENO reconstructed cell-avg. x-dir derivative of the velocity
        !!  @param dvelR_dy_vf Right, WENO reconstructed cell-avg. y-dir derivative of the velocity
        !!  @param flux_src_vf Intercell flux
        !!  @param norm_dir Dimensional splitting coordinate direction
        !!  @param ix Index bounds in  first coordinate direction
        !!  @param iy Index bounds in second coordinate direction
    subroutine s_compute_cartesian_viscous_source_flux(velL_vf, & ! -------------
                                                       dvelL_dx_vf, &
                                                       dvelL_dy_vf, &
                                                       velR_vf, &
                                                       dvelR_dx_vf, &
                                                       dvelR_dy_vf, &
                                                       flux_src_vf, &
                                                       norm_dir, &
                                                       ix, iy)

        type(scalar_field), &
            dimension(num_dims), &
            intent(IN) :: velL_vf, velR_vf, &
                          dvelL_dx_vf, dvelR_dx_vf, &
                          dvelL_dy_vf, dvelR_dy_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_src_vf

        integer, intent(IN) :: norm_dir

        type(int_bounds_info), intent(IN) :: ix, iy

        ! Arithmetic mean of the left and right, WENO-reconstructed, cell-
        ! boundary values of cell-average first-order spatial derivatives
        ! of velocity
        real(kind(0d0)), dimension(num_dims) :: dvel_avg_dx
        real(kind(0d0)), dimension(num_dims) :: dvel_avg_dy

        real(kind(0d0)), dimension(num_dims, num_dims) :: tau_Re !< Viscous stress tensor

        integer :: i, j, k, l !< Generic loop iterators

        ! Viscous Stresses in x-direction ==================================
        if (norm_dir == 1) then

            if (Re_size(1) > 0) then              ! Shear stresses
!$acc parallel loop collapse(2) gang vector default(present) private( dvel_avg_dx, tau_Re)
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dx(1) = 5d-1*(dvelL_dx_vf(1)%sf(j, k) &
                                                   + dvelR_dx_vf(1)%sf(j + 1, k))

                            tau_Re(1, 1) = (4d0/3d0)*dvel_avg_dx(1)/ &
                                           Re_avg_rsx_vf(j, k,1)

                            flux_src_vf(momxb)%sf(j, k) = &
                                flux_src_vf(momxb)%sf(j, k) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k) = &
                                flux_src_vf(E_idx)%sf(j, k) - &
                                vel_src_rsx_vf(j, k,1)* &
                                tau_Re(1, 1)

                        end do
                end do
            end if

            if (Re_size(2) > 0) then              ! Bulk stresses
!$acc parallel loop collapse(2) gang vector default(present) private( dvel_avg_dx, tau_Re)
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dx(1) = 5d-1*(dvelL_dx_vf(1)%sf(j, k) &
                                                   + dvelR_dx_vf(1)%sf(j + 1, k))

                            tau_Re(1, 1) = dvel_avg_dx(1)/ &
                                           Re_avg_rsx_vf(j, k,2)

                            flux_src_vf(momxb)%sf(j, k) = &
                                flux_src_vf(momxb)%sf(j, k) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k) = &
                                flux_src_vf(E_idx)%sf(j, k) - &
                                vel_src_rsx_vf(j, k,1)* &
                                tau_Re(1, 1)

                        end do
                end do
            end if

            if (n == 0) return

            if (Re_size(1) > 0) then              ! Shear stresses
!$acc parallel loop collapse(2) gang vector default(present) private(dvel_avg_dx, dvel_avg_dy, tau_Re)
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

!$acc loop seq
                            do i = 1, 2
                                dvel_avg_dy(i) = &
                                    5d-1*(dvelL_dy_vf(i)%sf(j, k) &
                                          + dvelR_dy_vf(i)%sf(j + 1, k))
                            end do

                            dvel_avg_dx(2) = 5d-1*(dvelL_dx_vf(2)%sf(j, k) &
                                                   + dvelR_dx_vf(2)%sf(j + 1, k))

                            tau_Re(1, 1) = -(2d0/3d0)*dvel_avg_dy(2)/ &
                                           Re_avg_rsx_vf(j, k,1)

                            tau_Re(1, 2) = (dvel_avg_dy(1) + dvel_avg_dx(2))/ &
                                           Re_avg_rsx_vf(j, k,1)

!$acc loop seq
                            do i = 1, 2

                                flux_src_vf(contxe + i)%sf(j, k) = &
                                    flux_src_vf(contxe + i)%sf(j, k) - &
                                    tau_Re(1, i)

                                flux_src_vf(E_idx)%sf(j, k) = &
                                    flux_src_vf(E_idx)%sf(j, k) - &
                                    vel_src_rsx_vf(j, k,i)* &
                                    tau_Re(1, i)

                            end do

                    end do
                end do
            end if

            if (Re_size(2) > 0) then              ! Bulk stresses
!$acc parallel loop collapse(2) gang vector default(present) private( dvel_avg_dy, tau_Re)
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dy(2) = 5d-1*(dvelL_dy_vf(2)%sf(j, k) &
                                                   + dvelR_dy_vf(2)%sf(j + 1, k))

                            tau_Re(1, 1) = dvel_avg_dy(2)/ &
                                           Re_avg_rsx_vf(j, k,2)

                            flux_src_vf(momxb)%sf(j, k) = &
                                flux_src_vf(momxb)%sf(j, k) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k) = &
                                flux_src_vf(E_idx)%sf(j, k) - &
                                vel_src_rsx_vf(j, k,1)* &
                                tau_Re(1, 1)

                        end do
                    end do
            end if
            ! END: Viscous Stresses in x-direction =============================

            ! Viscous Stresses in y-direction ==================================
        elseif (norm_dir == 2) then

            if (Re_size(1) > 0) then              ! Shear stresses
!$acc parallel loop collapse(2) gang vector default(present) private( dvel_avg_dx, dvel_avg_dy, tau_Re)
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

!$acc loop seq
                            do i = 1, 2

                                dvel_avg_dx(i) = &
                                    5d-1*(dvelL_dx_vf(i)%sf(j, k) &
                                          + dvelR_dx_vf(i)%sf(j, k + 1))

                                dvel_avg_dy(i) = &
                                    5d-1*(dvelL_dy_vf(i)%sf(j, k) &
                                          + dvelR_dy_vf(i)%sf(j, k + 1))

                            end do

                            tau_Re(2, 1) = (dvel_avg_dy(1) + dvel_avg_dx(2))/ &
                                           Re_avg_rsy_vf(k, j,1)

                            tau_Re(2, 2) = (4d0*dvel_avg_dy(2) &
                                            - 2d0*dvel_avg_dx(1))/ &
                                           (3d0*Re_avg_rsy_vf(k, j,1))

!$acc loop seq
                            do i = 1, 2

                                flux_src_vf(contxe + i)%sf(j, k) = &
                                    flux_src_vf(contxe + i)%sf(j, k) - &
                                    tau_Re(2, i)

                                flux_src_vf(E_idx)%sf(j, k) = &
                                    flux_src_vf(E_idx)%sf(j, k) - &
                                    vel_src_rsy_vf(k, j,i)* &
                                    tau_Re(2, i)

                            end do

                        end do
                    end do
            end if

            if (Re_size(2) > 0) then              ! Bulk stresses
!$acc parallel loop collapse(2) gang vector default(present) private( dvel_avg_dx, dvel_avg_dy, tau_Re)
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dx(1) = 5d-1*(dvelL_dx_vf(1)%sf(j, k) &
                                                   + dvelR_dx_vf(1)%sf(j, k + 1))

                            dvel_avg_dy(2) = 5d-1*(dvelL_dy_vf(2)%sf(j, k) &
                                                   + dvelR_dy_vf(2)%sf(j, k + 1))

                            tau_Re(2, 2) = (dvel_avg_dx(1) + dvel_avg_dy(2))/ &
                                           Re_avg_rsy_vf(k, j,2)

                            flux_src_vf(momxb + 1)%sf(j, k) = &
                                flux_src_vf(momxb + 1)%sf(j, k) - &
                                tau_Re(2, 2)

                            flux_src_vf(E_idx)%sf(j, k) = &
                                flux_src_vf(E_idx)%sf(j, k) - &
                                vel_src_rsy_vf(k, j,2)* &
                                tau_Re(2, 2)

                        end do
                    end do
            end if
        end if

    end subroutine s_compute_cartesian_viscous_source_flux ! -------------------------

    !>  Deallocation and/or disassociation procedures that are
        !!      needed to finalize the selected Riemann problem solver
        !!  @param flux_vf       Intercell fluxes
        !!  @param flux_src_vf   Intercell source fluxes
        !!  @param norm_dir Dimensional splitting coordinate direction
        !!  @param ix   Index bounds in  first coordinate direction
        !!  @param iy   Index bounds in second coordinate direction
    subroutine s_finalize_riemann_solver(flux_vf, flux_src_vf, & ! --------
                                         norm_dir, ix, iy)

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_vf, flux_src_vf

        integer, intent(IN) :: norm_dir

        type(int_bounds_info), intent(IN) :: ix, iy

        integer :: i, j, k, l !< Generic loop iterators

        ! Reshaping Outputted Data in y-direction ==========================
        if (norm_dir == 2) then
!$acc parallel loop collapse(3) gang vector default(present)
            do i = 1, sys_size
                    do j = is1%beg, is1%end
                        do k = is2%beg, is2%end
                            flux_vf(i)%sf(k, j) = &
                                flux_rsy_vf(j, k,i)
                        end do
                    end do
            end do

            !$acc parallel loop collapse(2) gang vector default(present)
                do j = is1%beg, is1%end
                    do k = is2%beg, is2%end
                        flux_src_vf(advxb)%sf(k, j) = &
                            flux_src_rsy_vf(j, k,advxb)
                    end do
                end do
        elseif (norm_dir == 1) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do i = 1, sys_size
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end
                            flux_vf(i)%sf(j, k) = &
                                flux_rsx_vf(j, k,i)
                        end do
                    end do
            end do

            !$acc parallel loop collapse(2) gang vector default(present)
                do k = is2%beg, is2%end
                    do j = is1%beg, is1%end
                        flux_src_vf(advxb)%sf(j, k) = &
                            flux_src_rsx_vf(j, k,advxb)
                    end do
                end do
        end if

        ! ==================================================================

        ! ==================================================================

    end subroutine s_finalize_riemann_solver ! -----------------------------

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_riemann_solvers_module() ! -----------------------

        ! Deallocating the variables that were utilized to formulate the
        ! left, right and average states of the Riemann problem, as well
        ! the Riemann problem solution

        integer :: i


        ! Disassociating procedural pointer to the subroutine which was
        ! utilized to calculate the solution of a given Riemann problem
        s_riemann_solver => null()

        ! Disassociating the procedural pointers to the procedures that were
        ! utilized to compute the average state and estimate the wave speeds
        s_compute_wave_speeds => null()

        ! Disassociating procedural pointer to the subroutine which was
        ! utilized to calculate the viscous source flux
        s_compute_viscous_source_flux => null()

        ! Disassociating the pointer to the procedure that was utilized to
        ! to convert mixture or species variables to the mixture variables
        s_convert_to_mixture_variables => null()

        if (Re_size(1) > 0) then
            deallocate (Re_avg_rsx_vf)
        end if
        deallocate (vel_src_rsx_vf)
        deallocate (flux_rsx_vf)
        deallocate (flux_src_rsx_vf)

        if (n == 0) return

        if (Re_size(1) > 0) then
            deallocate (Re_avg_rsy_vf)
        end if
        deallocate (vel_src_rsy_vf)
        deallocate (flux_rsy_vf)
        deallocate (flux_src_rsy_vf)

    end subroutine s_finalize_riemann_solvers_module ! ---------------------

end module m_riemann_solvers
