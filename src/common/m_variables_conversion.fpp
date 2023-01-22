!>
!! @file m_variables_conversion.f90
!! @brief Contains module m_variables_conversion

#:include 'macros.fpp'

!> @brief This module consists of subroutines used in the conversion of the
!!              conservative variables into the primitive ones and vice versa. In
!!              addition, the module also contains the subroutines used to obtain
!!              the mixture variables and the subroutines used to compute pressure.
module m_variables_conversion

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_variables_conversion_module, &
            s_convert_to_mixture_variables, &
            s_convert_species_to_mixture_variables, &
            s_convert_species_to_mixture_variables_acc, &
            s_convert_conservative_to_primitive_variables, &
            s_convert_primitive_to_conservative_variables, &
            s_convert_primitive_to_flux_variables, &
            s_compute_pressure, &
            s_finalize_variables_conversion_module

    !> Abstract interface to two subroutines designed for the transfer/conversion
    !! of the mixture/species variables to the mixture variables

    abstract interface ! =======================================================

        !> Structure of the s_convert_mixture_to_mixture_variables
        !!      and s_convert_species_to_mixture_variables subroutines
        !!  @param q_vf Conservative or primitive variables
        !!  @param i First-coordinate cell index
        !!  @param j First-coordinate cell index
        !!  @param k First-coordinate cell index
        !!  @param rho Density
        !!  @param gamma Specific heat ratio function
        !!  @param pi_inf Liquid stiffness function
        subroutine s_convert_xxxxx_to_mixture_variables(q_vf, i, j, &
                                                        rho, gamma, pi_inf, Re_K)

            ! Importing the derived type scalar_field from m_derived_types.f90
            ! and global variable sys_size, from m_global_variables.f90, as
            ! the abstract interface does not inherently have access to them
            import :: scalar_field, sys_size, num_fluids

            type(scalar_field), dimension(sys_size), intent(IN) :: q_vf

            integer, intent(IN) :: i, j

            real(kind(0d0)), intent(OUT), target :: rho
            real(kind(0d0)), intent(OUT), target :: gamma
            real(kind(0d0)), intent(OUT), target :: pi_inf

            real(kind(0d0)), optional, dimension(2), intent(OUT) :: Re_K

        end subroutine s_convert_xxxxx_to_mixture_variables

    end interface ! ============================================================

    integer, public :: ixb, ixe, iyb, iye
    !$acc declare create(ixb, ixe, iyb, iye)

    !! In simulation, gammas and pi_infs is already declared in m_global_variables
#ifndef MFC_SIMULATION
    real(kind(0d0)), allocatable, dimension(:) :: gammas, pi_infs
    !$acc declare create(gammas, pi_infs)
#endif

    real(kind(0d0)), allocatable, dimension(:, :) :: Res
    !$acc declare create(Res)

    integer :: is1b, is2b, is1e, is2e
    !$acc declare create(is1b, is2b, is1e, is2e)

    real(kind(0d0)), allocatable, dimension(:, :), target, public :: rho_sf !< Scalar density function
    real(kind(0d0)), allocatable, dimension(:, :), target, public :: gamma_sf !< Scalar sp. heat ratio function
    real(kind(0d0)), allocatable, dimension(:, :), target, public :: pi_inf_sf !< Scalar liquid stiffness function   

    procedure(s_convert_xxxxx_to_mixture_variables), pointer :: s_convert_to_mixture_variables => null() 

contains

    !>  This procedure conditionally calculates the appropriate pressure
        !! @param energy Energy
        !! @param dyn_p Dynamic Pressure
        !! @param pi_inf Liquid Stiffness
        !! @param gamma Specific Heat Ratio
        !! @param pres Pressure to calculate
    subroutine s_compute_pressure(energy, dyn_p, pi_inf, gamma, pres)      
!$acc routine seq

        real(kind(0d0)), intent(IN) :: energy, dyn_p, pi_inf, gamma
        real(kind(0d0)), intent(OUT) :: pres

        pres = (energy - dyn_p - pi_inf)/gamma

    end subroutine s_compute_pressure



    subroutine s_convert_species_to_mixture_variables(q_vf, k, l, rho, gamma, pi_inf, Re_K)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_vf

        integer, intent(IN) :: k, l

        real(kind(0d0)), intent(OUT), target :: rho
        real(kind(0d0)), intent(OUT), target :: gamma
        real(kind(0d0)), intent(OUT), target :: pi_inf

        real(kind(0d0)), optional, dimension(2), intent(OUT) :: Re_K

        real(kind(0d0)), dimension(num_fluids) :: alpha_rho_K, alpha_K !<
            !! Partial densities and volume fractions

        integer :: i, j !< Generic loop iterator

        real(kind(0d0)), pointer :: rho_K, gamma_K, pi_inf_K

        !> Post process requires rho_sf/gamma_sf/pi_inf_sf to be 
            !! updated alongside of rho/gamma/pi_inf. Therefore, the
            !! versions of these variables appended with '_K' represent
            !! pointers that target the correct variable. At the end, 
            !! rho/gamma/pi_inf are updated for post process.
#ifdef MFC_POST_PROCESS
        rho_K => rho_sf(k, l)
        gamma_K =>  gamma_sf(k, l)
        pi_inf_K => pi_inf_sf(k, l)
#else
        rho_K  => rho
        gamma_K => gamma
        pi_inf_K => pi_inf
#endif

        ! Computing the density, the specific heat ratio function and the
        ! liquid stiffness function, respectively

        do i = 1, num_fluids
            alpha_rho_K(i) = q_vf(i)%sf(k, l)
            alpha_K(i) = q_vf(advxb + i - 1)%sf(k, l)
        end do


        ! Calculating the density, the specific heat ratio function and the
        ! liquid stiffness function, respectively, from the species analogs
        rho_K = 0d0; gamma_K = 0d0; pi_inf_K = 0d0

        do i = 1, num_fluids
            rho_K = rho_K + alpha_rho_K(i)
            gamma_K = gamma_K + alpha_K(i)*gammas(i)
            pi_inf_K = pi_inf_K + alpha_K(i)*pi_infs(i)
        end do

#ifdef MFC_SIMULATION
        ! Computing the shear and bulk Reynolds numbers from species analogs
        do i = 1, 2

            Re_K(i) = dflt_real; if (Re_size(i) > 0) Re_K(i) = 0d0

            do j = 1, Re_size(i)
                Re_K(i) = alpha_K(Re_idx(i, j))/fluid_pp(Re_idx(i, j))%Re(i) &
                          + Re_K(i)
            end do

            Re_K(i) = 1d0/max(Re_K(i), sgm_eps)

        end do
#endif

#ifdef MFC_POST_PROCESS
        rho = rho_K
        gamma = gamma_K
        pi_inf = pi_inf_K
#endif

    end subroutine s_convert_species_to_mixture_variables ! ----------------

    subroutine s_convert_species_to_mixture_variables_acc(rho_K, &
                                                          gamma_K, pi_inf_K, &
                                                          alpha_K, alpha_rho_K, Re_K )
!$acc routine seq

        real(kind(0d0)), intent(OUT) :: rho_K, gamma_K, pi_inf_K
        real(kind(0d0)), dimension(num_fluids), intent(INOUT) :: alpha_rho_K, alpha_K 
        real(kind(0d0)), dimension(2), intent(OUT) :: Re_K

        integer :: i, j

#ifdef MFC_SIMULATION
        rho_K = 0d0
        gamma_K = 0d0
        pi_inf_K = 0d0

        do i = 1, num_fluids
            rho_K = rho_K + alpha_rho_K(i)
            gamma_K = gamma_K + alpha_K(i)*gammas(i)
            pi_inf_K = pi_inf_K + alpha_K(i)*pi_infs(i)
        end do

        if (any(Re_size > 0)) then

            do i = 1, 2
                Re_K(i) = dflt_real

                if (Re_size(i) > 0) Re_K(i) = 0d0

                do j = 1, Re_size(i)
                    Re_K(i) = alpha_K(Re_idx(i, j))/Res(i, j) &
                              + Re_K(i)
                end do

                Re_K(i) = 1d0/max(Re_K(i), sgm_eps)

            end do
        end if
#endif

    end subroutine s_convert_species_to_mixture_variables_acc ! ----------------


    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_variables_conversion_module() ! ----------------

        integer :: i, j
!$acc update device(momxb, momxe, advxb, advxe, contxb, contxe)

#ifdef MFC_PRE_PROCESS
        ixb = 0; iyb = 0
        ixe = m; iye = n
#else
        ixb = -buff_size
        ixe = m - ixb

        iyb = 0; iye = 0
        if (n > 0) then
            iyb = -buff_size; iye = n - iyb
       end if
#endif

        !$acc update device(ixb, ixe, iyb, iye)

        @:ALLOCATE(gammas (1:num_fluids))
        @:ALLOCATE(pi_infs(1:num_fluids))

        do i = 1, num_fluids
            gammas(i)  = fluid_pp(i)%gamma
            pi_infs(i) = fluid_pp(i)%pi_inf
        end do
        !$acc update device(gammas, pi_infs)

#ifdef MFC_SIMULATION

        if (any(Re_size > 0)) then
            @:ALLOCATE(Res(1:2, 1:maxval(Re_size)))
            
            do i = 1, 2
                do j = 1, Re_size(i)
                    Res(i, j) = fluid_pp(Re_idx(i, j))%Re(i)
                end do
            end do
            
            !$acc update device(Res, Re_idx, Re_size)
        end if
#endif


!$acc update device(dt, sys_size, gamma_idx, pi_inf_idx, E_idx, num_fluids, num_dims, weno_eps)

#ifdef MFC_POST_PROCESS
        ! Allocating the density, the specific heat ratio function and the
        ! liquid stiffness function, respectively

        ! Simulation is at least 2D
        if (n > 0) then
        ! Simulation is 2D
            allocate (rho_sf(-buff_size:m + buff_size, &
                             -buff_size:n + buff_size))
            allocate (gamma_sf(-buff_size:m + buff_size, &
                               -buff_size:n + buff_size ))
            allocate (pi_inf_sf(-buff_size:m + buff_size, &
                                -buff_size:n + buff_size ))
        ! Simulation is 1D
        else
            allocate (rho_sf(-buff_size:m + buff_size, &
                             0:0 ))
            allocate (gamma_sf(-buff_size:m + buff_size, &
                               0:0 ))
            allocate (pi_inf_sf(-buff_size:m + buff_size, &
                                0:0 ))
        end if
#endif

        ! Volume fraction model
        s_convert_to_mixture_variables => s_convert_species_to_mixture_variables

    end subroutine s_initialize_variables_conversion_module ! --------------

    !> The following procedure handles the conversion between
        !!      the conservative variables and the primitive variables.
        !! @param qK_cons_vf Conservative variables
        !! @param qK_prim_vf Primitive variables
    subroutine s_convert_conservative_to_primitive_variables(qK_cons_vf, &
                                                             qK_prim_vf )

        type(scalar_field), dimension(sys_size), intent(IN) :: qK_cons_vf
        type(scalar_field), dimension(sys_size), intent(INOUT) :: qK_prim_vf

        real(kind(0d0)), dimension(num_fluids) :: alpha_K, alpha_rho_K
        real(kind(0d0)), dimension(2) :: Re_K
        real(kind(0d0)) :: rho_K, gamma_K, pi_inf_K, dyn_pres_K
        real(kind(0d0)) :: pres

        integer :: i, j, k, l !< Generic loop iterators
        
        !$acc parallel loop collapse(3) gang vector default(present) private(alpha_K, alpha_rho_K, Re_K, rho_K, gamma_K, pi_inf_K, dyn_pres_K)
        do k = iyb, iye
            do j = ixb, ixe
                dyn_pres_K = 0d0
                
                !$acc loop seq
                do i = 1, num_fluids
                    alpha_rho_K(i) = qK_cons_vf(i)%sf(j, k)
                    alpha_K(i) = qK_cons_vf(advxb + i - 1)%sf(j, k)
                end do

                do i = 1, contxe
                    qK_prim_vf(i)%sf(j, k) = qK_cons_vf(i)%sf(j, k)
                end do

#ifdef MFC_SIMULATION
                ! If in simulation, use acc mixture subroutines
                call s_convert_species_to_mixture_variables_acc(rho_K, gamma_K, pi_inf_K, &
                                                                    alpha_K, alpha_rho_K, Re_K)
                                                                ! , k, l)
#else
                ! If pre-processing, use non acc mixture subroutines
                call s_convert_to_mixture_variables(qK_cons_vf, j, k, &
                                                    rho_K, gamma_K, pi_inf_K)
#endif

#ifdef MFC_SIMULATION
                rho_K = max(rho_K, sgm_eps)
#endif

                !$acc loop seq
                do i = momxb, momxe
                    qK_prim_vf(i)%sf(j, k) = qK_cons_vf(i)%sf(j, k) &
                                                /rho_K
                    dyn_pres_K = dyn_pres_K + 5d-1*qK_cons_vf(i)%sf(j, k) &
                                 *qK_prim_vf(i)%sf(j, k)
                end do
                call s_compute_pressure(qK_cons_vf(E_idx)%sf(j, k), &
                                        dyn_pres_K, pi_inf_K, gamma_K, pres)

                qK_prim_vf(E_idx)%sf(j, k) = pres

                do i = advxb, advxe
                    qK_prim_vf(i)%sf(j, k) = qK_cons_vf(i)%sf(j, k)
                end do
            end do
        end do
        !$acc end parallel loop

    end subroutine s_convert_conservative_to_primitive_variables ! ---------

    !>  The following procedure handles the conversion between
        !!      the primitive variables and the conservative variables.
        !!  @param qK_prim_vf Primitive variables
        !!  @param qK_cons_vf Conservative variables
        !!  @param gm_alphaK_vf Gradient magnitude of the volume fractions
        !!  @param ix Index bounds in the first coordinate direction
        !!  @param iy Index bounds in the second coordinate direction
        !!  @param iz Index bounds in the third coordinate direction
    subroutine s_convert_primitive_to_conservative_variables(q_prim_vf, &
                                                             q_cons_vf)

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_prim_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: q_cons_vf

        ! Density, specific heat ratio function, liquid stiffness function
        ! and dynamic pressure, as defined in the incompressible flow sense,
        ! respectively
        real(kind(0d0)) :: rho
        real(kind(0d0)) :: gamma
        real(kind(0d0)) :: pi_inf
        real(kind(0d0)) :: dyn_pres

        integer :: i, j, k

#ifndef MFC_SIMULATION
        ! Converting the primitive variables to the conservative variables
        do k = 0, n
            do j = 0, m

                ! Obtaining the density, specific heat ratio function
                ! and the liquid stiffness function, respectively
                call s_convert_to_mixture_variables(q_prim_vf, j, k,  &
                                                    rho, gamma, pi_inf)

                ! Transferring the continuity equation(s) variable(s)
                do i = 1, contxe
                    q_cons_vf(i)%sf(j, k) = q_prim_vf(i)%sf(j, k)
                end do

                ! Zeroing out the dynamic pressure since it is computed
                ! iteratively by cycling through the velocity equations
                dyn_pres = 0d0

                ! Computing momenta and dynamic pressure from velocity
                do i = momxb, momxe
                    q_cons_vf(i)%sf(j, k) = rho*q_prim_vf(i)%sf(j, k)
                    dyn_pres = dyn_pres + q_cons_vf(i)%sf(j, k)* &
                               q_prim_vf(i)%sf(j, k)/2d0
                end do

                ! Computing the energy from the pressure
                ! E = Gamma*P + \rho u u /2 + \pi_inf
                q_cons_vf(E_idx)%sf(j, k) = &
                    gamma*q_prim_vf(E_idx)%sf(j, k) + dyn_pres + pi_inf


                ! Transferring the advection equation(s) variable(s)
                do i = adv_idx%beg, adv_idx%end
                    q_cons_vf(i)%sf(j, k) = q_prim_vf(i)%sf(j, k)
                end do
            end do
        end do

#else
        if (proc_rank == 0) then
            print '(A)', 'Conversion from primitive to '// &
                'conservative variables not '// &
                'implemented. Exiting ...'
            call s_mpi_abort()
        end if
#endif

    end subroutine s_convert_primitive_to_conservative_variables ! ---------

    !>  The following subroutine handles the conversion between
        !!      the primitive variables and the Eulerian flux variables.
        !!  @param qK_prim_vf Primitive variables
        !!  @param FK_vf Flux variables
        !!  @param FK_src_vf Flux source variables
        !!  @param ix Index bounds in the first coordinate direction
        !!  @param iy Index bounds in the second coordinate direction
        !!  @param iz Index bounds in the third coordinate direction
    subroutine s_convert_primitive_to_flux_variables(qK_prim_vf, & ! ------
                                                     FK_vf, &
                                                     FK_src_vf, &
                                                     is1, is2, s2b)

        integer :: s2b
        real(kind(0d0)), dimension(0:, s2b:, 1:), intent(IN) :: qK_prim_vf
        real(kind(0d0)), dimension(0:, s2b:, 1:), intent(INOUT) :: FK_vf
        real(kind(0d0)), dimension(0:, s2b:, advxb:), intent(INOUT) :: FK_src_vf

        type(int_bounds_info), intent(IN) :: is1, is2

        ! Partial densities, density, velocity, pressure, energy, advection
        ! variables, the specific heat ratio and liquid stiffness functions,
        ! the shear and volume Reynolds numbers and the Weber numbers
        real(kind(0d0)), dimension(num_fluids) :: alpha_rho_K
        real(kind(0d0)), dimension(num_fluids) :: alpha_K
        real(kind(0d0)) :: rho_K
        real(kind(0d0)), dimension(num_dims) :: vel_K
        real(kind(0d0)) :: vel_K_sum
        real(kind(0d0)) :: pres_K
        real(kind(0d0)) :: E_K
        real(kind(0d0)) :: gamma_K
        real(kind(0d0)) :: pi_inf_K
        real(kind(0d0)), dimension(2) :: Re_K

        integer :: i, j, k !< Generic loop iterators

        is1b = is1%beg; is1e = is1%end
        is2b = is2%beg; is2e = is2%end

        !$acc update device(is1b, is2b, is1e, is2e)

        ! Computing the flux variables from the primitive variables, without
        ! accounting for the contribution of either viscosity or capillarity
#ifdef MFC_SIMULATION
!$acc parallel loop collapse(3) gang vector default(present) private(alpha_rho_K, vel_K, alpha_K, Re_K)
        do k = is2b, is2e
            do j = is1b, is1e

!$acc loop seq
                do i = 1, contxe
                    alpha_rho_K(i) = qK_prim_vf(j, k, i)
                end do

!$acc loop seq
                do i = advxb, advxe
                    alpha_K(i - E_idx) = qK_prim_vf(j, k, i)
                end do
!$acc loop seq
                do i = 1, num_dims
                    vel_K(i) = qK_prim_vf(j, k, contxe + i)
                end do

                vel_K_sum = 0d0
!$acc loop seq
                do i = 1, num_dims
                    vel_K_sum = vel_K_sum + vel_K(i)**2d0
                end do

                pres_K = qK_prim_vf(j, k, E_idx)
                call s_convert_species_to_mixture_variables_acc(rho_K, gamma_K, pi_inf_K, &
                                                                alpha_K, alpha_rho_K, Re_K)

                ! Computing the energy from the pressure
                E_K = gamma_K*pres_K + pi_inf_K + 5d-1*rho_K*vel_K_sum

                ! mass flux, this should be \alpha_i \rho_i u_i
!$acc loop seq
                do i = 1, contxe
                    FK_vf(j, k, i) = alpha_rho_K(i)*vel_K(dir_idx(1))
                end do

!$acc loop seq
                do i = 1, num_dims
                    FK_vf(j, k, contxe + dir_idx(i)) = &
                        rho_K*vel_K(dir_idx(1)) &
                        *vel_K(dir_idx(i)) &
                        + pres_K*dir_flg(dir_idx(i))
                end do

                ! energy flux, u(E+p)
                FK_vf(j, k, E_idx) = vel_K(dir_idx(1))*(E_K + pres_K)

                !$acc loop seq
                do i = advxb, advxe
                    FK_src_vf(j, k, i) = vel_K(dir_idx(1))
                end do
            end do
        end do
#endif

    end subroutine s_convert_primitive_to_flux_variables ! -----------------

    subroutine s_finalize_variables_conversion_module() ! ------------------

        ! Deallocating the density, the specific heat ratio function and the
        ! liquid stiffness function
#ifdef MFC_POST_PROCESS
        deallocate(rho_sf, gamma_sf, pi_inf_sf)
#endif

        @:DEALLOCATE(gammas, pi_infs)
        
        ! Nullifying the procedure pointer to the subroutine transfering/
        ! computing the mixture/species variables to the mixture variables
        s_convert_to_mixture_variables => null()

    end subroutine s_finalize_variables_conversion_module ! ----------------

end module m_variables_conversion
