module m_assign_patches

    ! Dependencies =============================================================
    use m_derived_types         ! Definitions of the derived types

    use m_global_parameters     ! Global parameters for the code

    use m_variables_conversion  ! Subroutines to change the state variables from
    ! one form to another

    ! ==========================================================================

    implicit none

    ! NOTE: The abstract interface allows for the declaration of a pointer to
    ! a procedure such that the choice of the model equations does not have to
    ! be queried every time the patch primitive variables are to be assigned in
    ! a cell in the computational domain.
    type(scalar_field), allocatable, dimension(:) :: q_prim_vf !< primitive variables
    type(scalar_field), allocatable, dimension(:) :: q_cons_vf !< conservative variables
    type(scalar_field) :: alf_sum

    real(kind(0d0)) :: x_centroid, y_centroid, z_centroid
    real(kind(0d0)) :: epsilon, beta
    integer :: smooth_patch_id

    real(kind(0d0)) :: eta !<
    !! In the case that smoothing of patch boundaries is enabled and the boundary
    !! between two adjacent patches is to be smeared out, this variable's purpose
    !! is to act as a pseudo volume fraction to indicate the contribution of each
    !! patch toward the composition of a cell's fluid state.

    integer, allocatable, dimension(:, :) :: patch_id_fp !<
    !! Bookkepping variable used to track the patch identities (id) associated
    !! with each of the cells in the computational domain. Note that only one
    !! patch identity may be associated with any one cell.

contains

    !>      This subroutine assigns the species primitive variables
        !!              of the patch designated by the patch_id, to the cell that
        !!              is designated by the indexes (j,k,l). In addition, the
        !!              variable bookkeeping the patch identities in the entire
        !!              domain is updated with the new assignment. Note that if
        !!              the smoothing of the patch's boundaries is employed, the
        !!              ensuing primitive variables in the cell will be a type of
        !!              combination of the current patch's primitive variables
        !!              with those of the smoothing patch. The specific details
        !!              of the combination may be found in Shyue's work (1998).
        !! @param patch_id the patch identifier
        !! @param j  the x-dir node index
        !! @param k  the y-dir node index
        !! @param l  the z-dir node index
    subroutine s_assign_patch_species_primitive_variables(patch_id, j, k)

        integer, intent(IN) :: patch_id
        integer, intent(IN) :: j, k

        real(kind(0d0)) :: rho
        real(kind(0d0)) :: gamma
        real(kind(0d0)) :: pi_inf
        real(kind(0d0)) :: orig_rho
        real(kind(0d0)) :: orig_gamma
        real(kind(0d0)) :: orig_pi_inf !<
            !! Density, the specific heat ratio function and the liquid stiffness
            !! function, respectively, obtained from the combination of primitive
            !! variables of the current and smoothing patches

        real(kind(0d0)), dimension(sys_size) :: orig_prim_vf !<
        ! Vector to hold original values of cell for smoothing purposes

        integer :: i !< generic loop iterator

        ! Transferring the identity of the smoothing patch
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id

        ! Transferring original primitive variables
        do i = 1, sys_size
            orig_prim_vf(i) = q_prim_vf(i)%sf(j, k)
        end do

        ! Computing Mixture Variables from Original Primitive Variables
        call s_convert_species_to_mixture_variables( &
            q_prim_vf, j, k, &
            orig_rho, &
            orig_gamma, &
            orig_pi_inf)

        ! Computing Mixture Variables of Current Patch =====================

        ! Partial densities
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf(j, k) = patch_icpp(patch_id)%alpha_rho(i)
        end do

        ! Volume fraction(s)
        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf(j, k) = patch_icpp(patch_id)%alpha(i - E_idx)
        end do

        ! Density and the specific heat ratio and liquid stiffness functions
        call s_convert_species_to_mixture_variables( &
            q_prim_vf, j, k, &
            patch_icpp(patch_id)%rho, &
            patch_icpp(patch_id)%gamma, &
            patch_icpp(patch_id)%pi_inf)

        ! ==================================================================

        ! Computing Mixture Variables of Smoothing Patch ===================

        ! Partial densities
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf(j, k) = &
                patch_icpp(smooth_patch_id)%alpha_rho(i)
        end do

        ! Volume fraction(s)
        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf(j, k) = &
                patch_icpp(smooth_patch_id)%alpha(i - E_idx)
        end do

        ! Density and the specific heat ratio and liquid stiffness functions
        call s_convert_species_to_mixture_variables( &
            q_prim_vf, j, k, &
            patch_icpp(smooth_patch_id)%rho, &
            patch_icpp(smooth_patch_id)%gamma, &
            patch_icpp(smooth_patch_id)%pi_inf)

        ! ==================================================================

        ! Partial densities
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf(j, k) = &
                eta*patch_icpp(patch_id)%alpha_rho(i) &
                + (1d0 - eta)*orig_prim_vf(i)
        end do
        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf(j, k) = &
                eta*patch_icpp(patch_id)%alpha(i - E_idx) &
                + (1d0 - eta)*orig_prim_vf(i)
        end do

        ! Density and the specific heat ratio and liquid stiffness functions
        call s_convert_species_to_mixture_variables(q_prim_vf, j, k, &
                                                    rho, gamma, pi_inf)

        ! Velocity
        do i = 1, E_idx - mom_idx%beg
            q_prim_vf(i + cont_idx%end)%sf(j, k) = &
                (eta*patch_icpp(patch_id)%vel(i) &
                 + (1d0 - eta)*orig_prim_vf(i + cont_idx%end))
        end do

        ! Pressure
        q_prim_vf(E_idx)%sf(j, k) = &
            (eta*patch_icpp(patch_id)%pres &
             + (1d0 - eta)*orig_prim_vf(E_idx))


        ! Updating the patch identities bookkeeping variable
        if (1d0 - eta < 1d-16) patch_id_fp(j, k) = patch_id

    end subroutine s_assign_patch_species_primitive_variables ! ------------

    !> Computation of parameters, allocation procedures, and/or
        !!              any other tasks needed to properly setup the module
    subroutine s_initialize_assign_patches_module() ! -------------------

        integer :: i !< generic loop iterator

        ! Allocating the primitive and conservative variables
        allocate (q_prim_vf(1:sys_size))
        allocate (q_cons_vf(1:sys_size))

        do i = 1, sys_size
            allocate (q_prim_vf(i)%sf(0:m, 0:n))
            allocate (q_cons_vf(i)%sf(0:m, 0:n))
        end do
        allocate (alf_sum%sf(0:m, 0:n))

        ! Allocating the patch identities bookkeeping variable
        allocate (patch_id_fp(0:m, 0:n))

        ! Setting default values for conservative and primitive variables so
        ! that in the case that the initial condition is wrongly laid out on
        ! the grid the simulation component will catch the problem on start-
        ! up. The conservative variables do not need to be similarly treated
        ! since they are computed directly from the primitive variables.
        do i = 1, sys_size
            q_cons_vf(i)%sf = dflt_real
            q_prim_vf(i)%sf = dflt_real
        end do

        ! Setting default values for patch identities bookkeeping variable.
        ! This is necessary to avoid any confusion in the assessment of the
        ! extent of application that the overwrite permissions give a patch
        ! when it is being applied in the domain.
        patch_id_fp = 0

    end subroutine s_initialize_assign_patches_module ! -----------------


    !>  Deallocation procedures for the module
    subroutine s_finalize_assign_patches_module() ! ---------------------

        integer :: i !< Generic loop iterator

        ! Dellocating the primitive and conservative variables
        do i = 1, sys_size
            deallocate (q_prim_vf(i)%sf)
            deallocate (q_cons_vf(i)%sf)
        end do

        deallocate (q_prim_vf)
        deallocate (q_cons_vf)

        ! Deallocating the patch identities bookkeeping variable
        deallocate (patch_id_fp)

    end subroutine s_finalize_assign_patches_module ! -------------------

end module m_assign_patches
