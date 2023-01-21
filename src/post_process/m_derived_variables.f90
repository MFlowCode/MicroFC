!>
!! @file m_derived_variables.f90
!! @brief Contains module m_derived_variables

!> @brief This module features subroutines that allow for the derivation of
!!      numerous flow variables from the conservative and primitive ones.
!!      Currently, the available derived variables include the unadvected
!!      volume fraction, specific heat ratio, liquid stiffness, speed of
!!      sound, vorticity and the numerical Schlieren function.
module m_derived_variables

    ! Dependencies =============================================================
    use m_derived_types         !< Definitions of the derived types

    use m_global_parameters     !< Global parameters for the code

    use m_mpi_proxy             !< Message passing interface (MPI) module proxy
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_derived_variables_module, &
                         s_compute_finite_difference_coefficients, &
                         s_derive_vorticity_component, &
                         s_derive_numerical_schlieren_function, &
                         s_finalize_derived_variables_module

    real(kind(0d0)), allocatable, dimension(:, :) :: gm_rho_sf !<
    !! Gradient magnitude (gm) of the density for each cell of the computational
    !! sub-domain. This variable is employed in the calculation of the numerical
    !! Schlieren function.

    !> @name Finite-difference (fd) coefficients in x-, y- and z-coordinate directions.
    !! Note that because sufficient boundary information is available for all the
    !! active coordinate directions, the centered family of the finite-difference
    !! schemes is used.
    !> @{
    real(kind(0d0)), allocatable, dimension(:, :), public :: fd_coeff_x
    real(kind(0d0)), allocatable, dimension(:, :), public :: fd_coeff_y
    !> @}

    integer, private :: flg  !<
    !! Flagging (flg) variable used to annotate the dimensionality of the dataset
    !! that is undergoing the post-process. A flag value of 1 indicates that the
    !! dataset is 3D, while a flag value of 0 indicates that it is not. This flg
    !! variable is necessary to avoid cycling through the third dimension of the
    !! flow variable(s) when the simulation is not 3D and the size of the buffer
    !! is non-zero. Note that a similar procedure does not have to be applied to
    !! the second dimension since in 1D, the buffer size is always zero.

contains

    !>  Computation of parameters, allocation procedures, and/or
        !!      any other tasks needed to properly setup the module
    subroutine s_initialize_derived_variables_module() ! ----------------------

        ! Allocating the gradient magnitude of the density variable provided
        ! that numerical Schlieren function is outputted during post-process
        if (schlieren_wrt) then
            allocate (gm_rho_sf(-offset_x%beg:m + offset_x%end, &
                                -offset_y%beg:n + offset_y%end ))
        end if

        ! Allocating the variables which will store the coefficients of the
        ! centered family of finite-difference schemes. Note that sufficient
        ! space is allocated so that the coefficients up to any chosen order
        ! of accuracy may be bookkept. However, if higher than fourth-order
        ! accuracy coefficients are wanted, the formulae required to compute
        ! these coefficients will have to be implemented in the subroutine
        ! s_compute_finite_difference_coefficients.

        ! Allocating centered finite-difference coefficients in x-direction
        if (omega_wrt(2) .or. schlieren_wrt) then
            allocate (fd_coeff_x(-fd_number:fd_number, &
                                 -offset_x%beg:m + offset_x%end))
        end if

        ! Allocating centered finite-difference coefficients in y-direction
        if (omega_wrt(1) .or. (n > 0 .and. schlieren_wrt)) then
            allocate (fd_coeff_y(-fd_number:fd_number, &
                                 -offset_y%beg:n + offset_y%end))
        end if

        ! Annotating the dimensionality of the dataset undergoing the post-
        ! process. A flag value of 1 indicates that the dataset is 3D, while
        ! a flag value of 0 indicates that it is not.
        flg = 0

    end subroutine s_initialize_derived_variables_module ! --------------------

    !> @name The purpose of this subroutine is to compute the finite-
        !!      difference coefficients for the centered schemes utilized
        !!      in computations of first order spatial derivatives in the
        !!      s-coordinate direction. The s-coordinate direction refers
        !!      to the x-, y- or z-coordinate direction, depending on the
        !!      subroutine's inputs. Note that coefficients of up to 4th
        !!      order accuracy are available.
        !!  @param q Number of cells in the s-coordinate direction
        !!  @param offset_s  Size of the ghost zone layer in the s-coordinate direction
        !!  @param s_cc Locations of the cell-centers in the s-coordinate direction
        !!  @param fd_coeff_s Finite-diff. coefficients in the s-coordinate direction
    subroutine s_compute_finite_difference_coefficients(q, offset_s, &
                                                        s_cc, fd_coeff_s)

        integer, intent(IN) :: q
        type(int_bounds_info), intent(IN) :: offset_s

        real(kind(0d0)), &
            dimension(-buff_size:q + buff_size), &
            intent(IN) :: s_cc

        real(kind(0d0)), &
            dimension(-fd_number:fd_number, -offset_s%beg:q + offset_s%end), &
            intent(INOUT) :: fd_coeff_s

        integer :: i !< Generic loop iterator

        ! Computing the 1st order finite-difference coefficients
        if (fd_order == 1) then
            do i = -offset_s%beg, q + offset_s%end
                fd_coeff_s(-1, i) = 0d0
                fd_coeff_s(0, i) = -1d0/(s_cc(i + 1) - s_cc(i))
                fd_coeff_s(1, i) = -fd_coeff_s(0, i)
            end do

            ! Computing the 2nd order finite-difference coefficients
        elseif (fd_order == 2) then
            do i = -offset_s%beg, q + offset_s%end
                fd_coeff_s(-1, i) = -1d0/(s_cc(i + 1) - s_cc(i - 1))
                fd_coeff_s(0, i) = 0d0
                fd_coeff_s(1, i) = -fd_coeff_s(-1, i)
            end do

            ! Computing the 4th order finite-difference coefficients
        else
            do i = -offset_s%beg, q + offset_s%end
                fd_coeff_s(-2, i) = 1d0/(s_cc(i - 2) - 8d0*s_cc(i - 1) &
                                         - s_cc(i + 2) + 8d0*s_cc(i + 1))
                fd_coeff_s(-1, i) = -8d0*fd_coeff_s(-2, i)
                fd_coeff_s(0, i) = 0d0
                fd_coeff_s(1, i) = -fd_coeff_s(-1, i)
                fd_coeff_s(2, i) = -fd_coeff_s(-2, i)
            end do

        end if

    end subroutine s_compute_finite_difference_coefficients ! --------------



    !>  This subroutine receives as inputs the indicator of the
        !!      component of the vorticity that should be outputted and
        !!      the primitive variables. From those inputs, it proceeds
        !!      to calculate values of the desired vorticity component,
        !!      which are subsequently stored in derived flow quantity
        !!      storage variable, q_sf.
        !!  @param i Vorticity component indicator
        !!  @param q_prim_vf Primitive variables
        !!  @param q_sf Vorticity component
    subroutine s_derive_vorticity_component(i, q_prim_vf, q_sf) ! ----------

        integer, intent(IN) :: i

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_prim_vf

        real(kind(0d0)), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end), &
                      intent(INOUT) :: q_sf

        integer :: j, k, l, r !< Generic loop iterators

        ! Computing the vorticity component in the x-coordinate direction
        if (i == 1) then
                do k = -offset_y%beg, n + offset_y%end
                    do j = -offset_x%beg, m + offset_x%end

                        q_sf(j, k) = 0d0

                        do r = -fd_number, fd_number
                            q_sf(j, k) = q_sf(j, k) + fd_coeff_y(r, k)* &
                                q_prim_vf(mom_idx%end)%sf(j, r + k) 
                        end do

                    end do
                end do

        ! Computing the vorticity component in the y-coordinate direction
        elseif (i == 2) then
            do k = -offset_y%beg, n + offset_y%end
                do j = -offset_x%beg, m + offset_x%end

                    q_sf(j, k) = 0d0

                    do r = -fd_number, fd_number
                        q_sf(j, k) = &
                            q_sf(j, k) - fd_coeff_x(r, j)* &
                            q_prim_vf(mom_idx%end)%sf(r + j, k)
                    end do

                end do
            end do
        end if

    end subroutine s_derive_vorticity_component ! --------------------------

    !>  This subroutine gets as inputs the conservative variables
        !!      and density. From those inputs, it proceeds to calculate
        !!      the values of the numerical Schlieren function, which are
        !!      subsequently stored in the derived flow quantity storage
        !!      variable, q_sf.
        !!  @param q_cons_vf Conservative variables
        !!  @param rho_sf Density
        !!  @param q_sf Numerical Schlieren function
    subroutine s_derive_numerical_schlieren_function(q_cons_vf, rho_sf, q_sf)

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_cons_vf

        real(kind(0d0)), &
            dimension(-buff_size:m + buff_size, &
                      -buff_size:n + buff_size), &
            intent(IN) :: rho_sf

        real(kind(0d0)), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end), &
            intent(INOUT) :: q_sf

        real(kind(0d0)) :: drho_dx, drho_dy, drho_dz !<
            !! Spatial derivatives of the density in the x-, y- and z-directions

        real(kind(0d0)), dimension(2) :: gm_rho_max !<
            !! Maximum value of the gradient magnitude (gm) of the density field
            !! in entire computational domain and not just the local sub-domain.
            !! The first position in the variable contains the maximum value and
            !! the second contains the rank of the processor on which it occured.

        real(kind(0d0)) :: alpha_unadv !< Unadvected volume fraction

        integer :: i, j, k, l !< Generic loop iterators

        ! Computing Gradient Magnitude of Density ==========================

        ! Contributions from the x- and y-coordinate directions
            do k = -offset_y%beg, n + offset_y%end
                do j = -offset_x%beg, m + offset_x%end

                    drho_dx = 0d0
                    drho_dy = 0d0

                    do i = -fd_number, fd_number
                        drho_dx = drho_dx + fd_coeff_x(i, j)*rho_sf(i + j, k)
                        drho_dy = drho_dy + fd_coeff_y(i, k)*rho_sf(j, i + k)
                    end do

                    gm_rho_sf(j, k) = drho_dx*drho_dx + drho_dy*drho_dy

                end do
            end do

        ! Up until now, only the dot product of the gradient of the density
        ! field has been calculated and stored in the gradient magnitude of
        ! density variable. So now we proceed to take the square-root as to
        ! complete the desired calculation.
        gm_rho_sf = sqrt(gm_rho_sf)

        ! ==================================================================

        ! Determining the local maximum of the gradient magnitude of density
        ! and bookkeeping the result, along with rank of the local processor
        gm_rho_max = (/maxval(gm_rho_sf), real(proc_rank, kind(0d0))/)

        ! Comparing the local maximum gradient magnitude of the density on
        ! this processor to the those computed on the remaining processors.
        ! This allows for the global maximum to be computed and the rank of
        ! the processor on which it has occured to be recorded.
        if (num_procs > 1) call s_mpi_reduce_maxloc(gm_rho_max)

        ! Computing Numerical Schlieren Function ===========================

        ! The form of the numerical Schlieren function depends on the choice
        ! of the multicomponent flow model. For the gamma/pi_inf model, the
        ! exponential of the negative, normalized, gradient magnitude of the
        ! density is computed. For the volume fraction model, the amplitude
        ! of the exponential's inside is also modulated with respect to the
        ! identity of the fluid in which the function is evaluated. For more
        ! information, refer to Marquina and Mulet (2003).
        do k = -offset_y%beg, n + offset_y%end
            do j = -offset_x%beg, m + offset_x%end

                q_sf(j, k) = 0d0

                do i = 1, adv_idx%end - E_idx
                    q_sf(j, k) = &
                        q_sf(j, k) - schlieren_alpha(i)* &
                        q_cons_vf(i + E_idx)%sf(j, k)* &
                        gm_rho_sf(j, k)/gm_rho_max(1)
                end do
            end do
        end do

        ! Up until now, only the inside of the exponential of the numerical
        ! Schlieren function has been evaluated and stored. Then, to finish
        ! the computation, the exponential of the inside quantity is taken.
        q_sf = exp(q_sf)

        ! ==================================================================

    end subroutine s_derive_numerical_schlieren_function ! -----------------

    !>  Deallocation procedures for the module
    subroutine s_finalize_derived_variables_module() ! -------------------

        ! Deallocating the variable containing the gradient magnitude of the
        ! density field provided that the numerical Schlieren function was
        ! was outputted during the post-process
        if (schlieren_wrt) deallocate (gm_rho_sf)

        ! Deallocating the variables that might have been used to bookkeep
        ! the finite-difference coefficients in the x-, y- and z-directions
        if (allocated(fd_coeff_x)) deallocate (fd_coeff_x)
        if (allocated(fd_coeff_y)) deallocate (fd_coeff_y)

    end subroutine s_finalize_derived_variables_module ! -----------------

end module m_derived_variables
