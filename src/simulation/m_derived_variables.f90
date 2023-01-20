!>
!! @file m_derived_variables.f90
!! @brief Contains module m_derived_variables

!> @brief This module features subroutines that allow for the derivation of
!!              numerous flow variables from the conservative and primitive ones.
!!              Currently, the available derived variables include the unadvected
!!              volume fraction, specific heat ratio, liquid stiffness, speed of
!!              sound, vorticity and the numerical Schlieren function.
module m_derived_variables

    ! Dependencies =============================================================
    use m_derived_types         !< Definitions of the derived types

    use m_global_parameters     !< Global parameters for the code

    use m_mpi_proxy             !< Message passing interface (MPI) module proxy

    use m_data_output           !< Data output module

    use m_time_steppers         !< Time-stepping algorithms
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_derived_variables_module, &
 s_initialize_derived_variables, &
 s_compute_derived_variables, &
 s_finalize_derived_variables_module

    !> @name Finite-difference coefficients
    !! Finite-difference (fd) coefficients in x-, y- and z-coordinate directions.
    !! Note that because sufficient boundary information is available for all the
    !! active coordinate directions, the centered family of the finite-difference
    !! schemes is used.
    !> @{
    real(kind(0d0)), public, allocatable, dimension(:, :) :: fd_coeff_x
    real(kind(0d0)), public, allocatable, dimension(:, :) :: fd_coeff_y
    !> @}

contains

    !>  Computation of parameters, allocation procedures, and/or
        !!      any other tasks needed to properly setup the module
    subroutine s_initialize_derived_variables_module() ! ----------------------

        ! Allocating the variables which will store the coefficients of the
        ! centered family of finite-difference schemes. Note that sufficient
        ! space is allocated so that the coefficients up to any chosen order
        ! of accuracy may be bookkept. However, if higher than fourth-order
        ! accuracy coefficients are wanted, the formulae required to compute
        ! these coefficients will have to be implemented in the subroutine
        ! s_compute_finite_difference_coefficients.

        ! Allocating centered finite-difference coefficients
        if (probe_wrt) then
            allocate (fd_coeff_x(-fd_number:fd_number, 0:m))
            if (n > 0) then
                allocate (fd_coeff_y(-fd_number:fd_number, 0:n))
            end if
        end if

    end subroutine s_initialize_derived_variables_module ! --------------------

    !> Allocate and open derived variables. Computing FD coefficients.
    subroutine s_initialize_derived_variables() ! -----------------------------

        ! Opening and writing header of CoM and flow probe files
        if (proc_rank == 0) then
            if (probe_wrt) then
                call s_open_probe_files()
            end if
        end if

        ! Computing centered finite difference coefficients
        if (probe_wrt) then
            call s_compute_finite_difference_coefficients(m, x_cc, fd_coeff_x)
            if (n > 0) then
                call s_compute_finite_difference_coefficients(n, y_cc, fd_coeff_y)
            end if
        end if

    end subroutine s_initialize_derived_variables ! -----------------------------

    !> Writes coherent body information, communication files, and probes.
        !!  @param t_step Current time-step
    subroutine s_compute_derived_variables(t_step) ! -----------------------

        integer, intent(IN) :: t_step

        integer :: i, j !< Generic loop iterators

        if (probe_wrt) then
            call s_write_probe_files(t_step, q_cons_ts(1)%vf)
        end if

    end subroutine s_compute_derived_variables ! ---------------------------

    !>  The purpose of this subroutine is to compute the finite-
        !!      difference coefficients for the centered schemes utilized
        !!      in computations of first order spatial derivatives in the
        !!      s-coordinate direction. The s-coordinate direction refers
        !!      to the x-, y- or z-coordinate direction, depending on the
        !!      subroutine's inputs. Note that coefficients of up to 4th
        !!      order accuracy are available.
        !!  @param q Number of cells in the s-coordinate direction
        !!  @param s_cc Locations of the cell-centers in the s-coordinate direction
        !!  @param fd_coeff_s Finite-diff. coefficients in the s-coordinate direction
    subroutine s_compute_finite_difference_coefficients(q, s_cc, fd_coeff_s)

        integer, intent(IN) :: q

        real(kind(0d0)), &
            dimension(-buff_size:q + buff_size), &
            intent(IN) :: s_cc

        real(kind(0d0)), &
            dimension(-fd_number:fd_number, 0:q), &
            intent(INOUT) :: fd_coeff_s

        integer :: i !< Generic loop iterator

        ! Computing the 1st order finite-difference coefficients
        if (fd_order == 1) then
            do i = 0, q
                fd_coeff_s(-1, i) = 0d0
                fd_coeff_s(0, i) = -1d0/(s_cc(i + 1) - s_cc(i))
                fd_coeff_s(1, i) = -fd_coeff_s(0, i)
            end do

            ! Computing the 2nd order finite-difference coefficients
        elseif (fd_order == 2) then
            do i = 0, q
                fd_coeff_s(-1, i) = -1d0/(s_cc(i + 1) - s_cc(i - 1))
                fd_coeff_s(0, i) = 0d0
                fd_coeff_s(1, i) = -fd_coeff_s(-1, i)
            end do

            ! Computing the 4th order finite-difference coefficients
        else
            do i = 0, q
                fd_coeff_s(-2, i) = 1d0/(s_cc(i - 2) - 8d0*s_cc(i - 1) - s_cc(i + 2) + 8d0*s_cc(i + 1))
                fd_coeff_s(-1, i) = -8d0*fd_coeff_s(-2, i)
                fd_coeff_s(0, i) = 0d0
                fd_coeff_s(1, i) = -fd_coeff_s(-1, i)
                fd_coeff_s(2, i) = -fd_coeff_s(-2, i)
            end do

        end if

    end subroutine s_compute_finite_difference_coefficients ! --------------



    !> Deallocation procedures for the module
    subroutine s_finalize_derived_variables_module() ! -------------------

        ! Closing CoM and flow probe files
        if (proc_rank == 0) then
            if (probe_wrt) then
                call s_close_probe_files()
            end if
        end if

        ! Deallocating the variables that might have been used to bookkeep
        ! the finite-difference coefficients in the x-, y- and z-directions
        if (allocated(fd_coeff_x)) deallocate (fd_coeff_x)
        if (allocated(fd_coeff_y)) deallocate (fd_coeff_y)

    end subroutine s_finalize_derived_variables_module ! -----------------

end module m_derived_variables
