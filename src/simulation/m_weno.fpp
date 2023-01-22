!>
!! @file m_weno.f90
!! @brief Contains module m_weno

#:include 'macros.fpp'

!> @brief  Weighted essentially non-oscillatory (WENO) reconstruction scheme
!!              that is supplemented with monotonicity preserving bounds (MPWENO)
!!              and a mapping function that boosts the accuracy of the non-linear
!!              weights (WENOM). MPWENO, see Balsara and Shu (2000), prevents the
!!              reconstructed values to lay outside the range set by the stencil,
!!              while WENOM, see Henrick et al. (2005), recovers the formal order
!!              of accuracy of the reconstruction at critical points. Please note
!!              that the basic WENO approach is implemented according to the work
!!              of Jiang and Shu (1996).
module m_weno
    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_variables_conversion !< State variables type conversion procedures

#ifdef _OPENACC
    use openacc
#endif

    use m_mpi_proxy
    ! ==========================================================================

    !implicit none

    private; public :: s_initialize_weno_module, s_initialize_weno, s_finalize_weno_module, s_weno

    !> @name The cell-average variables that will be WENO-reconstructed. Formerly, they
    !! are stored in v_vf. However, they are transferred to v_rs_wsL and v_rs_wsR
    !! as to be reshaped (RS) and/or characteristically decomposed. The reshaping
    !! allows the WENO procedure to be independent of the coordinate direction of
    !! the reconstruction. Lastly, notice that the left (L) and right (R) results
    !! of the characteristic decomposition are stored in custom-constructed WENO-
    !! stencils (WS) that are annexed to each position of a given scalar field.
    !> @{
        real(kind(0d0)), allocatable, dimension(:, :, :) :: v_rs_ws_x, v_rs_ws_y
    !> @}

    ! WENO Coefficients ========================================================

    !> @name Polynomial coefficients at the left and right cell-boundaries (CB) and at
    !! the left and right quadrature points (QP), in the x-, y- and z-directions.
    !! Note that the first dimension of the array identifies the polynomial, the
    !! second dimension identifies the position of its coefficients and the last
    !! dimension denotes the cell-location in the relevant coordinate direction.
    !> @{
    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: poly_coef_cbL_x
    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: poly_coef_cbL_y

    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: poly_coef_cbR_x
    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: poly_coef_cbR_y
    !> @}

    !> @name The ideal weights at the left and the right cell-boundaries and at the
    !! left and the right quadrature points, in x-, y- and z-directions. Note
    !! that the first dimension of the array identifies the weight, while the
    !! last denotes the cell-location in the relevant coordinate direction.
    !> @{
    real(kind(0d0)), target, allocatable, dimension(:, :) :: d_cbL_x
    real(kind(0d0)), target, allocatable, dimension(:, :) :: d_cbL_y

    real(kind(0d0)), target, allocatable, dimension(:, :) :: d_cbR_x
    real(kind(0d0)), target, allocatable, dimension(:, :) :: d_cbR_y
    !> @}

    !> @name Smoothness indicator coefficients in the x-, y-, and z-directions. Note
    !! that the first array dimension identifies the smoothness indicator, the
    !! second identifies the position of its coefficients and the last denotes
    !! the cell-location in the relevant coordinate direction.
    !> @{
    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: beta_coef_x
    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: beta_coef_y
    !> @}

    ! END: WENO Coefficients ===================================================

    integer :: v_size !< Number of WENO-reconstructed cell-average variables

    !> @name Indical bounds in the s1-, s2- and s3-directions
    !> @{
    type(int_bounds_info) :: is1, is2
    !> @}

!$acc declare create( &
!$acc                v_rs_ws_x, v_rs_ws_y, &
!$acc                poly_coef_cbL_x,poly_coef_cbL_y, &
!$acc                poly_coef_cbR_x,poly_coef_cbR_y,d_cbL_x,       &
!$acc                d_cbL_y,d_cbR_x,d_cbR_y,beta_coef_x,beta_coef_y,   &
!$acc                v_size, is1, is2)

contains

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_weno_module() ! --------------------------------

        if (weno_order == 1) return

        ! Allocating/Computing WENO Coefficients in x-direction ============
        is1%beg = -buff_size; is1%end = m - is1%beg
        if (n == 0) then
            is2%beg = 0
        else
            is2%beg = -buff_size; 
        end if

        is2%end = n - is2%beg


        @:ALLOCATE(poly_coef_cbL_x(is1%beg + weno_polyn:is1%end - weno_polyn, 0:weno_polyn, &
                                  0:weno_polyn - 1))
        @:ALLOCATE(poly_coef_cbR_x(is1%beg + weno_polyn:is1%end - weno_polyn, 0:weno_polyn, &
                                  0:weno_polyn - 1))

        @:ALLOCATE(d_cbL_x(0:weno_polyn, is1%beg + weno_polyn:is1%end - weno_polyn))
        @:ALLOCATE(d_cbR_x(0:weno_polyn, is1%beg + weno_polyn:is1%end - weno_polyn))

        @:ALLOCATE(beta_coef_x(is1%beg + weno_polyn:is1%end - weno_polyn, 0:weno_polyn, &
                              0:2*(weno_polyn - 1)))

        call s_compute_weno_coefficients(1, is1)

        @:ALLOCATE(v_rs_ws_x(is1%beg:is1%end, &
                                 is2%beg:is2%end,  1:sys_size))

        ! ==================================================================

        ! Allocating/Computing WENO Coefficients in y-direction ============
        if (n == 0) return

        is2%beg = -buff_size; is2%end = n - is2%beg
        is1%beg = -buff_size; is1%end = m - is1%beg

        @:ALLOCATE(poly_coef_cbL_y(is2%beg + weno_polyn:is2%end - weno_polyn, 0:weno_polyn, &
                                  0:weno_polyn - 1))
        @:ALLOCATE(poly_coef_cbR_y(is2%beg + weno_polyn:is2%end - weno_polyn, 0:weno_polyn, &
                                  0:weno_polyn - 1))

        @:ALLOCATE(d_cbL_y(0:weno_polyn, is2%beg + weno_polyn:is2%end - weno_polyn))
        @:ALLOCATE(d_cbR_y(0:weno_polyn, is2%beg + weno_polyn:is2%end - weno_polyn))

        @:ALLOCATE(beta_coef_y(is2%beg + weno_polyn:is2%end - weno_polyn, 0:weno_polyn, &
                              0:2*(weno_polyn - 1)))

        call s_compute_weno_coefficients(2, is2)

        @:ALLOCATE(v_rs_ws_y(is2%beg:is2%end, &
                                 is1%beg:is1%end, 1:sys_size))


    end subroutine s_initialize_weno_module ! ------------------------------

    !>  The purpose of this subroutine is to compute the grid
        !!      dependent coefficients of the WENO polynomials, ideal
        !!      weights and smoothness indicators, provided the order,
        !!      the coordinate direction and the location of the WENO
        !!      reconstruction.
        !! @param weno_dir Coordinate direction of the WENO reconstruction
        !! @param is Index bounds in the s-direction
    subroutine s_compute_weno_coefficients(weno_dir, is) ! -------

        integer, intent(IN) :: weno_dir
        type(int_bounds_info), intent(IN) :: is
        integer :: s

        real(kind(0d0)), pointer, dimension(:) :: s_cb => null() !<
            !! Cell-boundary locations in the s-direction

        type(int_bounds_info) :: bc_s !< Boundary conditions (BC) in the s-direction

        integer :: i !< Generic loop iterator

        ! Determining the number of cells, the cell-boundary locations and
        ! the boundary conditions in the coordinate direction selected for
        ! the WENO reconstruction
        if (weno_dir == 1) then
            s = m; s_cb => x_cb; bc_s = bc_x
        elseif (weno_dir == 2) then
            s = n; s_cb => y_cb; bc_s = bc_y
        end if

        #:for WENO_DIR, XYZ in [(1, 'x'), (2, 'y')]
            ! Computing WENO3 Coefficients =====================================
            if (weno_dir == ${WENO_DIR}$) then
                if (weno_order == 3) then
                    do i = is%beg - 1 + weno_polyn, is%end - 1 - weno_polyn

                        poly_coef_cbR_${XYZ}$ (i + 1, 0, 0) = (s_cb(i) - s_cb(i + 1))/ &
                                                              (s_cb(i) - s_cb(i + 2))
                        poly_coef_cbR_${XYZ}$ (i + 1, 1, 0) = (s_cb(i) - s_cb(i + 1))/ &
                                                              (s_cb(i - 1) - s_cb(i + 1))

                        poly_coef_cbL_${XYZ}$ (i + 1, 0, 0) = -poly_coef_cbR_${XYZ}$ (i + 1, 0, 0)
                        poly_coef_cbL_${XYZ}$ (i + 1, 1, 0) = -poly_coef_cbR_${XYZ}$ (i + 1, 1, 0)

                        d_cbR_${XYZ}$ (0, i + 1) = (s_cb(i - 1) - s_cb(i + 1))/ &
                                                   (s_cb(i - 1) - s_cb(i + 2))
                        d_cbL_${XYZ}$ (0, i + 1) = (s_cb(i - 1) - s_cb(i))/ &
                                                   (s_cb(i - 1) - s_cb(i + 2))

                        d_cbR_${XYZ}$ (1, i + 1) = 1d0 - d_cbR_${XYZ}$ (0, i + 1)
                        d_cbL_${XYZ}$ (1, i + 1) = 1d0 - d_cbL_${XYZ}$ (0, i + 1)

                        beta_coef_${XYZ}$ (i + 1, 0, 0) = 4d0*(s_cb(i) - s_cb(i + 1))**2d0/ &
                                                          (s_cb(i) - s_cb(i + 2))**2d0
                        beta_coef_${XYZ}$ (i + 1, 1, 0) = 4d0*(s_cb(i) - s_cb(i + 1))**2d0/ &
                                                          (s_cb(i - 1) - s_cb(i + 1))**2d0

                    end do
                    ! END: Computing WENO3 Coefficients ================================

                    ! Computing WENO5 Coefficients =====================================
                else

                    do i = is%beg - 1 + weno_polyn, is%end - 1 - weno_polyn

                        poly_coef_cbR_${XYZ}$ (i + 1, 0, 0) = &
                            ((s_cb(i) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i + 2)))/ &
                            ((s_cb(i) - s_cb(i + 3))*(s_cb(i + 3) - s_cb(i + 1)))
                        poly_coef_cbR_${XYZ}$ (i + 1, 1, 0) = &
                            ((s_cb(i - 1) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i)))/ &
                            ((s_cb(i - 1) - s_cb(i + 2))*(s_cb(i + 2) - s_cb(i)))
                        poly_coef_cbR_${XYZ}$ (i + 1, 1, 1) = &
                            ((s_cb(i) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i + 2)))/ &
                            ((s_cb(i - 1) - s_cb(i + 1))*(s_cb(i - 1) - s_cb(i + 2)))
                        poly_coef_cbR_${XYZ}$ (i + 1, 2, 1) = &
                            ((s_cb(i) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i - 1)))/ &
                            ((s_cb(i - 2) - s_cb(i))*(s_cb(i - 2) - s_cb(i + 1)))
                        poly_coef_cbL_${XYZ}$ (i + 1, 0, 0) = &
                            ((s_cb(i + 1) - s_cb(i))*(s_cb(i) - s_cb(i + 2)))/ &
                            ((s_cb(i) - s_cb(i + 3))*(s_cb(i + 3) - s_cb(i + 1)))
                        poly_coef_cbL_${XYZ}$ (i + 1, 1, 0) = &
                            ((s_cb(i) - s_cb(i - 1))*(s_cb(i) - s_cb(i + 1)))/ &
                            ((s_cb(i - 1) - s_cb(i + 2))*(s_cb(i) - s_cb(i + 2)))
                        poly_coef_cbL_${XYZ}$ (i + 1, 1, 1) = &
                            ((s_cb(i + 1) - s_cb(i))*(s_cb(i) - s_cb(i + 2)))/ &
                            ((s_cb(i - 1) - s_cb(i + 1))*(s_cb(i - 1) - s_cb(i + 2)))
                        poly_coef_cbL_${XYZ}$ (i + 1, 2, 1) = &
                            ((s_cb(i - 1) - s_cb(i))*(s_cb(i) - s_cb(i + 1)))/ &
                            ((s_cb(i - 2) - s_cb(i))*(s_cb(i - 2) - s_cb(i + 1)))

                        poly_coef_cbR_${XYZ}$ (i + 1, 0, 1) = &
                            ((s_cb(i) - s_cb(i + 2)) + (s_cb(i + 1) - s_cb(i + 3)))/ &
                            ((s_cb(i) - s_cb(i + 2))*(s_cb(i) - s_cb(i + 3)))* &
                            ((s_cb(i) - s_cb(i + 1)))
                        poly_coef_cbR_${XYZ}$ (i + 1, 2, 0) = &
                            ((s_cb(i - 2) - s_cb(i + 1)) + (s_cb(i - 1) - s_cb(i + 1)))/ &
                            ((s_cb(i - 1) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i - 2)))* &
                            ((s_cb(i + 1) - s_cb(i)))
                        poly_coef_cbL_${XYZ}$ (i + 1, 0, 1) = &
                            ((s_cb(i) - s_cb(i + 2)) + (s_cb(i) - s_cb(i + 3)))/ &
                            ((s_cb(i) - s_cb(i + 2))*(s_cb(i) - s_cb(i + 3)))* &
                            ((s_cb(i + 1) - s_cb(i)))
                        poly_coef_cbL_${XYZ}$ (i + 1, 2, 0) = &
                            ((s_cb(i - 2) - s_cb(i)) + (s_cb(i - 1) - s_cb(i + 1)))/ &
                            ((s_cb(i - 2) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i - 1)))* &
                            ((s_cb(i) - s_cb(i + 1)))

                        d_cbR_${XYZ}$ (0, i + 1) = &
                            ((s_cb(i - 2) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i - 1)))/ &
                            ((s_cb(i - 2) - s_cb(i + 3))*(s_cb(i + 3) - s_cb(i - 1)))
                        d_cbR_${XYZ}$ (2, i + 1) = &
                            ((s_cb(i + 1) - s_cb(i + 2))*(s_cb(i + 1) - s_cb(i + 3)))/ &
                            ((s_cb(i - 2) - s_cb(i + 2))*(s_cb(i - 2) - s_cb(i + 3)))
                        d_cbL_${XYZ}$ (0, i + 1) = &
                            ((s_cb(i - 2) - s_cb(i))*(s_cb(i) - s_cb(i - 1)))/ &
                            ((s_cb(i - 2) - s_cb(i + 3))*(s_cb(i + 3) - s_cb(i - 1)))
                        d_cbL_${XYZ}$ (2, i + 1) = &
                            ((s_cb(i) - s_cb(i + 2))*(s_cb(i) - s_cb(i + 3)))/ &
                            ((s_cb(i - 2) - s_cb(i + 2))*(s_cb(i - 2) - s_cb(i + 3)))

                        d_cbR_${XYZ}$ (1, i + 1) = 1d0 - d_cbR_${XYZ}$ (0, i + 1) - d_cbR_${XYZ}$ (2, i + 1)
                        d_cbL_${XYZ}$ (1, i + 1) = 1d0 - d_cbL_${XYZ}$ (0, i + 1) - d_cbL_${XYZ}$ (2, i + 1)

                        beta_coef_${XYZ}$ (i + 1, 0, 0) = &
                            4d0*(s_cb(i) - s_cb(i + 1))**2d0*(10d0*(s_cb(i + 1) - &
                                                                    s_cb(i))**2d0 + (s_cb(i + 1) - s_cb(i))*(s_cb(i + 2) - &
                                                                      s_cb(i + 1)) + (s_cb(i + 2) - s_cb(i + 1))**2d0)/((s_cb(i) - &
                                                                                 s_cb(i + 3))**2d0*(s_cb(i + 1) - s_cb(i + 3))**2d0)

                        beta_coef_${XYZ}$ (i + 1, 0, 1) = &
                            4d0*(s_cb(i) - s_cb(i + 1))**2d0*(19d0*(s_cb(i + 1) - &
                                                                    s_cb(i))**2d0 - (s_cb(i + 1) - s_cb(i))*(s_cb(i + 3) - &
                                                                        s_cb(i + 1)) + 2d0*(s_cb(i + 2) - s_cb(i))*((s_cb(i + 2) - &
                                                                              s_cb(i)) + (s_cb(i + 3) - s_cb(i + 1))))/((s_cb(i) - &
                                                                          s_cb(i + 2))*(s_cb(i) - s_cb(i + 3))**2d0*(s_cb(i + 3) - &
                                                                                                                       s_cb(i + 1)))

                        beta_coef_${XYZ}$ (i + 1, 0, 2) = &
                            4d0*(s_cb(i) - s_cb(i + 1))**2d0*(10d0*(s_cb(i + 1) - &
                                                                    s_cb(i))**2d0 + (s_cb(i + 1) - s_cb(i))*((s_cb(i + 2) - &
                                                                         s_cb(i)) + (s_cb(i + 3) - s_cb(i + 1))) + ((s_cb(i + 2) - &
                                                                         s_cb(i)) + (s_cb(i + 3) - s_cb(i + 1)))**2d0)/((s_cb(i) - &
                                                                                     s_cb(i + 2))**2d0*(s_cb(i) - s_cb(i + 3))**2d0)

                        beta_coef_${XYZ}$ (i + 1, 1, 0) = &
                            4d0*(s_cb(i) - s_cb(i + 1))**2d0*(10d0*(s_cb(i + 1) - &
                                                                    s_cb(i))**2d0 + (s_cb(i) - s_cb(i - 1))**2d0 + (s_cb(i) - &
                                                                             s_cb(i - 1))*(s_cb(i + 1) - s_cb(i)))/((s_cb(i - 1) - &
                                                                                     s_cb(i + 2))**2d0*(s_cb(i) - s_cb(i + 2))**2d0)

                        beta_coef_${XYZ}$ (i + 1, 1, 1) = &
                            4d0*(s_cb(i) - s_cb(i + 1))**2d0*((s_cb(i) - &
                                                               s_cb(i + 1))*((s_cb(i) - s_cb(i - 1)) + 20d0*(s_cb(i + 1) - &
                                                                         s_cb(i))) + (2d0*(s_cb(i) - s_cb(i - 1)) + (s_cb(i + 1) - &
                                                                                s_cb(i)))*(s_cb(i + 2) - s_cb(i)))/((s_cb(i + 1) - &
                                                                      s_cb(i - 1))*(s_cb(i - 1) - s_cb(i + 2))**2d0*(s_cb(i + 2) - &
                                                                                                                           s_cb(i)))

                        beta_coef_${XYZ}$ (i + 1, 1, 2) = &
                            4d0*(s_cb(i) - s_cb(i + 1))**2d0*(10d0*(s_cb(i + 1) - &
                                                                    s_cb(i))**2d0 + (s_cb(i + 1) - s_cb(i))*(s_cb(i + 2) - &
                                                                                 s_cb(i + 1)) + (s_cb(i + 2) - s_cb(i + 1))**2d0)/ &
                            ((s_cb(i - 1) - s_cb(i + 1))**2d0*(s_cb(i - 1) - &
                                                               s_cb(i + 2))**2d0)

                        beta_coef_${XYZ}$ (i + 1, 2, 0) = &
                            4d0*(s_cb(i) - s_cb(i + 1))**2d0*(12d0*(s_cb(i + 1) - &
                                                                    s_cb(i))**2d0 + ((s_cb(i) - s_cb(i - 2)) + (s_cb(i) - &
                                                                               s_cb(i - 1)))**2d0 + 3d0*((s_cb(i) - s_cb(i - 2)) + &
                                                                                (s_cb(i) - s_cb(i - 1)))*(s_cb(i + 1) - s_cb(i)))/ &
                            ((s_cb(i - 2) - s_cb(i + 1))**2d0*(s_cb(i - 1) - &
                                                               s_cb(i + 1))**2d0)

                        beta_coef_${XYZ}$ (i + 1, 2, 1) = &
                            4d0*(s_cb(i) - s_cb(i + 1))**2d0*(19d0*(s_cb(i + 1) - &
                                                                    s_cb(i))**2d0 + ((s_cb(i) - s_cb(i - 2))*(s_cb(i) - &
                                                                       s_cb(i + 1))) + 2d0*(s_cb(i + 1) - s_cb(i - 1))*((s_cb(i) - &
                                                                      s_cb(i - 2)) + (s_cb(i + 1) - s_cb(i - 1))))/((s_cb(i - 2) - &
                                                                          s_cb(i))*(s_cb(i - 2) - s_cb(i + 1))**2d0*(s_cb(i + 1) - &
                                                                                                                       s_cb(i - 1)))

                        beta_coef_${XYZ}$ (i + 1, 2, 2) = &
                            4d0*(s_cb(i) - s_cb(i + 1))**2d0*(10d0*(s_cb(i + 1) - &
                                                                    s_cb(i))**2d0 + (s_cb(i) - s_cb(i - 1))**2d0 + (s_cb(i) - &
                                                                             s_cb(i - 1))*(s_cb(i + 1) - s_cb(i)))/((s_cb(i - 2) - &
                                                                                     s_cb(i))**2d0*(s_cb(i - 2) - s_cb(i + 1))**2d0)

                    end do
                end if
            end if
        #:endfor

! END: Computing WENO5 Coefficients ================================
        if (weno_dir == 1) then
!$acc update device(poly_coef_cbL_x, poly_coef_cbR_x, d_cbL_x, d_cbR_x, beta_coef_x)
        elseif (weno_dir == 2) then
!$acc update device(poly_coef_cbL_y, poly_coef_cbR_y, d_cbL_y, d_cbR_y, beta_coef_y)
        end if

        ! Nullifying WENO coefficients and cell-boundary locations pointers

        nullify (s_cb)

    end subroutine s_compute_weno_coefficients ! ---------------------------

    subroutine s_weno(v_vf, vL_rs_vf_x, vL_rs_vf_y, vR_rs_vf_x, vR_rs_vf_y, & ! -------------------
                          weno_dir, &
                          is1_d, is2_d)

        type(scalar_field), dimension(1:), intent(IN) :: v_vf
        real(kind(0d0)), dimension(startx:, starty:, 1:), intent(INOUT) ::  vL_rs_vf_x, vL_rs_vf_y, vR_rs_vf_x, vR_rs_vf_y
        integer, intent(IN) :: weno_dir
        type(int_bounds_info), intent(IN) :: is1_d, is2_d

        real(kind(0d0)), dimension(-weno_polyn:weno_polyn - 1) :: dvd
        real(kind(0d0)), dimension(0:weno_polyn) :: poly
        real(kind(0d0)), dimension(0:weno_polyn) :: alpha
        real(kind(0d0)), dimension(0:weno_polyn) :: omega
        real(kind(0d0)), dimension(0:weno_polyn) :: beta

        integer :: i, j, k

        is1 = is1_d
        is2 = is2_d

!$acc update device(is1, is2)

        if (weno_order /= 1) then
            call s_initialize_weno(v_vf, weno_dir)
        end if

        if (weno_order == 1) then
            if (weno_dir == 1) then
!$acc parallel loop collapse(4) default(present)
                do i = 1, ubound(v_vf, 1)
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                vL_rs_vf_x(j, k, i) = v_vf(i)%sf(j, k)
                                vR_rs_vf_x(j, k, i) = v_vf(i)%sf(j, k)
                            end do
                        end do
                end do
!$acc end parallel loop
            else if (weno_dir == 2) then
!$acc parallel loop collapse(4) default(present)
                do i = 1, ubound(v_vf, 1)
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                vL_rs_vf_y(j, k, i) = v_vf(i)%sf(k, j)
                                vR_rs_vf_y(j, k, i) = v_vf(i)%sf(k, j)
                            end do
                        end do
                end do
!$acc end parallel loop
            end if

        elseif (weno_order == 3) then
            #:for WENO_DIR, XYZ in [(1, 'x'), (2, 'y')]
            if (weno_dir == ${WENO_DIR}$) then
!$acc parallel loop collapse(4) gang vector default(present) private(beta,dvd,poly,omega,alpha)
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                do i = 1, v_size
                                    ! reconstruct from left side

                                    dvd(0) = v_rs_ws_${XYZ}$(j + 1, k, i) &
                                             - v_rs_ws_${XYZ}$(j, k, i)
                                    dvd(-1) = v_rs_ws_${XYZ}$(j, k, i) &
                                              - v_rs_ws_${XYZ}$(j - 1, k, i)

                                    poly(0) = v_rs_ws_${XYZ}$(j, k, i) &
                                              + poly_coef_cbL_${XYZ}$(j, 0, 0)*dvd(0)
                                    poly(1) = v_rs_ws_${XYZ}$(j, k, i) &
                                              + poly_coef_cbL_${XYZ}$(j, 1, 0)*dvd(-1)

                                    beta(0) = beta_coef_${XYZ}$(j, 0, 0)*dvd(0)*dvd(0) &
                                              + weno_eps
                                    beta(1) = beta_coef_${XYZ}$(j, 1, 0)*dvd(-1)*dvd(-1) &
                                              + weno_eps

                                    alpha = d_cbL_${XYZ}$(:, j)/(beta*beta)

                                    omega = alpha/sum(alpha)


                                    vL_rs_vf_${XYZ}$(j, k, i) = omega(0)*poly(0) + omega(1)*poly(1)

                                    ! reconstruct from right side

                                    poly(0) = v_rs_ws_${XYZ}$(j, k, i) &
                                              + poly_coef_cbR_${XYZ}$(j, 0, 0)*dvd(0)
                                    poly(1) = v_rs_ws_${XYZ}$(j, k, i) &
                                              + poly_coef_cbR_${XYZ}$(j, 1, 0)*dvd(-1)

                                    alpha = d_cbR_${XYZ}$(:, j)/(beta*beta)

                                    omega = alpha/sum(alpha)


                                    vR_rs_vf_${XYZ}$(j, k, i) = omega(0)*poly(0) + omega(1)*poly(1)

                                end do
                            end do
                    end do
!$acc end parallel loop
            end if
            #:endfor
        else
            #:for WENO_DIR, XYZ in [(1, 'x'), (2, 'y')]
            if (weno_dir == ${WENO_DIR}$) then
!$acc parallel loop gang vector collapse (3)  default(present) private(dvd, poly, beta, alpha, omega)
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end
!$acc loop seq
                            do i = 1, v_size

                                dvd(1) = v_rs_ws_${XYZ}$(j + 2, k, i) &
                                         - v_rs_ws_${XYZ}$(j + 1, k, i)
                                dvd(0) = v_rs_ws_${XYZ}$(j + 1, k, i) &
                                         - v_rs_ws_${XYZ}$(j, k, i)
                                dvd(-1) = v_rs_ws_${XYZ}$(j, k, i) &
                                          - v_rs_ws_${XYZ}$(j - 1, k, i)
                                dvd(-2) = v_rs_ws_${XYZ}$(j - 1, k, i) &
                                          - v_rs_ws_${XYZ}$(j - 2, k, i)

                                poly(0) = v_rs_ws_${XYZ}$(j, k, i) &
                                          + poly_coef_cbL_${XYZ}$(j, 0, 0)*dvd(1) &
                                          + poly_coef_cbL_${XYZ}$(j, 0, 1)*dvd(0)
                                poly(1) = v_rs_ws_${XYZ}$(j, k, i) &
                                          + poly_coef_cbL_${XYZ}$(j, 1, 0)*dvd(0) &
                                          + poly_coef_cbL_${XYZ}$(j, 1, 1)*dvd(-1)
                                poly(2) = v_rs_ws_${XYZ}$(j, k, i) &
                                          + poly_coef_cbL_${XYZ}$(j, 2, 0)*dvd(-1) &
                                          + poly_coef_cbL_${XYZ}$(j, 2, 1)*dvd(-2)

                                beta(0) = beta_coef_${XYZ}$(j, 0, 0)*dvd(1)*dvd(1) &
                                          + beta_coef_${XYZ}$(j, 0, 1)*dvd(1)*dvd(0) &
                                          + beta_coef_${XYZ}$(j, 0, 2)*dvd(0)*dvd(0) &
                                          + weno_eps
                                beta(1) = beta_coef_${XYZ}$(j, 1, 0)*dvd(0)*dvd(0) &
                                          + beta_coef_${XYZ}$(j, 1, 1)*dvd(0)*dvd(-1) &
                                          + beta_coef_${XYZ}$(j, 1, 2)*dvd(-1)*dvd(-1) &
                                          + weno_eps
                                beta(2) = beta_coef_${XYZ}$(j, 2, 0)*dvd(-1)*dvd(-1) &
                                          + beta_coef_${XYZ}$(j, 2, 1)*dvd(-1)*dvd(-2) &
                                          + beta_coef_${XYZ}$(j, 2, 2)*dvd(-2)*dvd(-2) &
                                          + weno_eps

                                alpha = d_cbL_${XYZ}$(:, j)/(beta*beta)

                                omega = alpha/sum(alpha)



                                vL_rs_vf_${XYZ}$(j, k, i) = sum(omega*poly)

                                poly(0) = v_rs_ws_${XYZ}$(j, k, i) &
                                          + poly_coef_cbR_${XYZ}$(j, 0, 0)*dvd(1) &
                                          + poly_coef_cbR_${XYZ}$(j, 0, 1)*dvd(0)
                                poly(1) = v_rs_ws_${XYZ}$(j, k, i) &
                                          + poly_coef_cbR_${XYZ}$(j, 1, 0)*dvd(0) &
                                          + poly_coef_cbR_${XYZ}$(j, 1, 1)*dvd(-1)
                                poly(2) = v_rs_ws_${XYZ}$(j, k, i) &
                                          + poly_coef_cbR_${XYZ}$(j, 2, 0)*dvd(-1) &
                                          + poly_coef_cbR_${XYZ}$(j, 2, 1)*dvd(-2)

                                alpha = d_cbR_${XYZ}$(:, j)/(beta*beta)

                                omega = alpha/sum(alpha)


                                vR_rs_vf_${XYZ}$(j, k, i) = sum(omega*poly)

                            end do
                        end do
                end do
!$acc end parallel loop

            end if
            #:endfor
        end if

    end subroutine s_weno

    !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are required for the setup of the
        !!      WENO reconstruction.
        !! @param v_vf Cell-averaged variables
        !! @param vL_vf Left WENO reconstructed cell-boundary values
        !! @param vR_vf Right WENO reconstructed cell-boundary values
        !! @param weno_dir Coordinate direction of the WENO reconstruction
        !! @param is1 Index bounds in first coordinate direction
        !! @param is2 Index bounds in second coordinate direction
    subroutine s_initialize_weno(v_vf, weno_dir)

        type(scalar_field), dimension(:), intent(IN) :: v_vf
        integer, intent(IN) :: weno_dir
        integer :: j, k, l

        ! Determining the number of cell-average variables which will be
        ! WENO-reconstructed and mapping their indical bounds in the x-,
        ! y- and z-directions to those in the s1-, s2- and s3-directions
        ! as to reshape the inputted data in the coordinate direction of
        ! the WENO reconstruction
        v_size = ubound(v_vf, 1)

        !$acc update device(v_size)

        if (weno_dir == 1) then
!$acc parallel loop collapse(4) gang vector default(present)
            do j = 1, v_size
                do l = is2%beg, is2%end
                    do k = is1%beg - weno_polyn, is1%end + weno_polyn
                        v_rs_ws_x(k, l, j) = v_vf(j)%sf(k, l)
                    end do
                end do
            end do
!$acc end parallel loop
        end if

        ! ==================================================================

        ! Reshaping/Projecting onto Characteristic Fields in y-direction ===
        if (n == 0) return

        if (weno_dir == 2) then
!$acc parallel loop collapse(4) gang vector default(present)
            do j = 1, v_size
                do l = is2%beg, is2%end
                    do k = is1%beg - weno_polyn, is1%end + weno_polyn
                        v_rs_ws_y(k, l, j) = v_vf(j)%sf(l, k)
                    end do
                end do
            end do
!$acc end parallel loop
        end if

    end subroutine s_initialize_weno ! -------------------------------------


    !>  Module deallocation and/or disassociation procedures
    subroutine s_finalize_weno_module() ! ----------------------------------

        if (weno_order == 1) return

        ! Deallocating the WENO-stencil of the WENO-reconstructed variables

        !deallocate(vL_rs_vf_x, vR_rs_vf_x)
        @:DEALLOCATE(v_rs_ws_x)

        ! Deallocating WENO coefficients in x-direction ====================
        @:DEALLOCATE(poly_coef_cbL_x, poly_coef_cbR_x)
        @:DEALLOCATE(d_cbL_x, d_cbR_x)
        @:DEALLOCATE(beta_coef_x)
        ! ==================================================================

        ! Deallocating WENO coefficients in y-direction ====================
        if (n == 0) return

        !deallocate(vL_rs_vf_y, vR_rs_vf_y)
        @:DEALLOCATE(v_rs_ws_y)

        @:DEALLOCATE(poly_coef_cbL_y, poly_coef_cbR_y)
        @:DEALLOCATE(d_cbL_y, d_cbR_y)
        @:DEALLOCATE(beta_coef_y)
        ! ==================================================================

    end subroutine s_finalize_weno_module ! --------------------------------

end module m_weno
