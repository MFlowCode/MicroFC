!>
!! @file m_viscous.f90
!! @brief Contains module m_viscous

#:include 'macros.fpp'

!> @brief The module contains the subroutines used to compute viscous terms.
module m_viscous

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_weno
    ! ==========================================================================

    private; public  s_get_viscous, &
    s_initialize_viscous_module, &
    s_reconstruct_cell_boundary_values_visc_deriv, &
    s_finalize_viscous_module

    type(int_bounds_info) :: iv
    type(int_bounds_info) :: is1, is2
 !$acc declare create(is1, is2, iv)   

    real(kind(0d0)), allocatable, dimension(:, :) :: Res
!$acc declare create(Res)


    contains

    subroutine s_initialize_viscous_module()
        integer :: i, j !< generic loop iterators

        @:ALLOCATE(Res(1:2, 1:maxval(Re_size)))

        do i = 1, 2
            do j = 1, Re_size(i)
                Res(i, j) = fluid_pp(Re_idx(i, j))%Re(i)
            end do
        end do
        !$acc update device(Res, Re_idx, Re_size)

    end subroutine s_initialize_viscous_module


    !>  Computes the scalar gradient fields via finite differences
        !!  @param var Variable to compute derivative of
        !!  @param grad_x First coordinate direction component of the derivative
        !!  @param grad_y Second coordinate direction component of the derivative
    subroutine s_compute_fd_gradient(var, grad_x, grad_y, ix, iy)

        type(scalar_field), intent(IN) :: var
        type(scalar_field), intent(INOUT) :: grad_x
        type(scalar_field), intent(INOUT) :: grad_y

        integer :: j, k !< Generic loop iterators

        type(int_bounds_info) :: ix, iy

        ix%beg = -buff_size; ix%end = m + buff_size; 
        if (n > 0) then
            iy%beg = -buff_size; iy%end = n + buff_size
        else
            iy%beg = -1; iy%end = 1
        end if

        !$acc update device(ix, iy)

    !$acc parallel loop collapse(2) gang vector default(present)
            do k = iy%beg + 1, iy%end - 1
                do j = ix%beg + 1, ix%end - 1
                    grad_x%sf(j, k) = &
                        (var%sf(j + 1, k) - var%sf(j - 1, k))/ &
                        (x_cc(j + 1) - x_cc(j - 1))
                end do
            end do

        if (n > 0) then
    !$acc parallel loop collapse(2) gang vector
                do k = iy%beg + 1, iy%end - 1
                    do j = ix%beg + 1, ix%end - 1
                        grad_y%sf(j, k) = &
                            (var%sf(j, k + 1) - var%sf(j, k - 1))/ &
                            (y_cc(k + 1) - y_cc(k - 1))
                    end do
                end do
        end if

        ix%beg = -buff_size; ix%end = m + buff_size; 
        if (n > 0) then
            iy%beg = -buff_size; iy%end = n + buff_size
        else
            iy%beg = 0; iy%end = 0
        end if


        !$acc update device(ix, iy)

    !$acc parallel loop collapse(1) gang vector default(present)
            do k = iy%beg, iy%end
                grad_x%sf(ix%beg, k) = &
                    (-3d0*var%sf(ix%beg, k) + 4d0*var%sf(ix%beg + 1, k) - var%sf(ix%beg + 2, k))/ &
                    (x_cc(ix%beg + 2) - x_cc(ix%beg))
                grad_x%sf(ix%end, k) = &
                    (3d0*var%sf(ix%end, k) - 4d0*var%sf(ix%end - 1, k) + var%sf(ix%end - 2, k))/ &
                    (x_cc(ix%end) - x_cc(ix%end - 2))
            end do
        if (n > 0) then
    !$acc parallel loop collapse(1) gang vector default(present)
                do j = ix%beg, ix%end
                    grad_y%sf(j, iy%beg) = &
                        (-3d0*var%sf(j, iy%beg) + 4d0*var%sf(j, iy%beg + 1) - var%sf(j, iy%beg + 2))/ &
                        (y_cc(iy%beg + 2) - y_cc(iy%beg))
                    grad_y%sf(j, iy%end) = &
                        (3d0*var%sf(j, iy%end) - 4d0*var%sf(j, iy%end - 1) + var%sf(j, iy%end - 2))/ &
                        (y_cc(iy%end) - y_cc(iy%end - 2))
                end do
        end if

        if (bc_x%beg <= -3) then
    !$acc parallel loop collapse(1) gang vector default(present)
                do k = iy%beg, iy%end
                    grad_x%sf(0, k) = (-3d0*var%sf(0, k) + 4d0*var%sf(1, k) - var%sf(2, k))/ &
                                        (x_cc(2) - x_cc(0))
                end do
        end if
        if (bc_x%end <= -3) then
    !$acc parallel loop collapse(1) gang vector default(present)
                do k = iy%beg, iy%end
                    grad_x%sf(m, k) = (3d0*var%sf(m, k) - 4d0*var%sf(m - 1, k) + var%sf(m - 2, k))/ &
                                        (x_cc(m) - x_cc(m - 2))
                end do
        end if
        if (n > 0) then
            if (bc_y%beg <= -3 .and. bc_y%beg /= -13) then
    !$acc parallel loop collapse(1) gang vector default(present)
                    do j = ix%beg, ix%end
                        grad_y%sf(j, 0) = (-3d0*var%sf(j, 0) + 4d0*var%sf(j, 1) - var%sf(j, 2))/ &
                                            (y_cc(2) - y_cc(0))
                    end do
            end if
            if (bc_y%end <= -3) then
    !$acc parallel loop collapse(1) gang vector default(present)
                    do j = ix%beg, ix%end
                        grad_y%sf(j, n) = (3d0*var%sf(j, n) - 4d0*var%sf(j, n - 1) + var%sf(j, n - 2))/ &
                                            (y_cc(n) - y_cc(n - 2))
                    end do
            end if
        end if

    end subroutine s_compute_fd_gradient ! --------------------------------------

!>  Computes viscous terms
    !!  @param q_cons_vf Cell-averaged conservative variables
    !!  @param q_prim_vf Cell-averaged primitive variables
    !!  @param rhs_vf Cell-averaged RHS variables
    subroutine s_get_viscous(qL_prim_rsx_vf, qL_prim_rsy_vf, &
                             dqL_prim_dx_n, dqL_prim_dy_n, &
                             qL_prim, & 
                             qR_prim_rsx_vf, qR_prim_rsy_vf, &
                             dqR_prim_dx_n, dqR_prim_dy_n, &
                             qR_prim, &
                             q_prim_qp, &
                             dq_prim_dx_qp, dq_prim_dy_qp,  &
                             ix, iy)

        real(kind(0d0)), dimension(startx:, starty:, 1:), &
             intent(INOUT) :: qL_prim_rsx_vf, qR_prim_rsx_vf, &
                              qL_prim_rsy_vf, qR_prim_rsy_vf

        type(vector_field), dimension(sys_size) :: qL_prim, qR_prim

        type(vector_field) :: q_prim_qp

        type(vector_field), dimension(1:num_dims), &
            intent(INOUT) :: dqL_prim_dx_n, dqR_prim_dx_n, &
                             dqL_prim_dy_n, dqR_prim_dy_n

        type(vector_field) :: dq_prim_dx_qp, dq_prim_dy_qp

        integer :: i, j, k
        type(int_bounds_info), intent(IN) :: ix, iy

        do i = 1, num_dims

            iv%beg = mom_idx%beg; iv%end = mom_idx%end

            !$acc update device(iv)

            call s_reconstruct_cell_boundary_values_visc( &
                q_prim_qp%vf(iv%beg:iv%end), &
                qL_prim_rsx_vf, qL_prim_rsy_vf, &
                qR_prim_rsx_vf, qR_prim_rsy_vf, &
                i, qL_prim(i)%vf(iv%beg:iv%end), qR_prim(i)%vf(iv%beg:iv%end), &
                ix, iy)
        end do

        if (weno_Re_flux) then
            ! Compute velocity gradient at cell centers using scalar
            ! divergence theorem
            do i = 1, num_dims
                if (i == 1) then
                    call s_apply_scalar_divergence_theorem( &
                        qL_prim(i)%vf(iv%beg:iv%end), &
                        qR_prim(i)%vf(iv%beg:iv%end), &
                        dq_prim_dx_qp%vf(iv%beg:iv%end), i, &
                        ix, iy)
                elseif (i == 2) then
                    call s_apply_scalar_divergence_theorem( &
                        qL_prim(i)%vf(iv%beg:iv%end), &
                        qR_prim(i)%vf(iv%beg:iv%end), &
                        dq_prim_dy_qp%vf(iv%beg:iv%end), i, &
                        ix, iy)
                end if
            end do

        else ! Compute velocity gradient at cell centers using finite differences

            iv%beg = mom_idx%beg; iv%end = mom_idx%end
            !$acc update device(iv)

    !$acc parallel loop collapse(2) gang vector default(present)
                do k = iy%beg, iy%end
                    do j = ix%beg + 1, ix%end
    !$acc loop seq
                        do i = iv%beg, iv%end
                            dqL_prim_dx_n(1)%vf(i)%sf(j, k) = &
                                (q_prim_qp%vf(i)%sf(j, k) - &
                                q_prim_qp%vf(i)%sf(j - 1, k))/ &
                                (x_cc(j) - x_cc(j - 1))
                        end do
                    end do
                end do

    !$acc parallel loop collapse(2) gang vector default(present)
                do k = iy%beg, iy%end
                    do j = ix%beg, ix%end - 1
    !$acc loop seq
                        do i = iv%beg, iv%end
                            dqR_prim_dx_n(1)%vf(i)%sf(j, k) = &
                                (q_prim_qp%vf(i)%sf(j + 1, k) - &
                                q_prim_qp%vf(i)%sf(j, k))/ &
                                (x_cc(j + 1) - x_cc(j))
                        end do
                    end do
                end do

            if (n > 0) then

    !$acc parallel loop collapse(2) gang vector default(present)
                    do j = iy%beg + 1, iy%end
                        do k = ix%beg, ix%end
                            !$acc loop seq
                            do i = iv%beg, iv%end
                                dqL_prim_dy_n(2)%vf(i)%sf(k, j) = &
                                    (q_prim_qp%vf(i)%sf(k, j) - &
                                    q_prim_qp%vf(i)%sf(k, j - 1))/ &
                                    (y_cc(j) - y_cc(j - 1))
                            end do
                        end do
                    end do

    !$acc parallel loop collapse(2) gang vector default(present)
                    do j = iy%beg, iy%end - 1
                        do k = ix%beg, ix%end
                            !$acc loop seq
                            do i = iv%beg, iv%end
                                dqR_prim_dy_n(2)%vf(i)%sf(k, j) = &
                                    (q_prim_qp%vf(i)%sf(k, j + 1) - &
                                    q_prim_qp%vf(i)%sf(k, j))/ &
                                    (y_cc(j + 1) - y_cc(j))
                            end do
                        end do
                    end do

    !$acc parallel loop collapse(2) gang vector default(present)
                    do j = iy%beg + 1, iy%end
                        do k = ix%beg + 1, ix%end - 1
    !$acc loop seq
                            do i = iv%beg, iv%end
                                dqL_prim_dx_n(2)%vf(i)%sf(k, j) = &
                                    (dqL_prim_dx_n(1)%vf(i)%sf(k, j) + &
                                    dqR_prim_dx_n(1)%vf(i)%sf(k, j) + &
                                    dqL_prim_dx_n(1)%vf(i)%sf(k, j - 1) + &
                                    dqR_prim_dx_n(1)%vf(i)%sf(k, j - 1))

                                dqL_prim_dx_n(2)%vf(i)%sf(k, j) = 25d-2* &
                                                                    dqL_prim_dx_n(2)%vf(i)%sf(k, j)
                            end do
                        end do
                    end do

    !$acc parallel loop collapse(2) gang vector default(present)
                    do j = iy%beg, iy%end - 1
                        do k = ix%beg + 1, ix%end - 1
    !$acc loop seq
                            do i = iv%beg, iv%end
                                dqR_prim_dx_n(2)%vf(i)%sf(k, j) = &
                                    (dqL_prim_dx_n(1)%vf(i)%sf(k, j + 1) + &
                                    dqR_prim_dx_n(1)%vf(i)%sf(k, j + 1) + &
                                    dqL_prim_dx_n(1)%vf(i)%sf(k, j) + &
                                    dqR_prim_dx_n(1)%vf(i)%sf(k, j))

                                dqR_prim_dx_n(2)%vf(i)%sf(k, j) = 25d-2* &
                                                                    dqR_prim_dx_n(2)%vf(i)%sf(k, j)

                            end do
                        end do
                    end do

    !$acc parallel loop collapse(2) gang vector default(present)
                    do k = iy%beg + 1, iy%end - 1
                        do j = ix%beg + 1, ix%end
    !$acc loop seq
                            do i = iv%beg, iv%end
                                dqL_prim_dy_n(1)%vf(i)%sf(j, k) = &
                                    (dqL_prim_dy_n(2)%vf(i)%sf(j, k) + &
                                    dqR_prim_dy_n(2)%vf(i)%sf(j, k) + &
                                    dqL_prim_dy_n(2)%vf(i)%sf(j - 1, k) + &
                                    dqR_prim_dy_n(2)%vf(i)%sf(j - 1, k))

                                dqL_prim_dy_n(1)%vf(i)%sf(j, k) = 25d-2* &
                                                                    dqL_prim_dy_n(1)%vf(i)%sf(j, k)

                            end do
                        end do
                    end do

    !$acc parallel loop collapse(2) gang vector default(present)
                    do k = iy%beg + 1, iy%end - 1
                        do j = ix%beg, ix%end - 1
    !$acc loop seq
                            do i = iv%beg, iv%end
                                dqR_prim_dy_n(1)%vf(i)%sf(j, k) = &
                                    (dqL_prim_dy_n(2)%vf(i)%sf(j + 1, k) + &
                                    dqR_prim_dy_n(2)%vf(i)%sf(j + 1, k) + &
                                    dqL_prim_dy_n(2)%vf(i)%sf(j, k) + &
                                    dqR_prim_dy_n(2)%vf(i)%sf(j, k))

                                dqR_prim_dy_n(1)%vf(i)%sf(j, k) = 25d-2* &
                                                                    dqR_prim_dy_n(1)%vf(i)%sf(j, k)

                            end do
                        end do
                    end do


                    do i = iv%beg, iv%end
                        call s_compute_fd_gradient(q_prim_qp%vf(i), &
                                                dq_prim_dx_qp%vf(i), &
                                                dq_prim_dy_qp%vf(i), &
                                                ix, iy)
                    end do


            else
                do i = iv%beg, iv%end
                    call s_compute_fd_gradient(q_prim_qp%vf(i), &
                                            dq_prim_dx_qp%vf(i), &
                                            dq_prim_dx_qp%vf(i), &
                                            ix, iy)
                end do

            end if

        end if

    end subroutine s_get_viscous

        !>  The purpose of this subroutine is to employ the inputted
        !!      left and right cell-boundary integral-averaged variables
        !!      to compute the relevant cell-average first-order spatial
        !!      derivatives in the x-, y- or z-direction by means of the
        !!      scalar divergence theorem.
        !!  @param vL_vf Left cell-boundary integral averages
        !!  @param vR_vf Right cell-boundary integral averages
        !!  @param dv_ds_vf Cell-average first-order spatial derivatives
        !!  @param norm_dir Splitting coordinate direction
    subroutine s_apply_scalar_divergence_theorem(vL_vf, vR_vf, & ! --------
                                                 dv_ds_vf, &
                                                 norm_dir, &
                                                 ix, iy)

        type(scalar_field), &
            dimension(iv%beg:iv%end), &
            intent(IN) :: vL_vf, vR_vf

        type(scalar_field), &
            dimension(iv%beg:iv%end), &
            intent(INOUT) :: dv_ds_vf

        integer, intent(IN) :: norm_dir

        integer :: i, j, k !< Generic loop iterators

        type(int_bounds_info) :: ix, iy

        !$acc update device(ix, iy, iv)

        ! First-Order Spatial Derivatives in x-direction ===================
        if (norm_dir == 1) then

            ! A general application of the scalar divergence theorem that
            ! utilizes the left and right cell-boundary integral-averages,
            ! inside each cell, or an arithmetic mean of these two at the
            ! cell-boundaries, to calculate the cell-averaged first-order
            ! spatial derivatives inside the cell.

!$acc parallel loop collapse(2) gang vector default(present)
                do k = iy%beg, iy%end
                    do j = ix%beg + 1, ix%end - 1
!$acc loop seq
                        do i = iv%beg, iv%end

                            dv_ds_vf(i)%sf(j, k) = &
                                1d0/dx(j) &
                                *( &
                                vR_vf(i)%sf(j, k) &
                                - vL_vf(i)%sf(j, k) &
                                )
                        end do
                    end do
                end do

            ! END: First-Order Spatial Derivatives in x-direction ==============

            ! First-Order Spatial Derivatives in y-direction ===================
        elseif (norm_dir == 2) then

            ! A general application of the scalar divergence theorem that
            ! utilizes the left and right cell-boundary integral-averages,
            ! inside each cell, or an arithmetic mean of these two at the
            ! cell-boundaries, to calculate the cell-averaged first-order
            ! spatial derivatives inside the cell.

!$acc parallel loop collapse(2) gang vector default(present)

                do k = iy%beg + 1, iy%end - 1
                    do j = ix%beg, ix%end
!$acc loop seq
                        do i = iv%beg, iv%end
                            dv_ds_vf(i)%sf(j, k) = &
                                1d0/dy(k) &
                                *( &
                                vR_vf(i)%sf(j, k) &
                                - vL_vf(i)%sf(j, k) &
                                )
                        end do
                    end do
                end do
        end if

    end subroutine s_apply_scalar_divergence_theorem ! ---------------------

    subroutine s_reconstruct_cell_boundary_values_visc(v_vf, vL_x, vL_y, vR_x, vR_y, & ! -
                                                       norm_dir, vL_prim_vf, vR_prim_vf, ix, iy)

        type(scalar_field), dimension(iv%beg:iv%end), intent(IN) :: v_vf
        type(scalar_field), dimension(iv%beg:iv%end), intent(INOUT) :: vL_prim_vf, vR_prim_vf

        real(kind(0d0)), dimension(startx:, starty:, 1:), intent(INOUT) :: vL_x, vL_y, vR_x, vR_y

        integer, intent(IN) :: norm_dir

        integer :: weno_dir !< Coordinate direction of the WENO reconstruction

        integer :: i, j, k

        type(int_bounds_info) :: ix, iy
        ! Reconstruction in s1-direction ===================================

        if (norm_dir == 1) then
            is1 = ix; is2 = iy
            weno_dir = 1; is1%beg = is1%beg + weno_polyn
            is1%end = is1%end - weno_polyn

        elseif (norm_dir == 2) then
            is1 = iy; is2 = ix
            weno_dir = 2; is1%beg = is1%beg + weno_polyn
            is1%end = is1%end - weno_polyn
        end if

        !$acc update device(is1, is2, is3, iv)

        if (n > 0) then
            call s_weno(v_vf(iv%beg:iv%end), &
                vL_x(:, :, iv%beg:iv%end), vL_y(:, :, iv%beg:iv%end),  vR_x(:, :, iv%beg:iv%end), vR_y(:, :, iv%beg:iv%end), weno_dir, is1, is2)
        else
            call s_weno(v_vf(iv%beg:iv%end), &
                        vL_x(:, :, iv%beg:iv%end), vL_y(:, :, :), vR_x(:, :, iv%beg:iv%end), vR_y(:, :, :), weno_dir, is1, is2)
        end if

        if (any(Re_size > 0)) then
            if (weno_Re_flux) then
                if (norm_dir == 2) then
!$acc parallel loop collapse(3) gang vector default(present)
                    do i = iv%beg, iv%end
                            do j = is1%beg, is1%end
                                do k = is2%beg, is2%end
                                    vL_prim_vf(i)%sf(k, j) = vL_y(j, k, i)
                                    vR_prim_vf(i)%sf(k, j) = vR_y(j, k, i)
                                end do
                        end do
                    end do
                elseif (norm_dir == 1) then
!$acc parallel loop collapse(3) gang vector default(present)
                    do i = iv%beg, iv%end
                            do k = is2%beg, is2%end
                                do j = is1%beg, is1%end
                                    vL_prim_vf(i)%sf(j, k) = vL_x(j, k, i)
                                    vR_prim_vf(i)%sf(j, k) = vR_x(j, k, i)
                                end do
                            end do
                    end do
                end if
            end if
        end if

        ! ==================================================================

    end subroutine s_reconstruct_cell_boundary_values_visc ! --------------------

    subroutine s_reconstruct_cell_boundary_values_visc_deriv(v_vf, vL_x, vL_y, vR_x, vR_y, & ! -
                                                             norm_dir, vL_prim_vf, vR_prim_vf, ix, iy)

        type(scalar_field), dimension(iv%beg:iv%end), intent(IN) :: v_vf
        type(scalar_field), dimension(iv%beg:iv%end), intent(INOUT) :: vL_prim_vf, vR_prim_vf

        type(int_bounds_info) :: ix, iy

        real(kind(0d0)), dimension(startx:, starty:, iv%beg:), intent(INOUT) :: vL_x, vL_y, vR_x, vR_y

        integer, intent(IN) :: norm_dir

        integer :: weno_dir !< Coordinate direction of the WENO reconstruction

        integer :: i, j, k
        ! Reconstruction in s1-direction ===================================

        if (norm_dir == 1) then
            is1 = ix; is2 = iy
            weno_dir = 1; is1%beg = is1%beg + weno_polyn
            is1%end = is1%end - weno_polyn

        elseif (norm_dir == 2) then
            is1 = iy; is2 = ix
            weno_dir = 2; is1%beg = is1%beg + weno_polyn
            is1%end = is1%end - weno_polyn
        end if

        !$acc update device(is1, is2, iv)

        if (n > 0) then
            call s_weno(v_vf(iv%beg:iv%end), &
                vL_x(:, :, iv%beg:iv%end), vL_y(:, :, iv%beg:iv%end), vR_x(:, :, iv%beg:iv%end), vR_y(:, :, iv%beg:iv%end), weno_dir, is1, is2)
        else
            call s_weno(v_vf(iv%beg:iv%end), &
                        vL_x(:, :, iv%beg:iv%end), vL_y(:, :, :), vR_x(:, :, iv%beg:iv%end), vR_y(:, :, :), weno_dir, is1, is2)
        end if

        if (any(Re_size > 0)) then
            if (weno_Re_flux) then
                if (norm_dir == 2) then
!$acc parallel loop collapse(3) gang vector default(present)
                    do i = iv%beg, iv%end
                        do j = is1%beg, is1%end
                            do k = is2%beg, is2%end
                                vL_prim_vf(i)%sf(k, j) = vL_y(j, k, i)
                                vR_prim_vf(i)%sf(k, j) = vR_y(j, k, i)
                            end do
                        end do
                    end do
                elseif (norm_dir == 1) then
!$acc parallel loop collapse(3) gang vector default(present)
                    do i = iv%beg, iv%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                vL_prim_vf(i)%sf(j, k) = vL_x(j, k, i)
                                vR_prim_vf(i)%sf(j, k) = vR_x(j, k, i)
                            end do
                        end do
                    end do
                end if
            end if
        end if

        ! ==================================================================

    end subroutine s_reconstruct_cell_boundary_values_visc_deriv ! --------------------

    subroutine s_finalize_viscous_module()
        @:DEALLOCATE(Res)
    end subroutine s_finalize_viscous_module

end module m_viscous
