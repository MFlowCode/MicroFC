module m_create_patches

    ! Dependencies =============================================================
    use m_derived_types         ! Definitions of the derived types

    use m_global_parameters     ! Global parameters for the code

    use m_variables_conversion  ! Subroutines to change the state variables from
    ! one form to another

    use m_assign_patches
    ! ==========================================================================

    implicit none


    real(kind(0d0)) :: radius
    real(kind(0d0)) :: smooth_coeff !<
    !! These variables are analogous in both meaning and use to the similarly
    !! named components in the ic_patch_parameters type (see m_derived_types.f90
    !! for additional details). They are employed as a means to more concisely
    !! perform the actions necessary to lay out a particular patch on the grid.

    type(bounds_info) :: x_boundary, y_boundary !<
    !! These variables combine the centroid and length parameters associated with
    !! a particular patch to yield the locations of the patch boundaries in the
    !! x-, y- and z-coordinate directions. They are used as a means to concisely
    !! perform the actions necessary to lay out a particular patch on the grid.

    real(kind(0d0)) :: length_x, length_y

    real(kind(0d0)) :: a, b, c, d !<
    !! When a line or a plane sweep patch geometry is employed, these variables
    !! represent the coefficients associated with the equation describing the
    !! said line or plane.

    real(kind(0d0)) :: cart_y
    
contains

    !>          The line segment patch is a 1D geometry that may be used,
        !!              for example, in creating a Riemann problem. The geometry
        !!              of the patch is well-defined when its centroid and length
        !!              in the x-coordinate direction are provided. Note that the
        !!              line segment patch DOES NOT allow for the smearing of its
        !!              boundaries.
        !! @param patch_id patch identifier
    subroutine s_line_segment(patch_id) ! ----------------------------------

        integer, intent(IN) :: patch_id

        real(kind(0d0)) :: pi_inf, gamma, lit_gamma

        integer :: i, j  !< Generic loop operators

        pi_inf = fluid_pp(1)%pi_inf
        gamma = fluid_pp(1)%gamma
        lit_gamma = (1d0 + gamma)/gamma

        ! Transferring the line segment's centroid and length information
        x_centroid = patch_icpp(patch_id)%x_centroid
        length_x = patch_icpp(patch_id)%length_x

        ! Computing the beginning and end x-coordinates of the line segment
        ! based on its centroid and length
        x_boundary%beg = x_centroid - 0.5d0*length_x
        x_boundary%end = x_centroid + 0.5d0*length_x

        ! Since the line segment patch does not allow for its boundaries to
        ! be smoothed out, the pseudo volume fraction is set to 1 to ensure
        ! that only the current patch contributes to the fluid state in the
        ! cells that this patch covers.
        eta = 1d0

        ! Checking whether the line segment covers a particular cell in the
        ! domain and verifying whether the current patch has the permission
        ! to write to that cell. If both queries check out, the primitive
        ! variables of the current patch are assigned to this cell.
        do i = 0, m
            if (x_boundary%beg <= x_cc(i) .and. &
                x_boundary%end >= x_cc(i) .and. &
                patch_icpp(patch_id)%alter_patch(patch_id_fp(i, 0))) then

                call s_assign_patch_species_primitive_variables(patch_id, i, 0)

            end if
        end do

    end subroutine s_line_segment ! ----------------------------------------

    !> The circular patch is a 2D geometry that may be used, for
        !!              example, in creating a bubble or a droplet. The geometry
        !!              of the patch is well-defined when its centroid and radius
        !!              are provided. Note that the circular patch DOES allow for
        !!              the smoothing of its boundary.
        !! @param patch_id is the patch identifier
    subroutine s_circle(patch_id) ! ----------------------------------------

        integer, intent(IN) :: patch_id

        integer :: i, j !< Generic loop iterators

        ! Transferring the circular patch's radius, centroid, smearing patch
        ! identity and smearing coefficient information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        radius = patch_icpp(patch_id)%radius
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
        smooth_coeff = patch_icpp(patch_id)%smooth_coeff

        ! Initializing the pseudo volume fraction value to 1. The value will
        ! be modified as the patch is laid out on the grid, but only in the
        ! case that smoothing of the circular patch's boundary is enabled.
        eta = 1d0

        ! Checking whether the circle covers a particular cell in the domain
        ! and verifying whether the current patch has permission to write to
        ! that cell. If both queries check out, the primitive variables of
        ! the current patch are assigned to this cell.
        do j = 0, n
            do i = 0, m

                if (patch_icpp(patch_id)%smoothen) then

                    eta = tanh(smooth_coeff/min(dx, dy)* &
                               (sqrt((x_cc(i) - x_centroid)**2 &
                                     + (y_cc(j) - y_centroid)**2) &
                                - radius))*(-0.5d0) + 0.5d0

                end if

                if (((x_cc(i) - x_centroid)**2 &
                     + (y_cc(j) - y_centroid)**2 <= radius**2 &
                     .and. &
                     patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j))) &
                    .or. &
                    patch_id_fp(i, j) == smooth_patch_id) &
                    then

                    call s_assign_patch_species_primitive_variables(patch_id, i, j)

                end if

            end do
        end do

    end subroutine s_circle ! ----------------------------------------------

    !>             The varcircle patch is a 2D geometry that may be used
        !!             . It  generatres an annulus
        !! @param patch_id is the patch identifier
    subroutine s_varcircle(patch_id) ! ----------------------------------------

        ! Patch identifier
        integer, intent(IN) :: patch_id

        ! Generic loop iterators
        integer :: i, j

        real(kind(0d0)) :: myr, thickness

        ! Transferring the circular patch's radius, centroid, smearing patch
        ! identity and smearing coefficient information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        radius = patch_icpp(patch_id)%radius
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
        smooth_coeff = patch_icpp(patch_id)%smooth_coeff
        thickness = patch_icpp(patch_id)%epsilon

        ! Initializing the pseudo volume fraction value to 1. The value will
        ! be modified as the patch is laid out on the grid, but only in the
        ! case that smoothing of the circular patch's boundary is enabled.
        eta = 1d0

        ! Checking whether the circle covers a particular cell in the domain
        ! and verifying whether the current patch has permission to write to
        ! that cell. If both queries check out, the primitive variables of
        ! the current patch are assigned to this cell.
        do j = 0, n
            do i = 0, m
                myr = dsqrt((x_cc(i) - x_centroid)**2 &
                            + (y_cc(j) - y_centroid)**2)

                if (myr <= radius + thickness/2.d0 .and. &
                    myr >= radius - thickness/2.d0 .and. &
                    patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j))) then

                    call s_assign_patch_species_primitive_variables(patch_id, i, j)
                end if

            end do
        end do

    end subroutine s_varcircle ! ----------------------------------------------


    !>      The elliptical patch is a 2D geometry. The geometry of
        !!      the patch is well-defined when its centroid and radii
        !!      are provided. Note that the elliptical patch DOES allow
        !!      for the smoothing of its boundary
        !! @param patch_id is the patch identifier
    subroutine s_ellipse(patch_id) ! ---------------------------------------

        integer, intent(IN) :: patch_id

        integer :: i, j !< Generic loop operators

        ! Transferring the elliptical patch's radii, centroid, smearing
        ! patch identity, and smearing coefficient information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        a = patch_icpp(patch_id)%radii(1)
        b = patch_icpp(patch_id)%radii(2)
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
        smooth_coeff = patch_icpp(patch_id)%smooth_coeff

        ! Initializing the pseudo volume fraction value to 1. The value
        ! be modified as the patch is laid out on the grid, but only in
        ! the case that smoothing of the elliptical patch's boundary is
        ! enabled.
        eta = 1d0

        ! Checking whether the ellipse covers a particular cell in the
        ! domain and verifying whether the current patch has permission
        ! to write to that cell. If both queries check out, the primitive
        ! variables of the current patch are assigned to this cell.
        do j = 0, n
            do i = 0, m

                if (patch_icpp(patch_id)%smoothen) then
                    eta = tanh(smooth_coeff/min(dx, dy)* &
                               (sqrt(((x_cc(i) - x_centroid)/a)**2 + &
                                     ((y_cc(j) - y_centroid)/b)**2) &
                                - 1d0))*(-0.5d0) + 0.5d0
                end if

                if ((((x_cc(i) - x_centroid)/a)**2 + &
                     ((y_cc(j) - y_centroid)/b)**2 <= 1d0 &
                     .and. &
                     patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j))) &
                    .or. &
                    patch_id_fp(i, j) == smooth_patch_id) &
                    then

                    call s_assign_patch_species_primitive_variables(patch_id, i, j)
                end if
            end do
        end do

    end subroutine s_ellipse ! ---------------------------------------------


    !>      The rectangular patch is a 2D geometry that may be used,
        !!              for example, in creating a solid boundary, or pre-/post-
        !!              shock region, in alignment with the axes of the Cartesian
        !!              coordinate system. The geometry of such a patch is well-
        !!              defined when its centroid and lengths in the x- and y-
        !!              coordinate directions are provided. Please note that the
        !!              rectangular patch DOES NOT allow for the smoothing of its
        !!              boundaries.
        !! @param patch_id is the patch identifier
    subroutine s_rectangle(patch_id) ! -------------------------------------

        integer, intent(IN) :: patch_id

        real(kind(0d0)) :: pi_inf, gamma, lit_gamma !< Equation of state parameters

        integer :: i, j !< generic loop iterators

        pi_inf = fluid_pp(1)%pi_inf
        gamma = fluid_pp(1)%gamma
        lit_gamma = (1d0 + gamma)/gamma

        ! Transferring the rectangle's centroid and length information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        length_x = patch_icpp(patch_id)%length_x
        length_y = patch_icpp(patch_id)%length_y

        ! Computing the beginning and the end x- and y-coordinates of the
        ! rectangle based on its centroid and lengths
        x_boundary%beg = x_centroid - 0.5d0*length_x
        x_boundary%end = x_centroid + 0.5d0*length_x
        y_boundary%beg = y_centroid - 0.5d0*length_y
        y_boundary%end = y_centroid + 0.5d0*length_y

        ! Since the rectangular patch does not allow for its boundaries to
        ! be smoothed out, the pseudo volume fraction is set to 1 to ensure
        ! that only the current patch contributes to the fluid state in the
        ! cells that this patch covers.
        eta = 1d0

        ! Checking whether the rectangle covers a particular cell in the
        ! domain and verifying whether the current patch has the permission
        ! to write to that cell. If both queries check out, the primitive
        ! variables of the current patch are assigned to this cell.
        do j = 0, n
            do i = 0, m
                if (x_boundary%beg <= x_cc(i) .and. &
                    x_boundary%end >= x_cc(i) .and. &
                    y_boundary%beg <= y_cc(j) .and. &
                    y_boundary%end >= y_cc(j) &
                    .and. &
                    patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j))) &
                    then

                    call s_assign_patch_species_primitive_variables(patch_id, i, j)

                end if
            end do
        end do

    end subroutine s_rectangle ! -------------------------------------------

    !>  The swept line patch is a 2D geometry that may be used,
        !!      for example, in creating a solid boundary, or pre-/post-
        !!      shock region, at an angle with respect to the axes of the
        !!      Cartesian coordinate system. The geometry of the patch is
        !!      well-defined when its centroid and normal vector, aimed
        !!      in the sweep direction, are provided. Note that the sweep
        !!      line patch DOES allow the smoothing of its boundary.
        !! @param patch_id is the patch identifier
    subroutine s_sweep_line(patch_id) ! ------------------------------------

        integer, intent(IN) :: patch_id

        integer :: i, j !< Generic loop operators

        ! Transferring the centroid information of the line to be swept
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
        smooth_coeff = patch_icpp(patch_id)%smooth_coeff

        ! Obtaining coefficients of the equation describing the sweep line
        a = patch_icpp(patch_id)%normal(1)
        b = patch_icpp(patch_id)%normal(2)
        c = -a*x_centroid - b*y_centroid

        ! Initializing the pseudo volume fraction value to 1. The value will
        ! be modified as the patch is laid out on the grid, but only in the
        ! case that smoothing of the sweep line patch's boundary is enabled.
        eta = 1d0

        ! Checking whether the region swept by the line covers a particular
        ! cell in the domain and verifying whether the current patch has the
        ! permission to write to that cell. If both queries check out, the
        ! primitive variables of the current patch are written to this cell.
        do j = 0, n
            do i = 0, m

                if (patch_icpp(patch_id)%smoothen) then
                    eta = 5d-1 + 5d-1*tanh(smooth_coeff/min(dx, dy) &
                                           *(a*x_cc(i) + b*y_cc(j) + c) &
                                           /sqrt(a**2 + b**2))
                end if

                if ((a*x_cc(i) + b*y_cc(j) + c >= 0d0 &
                     .and. &
                     patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j))) &
                    .or. &
                    patch_id_fp(i, j) == smooth_patch_id) &
                    then
                    call s_assign_patch_species_primitive_variables(patch_id, i, j)

                end if

            end do
        end do

    end subroutine s_sweep_line ! ------------------------------------------

    !> The isentropic vortex is a 2D geometry that may be used,
        !!              for example, to generate an isentropic flow disturbance.
        !!              Geometry of the patch is well-defined when its centroid
        !!              and radius are provided. Notice that the patch DOES NOT
        !!              allow for the smoothing of its boundary.
        !! @param patch_id is the patch identifier
    subroutine s_isentropic_vortex(patch_id) ! ----------------------------

        ! Patch identifier
        integer, intent(IN) :: patch_id

        ! Generic loop iterators
        integer :: i, j

        ! Transferring isentropic vortex patch's centroid and radius info
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        radius = patch_icpp(patch_id)%radius

        ! Since the isentropic vortex patch does not allow for its boundary
        ! to get smoothed, the pseudo volume fraction is set to 1 to ensure
        ! that only the current patch contributes to the fluid state in the
        ! cells that this patch covers.
        eta = 1d0

        ! Verifying whether the isentropic vortex includes a particular cell
        ! and verifying whether the current patch has permission to write to
        ! that cell. If both queries work out the primitive variables of the
        ! the current patch are assigned to this cell.
        do j = 0, n
            do i = 0, m

                if ((x_cc(i) - x_centroid)**2 &
                    + (y_cc(j) - y_centroid)**2 <= radius**2 &
                    .and. &
                    patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j))) &
                    then

                    call s_assign_patch_species_primitive_variables(patch_id, &
                                                            i, j)

                end if

            end do
        end do

    end subroutine s_isentropic_vortex ! -----------------------------------

    !>  This patch assigns the primitive variables as analytical
        !!  functions such that the code can be verified.
        !!  @param patch_id is the patch identifier
    subroutine s_1D_analytical(patch_id) ! ---------------------------------

        ! Patch identifier
        integer, intent(IN) :: patch_id

        ! Placeholders for the cell boundary values
        real(kind(0d0)) :: a, b, c, d, pi_inf, gamma, lit_gamma

        ! Generic loop iterators
        integer :: i, j

        pi_inf = fluid_pp(1)%pi_inf
        gamma = fluid_pp(1)%gamma
        lit_gamma = (1d0 + gamma)/gamma

        ! Transferring the patch's centroid and length information
        x_centroid = patch_icpp(patch_id)%x_centroid
        length_x = patch_icpp(patch_id)%length_x

        ! Computing the beginning and the end x- and y-coordinates
        ! of the patch based on its centroid and lengths
        x_boundary%beg = x_centroid - 0.5d0*length_x
        x_boundary%end = x_centroid + 0.5d0*length_x

        ! Since the patch doesn't allow for its boundaries to be
        ! smoothed out, the pseudo volume fraction is set to 1 to
        ! ensure that only the current patch contributes to the fluid
        ! state in the cells that this patch covers.
        eta = 1d0

        ! Checking whether the line segment covers a particular cell in the
        ! domain and verifying whether the current patch has the permission
        ! to write to that cell. If both queries check out, the primitive
        ! variables of the current patch are assigned to this cell.
        do i = 0, m
            if (x_boundary%beg <= x_cc(i) .and. &
                x_boundary%end >= x_cc(i) .and. &
                patch_icpp(patch_id)%alter_patch(patch_id_fp(i, 0))) then

                call s_assign_patch_species_primitive_variables(patch_id, i, 0)

                !what variables to alter
                !bump in pressure
                q_prim_vf(E_idx)%sf(i, 0) = q_prim_vf(E_idx)%sf(i, 0)* &
                                               (1d0 + 0.2d0*dexp(-1d0*((x_cb(i) - x_centroid)**2.d0)/(2.d0*0.005d0)))

            end if
        end do

    end subroutine s_1D_analytical ! ---------------------------------------


    !>  This patch assigns the primitive variables as analytical
        !!  functions such that the code can be verified.
        !!  @param patch_id is the patch identifier
    subroutine s_2D_analytical(patch_id) ! ---------------------------------

        integer, intent(IN) :: patch_id

        real(kind(0d0)) :: a, b, c, d !< placeholderrs for the cell boundary values
        real(kind(0d0)) :: pi_inf, gamma, lit_gamma !< equation of state parameters

        integer :: i, j !< generic loop iterators

        pi_inf = fluid_pp(1)%pi_inf
        gamma = fluid_pp(1)%gamma
        lit_gamma = (1d0 + gamma)/gamma

        ! Transferring the patch's centroid and length information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        length_x = patch_icpp(patch_id)%length_x
        length_y = patch_icpp(patch_id)%length_y

        ! Computing the beginning and the end x- and y-coordinates
        ! of the patch based on its centroid and lengths
        x_boundary%beg = x_centroid - 0.5d0*length_x
        x_boundary%end = x_centroid + 0.5d0*length_x
        y_boundary%beg = y_centroid - 0.5d0*length_y
        y_boundary%end = y_centroid + 0.5d0*length_y

        ! Since the patch doesn't allow for its boundaries to be
        ! smoothed out, the pseudo volume fraction is set to 1 to
        ! ensure that only the current patch contributes to the fluid
        ! state in the cells that this patch covers.
        eta = 1d0

        ! Checking whether the patch covers a particular cell in the
        ! domain and verifying whether the current patch has the
        ! permission to write to that cell. If both queries check out,
        ! the primitive variables of the current patch are assigned
        ! to this cell.
        do j = 0, n
            do i = 0, m
                if (x_boundary%beg <= x_cc(i) .and. &
                    x_boundary%end >= x_cc(i) .and. &
                    y_boundary%beg <= y_cc(j) .and. &
                    y_boundary%end >= y_cc(j) .and. &
                    patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j))) then

                    call s_assign_patch_species_primitive_variables(patch_id, i, j)

                    !what variables to alter
                    !x-y bump in pressure
                    q_prim_vf(E_idx)%sf(i, j) = q_prim_vf(E_idx)%sf(i, j)* &
                               (1d0 + 0.2d0*dexp(-1d0*((x_cb(i) - x_centroid)**2.d0 + (y_cb(j) - y_centroid)**2.d0)/(2.d0*0.005d0)))

                end if
            end do
        end do

    end subroutine s_2D_analytical ! ---------------------------------------

end module m_create_patches
