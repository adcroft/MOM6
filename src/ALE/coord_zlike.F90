! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Regrid columns for a z-like coordinate (z-star, z-level)
module coord_zlike

use MOM_error_handler, only : MOM_error, FATAL
use numerical_testing_type, only : testing

implicit none ; private

!> Control structure containing required parameters for a z-like coordinate
type, public :: zlike_CS ; private

  !> Number of levels to be generated
  integer :: nk

  !> Minimum thickness allowed for layers, in the same thickness units (perhaps [H ~> m or kg m-2])
  !! that will be used in all subsequent calls to build_zstar_column with this structure.
  real :: min_thickness

  !> Target coordinate resolution, usually in [Z ~> m]
  real, allocatable, dimension(:) :: coordinateResolution
end type zlike_CS

public init_coord_zlike, set_zlike_params, end_coord_zlike
public build_zstar_open_ocean_column, build_zstar_column
public coord_zlike_unit_tests

contains

!> Initialise a zlike_CS with pointers to parameters
subroutine init_coord_zlike(CS, nk, coordinateResolution, z_scale)
  type(zlike_CS),     pointer    :: CS !< Unassociated pointer to hold the control structure
  integer,            intent(in) :: nk !< Number of levels in the grid
  real, dimension(:), intent(in) :: coordinateResolution !< Target coordinate resolution [Z ~> m]
  real, optional,     intent(in) :: z_scale !< Scaling factor from the target coordinate resolution
                                            !! in Z to desired units for zInterface, perhaps Z_to_H,
                                            !! often [nondim] or [H Z-1 ~> 1 or kg m-3]

  if (associated(CS)) call MOM_error(FATAL, "init_coord_zlike: CS already associated!")
  allocate(CS)
  allocate(CS%coordinateResolution(nk))

  CS%nk                   = nk
  CS%coordinateResolution = coordinateResolution
  if (present(z_scale)) CS%coordinateResolution = CS%coordinateResolution * z_scale
end subroutine init_coord_zlike

!> Deallocates the zlike control structure
subroutine end_coord_zlike(CS)
  type(zlike_CS), pointer :: CS !< Coordinate control structure

  ! Nothing to do
  if (.not. associated(CS)) return
  deallocate(CS%coordinateResolution)
  deallocate(CS)
end subroutine end_coord_zlike

!> Set parameters in the zlike structure
subroutine set_zlike_params(CS, min_thickness)
  type(zlike_CS), pointer    :: CS !< Coordinate control structure
  real, optional, intent(in) :: min_thickness !< Minimum allowed thickness [H ~> m or kg m-2]

  if (.not. associated(CS)) call MOM_error(FATAL, "set_zlike_params: CS not associated")

  if (present(min_thickness)) CS%min_thickness = min_thickness
end subroutine set_zlike_params

!> Builds a z* coordinate for the open ocean with a minimum thickness
!! for vanished layers
subroutine build_zstar_open_ocean_column(CS, depth, total_thickness, zInterface)
  type(zlike_CS),           intent(in)    :: CS !< Coordinate control structure
  real,                     intent(in)    :: depth !< Depth of ocean bottom (positive downward in the
                                                   !! output units), units may be [Z ~> m] or [H ~> m or kg m-2]
  real,                     intent(in)    :: total_thickness !< Column thickness (positive definite in the same
                                                   !! units as depth) [Z ~> m] or [H ~> m or kg m-2]
  real, dimension(CS%nk+1), intent(inout) :: zInterface !< Absolute positions of interfaces (in the same
                                                   !! units as depth) [Z ~> m] or [H ~> m or kg m-2]
                                                   !! Note: incoming values are NOT used.
  ! Local variables
  real :: stretching ! A stretching factor for the coordinate [nondim]
  real :: dh, min_thickness, z0_top, z_star ! Thicknesses or heights [Z ~> m] or [H ~> m or kg m-2]
  integer :: k

  min_thickness = min( CS%min_thickness, total_thickness/real(CS%nk) )

  ! Conventional z* coordinate:
  !   z* = (z-eta) / stretching   where stretching = (H+eta)/H
  !   z = eta + stretching * z*
  ! The above gives z*(z=eta) = 0, z*(z=-H) = -H.
  stretching = total_thickness / depth

  ! Integrate down from the top for a notional new grid, ignoring topography
  ! The starting position is offset by z0_top which, if z0_top<0, will place
  ! interfaces above the rigid boundary.
  zInterface(1) = total_thickness - depth
  do k = 1,CS%nk
    dh = stretching * CS%coordinateResolution(k) ! Notional grid spacing
    zInterface(k+1) = zInterface(k) - dh
  enddo

  ! Integrating up from the bottom adjusting interface position to accommodate
  ! inflating layers without disturbing the interface above
  zInterface(CS%nk+1) = -depth
  do k = CS%nk,1,-1
    if ( zInterface(k) < (zInterface(k+1) + min_thickness) ) then
      zInterface(k) = zInterface(k+1) + min_thickness
    endif
  enddo

end subroutine build_zstar_open_ocean_column

!> Builds a z* coordinate with a minimum thickness
subroutine build_zstar_column(CS, depth, total_thickness, zInterface, &
                              z_rigid_top, eta_orig)
  type(zlike_CS),           intent(in)    :: CS !< Coordinate control structure
  real,                     intent(in)    :: depth !< Depth of ocean bottom (positive downward in the
                                                   !! output units), units may be [Z ~> m] or [H ~> m or kg m-2]
  real,                     intent(in)    :: total_thickness !< Column thickness (positive definite in the same
                                                   !! units as depth) [Z ~> m] or [H ~> m or kg m-2]
  real, dimension(CS%nk+1), intent(inout) :: zInterface !< Absolute positions of interfaces (in the same
                                                   !! units as depth) [Z ~> m] or [H ~> m or kg m-2]
                                                   !! Note: incoming values are NOT used.
  real,                     intent(in)    :: z_rigid_top !< The height of a rigid top (positive upward in the same
                                                   !! units as depth) [Z ~> m] or [H ~> m or kg m-2]
  real,                     intent(in)    :: eta_orig !< The actual original height of the top (in the same
                                                   !! units as depth) [Z ~> m] or [H ~> m or kg m-2]
  ! Local variables
  real :: stretching ! A stretching factor for the coordinate [nondim]
  real :: dh, min_thickness, z0_top, z_star ! Thicknesses or heights [Z ~> m] or [H ~> m or kg m-2]
  integer :: k

  min_thickness = min( CS%min_thickness, total_thickness/real(CS%nk) )
  z0_top = z_rigid_top

  ! With a rigid top boundary at eta = z0_top then
  !   z* = z0 + (z-eta) / stretching   where stretching = (H+eta)/(H+z0)
  !   z = eta + stretching * (z*-z0) * stretching
  stretching = total_thickness / ( depth + z0_top )

  ! z_star is the notional z* coordinate in absence of upper/lower topography
  z_star = 0. ! z*=0 at the free-surface
  zInterface(1) = eta_orig ! The actual position of the top of the column
  do k = 2,CS%nk
    z_star = z_star - CS%coordinateResolution(k-1)
    ! This ensures that z is below a rigid upper surface (ice shelf bottom)
    zInterface(k) = min( eta_orig + stretching * ( z_star - z0_top ), z0_top )
    ! This ensures that the layer in inflated
    zInterface(k) = min( zInterface(k), zInterface(k-1) - min_thickness )
    ! This ensures that z is above or at the topography
    zInterface(k) = max( zInterface(k), -depth + real(CS%nk+1-k) * min_thickness )
  enddo
  zInterface(CS%nk+1) = -depth

end subroutine build_zstar_column

!> Runs unit tests and returns True for any fails, False otherwise
logical function coord_zlike_unit_tests(verbose)
  logical, intent(in)    :: verbose !< True, if verbose
  ! Local variables
  type(testing) :: test ! numerical_testing_type
  type(zlike_CS), pointer :: CS ! Coordinate control structure
  real :: zi(5) ! Interface positions [A]

  call test%set( verbose=verbose )

  call init_coord_zlike(CS, 4, [1., 3., 5., 7.]) ! Nominal z-grid, 16 deep

  ! Using just "depth" and "total thickness" to generate grid acts like z* in
  ! open ocean, and like top-down sigma-coordinates under an ice shelf.
  !
  ! In this mode, the nominal dz is clipped, or expanded, at the bottom depth
  ! and then the result is stretched to have the given total thickness.
  ! There is no clipping at the top by an ice shelf, which is why it resorts
  ! to a sigma-coordinate like behavior in ice shelf cavities.

  call build_zstar_open_ocean_column(CS, 16., & ! depth
                              16., & ! total thickness
                              zi )   ! current and future interface positions
  call test%real_arr(5, zi, (/0., -1., -4., -9., -16./), 'zi matches nominal grid')

  call build_zstar_open_ocean_column(CS, 16., & ! depth
                              32., & ! total thickness
                              zi )   ! current and future interface positions
  call test%real_arr(5, zi, (/16., 14., 8., -2., -16./), 'Double total thickness')

  call build_zstar_open_ocean_column(CS, 16., & ! depth
                               8., & ! total thickness
                              zi )   ! current and future interface positions
  call test%real_arr(5, zi, (/-8., -8.5, -10., -12.5, -16./), 'Half total thickness')

  call build_zstar_open_ocean_column(CS, 10., & ! depth
                              10., & ! total thickness
                              zi )   ! current and future interface positions
  call test%real_arr(5, zi, (/0., -1., -4., -9., -10./), 'Shallower bottom')

  call build_zstar_open_ocean_column(CS, 10., & ! depth
                              20., & ! total thickness
                              zi )   ! current and future interface positions
  call test%real_arr(5, zi, (/10., 8., 2., -8., -10./), 'Shallower bottom, double thickness')

  call build_zstar_open_ocean_column(CS, 18., & ! depth
                              18., & ! total thickness
                              zi )   ! current and future interface positions
  call test%real_arr(5, zi, (/0., -1., -4., -9., -18./), 'Deeper bottom')

  call build_zstar_open_ocean_column(CS, 18., & ! depth
                              36., & ! total thickness
                              zi )   ! current and future interface positions
  call test%real_arr(5, zi, (/18., 16., 10., 0., -18./), 'Deeper bottom, double thickness')

  ! Using all of "depth", "total thickness", "z_rigid_top" and "eta_orig" to generate grid
  ! then clipping is as expected when depressed
  call build_zstar_column(CS, 16., & ! depth
                              13., & ! total thickness
                              zi,  & ! current and future interface positions
                              -3., & ! z_rigid_top
                              -3.)
  call test%real_arr(5, zi, (/-3., -3., -4., -9., -16./), 'All args are consistent w. clipped')

  ! Using all of "depth", "total thickness", "z_rigid_top" and "eta_orig" to generate grid
  ! then surplus volume inflates only the first layer
  call build_zstar_column(CS, 16., & ! depth
                              24., & ! total thickness
                              zi,  & ! current and future interface positions
                              8.,  & ! z_rigid_top
                              8.)    ! eta_orig
  call test%real_arr(5, zi, (/8., -1., -4., -9., -16./), 'All args are consistent w. inflation')

  ! If z_rigid_top is consistent with total_thickness, eta_top changes the top
  ! position, with nominal thicknesses apparently applied unclipped.
  ! Yields thicknesses of 1, 3, 5, 8
  call build_zstar_column(CS, 16., & ! depth
                              16., & ! total thickness
                              zi,  & ! current and future interface positions
                              0.,  & ! z_rigid_top
                              1.)    ! eta_orig
  call test%real_arr(5, zi, (/1., 0., -3., -8., -16./), 'Eta is inconsistent > top')

  ! If z_rigid_top is consistent with total_thickness, eta_top changes the top
  ! position, with nominal thicknesses apparently applied unclipped.
  ! Yields thicknesses of 1, 3, 5, 6
  call build_zstar_column(CS, 16., & ! depth
                              16., & ! total thickness
                              zi,  & ! current and future interface positions
                              0.,  & ! z_rigid_top
                              -1.)   ! eta_orig
  call test%real_arr(5, zi, (/-1., -2., -5., -10., -16./), 'Eta is inconsistent < top')

  ! If eta_top is consistent with total_thickness, z_rigid_top changes the top
  call build_zstar_column(CS, 16., & ! depth
                              32., & ! total thickness
                              zi,  & ! current and future interface positions
                              0.,  & ! z_rigid_top
                              16.)   ! eta_orig
  ! Yields thicknesses of 16, 0, 2, 14
  call test%real_arr(5, zi, (/16., 0., 0., -2., -16./), 'rigid top is inconsistent < eta', tol=2.e-16)

  ! If eta_top is consistent with total_thickness, z_rigid_top changes the top
  call build_zstar_column(CS, 16., & ! depth
                              16., & ! total thickness
                              zi,  & ! current and future interface positions
                              16., & ! z_rigid_top
                              0.)    ! eta_orig
  ! Yields thicknesses of 8.5, 1.5, 2.5, 3.5
  call test%real_arr(5, zi, (/0., -8.5, -10., -12.5, -16./), 'rigid top is inconsistent > eta')

  coord_zlike_unit_tests = test%summarize('coord_zlike_unit_tests') ! Return true if a fail occurred

end function coord_zlike_unit_tests

end module coord_zlike
