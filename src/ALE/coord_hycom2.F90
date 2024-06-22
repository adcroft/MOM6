!> Regrid columns for the HyCOM coordinate
module coord_hycom2

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL
use MOM_remapping,     only : remapping_CS, remapping_core_h
use MOM_EOS,           only : EOS_type, calculate_density
use regrid_interp,     only : interp_CS_type, build_and_interpolate_grid, regridding_set_ppolys
use regrid_interp,     only : DEGREE_MAX

use iso_fortran_env,  only : stdout=>output_unit, stderr=>error_unit

implicit none ; private

!> Control structure containing required parameters for the HyCOM coordinate
type, public :: hycom2_CS ; private

  !> Number of layers/levels in generated grid
  integer :: nk

  !> Nominal near-surface resolution [Z ~> m]
  real, allocatable, dimension(:) :: coordinateResolution

  !> Nominal density of interfaces [R ~> kg m-3]
  real, allocatable, dimension(:) :: target_density

  !> Maximum depths of interfaces [H ~> m or kg m-2]
  real, allocatable, dimension(:) :: max_interface_depths

  !> Maximum thicknesses of layers [H ~> m or kg m-2]
  real, allocatable, dimension(:) :: max_layer_thickness

  !> If true, an interface only moves if it improves the density fit
  logical :: only_improves = .false.

  !> Interpolation control structure
  type(interp_CS_type) :: interp_CS
end type hycom2_CS

public init_coord_hycom2, set_hycom2_params, build_hycom2_column, end_coord_hycom2
public coord_hycom2_unit_tests

contains

!> Initialise a hycom2_CS with pointers to parameters
subroutine init_coord_hycom2(CS, nk, coordinateResolution, target_density, interp_CS)
  type(hycom2_CS),      pointer    :: CS !< Unassociated pointer to hold the control structure
  integer,              intent(in) :: nk !< Number of layers in generated grid
  real, dimension(nk),  intent(in) :: coordinateResolution !< Nominal near-surface resolution [Z ~> m]
  real, dimension(nk+1),intent(in) :: target_density !< Interface target densities [R ~> kg m-3]
  type(interp_CS_type), intent(in) :: interp_CS !< Controls for interpolation

  if (associated(CS)) call MOM_error(FATAL, "init_coord_hycom2: CS already associated!")
  allocate(CS)
  allocate(CS%coordinateResolution(nk))
  allocate(CS%target_density(nk+1))

  CS%nk                      = nk
  CS%coordinateResolution(:) = coordinateResolution(:)
  CS%target_density(:)       = target_density(:)
  CS%interp_CS               = interp_CS

end subroutine init_coord_hycom2

!> This subroutine deallocates memory in the control structure for the coord_hycom module
subroutine end_coord_hycom2(CS)
  type(hycom2_CS), pointer :: CS !< Coordinate control structure

  ! nothing to do
  if (.not. associated(CS)) return
  deallocate(CS%coordinateResolution)
  deallocate(CS%target_density)
  if (allocated(CS%max_layer_thickness)) deallocate(CS%max_layer_thickness)
  deallocate(CS)
end subroutine end_coord_hycom2

!> This subroutine can be used to set the parameters for the coord_hycom module
subroutine set_hycom2_params(CS, max_layer_thickness, only_improves, interp_CS)
  type(hycom2_CS),                 pointer    :: CS !< Coordinate control structure
  real, dimension(:),   optional, intent(in) :: max_layer_thickness  !< Maximum thicknesses of layers [H ~> m or kg m-2]
  logical, optional, intent(in) :: only_improves !< If true, an interface only moves if it improves the density fit
  type(interp_CS_type), optional, intent(in) :: interp_CS !< Controls for interpolation

  if (.not. associated(CS)) call MOM_error(FATAL, "set_hycom2_params: CS not associated")

  if (present(max_layer_thickness)) then
    if (size(max_layer_thickness) /= CS%nk) &
      call MOM_error(FATAL, "set_hycom2_params: max_layer_thickness inconsistent size")
    allocate(CS%max_layer_thickness(CS%nk))
    CS%max_layer_thickness(:) = max_layer_thickness(:)
  endif

  if (present(only_improves)) CS%only_improves = only_improves

  if (present(interp_CS)) CS%interp_CS = interp_CS
end subroutine set_hycom2_params

!> Build a HyCOM coordinate column
subroutine build_hycom2_column(CS, remapCS, eqn_of_state, nz, h, T, S, h_new)
  type(hycom2_CS),       intent(in)    :: CS    !< Coordinate control structure
  type(remapping_CS),    intent(in)    :: remapCS !< Remapping parameters and options
  type(EOS_type),        intent(in)    :: eqn_of_state !< Equation of state structure
  integer,               intent(in)    :: nz    !< Number of levels
  real, dimension(nz),   intent(in)    :: T     !< Temperature of column [C ~> degC]
  real, dimension(nz),   intent(in)    :: S     !< Salinity of column [S ~> ppt]
  real, dimension(nz),   intent(in)    :: h     !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(nz),   intent(in)    :: h_new !< New thicknesses [H ~> m or kg m-2]

  ! Local variables
  integer   :: k
  real, dimension(nz)      :: rho_col   ! Layer densities in a column [R ~> kg m-3]
  real, dimension(CS%nk)   :: h_col_new ! New layer thicknesses [H ~> m or kg m-2]
  real, dimension(CS%nk)   :: r_col_new ! New layer densities [R ~> kg m-3]
  real, dimension(CS%nk)   :: T_col_new ! New layer temperatures [C ~> degC]
  real, dimension(CS%nk)   :: S_col_new ! New layer salinities [S ~> ppt]
  real, dimension(CS%nk)   :: p_col_new ! New layer pressure [R L2 T-2 ~> Pa]
  real, dimension(CS%nk+1) :: RiA_ini   ! Initial nk+1 interface density anomaly w.r.t. the
                                        ! interface target densities [R ~> kg m-3]
  real, dimension(CS%nk+1) :: RiA_new   ! New interface density anomaly w.r.t. the
                                        ! interface target densities [R ~> kg m-3]
  real :: z_1, z_nz  ! mid point of 1st and last layers [H ~> m or kg m-2]
  real :: z_scale    ! A scaling factor from the input thicknesses to the target thicknesses,
                     ! perhaps 1 or a factor in [H Z-1 ~> 1 or kg m-3]
  real :: stretching ! z* stretching, converts z* to z [nondim].
  real :: nominal_z ! Nominal depth of interface when using z* [H ~> m or kg m-2]
  logical :: maximum_h_set      ! If true, the maximum layer thicknesses have been set.


! ! Reconstruct T/S for all layers
! TR = remapCS%fit( nz, h, T )
! SR = remapCS%fit( nz, h, S )
!
! ! Start at top
! z_old = 0.
! k_old = 1

! ! For each interface (always using the reference pressure for the interface)
! do K_new = 2, CS%nk+1
!   p_ref = CS%p_ref(K_new)
!
!   ! Search down the source column to find the layer that contains this density of this interface
!   in_layer_k_old = .false.
!   h_new(K_new-1) = 0.
!   scan_k_old: do while (.not. in_layer_k_old)

!     ! Calculate the density of the top of the current source layer
!     S_top = SR(k_old)%eval(0.)
!     T_top = TR(k_old)%eval(0.)
!     rho_top = rho( S_top, T_top, p_ref )
!     if (CS$target_density(K_new) < rho_top) then
!       ! The target is lighter than the top of this layer so we do not move z_new
!       in_layer_k_old = .true.
!     else
!       ! Calculate the density of the bottom of the current source layer
!       S_bot = SR(k_old)%eval(1.)
!       T_bot = TR(k_old)%eval(1.)
!       rho_bot = rho( S_bot, T_bot, p_ref )
!       if (CS$target_density(K_new) >= rho_bot) then
!         ! The target is denser than the bottom of this layer so we will move to the next source layer
!         ! and accumulate the thickness
!         h_new(K_new-1) = h_new(K_new-1) + h(k_old)
!       else
!         ! The top/bot densities of layer k_old span the target
!         in_layer_k_old = .true. ! Combined with h_new(K_new-1)=0 says point to the top of the layer
!         ! Interpolate for the position of the target within the layer.
!         ! We start with a linear interpolation.
!         drhol = rho_bot - rho_top ! Density difference across the layer
!         if (drhol<=0.) then
!           dhp = 1. ! For a PCM or inverted layer, we snap to the bottom
!         else
!           drhop = CS$target_density(K_new) - rho_top
!           dhp = min(drhop / drhol, 1.) ! The min() here might not be necessary; just being cautious. -aja

!           ! Now refine dhp
!           S_mid = SR(k_old)%eval(dhp)
!           T_mid = TR(k_old)%eval(dhp)
!           rho_mid = rho( S_mid, T_mid, p_ref )
!           a = ( rho_top - rho_mid ) + ( rho_bot - rho_mid )
!           b = 2. * ( rho_mid - rho_bot )
!           c = drhop
!           d = sqrt( 0.25*b*b - a*c )
!           if (a > 0.) then
!             ! we know the value is to the right of the mid point
!           elseif (a < 0.) then
!             ! we know the value is to the left of the mid point
!           else ! a=0
!           endif


!         endif


!         h_new(K_new-1) = h_new(K_new-1) + h(k_old) * dhp
!       endif


!     endif
!   enddo ! scan_k_old

! enddo

end subroutine build_hycom2_column

!> Return solution of f(xs) = ys where f(xs) is a parabola that
!! passes through (x1,y1), (0,0), and (x2,y2).
!!
!! x1,0,x2 and y1,0,y2 must be distinct and increasing.
!! ys should be in the range(y1..y2).
real function qroot(x1, x2, y1, y2, ys)
!real pure function qroot(x1, x2, y1, y2, ys)
  real, intent(in) :: x1 !< Position of y1 (required to be negative)
  real, intent(in) :: x2 !< Position of y2 (required to be positive)
  real, intent(in) :: y1 !< Value of f(x) at x1 (required to be negative)
  real, intent(in) :: y2 !< Value of f(x) at y1 (required to be positive)
  real, intent(in) :: ys !< Value of f(x) to solve for (must be in range y1..y2)
  ! Local variables
  real :: a,b ! Polynomial coefficients
  real :: xm ! Position of parabola extremum
  real :: delta ! Distance between extremum and roots
  real :: d ! determinant (or reciprocal) of 2x2 matrix for a and b
  real :: e ! common term in series
  real :: ra, rb ! 1/a and 1/b
  real, parameter :: tt = 1./8192. ! threshold of series truncation term

  ! Determine the coefficients of the parabola f(x) = a x^2 + b x
  d = 1. / ( ( x1 * x2 ) * ( x1 - x2 ) )
  a = d * ( x2 * y1 - x1 * y2 )
  b = d * ( x1**2 * y2 - x2**2 * y1 )
  ! We only ever need 1/b
  rb = 1. / b

  ! The roots of "a x^2 + b x - ys = 0" are "xm +/- delta" where
  !   xm = -b / 2 a
  !   delta^2 = xm^2 + ys/a
  ! if a>0 we need the right root, xm + delta
  ! if a<0 we need the left root, xm - delta
  ! if a<<1 there is an inaccuracy due to canceling large terms in
  ! the expression xm - sqrt( xm^2 + ys/a ) so we expand this as
  ! a series in the quantity "a ys / b^2".
  e = ( a * ys ) * b**2

  if (abs(e)**3 > epsilon(e)) then
    ! sqrt(1+e)-1 is accurate enough
    ra = 1. / a
    xm = -0.5 * ( b * ra )
    delta = sqrt( xm**2 + ys * ra )
    if (a>0.) then
      qroot = xm + delta
    else
      qroot = xm - delta
    endif
  else
    ! series for root is valid for -1<4e<1 and but truncated to
    ! four terms is precise enough only for 14 e^4 < epsilon:
    !   ys/b * ( 1 -e + 2e^2 -5e^3 +14e^4 - ...A )
    !qroot = (ys * rb) * ( 1. - e * ( 1. - e * (2. -  e * ( 5. - e * 14. ) ) ) )
    !qroot = (ys * rb) * ( 1. - e )
    qroot = (ys * rb) * ( 1. - e * ( 1. - e * (2. -  e * 5. ) ) )
  endif

  if (qroot<x1 .or. qroot>x2) then
    ! Fall back to linear interpolation
    if (ys>=0) then
      qroot = ( ys / y2 ) * x2
    else
      qroot = ( ys / y1 ) * x1
    endif
  endif
end function qroot

logical function coord_hycom2_unit_tests(verbose)
  logical, intent(in) :: verbose !< If true, display tests to stdout

  coord_hycom2_unit_tests = .false.

  ! Regardless of fit, 0 -> 0
  coord_hycom2_unit_tests = coord_hycom2_unit_tests .or. test_scalar(verbose, &
        qroot(-3., 5., -7., 11., 0.), 0., 'p(a,b)(0) = 0')

  ! The line "0 x^2 + x" passes through (-1,-1), (0,0), (3,3)
  coord_hycom2_unit_tests = coord_hycom2_unit_tests .or. test_scalar(verbose, &
        qroot(-1., 3., -1., 3., 1.), 1., 'p(0,1)(1) = 1')

  ! The parabola "x^2 + 2 x" passes through (-1,-1), (0,0), (2,8)
  ! Test left limit
  coord_hycom2_unit_tests = coord_hycom2_unit_tests .or. test_scalar(verbose, &
        qroot(-1., 2., -1., 8., -1.), -1., 'p(1,2)(-1) = -1')
  ! Test right limit
  coord_hycom2_unit_tests = coord_hycom2_unit_tests .or. test_scalar(verbose, &
        qroot(-1., 2., -1., 8., 8.), 2., 'p(1,2)(2) = 8')
  ! Test middle value
  coord_hycom2_unit_tests = coord_hycom2_unit_tests .or. test_scalar(verbose, &
        qroot(-1., 2., -1., 8., 3.), 1., 'p(1,2)(1) = 3')

  ! The parabola "-x^2 + 2 x" passes through (-1,-3), (0,0), (1,1)
  ! Test left limit
  coord_hycom2_unit_tests = coord_hycom2_unit_tests .or. test_scalar(verbose, &
        qroot(-2., 1., -8., 1., -8.), -2., 'p(-1,2)(-2) = -8')
  ! Test right limit
  coord_hycom2_unit_tests = coord_hycom2_unit_tests .or. test_scalar(verbose, &
        qroot(-2., 1., -8., 1., 1.), 1., 'p(-1,2)(1) = 1')
  ! Test middle value
  coord_hycom2_unit_tests = coord_hycom2_unit_tests .or. test_scalar(verbose, &
        qroot(-2., 1., -8., 1., -3.), -1., 'p(-1,2)(-1) = -3')

  ! A parabola that passes through (-1,-1), (0,0), (1,1+eps) has tiny curvature ~O(eps/2)
  ! eps ~ 1e-9 tests the series form
  coord_hycom2_unit_tests = coord_hycom2_unit_tests .or. test_scalar(verbose, &
        qroot(-1., 1., -1., 1.+2.e-9, -1.), -1., 'p(1e-10,1)(-1) = -1')

  ! A parabola that passes through (-1,-1), (0,0), (1,1-eps) has tiny curvature ~O(-eps/2)
  ! eps ~ 1e-9 tests the series form
  coord_hycom2_unit_tests = coord_hycom2_unit_tests .or. test_scalar(verbose, &
        qroot(-1., 1., -1., 1.-2.e-9, -1.), -1., 'p(-1e-10,1)(-1) = -1')

  ! A parabola that passes through (-1,-1), (0,0), (1,1+eps) has tiny curvature ~O(eps/2)
  ! eps ~ 1e-4 tests the series form that needs bounding
  coord_hycom2_unit_tests = coord_hycom2_unit_tests .or. test_scalar(verbose, &
        qroot(-1., 1., -1., 1.+1.e-5, -1.), -1., 'p(5e-6,1)(-1) = -1')

  ! A parabola that passes through (-1,-1-eps), (0,0), (1,1) has tiny curvature ~O(-eps/2)
  ! eps ~ 1e-4 tests the series form that needs bounding
  coord_hycom2_unit_tests = coord_hycom2_unit_tests .or. test_scalar(verbose, &
        qroot(-1., 1., -1.-1.e-5, 1., 1.), 1., 'p(-5e-6,1)(1) = 1')

end function coord_hycom2_unit_tests

!> Returns true if any cell of u and u_true are not identical. Returns false otherwise.
logical function test_scalar(verbose, u, u_true, label)
  logical,            intent(in) :: verbose !< If true, write results to stdout
  real,               intent(in) :: u      !< Values to test [A]
  real,               intent(in) :: u_true !< Values to test against (correct answer) [A]
  character(len=*),   intent(in) :: label  !< Message

  test_scalar = .false.
  if (abs(u - u_true) > 0.) test_scalar = .true.
  if (test_scalar .or. verbose) then
    if (abs(u - u_true) > 0.) then
      write(stdout,'("coord_hycom2:",2(x,a,1pe24.16),x,2a)') 'u=',u,'err=',u-u_true,label,' <<<<<<< WRONG'
      write(stderr,'("coord_hycom2:",2(x,a,1pe24.16),x,2a)') 'u=',u,'err=',u-u_true,label,' <<<<<<< WRONG'
      write(stderr,'(3a)') 'coord_hycom2: test "',label,' FAILED!'
    else
      write(stdout,'("coord_hycom2:",2(x,a,1pe24.16),x,a)') 'u=',u,'err=',u-u_true,label
    endif
  endif

end function test_scalar

end module coord_hycom2
