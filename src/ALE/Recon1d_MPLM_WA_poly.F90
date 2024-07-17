!> Monotonized Piecewise Linear Method 1D reconstruction using polynomial representation
!!
!! This implementation of PLM follows White and Adcroft, 2008. The PLM slopes are first limited following
!! Colella and Woodward, 1984, but are then further limited to ensure the edge values moving across cell
!! boundaries are monotone.  The first and last cells are always limited to PCM.
!!
!! This stores and evaluates the reconstruction using a polynomial representation which is not preferred
!! but was the form used in OM4.
module Recon1d_MPLM_WA_poly

! This file is part of MOM6. See LICENSE.md for the license.

use Recon1d_MPLM_WA, only : MPLM_WA, testing

implicit none ; private

public MPLM_WA_poly, testing

!> Limited Monotonic PLM reconstruction following White and Adcroft, 2008
!!
!! The following methods are defined in the PLM_CW parent class (inherited via MPLM_WA):
!!   %lr_edge()
!!   %inv_f()
!!   %destroy()
!! The following methods are defined in the Recon1d base class:
!!   %cell_mean()
!!   %remap_to_sub_grid()
!! All other methods are defined in this module.
type, extends (MPLM_WA) :: MPLM_WA_poly

  ! Legacy representation
  integer :: degree !< Degree of polynomial used in legacy representation
  real, allocatable, dimension(:,:) :: poly_coef !< Polynomial coefficients in legacy representation

contains
  !> Implementation of the PLM initialization
  procedure :: init => init
  !> Implementation of the PLM reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of the PLM average over an interval [A]
  procedure :: average => average
  !> Implementation of unit tests for the PLM reconstruction
  procedure :: unit_tests => unit_tests

  !> Duplicate interface to reconstruct()
  procedure :: init_parent => init
  !> Duplicate interface to reconstruct()
  procedure :: reconstruct_parent => reconstruct

end type MPLM_WA_poly

contains

!> Initialize a 1D PLM reconstruction for n cells
subroutine init(this, n, h_neglect)
  class(MPLM_WA_poly), intent(out) :: this !< This reconstruction
  integer,        intent(in)  :: n    !< Number of cells in this column
  real, optional, intent(in)  :: h_neglect !< A negligibly small width used in cell reconstructions [H]

  this%n = n

  allocate( this%u_mean(n) )
  allocate( this%ul(n) )
  allocate( this%ur(n) )

  this%h_neglect = tiny( this%u_mean(1) )
  if (present(h_neglect)) this%h_neglect = h_neglect

  this%degree = 2
  allocate( this%poly_coef(n,2) )

end subroutine init

!> Calculate a 1D PLM reconstructions based on h(:) and u(:)
subroutine reconstruct(this, h, u)
  class(MPLM_WA_poly), intent(inout) :: this !< This reconstruction
  real,           intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,           intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables
  real :: slp(this%n) ! The PLM slopes (difference across cell) [A]
  real :: mslp(this%n) ! The monotonized PLM slopes [A]
  real :: e_r, edge ! Edge values [A]
  real :: almost_one  ! A value that is slightly smaller than 1 [nondim]
  integer :: k, n

  n = this%n

  ! Loop over all cells
  do k = 1, n
    this%u_mean(k) = u(k)
  enddo

  ! Loop on interior cells
  do k = 2, n-1
    slp(k) = PLM_slope_wa(h(k-1), h(k), h(k+1), this%h_neglect, u(k-1), u(k), u(k+1))
  enddo ! end loop on interior cells

  ! Boundary cells use PCM. Extrapolation is handled after monotonization.
  slp(1) = 0.
  slp(n) = 0.

  ! This loop adjusts the slope so that edge values are monotonic.
  do k = 2, n-1
    mslp(k) = PLM_monotonized_slope( u(k-1), u(k), u(k+1), slp(k-1), slp(k), slp(k+1) )
  enddo ! end loop on interior cells
  mslp(1) = 0.
  mslp(n) = 0.

  ! Store and return edge values and polynomial coefficients.
  almost_one = 1. - epsilon(e_r)
  this%ul(1) = u(1)
  this%ur(1) = u(1)
  this%poly_coef(1,1) = u(1)
  this%poly_coef(1,2) = 0.
  do k = 2, n-1
    this%ul(k) = u(k) - 0.5 * mslp(k) ! Left edge value of cell k
    this%ur(k) = u(k) + 0.5 * mslp(k) ! Right edge value of cell k

    this%poly_coef(k,1) = this%ul(k)
    this%poly_coef(k,2) = this%ur(k) - this%ul(k)
    ! Check to see if this evaluation of the polynomial at x=1 would be
    ! monotonic w.r.t. the next cell's edge value. If not, scale back!
    edge = this%poly_coef(k,2) + this%poly_coef(k,1)
    e_r = u(k+1) - 0.5 * sign( mslp(k+1), slp(k+1) )
    if ( (edge-u(k))*(e_r-edge)<0.) then
      this%poly_coef(k,2) = this%poly_coef(k,2) * almost_one
    endif
  enddo
  this%ul(n) = u(n)
  this%ur(n) = u(n)
  this%poly_coef(n,1) = u(n)
  this%poly_coef(n,2) = 0.

end subroutine reconstruct

!> Returns a limited PLM slope following White and Adcroft, 2008, in the same arbitrary
!! units [A] as the input values.
!! Note that this is not the same as the Colella and Woodward method.
real elemental pure function PLM_slope_wa(h_l, h_c, h_r, h_neglect, u_l, u_c, u_r)
  real, intent(in) :: h_l !< Thickness of left cell in arbitrary grid thickness units [H]
  real, intent(in) :: h_c !< Thickness of center cell in arbitrary grid thickness units [H]
  real, intent(in) :: h_r !< Thickness of right cell in arbitrary grid thickness units [H]
  real, intent(in) :: h_neglect !< A negligible thickness [H]
  real, intent(in) :: u_l !< Value of left cell in arbitrary units [A]
  real, intent(in) :: u_c !< Value of center cell in arbitrary units [A]
  real, intent(in) :: u_r !< Value of right cell in arbitrary units [A]
  ! Local variables
  real :: sigma_l, sigma_c, sigma_r ! Left, central and right slope estimates as
                                    ! differences across the cell [A]
  real :: u_min, u_max ! Minimum and maximum value across cell [A]

  ! Side differences
  sigma_r = u_r - u_c
  sigma_l = u_c - u_l

  ! Quasi-second order difference
  sigma_c = 2.0 * ( u_r - u_l ) * ( h_c / ( h_l + 2.0*h_c + h_r + h_neglect) )

  ! Limit slope so that reconstructions are bounded by neighbors
  u_min = min( u_l, u_c, u_r )
  u_max = max( u_l, u_c, u_r )
  if ( (sigma_l * sigma_r) > 0.0 ) then
    ! This limits the slope so that the edge values are bounded by the
    ! two cell averages spanning the edge.
    PLM_slope_wa = sign( min( abs(sigma_c), 2.*min( u_c - u_min, u_max - u_c ) ), sigma_c )
  else
    ! Extrema in the mean values require a PCM reconstruction avoid generating
    ! larger extreme values.
    PLM_slope_wa = 0.0
  endif

  ! This block tests to see if roundoff causes edge values to be out of bounds
  if (u_c - 0.5*abs(PLM_slope_wa) < u_min .or.  u_c + 0.5*abs(PLM_slope_wa) > u_max) then
    PLM_slope_wa = PLM_slope_wa * ( 1. - epsilon(PLM_slope_wa) )
  endif

  ! An attempt to avoid inconsistency when the values become unrepresentable.
  ! ### The following 1.E-140 is dimensionally inconsistent. A newer version of
  ! PLM is progress that will avoid the need for such rounding.
  if (abs(PLM_slope_wa) < 1.E-140) PLM_slope_wa = 0.

end function PLM_slope_wa

!> Returns a limited PLM slope following Colella and Woodward 1984, in the same
!! arbitrary units as the input values [A].
real elemental pure function PLM_monotonized_slope(u_l, u_c, u_r, s_l, s_c, s_r)
  real, intent(in) :: u_l !< Value of left cell in arbitrary units [A]
  real, intent(in) :: u_c !< Value of center cell in arbitrary units [A]
  real, intent(in) :: u_r !< Value of right cell in arbitrary units [A]
  real, intent(in) :: s_l !< PLM slope of left cell [A]
  real, intent(in) :: s_c !< PLM slope of center cell [A]
  real, intent(in) :: s_r !< PLM slope of right cell [A]
  ! Local variables
  real :: e_r, e_l, edge ! Right, left and temporary edge values [A]
  real :: almost_two ! The number 2, almost [nondim]
  real :: slp ! Magnitude of PLM central slope [A]

  almost_two = 2. * ( 1. - epsilon(s_c) )

  ! Edge values of neighbors abutting this cell
  e_r = u_l + 0.5*s_l
  e_l = u_r - 0.5*s_r
  slp = abs(s_c)

  ! Check that left edge is between right edge of cell to the left and this cell mean
  edge = u_c - 0.5 * s_c
  if ( ( edge - e_r ) * ( u_c - edge ) < 0. ) then
    edge = 0.5 * ( edge + e_r )
    slp = min( slp, abs( edge - u_c ) * almost_two )
  endif

  ! Check that right edge is between left edge of cell to the right and this cell mean
  edge = u_c + 0.5 * s_c
  if ( ( edge - u_c ) * ( e_l - edge ) < 0. ) then
    edge = 0.5 * ( edge + e_l )
    slp = min( slp, abs( edge - u_c ) * almost_two )
  endif

  PLM_monotonized_slope = sign( slp, s_c )

end function PLM_monotonized_slope

!> Average between xa and xb for cell k of a 1D PLM reconstruction [A]
!! Note: this uses the simple polynomial form a + b * x  on x E (0,1)
!! which can overshoot at x=1
real function average(this, k, xa, xb)
  class(MPLM_WA_poly), intent(in) :: this !< This reconstruction
  integer,        intent(in) :: k    !< Cell number
  real,           intent(in) :: xa   !< Start of averaging interval on element (0 to 1)
  real,           intent(in) :: xb   !< End of averaging interval on element (0 to 1)

  average = this%poly_coef(k,1) &
          + this%poly_coef(k,2) * 0.5 * ( xb + xa )

end function average

!> Runs PLM reconstruction unit tests and returns True for any fails, False otherwise
logical function unit_tests(this, verbose, stdout, stderr)
  class(MPLM_WA_poly), intent(inout) :: this    !< This reconstruction
  logical,        intent(in)    :: verbose !< True, if verbose
  integer,        intent(in)    :: stdout  !< I/O channel for stdout
  integer,        intent(in)    :: stderr  !< I/O channel for stderr
  ! Local variables
  real, allocatable :: ul(:), ur(:), um(:) ! test values [A]
  real, allocatable :: ull(:), urr(:) ! test values [A]
  type(testing) :: test ! convenience functions
  integer :: k

  call test%set( stdout=stdout ) ! Sets the stdout channel in test
  call test%set( stderr=stderr ) ! Sets the stderr channel in test
  call test%set( verbose=verbose ) ! Sets the verbosity flag in test

  call this%init(3)
  call test%test( this%n /= 3, 'Setting number of levels')
  allocate( um(3), ul(3), ur(3), ull(3), urr(3) )

  call this%reconstruct( (/2.,2.,2./), (/1.,3.,5./) )
  call test%real_arr(3, this%u_mean, (/1.,3.,5./), 'Setting cell values')
  do k = 1, 3
    um(k) = this%cell_mean(k)
  enddo
  call test%real_arr(3, um, (/1.,3.,5./), 'Return cell mean')

  do k = 1, 3
    call this%lr_edge(k, ul(k), ur(k))
  enddo
  call test%real_arr(3, ul, (/1.,2.,5./), 'Return left edge')
  call test%real_arr(3, ur, (/1.,4.,5./), 'Return right edge')

  do k = 1, 3
    um(k) = this%average(k, 0.5, 0.75) ! Average from x=0.25 to 0.75 in each cell
  enddo
  call test%real_arr(3, um, (/1.,3.25,5./), 'Return interval average')

  do k = 1, 3
    ull(k) = this%inv_f(k, real(2*k-2))
    ul(k) = this%inv_f(k, real(2*k-1)-0.5)
    um(k) = this%inv_f(k, real(2*k-1))
    ur(k) = this%inv_f(k, real(2*k-1)+0.5)
    urr(k) = this%inv_f(k, real(2*k-0))
  enddo
  call test%real_arr(3, ull, (/0.,0.,0./), 'Return position of f<<')
  call test%real_arr(3, ul, (/0.,0.25,0./), 'Return position of f<')
  call test%real_arr(3, um, (/0.5,0.5,0.5/), 'Return position of f=')
  call test%real_arr(3, ur, (/1.,0.75,1./), 'Return position of f>')
  call test%real_arr(3, urr, (/1.,1.,1./), 'Return position of f>>')

  unit_tests = test%summarize('MPLM_WA_poly:unit_tests')

end function unit_tests

!> \namespace recon1d_mplm_wa_poly
!!

end module Recon1d_MPLM_WA_poly
