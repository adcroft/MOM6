!> Piecewise Linear Method 1D reconstruction
!!
!! This implementation of PLM follows Colella and Woodward, 1984, with cells resorting to PCM for
!! extrema including first and last cells in column. The cell-wise reconstructions are limited so
!! that the edge values (which are also the extrema in a cell) are bounded by the neighbors. The
!! limiter yields monotonicity for the CFL<1 transport problem where parts of a cell can only move
!! to a neighboring cell, but does not yield monotonic profiles for the general remapping problem.
!! The first and last cells are always limited to PCM.
module Recon1d_PLM_CW

! This file is part of MOM6. See LICENSE.md for the license.

use Recon1d_type, only : Recon1d, testing

implicit none ; private

public PLM_CW, testing

!> PLM reconstruction following Colella and Woodward, 1984
!!
!! The source for the methods ultimately used by this class are:
!!   init()                    *locally defined
!!   reconstruct()             *locally defined
!!   lr_edge()                 *locally defined
!!   average()                 *locally defined
!!   inv_f()                   *locally defined
!!   check_reconstruction()    *locally defined
!!   unit_tests()              *locally defined
!!   destroy()                 *locally defined
!!   cell_mean()            -> Recon1d%cell_mean()
!!   remap_to_sub_grid()    -> Recon1d%remap_to_sub_grid()
!!   init_parent()          -> init()
!!   reconstruct_parent()   -> reconstruct()
type, extends (Recon1d) :: PLM_CW

  real, allocatable :: ul(:) !< Left edge value [A]
  real, allocatable :: ur(:) !< Right edge value [A]

contains
  !> Implementation of the PLM_CW initialization
  procedure :: init => init
  !> Implementation of the PLM_CW reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of function returning the PLM_CW edge values
  procedure :: lr_edge => lr_edge
  !> Implementation of the PLM_CW average over an interval [A]
  procedure :: average => average
  !> Implementation of finding the PLM_CW position of a value
  procedure :: inv_f => inv_f
  !> Implementation of deallocation for PLM_CW
  procedure :: destroy => destroy
  !> Implementation of check reconstruction for the PLM_CW reconstruction
  procedure :: check_reconstruction => check_reconstruction
  !> Implementation of unit tests for the PLM_CW reconstruction
  procedure :: unit_tests => unit_tests

  !> Duplicate interface to init()
  procedure :: init_parent => init
  !> Duplicate interface to reconstruct()
  procedure :: reconstruct_parent => reconstruct

end type PLM_CW

contains

!> Initialize a 1D PLM reconstruction for n cells
subroutine init(this, n, h_neglect, check)
  class(PLM_CW),     intent(out) :: this      !< This reconstruction
  integer,           intent(in)  :: n         !< Number of cells in this column
  real, optional,    intent(in)  :: h_neglect !< A negligibly small width used in cell reconstructions [H]
  logical, optional, intent(in)  :: check     !< If true, enable some consistency checking

  this%n = n

  allocate( this%u_mean(n) )
  allocate( this%ul(n) )
  allocate( this%ur(n) )

  this%h_neglect = tiny( this%u_mean(1) )
  if (present(h_neglect)) this%h_neglect = h_neglect
  this%check = .false.
  if (present(check)) this%check = check

end subroutine init

!> Calculate a 1D PLM reconstructions based on h(:) and u(:)
subroutine reconstruct(this, h, u)
  class(PLM_CW), intent(inout) :: this !< This reconstruction
  real,          intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,          intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables
  real :: slp ! The PLM slopes (difference across cell) [A]
  real :: sigma_l, sigma_c, sigma_r ! Left, central and right slope estimates as
                                    ! differences across the cell [A]
  real :: u_min, u_max ! Minimum and maximum value across cell [A]
  real :: u_l, u_r, u_c ! Left, right, and center values [A]
  real :: h_l, h_c, h_r ! Thickness of left, center and right cells [H]
  real :: h_c0 ! Thickness of center with h_neglect added [H]
  integer :: k, n

  n = this%n

  ! Loop over all cells
  do k = 1, n
    this%u_mean(k) = u(k)
  enddo

  ! Boundary cells use PCM
  this%ul(1) = u(1)
  this%ur(1) = u(1)

  ! Loop over interior cells
  do k = 2, n-1
    u_l = u(k-1)
    u_c = u(k)
    u_r = u(k+1)

    ! Side differences
    sigma_r = u_r - u_c
    sigma_l = u_c - u_l

    h_l = h(k-1)
    h_c = h(k)
    h_r = h(k+1)
    ! Avoids division by zero
    h_c0 = h_c + this%h_neglect

    ! This is the second order slope given by equation 1.7 of
    ! Piecewise Parabolic Method, Colella and Woodward (1984),
    ! http://dx.doi.org/10.1016/0021-991(84)90143-8.
    ! For uniform resolution it simplifies to ( u_r - u_l )/2 .
    sigma_c = ( h_c / ( h_c0 + ( h_l + h_r ) ) ) * ( &
                  ( 2.*h_l + h_c ) / ( h_r + h_c0 ) * sigma_r &
                + ( 2.*h_r + h_c ) / ( h_l + h_c0 ) * sigma_l )

    ! Limit slope so that reconstructions are bounded by neighbors
    u_min = min( u_l, u_c, u_r )
    u_max = max( u_l, u_c, u_r )

    if ( (sigma_l * sigma_r) > 0.0 ) then
      ! This limits the slope so that the edge values are bounded by the two cell averages spanning the edge
      slp = sign( min( abs(sigma_c), 2.*min( u_c - u_min, u_max - u_c ) ), sigma_c )
    else
      ! Extrema in the mean values require a PCM reconstruction
      slp = 0.0
    endif

    ! Left edge
    u_min = min( u_c, u_l )
    u_max = max( u_c, u_l )
    u_l = u_c - 0.5 * slp
    this%ul(k) = max( min( u_l, u_max), u_min )

    ! Right edge
    u_min = min( u_c, u_r )
    u_max = max( u_c, u_r )
    u_r = u_c + 0.5 * slp
    this%ur(k) = max( min( u_r, u_max), u_min )
  enddo

  ! Boundary cells use PCM
  this%ul(n) = u(n)
  this%ur(n) = u(n)

end subroutine reconstruct

!> Returns the left/right edge values in cell k of a 1D PLM reconstruction
subroutine lr_edge(this, k, u_left, u_right)
  class(PLM_CW), intent(in)  :: this    !< This reconstruction
  integer,       intent(in)  :: k       !< Cell number
  real,          intent(out) :: u_left  !< Left edge value [A]
  real,          intent(out) :: u_right !< Right edge value [A]

  u_left = this%ul(k)
  u_right = this%ur(k)

end subroutine lr_edge

!> Position at which 1D PLM reconstruction has a value f in cell k [0..1]
real function inv_f(this, k, f)
  class(PLM_CW), intent(in) :: this !< This reconstruction
  integer,       intent(in) :: k    !< Cell number
  real,          intent(in) :: f    !< Value of reconstruction to solve for [A]

  inv_f = 0.5 ! Avoid compiler warnings

  if ( f < this%ul(k) ) then
    inv_f = 0.
  elseif ( f > this%ur(k) ) then
    inv_f = 1.
  elseif ( this%ur(k) - this%ul(k) == 0. ) then
    inv_f = 0.5 ! Same as for PCM
  else ! ul < f < ur
    ! Note that if ur=ul (i.e. ur-ul=0 then we would meet one of the previous
    ! conditions that avoid division by zero
    inv_f = ( f - this%ul(k) ) / ( this%ur(k) - this%ul(k) )
  endif

end function inv_f

!> Average between xa and xb for cell k of a 1D PLM reconstruction [A]
real function average(this, k, xa, xb)
  class(PLM_CW), intent(in) :: this !< This reconstruction
  integer,       intent(in) :: k    !< Cell number
  real,          intent(in) :: xa   !< Start of averaging interval on element (0 to 1)
  real,          intent(in) :: xb   !< End of averaging interval on element (0 to 1)
  real :: xmab ! Mid-point between xa and xb (0 to 1)
  !real :: u_a, u_b ! Values at xa and xb [A]

  !u_a = this%ul(k) * ( 1. - xa ) + this%ur(k) * xa
  !u_b = this%ul(k) * ( 1. - xb ) + this%ur(k) * xb
  !average = 0.5 * ( u_a + u_b )

  ! Mid-point between xa and xb
  xmab = 0.5 * ( xa + xb )

  ! The following expression is exact at xmab=0 and xmab=1,
  ! i.e. gives the numerically correct values.
  ! It is not obvious that the expression is monotonic but according to
  ! https://math.stackexchange.com/questions/907329/accurate-floating-point-linear-interpolation
  ! it will be for the default rounding behavior. Otherwise is it
  ! then possible this expression can be outside the range of ul and ur?
  average = this%ul(k) * ( 1. - xmab ) + this%ur(k) * xmab

! ! The following is more complicated but seems to ensure being within bounds.
! ! This expression for u_a can overshoot u_r but is good for xmab<<1
! u_a = this%ul(k) + ( this%ur(k)  - this%ul(k) ) * xmab
! ! This expression for u_b can overshoot u_l but is good for 1-xmab<<1
! u_b = this%ur(k) + ( this%ul(k)  - this%ur(k) ) * ( 1. - xmab )
! ! Replace xmab with -1 for xmab<0.5, 1 for xmab>=0.5
! xmab = sign(1., xmab-0.5)
! ! Select either u_a or u_b, depending whether mid-point of xa, xb is smaller/larger than 0.5
! average = xmab * u_b + ( 1. - xmab ) * u_a

end function average

!> Deallocate the PLM reconstruction
subroutine destroy(this)
  class(PLM_CW), intent(inout) :: this !< This reconstruction

  deallocate( this%u_mean, this%ul, this%ur )

end subroutine destroy

!> Checks the PLM_CW reconstruction for consistency
logical function check_reconstruction(this, h, u)
  class(PLM_CW), intent(in) :: this !< This reconstruction
  real,          intent(in) :: h(*) !< Grid spacing (thickness) [typically H]
  real,          intent(in) :: u(*) !< Cell mean values [A]
  ! Local variables
  integer :: k

  check_reconstruction = .false.

  do k = 1, this%n
    if ( abs( this%u_mean(k) - u(k) ) > 0. ) check_reconstruction = .true.
  enddo

  ! Check implied curvature
  do k = 1, this%n
    if ( ( this%u_mean(k) - this%ul(k) ) * ( this%ur(k) - this%u_mean(k) ) < 0. ) check_reconstruction = .true.
  enddo

  ! Check bounding of right edges
  do K = 1, this%n-1
    if ( ( this%ur(k) - this%u_mean(k) ) * ( this%u_mean(k+1) - this%ur(k) ) < 0. ) check_reconstruction = .true.
  enddo

  ! Check bounding of left edges
  do K = 2, this%n
    if ( ( this%u_mean(k) - this%ul(k) ) * ( this%ul(k) - this%u_mean(k-1) ) < 0. ) check_reconstruction = .true.
  enddo

end function check_reconstruction

!> Runs PLM reconstruction unit tests and returns True for any fails, False otherwise
logical function unit_tests(this, verbose, stdout, stderr)
  class(PLM_CW), intent(inout) :: this    !< This reconstruction
  logical,       intent(in)    :: verbose !< True, if verbose
  integer,       intent(in)    :: stdout  !< I/O channel for stdout
  integer,       intent(in)    :: stderr  !< I/O channel for stderr
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

  call this%destroy()
  deallocate( um, ul, ur, ull, urr )

  allocate( um(4), ul(4), ur(4) )
  call this%init(4)

  ! These values lead to non-monotonic reconstuctions which are
  ! valid for transport problems but not always appropriate for
  ! remapping to arbitrary resolution grids.
  ! The O(h^2) slopes are -, 2, 2, - and the limited
  ! slopes are 0, 1, 1, 0 so the everywhere the reconstructions
  ! are bounded by neighbors but ur(2) and ul(3) are out-of-order.
  call this%reconstruct( (/1.,1.,1.,1./), (/0.,3.,4.,7./) )
  do k = 1, 4
    call this%lr_edge(k, ul(k), ur(k))
  enddo
  call test%real_arr(4, ul, (/0.,2.,3.,7./), 'Return left edge')
  call test%real_arr(4, ur, (/0.,4.,5.,7./), 'Return right edge')

  deallocate( um, ul, ur )

  unit_tests = test%summarize('PLM_CW:unit_tests')

end function unit_tests

!> \namespace recon1d_plm_cw
!!

end module Recon1d_PLM_CW
