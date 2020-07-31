! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> A thin layer to give the MOM_intrinsics_function module access to the
!! intrinsics normally used in Fortran by prepending "f_"
module fortran_intrinsics

public :: f_sin
public :: f_atan
public :: f_tan
contains

!> sin(x)
real elemental function f_sin(x)
  real, intent(in) :: x !< Argument to sin(x)
  f_sin = sin(x)
end function f_sin

!> atan(x)
real elemental function f_atan(x)
  real, intent(in) :: x !< Argument to atan(x)
  f_atan = atan(x)
end function f_atan

!> tan(x)
real elemental function f_tan(x)
  real, intent(in) :: x !< Argument to tan(x)
  f_tan = tan(x)
end function f_tan

end module fortran_intrinsics

!> A module with intrinsic functions written for reproducible answers.
!! Originally this module existed because some intrinsic functions were
!! not supported by some compilers. We now also supply intrinsic functions
!! calculated from series or other means. These functions are not as efficient
!! as those in the math libraries but can be used during non-time critical
!! initialization of data.
module MOM_intrinsic_functions

use iso_fortran_env, only : stdout => output_unit, stderr => error_unit
use iso_fortran_env, only : int64, real64
use fortran_intrinsics, only : f_sin
use fortran_intrinsics, only : f_atan
use fortran_intrinsics, only : f_tan

implicit none ; private

public :: intrinsic_functions_unit_tests
public :: invcosh
public :: cuberoot
public :: sin_m6
public :: cos_m6
public :: sind_m6
public :: cosd_m6
public :: pi
public :: pi_180
public :: rootin

! Floating point model, if bit layout from high to low is (sign, exp, frac)

integer, parameter :: bias = maxexponent(1.) - 1
  !< The double precision exponent offset
integer, parameter :: signbit = storage_size(1.) - 1
  !< Position of sign bit
integer, parameter :: explen = 1 + ceiling(log(real(bias))/log(2.))
  !< Bit size of exponent
integer, parameter :: expbit = signbit - explen
  !< Position of lowest exponent bit
integer, parameter :: fraclen = expbit
  !< Length of fractional part

!> The ratio of a circle's circumference to its diameter, approximately
!! 22/7, or 355/113, ...
!! Some people can recite hundreds of digits (base 10). Here are just
!! 41 digits: 3.141592653589793238462643383279502884197169399...
!! In IEEE 754 single precision floating point is 3.141593
!! (24 mantissa bits or 7 digits).
!! In IEEE 754 double precision floating point representation pi is
!! 3.14159265358979323846 (52 mantissa bits or 21 digits) which is the
!! value found in the C library math.h. We provide more digits (40)
!! here (rounded) for no better reason than the compilers handle it.
real, parameter :: pi = 3.141592653589793238462643383279502884197
real, parameter :: pi_180 = 0.01745329251994329576923690768488612713443

logical :: use_fortran_intrinsics = .false.

contains

!> Evaluate the inverse cosh, either using a math library or an
!! equivalent expression
function invcosh(x)
  real, intent(in) :: x !< The argument of the inverse of cosh [nondim].  NaNs will
                        !! occur if x<1, but there is no error checking
  real :: invcosh  ! The inverse of cosh of x [nondim]

#ifdef __INTEL_COMPILER
  invcosh = acosh(x)
#else
  invcosh = log(x+sqrt(x*x-1))
#endif

end function invcosh

!> Returns the cube root of a real argument at roundoff accuracy, in a form that works properly with
!! rescaling of the argument by integer powers of 8.  If the argument is a NaN, a NaN is returned.
elemental function cuberoot(x) result(root)
  real, intent(in) :: x !< The argument of cuberoot in arbitrary units cubed [A3]
  real :: root !< The real cube root of x in arbitrary units [A]

  real :: asx ! The absolute value of x rescaled by an integer power of 8 to put it into
              ! the range from 0.125 < asx <= 1.0, in ambiguous units cubed [B3]
  real :: root_asx ! The cube root of asx [B]
  real :: ra_3 ! root_asx cubed [B3]
  real :: num ! The numerator of an expression for the evolving estimate of the cube root of asx
              ! in arbitrary units that can grow or shrink with each iteration [B C]
  real :: den ! The denominator of an expression for the evolving estimate of the cube root of asx
              ! in arbitrary units that can grow or shrink with each iteration [C]
  real :: num_prev ! The numerator of an expression for the previous iteration of the evolving estimate
              ! of the cube root of asx in arbitrary units that can grow or shrink with each iteration [B D]
  real :: np_3 ! num_prev cubed  [B3 D3]
  real :: den_prev ! The denominator of an expression for the previous iteration of the evolving estimate of
              ! the cube root of asx in arbitrary units that can grow or shrink with each iteration [D]
  real :: dp_3 ! den_prev cubed  [C3]
  real :: r0  ! Initial value of the iterative solver. [B C]
  real :: r0_3 ! r0 cubed [B3 C3]
  integer :: itt

  integer(kind=int64) :: e_x, s_x

  if ((x >= 0.0) .eqv. (x <= 0.0)) then
    ! Return 0 for an input of 0, or NaN for a NaN input.
    root = x
  else
    call rescale_cbrt(x, asx, e_x, s_x)

    !   Iteratively determine root_asx = asx**1/3 using Halley's method and then Newton's method,
    ! noting that Halley's method onverges monotonically and needs no bounding.  Halley's method is
    ! slightly more complicated that Newton's method, but converges in a third fewer iterations.
    !   Keeping the estimates in a fractional form Root = num / den allows this calculation with
    ! no real divisions during the iterations before doing a single real division at the end,
    ! and it is therefore more computationally efficient.

    ! This first estimate gives the same magnitude of errors for 0.125 and 1.0 after two iterations.
    ! The first iteration is applied explicitly.
    r0 = 0.707106
    r0_3 = r0 * r0 * r0
    num = r0 * (r0_3 + 2.0 * asx)
    den = 2.0 * r0_3 + asx

    do itt=1,2
      ! Halley's method iterates estimates as Root = Root * (Root**3 + 2.*asx) / (2.*Root**3 + asx).
      num_prev = num ; den_prev = den

      ! Pre-compute these as integer powers, to avoid `pow()`-like intrinsics.
      np_3 = num_prev * num_prev * num_prev
      dp_3 = den_prev * den_prev * den_prev

      num = num_prev * (np_3 + 2.0 * asx * dp_3)
      den = den_prev * (2.0 * np_3 + asx * dp_3)
      ! Equivalent to:  root_asx = root_asx * (root_asx**3 + 2.*asx) / (2.*root_asx**3 + asx)
    enddo
    ! At this point the error in root_asx is better than 1 part in 3e14.
    root_asx = num / den

    ! One final iteration with Newton's method polishes up the root and gives a solution
    ! that is within the last bit of the true solution.
    ra_3 = root_asx * root_asx * root_asx
    root_asx = root_asx - (ra_3 - asx) / (3.0 * (root_asx * root_asx))

    root = descale(root_asx, e_x, s_x)
  endif
end function cuberoot


!> Rescale `a` to the range [0.125, 1) and compute its cube-root exponent.
pure subroutine rescale_cbrt(a, x, e_r, s_a)
  real, intent(in) :: a
    !< The real parameter to be rescaled for cube root in arbitrary units cubed [A3]
  real, intent(out) :: x
    !< The rescaled value of a in the range from 0.125 < asx <= 1.0, in ambiguous units cubed [B3]
  integer(kind=int64), intent(out) :: e_r
    !< Cube root of the exponent of the rescaling of `a`
  integer(kind=int64), intent(out) :: s_a
    !< The sign bit of a

  integer(kind=int64) :: xb
    ! Floating point value of a, bit-packed as an integer
  integer(kind=int64) :: e_a
    ! Unscaled exponent of a
  integer(kind=int64) :: e_x
    ! Exponent of x
  integer(kind=int64) :: e_div, e_mod
    ! Quotient and remainder of e in e = 3*(e/3) + modulo(e,3).

  ! Pack bits of a into xb and extract its exponent and sign.
  xb = transfer(a, 1_int64)
  s_a = ibits(xb, signbit, 1)
  e_a = ibits(xb, expbit, explen) - bias

  ! Compute terms of exponent decomposition e = 3*(e/3) + modulo(e,3).
  ! (Fortran division is round-to-zero, so we must emulate floor division.)
  e_mod = modulo(e_a, 3_int64)
  e_div = (e_a - e_mod)/3

  ! Our scaling decomposes e_a into e = {3*(e/3) + 3} + {modulo(e,3) - 3}.

  ! The first term is a perfect cube, whose cube root is computed below.
  e_r = e_div + 1

  ! The second term ensures that x is shifted to [0.125, 1).
  e_x = e_mod - 3

  ! Insert the new 11-bit exponent into xb and write to x and extend the
  ! bitcount to 12, so that the sign bit is zero and x is always positive.
  call mvbits(e_x + bias, 0, explen + 1, xb, fraclen)
  x = transfer(xb, 1.)
end subroutine rescale_cbrt


!> Undo the rescaling of a real number back to its original base.
pure function descale(x, e_a, s_a) result(a)
  real, intent(in) :: x
    !< The rescaled value which is to be restored in ambiguous units [B]
  integer(kind=int64), intent(in) :: e_a
    !< Exponent of the unscaled value
  integer(kind=int64), intent(in) :: s_a
    !< Sign bit of the unscaled value
  real :: a
    !< Restored value with the corrected exponent and sign in arbitrary units [A]

  integer(kind=int64) :: xb
    ! Bit-packed real number into integer form
  integer(kind=int64) :: e_x
    ! Biased exponent of x

  ! Apply the corrected exponent and sign to x.
  xb = transfer(x, 1_int64)
  e_x = ibits(xb, expbit, explen)
  call mvbits(e_a + e_x, 0, explen, xb, expbit)
  call mvbits(s_a, 0, 1, xb, signbit)
  a = transfer(xb, 1.)
end function descale

!> Returns sin(x) where x is in radians
real elemental function sin_m6(x)
  real, intent(in) :: x !< Argument of sin [radians]
  real :: a ! abs(x) or multiples thereof
  real :: s ! 1 or -1 depending on quadrant

  s = 1.0
  a = abs(x)
  if (a>2.0*pi) then
    ! Reduce range to 0..2pi
    a = mod(a, 2.0*pi)
  endif
  if (a>pi) then
    ! Reduce range to 0..pi
    a = a - pi
    s = -1.0
  endif
  if (a<=0.5*pi) then
    sin_m6 = sign(sin_Taylor(a), s*x)
  elseif (a<=pi) then
    sin_m6 = sign(sin_Taylor(pi-a), s*x)
  endif

  if (use_fortran_intrinsics) sin_m6 = sin(x)

end function sin_m6

!> Returns sin(x) where x is in degrees
real elemental function sind_m6(x)
  real, intent(in) :: x !< Argument of sin [degrees]
  real :: a ! abs(x) or multiples thereof
  real :: s ! 1 or -1 depending on quadrant
  real :: d2r ! degrees to radians

  d2r = pi / 180.
  s = 1.0
  a = abs(x)
  if (a>360.) then
    ! Reduce range to 0..360
    a = mod(a, 360.)
  endif
  if (a>180.) then
    ! Reduce range to 0..180
    a = a - 180.
    s = -1.0
  endif
  if (a<=90.) then
    sind_m6 = sign(sin_Taylor(d2r*a), s*x)
  elseif (a<=180.) then
    sind_m6 = sign(sin_Taylor(d2r*(180.-a)), s*x)
  endif

  if (use_fortran_intrinsics) sind_m6 = sin(d2r*x)

end function sind_m6

!> Returns sin(x) if x is in range -pi/2..pi/2 calculated using Taylor
!! series. This approach adds Taylor series terms from smallest to largest
!! and is thus as accurate as the underlying f.p. representation can be.
real elemental function sin_Taylor(x)
  real, intent(in) :: x !< Argument of sin in range -pi/2..pi/2 [radians]
  ! Local variables
  integer, parameter :: n = 16 ! N-1 number of terms in series
  ! Coefficients in Taylor series
  ! https://en.wikipedia.org/wiki/Sine#Series_definition
  ! 15 terms of the Taylor series are needed for 64 bit precision when x=pi
  ! 12 terms of the Taylor series are needed for 64 bit precision when x=pi/2
  !  9 terms of the Taylor series are needed for 64 bit precision when x=pi/4
  real, parameter :: C(16) = (/0.16666666666666666, 8.3333333333333333E-003, &
                           1.9841269841269839E-004, 2.7557319223985884E-006, &
                           2.5052108385441710E-008, 1.6059043836821608E-010, &
                           7.6471637318198144E-013, 2.8114572543455198E-015, &
                           8.2206352466243264E-018, 1.9572941063391257E-020, &
                           3.8681701706306830E-023, 6.4469502843844724E-026, &
                           9.1836898637955449E-029, 1.1309962886447716E-031, &
                           1.2161250415535179E-034, 1.1516335620771951E-037/)
  real :: x2 ! x**2
  real :: xxx(n) ! computed powers of x**2
  real :: r ! accumulated terms
  integer :: j ! term number

  x2 = x*x
  xxx(1) = -x2
  do j = 2, n
    xxx(j) = -xxx(j-1) * x2
  enddo
  r = 0.0
  do j = n, 1, -1
    r = r + C(j) * xxx(j)
  enddo

  sin_Taylor = ( 1.0 + r ) * x

end function sin_Taylor

! !> Returns sin(x) if x is in range -pi/2..pi/2 calculated using Horner's
! !! method applied to the Taylor series polynomial.
! real elemental function sin_Horner(x)
!   real, intent(in) :: x !< Argument of sin in range -pi/2..pi/2 [radians]
!   ! Local variables
!   integer, parameter :: n = 16 ! N-1 number of terms in series
!   ! Coefficients in Taylor series, divided by coefficient of previous term
!   ! 15 terms of the Taylor series are needed for 64 bit precision when x=pi
!   ! 12 terms of the Taylor series are needed for 64 bit precision when x=pi/2
!   !  9 terms of the Taylor series are needed for 64 bit precision when x=pi/4
!   real, parameter :: C(19) = (/0.16666666666666667,0.05,0.0238095238095238081, &
!                    0.013888888888888889,0.00909090909090909,0.00641025641025641, &
!                    0.004761904761904762,0.003676470588235294,0.0029239766081871343, &
!                    0.002380952380952381,0.001976284584980237,0.0016666666666666667, &
!                    0.0014245014245014246,0.0012315270935960591,0.001075268817204301, &
!                    0.000946969696969697,0.0008403361344537816,0.0007507507507507508, &
!                    0.0006747638326585695/) ! https://en.wikipedia.org/wiki/Sine#Series_definition
!   real :: x2 ! x**2
!   real :: r ! accumulated terms
!   integer :: j ! term number

!   x2 = x*x
!   r = 1.0
!   do j = n, 1, -1
!     r = 1.0 - ( C(j) * x2 ) * r
!   enddo

!   sin_Horner = r * x

! end function sin_Horner

!> Returns cos(x) where x is in radians
real elemental function cos_m6(x)
  real, intent(in) :: x !< Argument of cos [radians]
  real :: a ! abs(x) or multiples thereof

  a = abs(x)
  if (a>=2.0*pi) then
    ! Reduce range to 0..2pi
    a = mod(a, 2.0*pi)
  endif
  if (a>pi) then
    ! Reduce range to 0..pi
    a = abs(pi - a)
  endif
  if (a<=0.5*pi) then
    cos_m6 = cos_Taylor(a)
  elseif (a<=pi) then
    cos_m6 = -cos_Taylor(pi-a)
  endif

  if (use_fortran_intrinsics) cos_m6 = cos(x)

end function cos_m6

!> Returns cos(x) where x is in degrees
real elemental function cosd_m6(x)
  real, intent(in) :: x !< Argument of cos [degrees]
  real :: a ! abs(x) or multiples thereof
  real :: d2r ! degrees to radians

  d2r = pi / 180.
  a = abs(x)
  if (a>=360.) then
    ! Reduce range to 0..2pi
    a = mod(a, 360.)
  endif
  if (a>180.) then
    ! Reduce range to 0..pi
    a = abs(180. - a)
  endif
  if (a<=90.) then
    cosd_m6 = cos_Taylor(d2r*a)
  elseif (a<=180.) then
    cosd_m6 = -cos_Taylor(d2r*(180.-a))
  endif

  if (use_fortran_intrinsics) cosd_m6 = cos(d2r*x)

end function cosd_m6

!> Returns cos(x) if x is in range -pi/2..pi/2 calculated using Taylor
!! series. This approach adds Taylor series terms from smallest to largest
!! and is thus as accurate as the underlying f.p. representation can be.
real elemental function cos_Taylor(x)
  real, intent(in) :: x !< Argument of sin in range -pi/2..pi/2 [radians]
  ! Local variables
  integer, parameter :: n = 16 ! N-1 number of terms in series
  ! Coefficients in Taylor series
  ! https://en.wikipedia.org/wiki/Sine#Series_definition
  real, parameter :: C(20) = (/ 0.50000000000000000, 4.1666666666666664E-002, &
                            1.3888888888888887E-003, 2.4801587301587298E-005, &
                            2.7557319223985888E-007, 2.0876756987868096E-009, &
                            1.1470745597729723E-011, 4.7794773323873846E-014, &
                            1.5619206968586225E-016, 4.1103176233121644E-019, &
                            8.8967913924505722E-022, 1.6117375710961182E-024, &
                            2.4795962632247972E-027, 3.2798892370698378E-030, &
                            3.7699876288159054E-033, 3.8003907548547434E-036, &
                            3.3871575355211618E-039, 2.6882202662866363E-042, &
                            1.9119632050402820E-045, 1.2256174391283858E-048/)
  real :: x2 ! x**2
  real :: xxx(n) ! computed powers of x**2
  real :: r ! accumulated terms
  integer :: j ! term number

  x2 = x*x
  xxx(1) = -x2
  do j = 2, n
    xxx(j) = -xxx(j-1) * x2
  enddo
  r = 0.0
  do j = n, 1, -1
    r = r + C(j) * xxx(j)
  enddo

  cos_Taylor = 1.0 + r

end function cos_Taylor

!> Returns x**(1/n), the integer n'th root of x
real elemental function rootin(x, n)
  real, intent(in) :: x !< Argument to be raised to (1/n)
  integer, intent(in) :: n !< Inverse power (assumed >1)
  ! Local variables
  real :: a ! x/n
  real :: nm1on ! (n-1)/n
  real :: recip ! 1/x**(n-1)
  real :: xn, xnm1 ! previous iteration results

  a = abs(x)
  if (a<=0.) then ! zero is a special case
    rootin = 0.
    return
  endif
  xn = -1.
  rootin = 1.0 ! Since n'th roots are closer to 1 than "x" we start at 1
  do while ( abs(rootin - xn)>0. ) ! rootin==xn means we found a unique solution
    xnm1 = xn
    xn = rootin ! previous iteration to allow testing for convergence
    recip = 1.0 / ( float(n) * xn**(n-1) )
    rootin = xn + recip * ( A - xn**n )
    if ( abs(rootin - xnm1)<=0. ) then
      ! If rootin==xnm1 then the solution is oscillating in the last bit
      ! so we choose the larger value and stop iterating
      rootin = max(xn, xnm1)
      exit
    endif
  enddo

! ! Heron's method is mathematically the same as Newton's method but
! ! arrives as a different result (last-bit difference)
! a = 1.0 / float(n)
! nm1on = float(n-1) * a ! (n-1)/n
! a = a * abs(x) ! x/n
!
! xnm1 = -1.
! xn = 0.
! rootin = 1.
! do while ( abs(rootin - xn)>0. .and. abs(rootin - xnm1)>0. )
!   xnm1 = xn
!   xn = rootin ! previous iteration to allow testing for convergence
!   recip = 1.0 / ( rootin**(n-1) ) ! 1/x
!   rootin = nm1on * rootin + recip * a ! x <- (n-1)/n * x + A / x**(n-1)
! enddo

  if (mod(n,2)==1) rootin = sign(rootin, x) ! Allow negative roots for odd powers
  if (use_fortran_intrinsics) rootin = x**(1.0/float(n))

end function rootin

!> Returns true if any unit test of intrinsic_functions fails, or false if they all pass.
function intrinsic_functions_unit_tests(verbose) result(fail)
  logical, intent(in) :: verbose !< If true, write results to stdout
  logical :: fail !< True if any of the unit tests fail

  ! Local variables
  real :: testval  ! A test value for self-consistency testing [nondim]
  real :: x ! Temporary argument
  integer :: n

  fail = .false.  ! Start with no fails
  if (verbose) write(stdout,*) '==== MOM_intrinsic_functions: intrinsics_unit_tests ==='

  ! Cube root tests
  fail = fail .or. Test_cuberoot(verbose, 1.2345678901234e9)
  fail = fail .or. Test_cuberoot(verbose, -9.8765432109876e-21)
  fail = fail .or. Test_cuberoot(verbose, 64.0)
  fail = fail .or. Test_cuberoot(verbose, -0.5000000000001)
  fail = fail .or. Test_cuberoot(verbose, 0.0)
  fail = fail .or. Test_cuberoot(verbose, 1.0)
  fail = fail .or. Test_cuberoot(verbose, 0.125)
  fail = fail .or. Test_cuberoot(verbose, 0.965)
  fail = fail .or. Test_cuberoot(verbose, 1.0 - epsilon(1.0))
  fail = fail .or. Test_cuberoot(verbose, 1.0 - 0.5*epsilon(1.0))

  testval = 1.0e-99
  do n=-160,160
    fail = fail .or. Test_cuberoot(verbose, testval)
    testval = (-2.908 * (1.414213562373 + 1.2345678901234e-5*n)) * testval
  enddo

  ! Trig tests
  if (verbose) write(stdout,'(a25,1pe24.16)') 'module pi:',pi
  if (verbose) write(stdout,'(a25,2a24,x,a)') '','result','correct result','err/epsilon'

  fail = test(fail, pi, 4.0 * atan( 1.0 ), 'module pi (v. library)')

  ! Sine tests
  if (verbose) write(stdout,*) 'Tests of sin()'
  fail = test(fail, sin_m6(0.0), 0., 'sin(0)')
  fail = test(fail, sin_m6(pi/12.), 0.25*(sqrt(6.)-sqrt(2.)), 'sin(pi/12)=0.2588...', inexact=1.)
  fail = test(fail, sin_m6(pi/6.), .5, 'sin(pi/6)=0.5', inexact=1.)
  fail = test(fail, sin_m6(0.25*pi), 0.5*sqrt(2.), 'sin(pi/4)=sqrt(0.5)', inexact=1.)
  fail = test(fail, sin_m6(pi/3.), 0.5*sqrt(3.), 'sin(pi/3)=sqrt(3/4)', inexact=1.)
  fail = test(fail, sin_m6(0.5*pi), 1.0, 'sin(pi/2)=1')
  fail = test(fail, sin_m6(pi), 0., 'sin(pi)')
  fail = test(fail, sin_m6(1.5*pi), -1.0, 'sin(3/2 pi)')
  fail = test(fail, sin_m6(2.5*pi), 1.0, 'sin(5/2 pi)')
  fail = test(fail, sin_m6(-2.5*pi), -1.0, 'sin(-5/2 pi)')

  ! Cosine tests
  if (verbose) write(stdout,*) 'Tests of cos()'
  fail = test(fail, cos_m6(0.), 1., 'cos(0)=1')
  fail = test(fail, cos_m6(0.25*pi), sqrt(0.5), 'cos(pi/4)=sqrt(0.5)',inexact=1.)
  fail = test(fail, cos_m6(0.5*pi), 0., 'cos(pi/2)=0')
  fail = test(fail, cos_m6(pi), -1., 'cos(pi)=-1')
  fail = test(fail, cos_m6(1.5*pi), 0., 'cos(3/2 pi)=0')
  fail = test(fail, cos_m6(2.0*pi), 1., 'cos(2pi)=-1')

  ! Tests that sin(x)**2 + cos(x)**2 = 1 (or less within a bit)
  if (verbose) write(stdout,*) 'Tests of sin(x)**2 + cos(x)**2'
  x = 0.5*pi
  fail = test(fail, cos_m6(x)**2+sin_m6(x)**2, 1., 'cos^2+sin^2, x=pi/2')
  x = pi/3.
  fail = test(fail, cos_m6(x)**2+sin_m6(x)**2, 1., 'cos^2+sin^2, x=pi/3')
  x = 0.25
  fail = test(fail, cos_m6(x)**2+sin_m6(x)**2, 1., 'cos^2+sin^2, x=1/4', inexact=1.)
  x = 0.5
  fail = test(fail, cos_m6(x)**2+sin_m6(x)**2, 1., 'cos^2+sin^2, x=1/2')

  ! Sine tests in degrees
  if (verbose) write(stdout,*) 'Tests of sind() in degrees'
  fail = test(fail, sind_m6(0.0), 0., 'sin(0)')
  fail = test(fail, sind_m6(30.), .5, 'sin(30)=0.5', inexact=1.)
  fail = test(fail, sind_m6(45.), 0.5*sqrt(2.), 'sin(45)=sqrt(0.5)', inexact=1.)
  fail = test(fail, sind_m6(60.), 0.5*sqrt(3.), 'sin(60)=sqrt(3/4)', inexact=1.)
  fail = test(fail, sind_m6(90.), 1.0, 'sin(90)=1')
  fail = test(fail, sind_m6(180.), 0., 'sin(180)')
  fail = test(fail, sind_m6(270.), -1.0, 'sin(270)')

  ! Cosine tests in degrees
  if (verbose) write(stdout,*) 'Tests of cosd() in degrees'
  fail = test(fail, cosd_m6(0.), 1., 'cos(0)=1')
  fail = test(fail, cosd_m6(45.), sqrt(0.5), 'cos(45)=sqrt(0.5)',inexact=1.)
  fail = test(fail, cosd_m6(90.), 0., 'cos(90)=0')
  fail = test(fail, cosd_m6(180.), -1., 'cos(180)=-1')
  fail = test(fail, cosd_m6(270.), 0., 'cos(270)=0')
  fail = test(fail, cosd_m6(360.), 1., 'cos(360)=-1')

  ! Test rootin()
  if (verbose) write(stdout,*) 'Tests of rootin()'
  fail = test(fail, rootin(0.,2), 0., 'rootin(0,2)=0')
  fail = test(fail, rootin(0.,5), 0., 'rootin(0,5)=0')
  fail = test(fail, rootin(1.,2), 1., 'rootin(1,2)=0')
  fail = test(fail, rootin(1.,3), 1., 'rootin(1,3)=0')
  fail = test(fail, rootin(1.,4), 1., 'rootin(1,4)=0')
  fail = test(fail, rootin(1.,5), 1., 'rootin(1,5)=0')
  fail = test(fail, rootin(0.25,2), 0.5, 'rootin(1/4,2)=1/2')
  fail = test(fail, rootin(0.125,3), 0.5, 'rootin(1/8,2)=1/2')
  fail = test(fail, rootin(-0.125,3), -0.5, 'rootin(-1/8,2)=-1/2')
  fail = test(fail, rootin(4.,2), 2., 'rootin(4,2)=2')
  fail = test(fail, rootin(9.,2), 3., 'rootin(9,2)=3')
  fail = test(fail, rootin(16.,2), 4., 'rootin(16,2)=4')
  fail = test(fail, rootin(25.,2), 5., 'rootin(25,2)=5')
  fail = test(fail, rootin(-125.,3), -5., 'rootin(-125,3)=-5')
  fail = test(fail, rootin(2.**30,30), 2., 'rootin(2**30,30)=2')
  fail = test(fail, rootin(-5.**13,13), -5., 'rootin(-5**13,13)=-5', inexact=1.)
  fail = test(fail, rootin(0.5,2), sqrt(0.5), 'rootin(1/2,2)=sqrt(1/2)', inexact=1.)
  fail = test(fail, rootin(2.,2), sqrt(2.0), 'rootin(2,2)=sqrt(2)', inexact=1.)
  x = rootin(2.,2)
  fail = test(fail, x**2, 2.0, 'rootin(2,2)**2=2', inexact=1.)

  if (verbose .and. .not. fail) write(stdout,*) 'Pass'

  contains

  !> Returns true if a==b or previous failed
  logical function test(previous_test, val, correctval, msg, inexact)
    logical, intent(in) :: previous_test !< True if an earlier test failed
    real, intent(in) :: val !< Value to test against correct value
    real, intent(in) :: correctval !< Correct value
    character(len=*), intent(in) :: msg !< Label
    real, optional :: inexact !< If present allow a relative error of inexact*epsilon
    ! Local variables
    integer :: chan ! Stream to print to
    real :: err ! error
    err = abs( val - correctval )
    test = err > 0. ! Comparison must be exact
    if (present(inexact).and.test) then
      test = err / ( abs(correctval) * epsilon(val) ) > inexact ! Allow a round off difference
    endif
    chan = stdout
    if (test.and..not.verbose) chan = stderr
    if (verbose.or.test) then
      if (test) then
        err = err / ( abs(val) * epsilon(correctval) )
        write(chan,'(a24,":",2(1pe24.16),g9.2,a)') msg,val,correctval,err,' <--- FAIL!'
      elseif (err<=0.) then
        write(chan,'(a24,":",2(1pe24.16))') msg,val,correctval
      else
        err = err / ( abs(val) * epsilon(correctval) )
        if (err>1.) then
          write(chan,'(a24,":",2(1pe24.16),g9.2,a)') msg,val,correctval,err,' <- differ in last bits'
        else
          write(chan,'(a24,":",2(1pe24.16),g9.2,a)') msg,val,correctval,err,' <- differs in last bit'
        endif
      endif
    endif
    test = test .or. previous_test
  end function test

  !> True if the cube of cuberoot(val) does not closely match val. False otherwise.
  logical function Test_cuberoot(verbose, val)
    logical, intent(in) :: verbose !< If true, write results to stdout
    real, intent(in) :: val  !< The real value to test, in arbitrary units [A]
    ! Local variables
    real :: diff ! The difference between val and the cube root of its cube [A].

    diff = val - cuberoot(val)**3
    Test_cuberoot = (abs(diff) > 2.0e-15*abs(val))

    if (Test_cuberoot) then
      write(stdout, '("For val = ",ES22.15,", (val - cuberoot(val**3))) = ",ES9.2," <-- FAIL")') val, diff
    elseif (verbose) then
      write(stdout, '("For val = ",ES22.15,", (val - cuberoot(val**3))) = ",ES9.2)') val, diff

    endif
  end function Test_cuberoot

end function intrinsic_functions_unit_tests

end module MOM_intrinsic_functions
