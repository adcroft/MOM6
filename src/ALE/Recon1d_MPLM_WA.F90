!> Monotonized Piecewise Linear Method 1D reconstruction
!!
!! This implementation of PLM follows White and Adcroft, 2008. The PLM slopes are first limited following
!! Colella and Woodward, 1984, but are then further limited to ensure the edge values moving across cell
!! boundaries are monotone.  The first and last cells are always limited to PCM.
module Recon1d_MPLM_WA

! This file is part of MOM6. See LICENSE.md for the license.

use Recon1d_PLM_CW, only : PLM_CW, testing

implicit none ; private

public MPLM_WA, testing

!> Limited Monotonic PLM reconstruction following White and Adcroft, 2008
!!
!! The following methods are defined in the PLM_CW parent class:
!!   %init()
!!   %lr_edge()
!!   %average()
!!   %inv_f()
!!   %destroy()
!!   %init_parent()
!!   %destroy_parent()
!! The following methods are defined in the Recon1d base class:
!!   %cell_mean()
!!   %remap_to_sub_grid()
!! All other methods are defined in this module.
type, extends (PLM_CW) :: MPLM_WA

contains
  !> Implementation of the PLM reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of unit tests for the PLM reconstruction
  procedure :: unit_tests => unit_tests

  !> Duplicate interface to reconstruct()
  procedure :: reconstruct_parent => reconstruct

#ifdef REMAP_INLINE
! This block is here to test whether the compiler can do better if we have local copies of
! the remapping functions.
  !> Remaps the column to subgrid h_sub
  procedure :: remap_to_sub_grid => remap_to_sub_grid
#endif

end type MPLM_WA

contains

!> Calculate a 1D PLM reconstructions based on h(:) and u(:)
subroutine reconstruct(this, h, u)
  class(MPLM_WA), intent(inout) :: this !< This reconstruction
  real,           intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,           intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables
  real :: slp(this%n) ! The PLM slopes (difference across cell) [A]
  real :: mslp(this%n) ! The monotonized PLM slopes [A]
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
  this%ul(1) = u(1)
  this%ur(1) = u(1)
  do k = 2, n-1
    this%ul(k) = u(k) - 0.5 * mslp(k) ! Left edge value of cell k
    this%ur(k) = u(k) + 0.5 * mslp(k) ! Right edge value of cell k
  enddo
  this%ul(n) = u(n)
  this%ur(n) = u(n)

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
  real :: neighbor_edge ! Edge value of nieghbor cell [A]
  real :: this_edge ! Edge value of this cell [A]
  real :: slp ! Magnitude of PLM central slope [A]

  ! Comparison are made assuming +ve slopes
  slp = abs(s_c)

  ! Check that left edge is between right edge of cell to the left and this cell mean
  neighbor_edge = u_l + 0.5 * s_l
  this_edge = u_c - 0.5 * s_c
  if ( ( this_edge - neighbor_edge ) * ( u_c - this_edge ) < 0. ) then
    ! Using the midpoint works because the neighbor is similarly adjusted
    this_edge = 0.5 * ( this_edge + neighbor_edge )
    slp = min( slp, abs( this_edge - u_c ) * 2. )
  endif

  ! Check that right edge is between left edge of cell to the right and this cell mean
  neighbor_edge = u_r - 0.5 * s_r
  this_edge = u_c + 0.5 * s_c
  if ( ( this_edge - u_c ) * ( neighbor_edge - this_edge ) < 0. ) then
    ! Using the midpoint works because the neighbor is similarly adjusted
    this_edge = 0.5 * ( this_edge + neighbor_edge )
    slp = min( slp, abs( this_edge - u_c ) * 2. )
  endif

  PLM_monotonized_slope = sign( slp, s_c )

end function PLM_monotonized_slope

#ifdef REMAP_INLINE
! This block is here to test whether the compiler can do better if we have local copies of
! the remapping functions.

!> Remaps the column to subgrid h_sub
!!
!! It is assumed that h_sub is a perfect sub-grid of h0, meaning each h0 cell
!! can be constructed by joining a contiguous set of h_sub cells. The integer
!! indices isrc_start, isrc_end, isub_src provide this mapping, and are
!! calculated in MOM_remapping
subroutine remap_to_sub_grid(this, h0, u0, n1, h_sub, &
                                   isrc_start, isrc_end, isrc_max, isub_src, &
                                   u_sub, uh_sub, u02_err)
  class(MPLM_WA), intent(in) :: this !< 1-D reconstruction type
  real,    intent(in)  :: h0(*)  !< Source grid widths (size n0) [H]
  real,    intent(in)  :: u0(*)  !< Source grid widths (size n0) [H]
  integer, intent(in)  :: n1      !< Number of cells in target grid
  real,    intent(in)  :: h_sub(*) !< Overlapping sub-cell thicknesses, h_sub [H]
  integer, intent(in)  :: isrc_start(*) !< Index of first sub-cell within each source cell
  integer, intent(in)  :: isrc_end(*) !< Index of last sub-cell within each source cell
  integer, intent(in)  :: isrc_max(*) !< Index of thickest sub-cell within each source cell
  integer, intent(in)  :: isub_src(*) !< Index of source cell for each sub-cell
  real,    intent(out) :: u_sub(*) !< Sub-cell cell averages (size n1) [A]
  real,    intent(out) :: uh_sub(*) !< Sub-cell cell integrals (size n1) [A H]
  real,    intent(out) :: u02_err !< Integrated reconstruction error estimates [A H]
  ! Local variables
  integer :: i_sub ! Index of sub-cell
  integer :: i0 ! Index into h0(1:n0), source column
  integer :: i_max ! Used to record which sub-cell is the largest contribution of a source cell
  real :: dh_max ! Used to record which sub-cell is the largest contribution of a source cell [H]
  real :: xa, xb ! Non-dimensional position within a source cell (0..1) [nondim]
  real :: dh ! The width of the sub-cell [H]
  real :: duh ! The total amount of accumulated stuff (u*h) [A H]
  real :: dh0_eff ! Running sum of source cell thickness [H]
  integer :: i0_last_thick_cell, n0

  n0 = this%n

  i0_last_thick_cell = 0
  do i0 = 1, n0
    if (h0(i0)>0.) i0_last_thick_cell = i0
  enddo

  ! Loop over each sub-cell to calculate average/integral values within each sub-cell.
  ! Uses: h_sub, isub_src, h0_eff
  ! Sets: u_sub, uh_sub
  xa = 0.
  dh0_eff = 0.
  u02_err = 0.
  do i_sub = 1, n0+n1

    ! Sub-cell thickness from loop above
    dh = h_sub(i_sub)

    ! Source cell
    i0 = isub_src(i_sub)

    ! Evaluate average and integral for sub-cell i_sub.
    ! Integral is over distance dh but expressed in terms of non-dimensional
    ! positions with source cell from xa to xb  (0 <= xa <= xb <= 1).
    dh0_eff = dh0_eff + dh ! Cumulative thickness within the source cell
    if (h0(i0)>0.) then
      xb = dh0_eff / h0(i0) ! This expression yields xa <= xb <= 1.0
      xb = min(1., xb) ! This is only needed when the total target column is wider than the source column
      u_sub(i_sub) = this%average( i0, xa, xb )
    else ! Vanished cell
      xb = 1.
      u_sub(i_sub) = u0(i0)
    endif
    uh_sub(i_sub) = dh * u_sub(i_sub)

    if (isub_src(i_sub+1) /= i0) then
      ! If the next sub-cell is in a different source cell, reset the position counters
      dh0_eff = 0.
      xa = 0.
    else
      xa = xb ! Next integral will start at end of last
    endif

  enddo
  i_sub = n0+n1+1
  ! Sub-cell thickness from loop above
  dh = h_sub(i_sub)
  ! Source cell
  i0 = isub_src(i_sub)

  ! Evaluate average and integral for sub-cell i_sub.
  ! Integral is over distance dh but expressed in terms of non-dimensional
  ! positions with source cell from xa to xb  (0 <= xa <= xb <= 1).
  dh0_eff = dh0_eff + dh ! Cumulative thickness within the source cell
  if (h0(i0)>0.) then
    xb = dh0_eff / h0(i0) ! This expression yields xa <= xb <= 1.0
    xb = min(1., xb) ! This is only needed when the total target column is wider than the source column
    u_sub(i_sub) = this%average( i0, xa, xb )
  else ! Vanished cell
    xb = 1.
    u_sub(i_sub) = u0(i0)
  endif
  uh_sub(i_sub) = dh * u_sub(i_sub)

  ! Loop over each source cell substituting the integral/average for the thickest sub-cell (within
  ! the source cell) with the residual of the source cell integral minus the other sub-cell integrals
  ! aka a genius algorithm for accurate conservation when remapping from Robert Hallberg (@Hallberg-NOAA).
  ! Uses: i0_last_thick_cell, isrc_max, h_sub, isrc_start, isrc_end, uh_sub, u0, h0
  ! Updates: uh_sub
  do i0 = 1, i0_last_thick_cell
    i_max = isrc_max(i0)
    dh_max = h_sub(i_max)
    if (dh_max > 0.) then
      ! duh will be the sum of sub-cell integrals within the source cell except for the thickest sub-cell.
      duh = 0.
      do i_sub = isrc_start(i0), isrc_end(i0)
        if (i_sub /= i_max) duh = duh + uh_sub(i_sub)
      enddo
      uh_sub(i_max) = u0(i0)*h0(i0) - duh
      u02_err = u02_err + max( abs(uh_sub(i_max)), abs(u0(i0)*h0(i0)), abs(duh) )
    endif
  enddo

end subroutine remap_to_sub_grid
#endif

!> Runs PLM reconstruction unit tests and returns True for any fails, False otherwise
logical function unit_tests(this, verbose, stdout, stderr)
  class(MPLM_WA), intent(inout) :: this    !< This reconstruction
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

  unit_tests = test%summarize('MPLM_WA:unit_tests')

end function unit_tests

!> \namespace recon1d_mplm_wa
!!

end module Recon1d_MPLM_WA
