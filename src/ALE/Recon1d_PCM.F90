!> 1D reconstructions using the Piecewise Constant Method (PCM)
module Recon1d_PCM

! This file is part of MOM6. See LICENSE.md for the license.

use Recon1d_type, only : Recon1d, testing

implicit none ; private

public PCM

!> PCM (piecewise constant) reconstruction
!!
!! The following methods are defined in the Recon1d base class:
!!   %cell_mean()
!!   %remap_to_sub_grid()
!! All other methods are defined in this module.
type, extends (Recon1d) :: PCM

contains
  !> Implementation of the PCM initialization
  procedure :: init => init
  !> Implementation of the PCM reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of function returning the PCM edge values
  procedure :: lr_edge => lr_edge
  !> Implementation of the PCM average over an interval [A]
  procedure :: average => average
  !> Implementation of finding the PCM position of a value
  procedure :: inv_f => inv_f
  !> Implementation of deallocation for PCM
  procedure :: destroy => destroy
  !> Implementation of check reconstruction for the PLM reconstruction
  procedure :: check_reconstruction => check_reconstruction
  !> Implementation of unit tests for the PCM reconstruction
  procedure :: unit_tests => unit_tests

  !> Duplicate interface to init()
  procedure :: init_parent => init
  !> Duplicate interface to reconstruct()
  procedure :: reconstruct_parent => reconstruct

end type PCM

contains

!> Initialize a 1D PCM reconstruction for n cells
subroutine init(this, n, h_neglect, check)
  class(PCM),        intent(out) :: this      !< This reconstruction
  integer,           intent(in)  :: n         !< Number of cells in this column
  real, optional,    intent(in)  :: h_neglect !< A negligibly small width used in cell reconstructions [H].
                                              !! Not used by PCM.
  logical, optional, intent(in)  :: check     !< If true, enable some consistency checking

  if (present(h_neglect)) this%n = n ! no-op to avoid compiler warning about unused dummy argument
  if (present(check)) this%check = check

  this%n = n

  allocate( this%u_mean(n) )

end subroutine init

!> Calculate a 1D PCM reconstructions based on h(:) and u(:)
subroutine reconstruct(this, h, u)
  class(PCM), intent(inout) :: this !< This reconstruction
  real,       intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,       intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables
  integer :: k

  this%u_mean(1) = h(1) ! no-op to avoid compiler warning about unused dummy argument

  do k = 1, this%n
    this%u_mean(k) = u(k)
  enddo

end subroutine reconstruct

!> Returns the left/right edge values in cell k of a 1D PCM reconstruction
subroutine lr_edge(this, k, u_left, u_right)
  class(PCM), intent(in)  :: this    !< This reconstruction
  integer,    intent(in)  :: k       !< Cell number
  real,       intent(out) :: u_left  !< Left edge value [A]
  real,       intent(out) :: u_right !< Right edge value [A]

  u_left = this%u_mean(k)
  u_right = this%u_mean(k)

end subroutine lr_edge

!> Position at which 1D PCM reconstruction has a value f in cell k [0..1]
real function inv_f(this, k, f)
  class(PCM), intent(in) :: this !< This reconstruction
  integer,    intent(in) :: k    !< Cell number
  real,       intent(in) :: f    !< Value of reconstruction to solve for [A]

  inv_f = 0.5 ! For PCM, every value between x=0 and x=1 has the same value
                    ! but x=0.5 is the average position
  if ( f < this%u_mean(k) ) then
    inv_f = 0.
  elseif ( f > this%u_mean(k) ) then
    inv_f = 1.
  endif

end function inv_f

!> Average between xa and xb for cell k of a 1D PCM reconstruction [A]
real function average(this, k, xa, xb)
  class(PCM), intent(in) :: this !< This reconstruction
  integer,    intent(in) :: k    !< Cell number
  real,       intent(in) :: xa   !< Start of averaging interval on element (0 to 1)
  real,       intent(in) :: xb   !< End of averaging interval on element (0 to 1)

  average = xb + xa ! no-op to avoid compiler warnings about unused dummy argument
  average = this%u_mean(k)

end function average

!> Deallocate the PCM reconstruction
subroutine destroy(this)
  class(PCM), intent(inout) :: this !< This reconstruction

  deallocate( this%u_mean )

end subroutine destroy

!> Checks the PCM reconstruction for consistency
logical function check_reconstruction(this, h, u)
  class(PCM), intent(in) :: this !< This reconstruction
  real,       intent(in) :: h(*) !< Grid spacing (thickness) [typically H]
  real,       intent(in) :: u(*) !< Cell mean values [A]
  ! Local variables
  integer :: k

  check_reconstruction = .false.

  do k = 1, this%n
    if ( abs( this%u_mean(k) - u(k) ) > 0. ) check_reconstruction = .true.
  enddo

end function check_reconstruction

!> Runs PCM reconstruction unit tests and returns True for any fails, False otherwise
logical function unit_tests(this, verbose, stdout, stderr)
  class(PCM), intent(inout) :: this    !< This reconstruction
  logical,    intent(in)    :: verbose !< True, if verbose
  integer,    intent(in)    :: stdout  !< I/O channel for stdout
  integer,    intent(in)    :: stderr  !< I/O channel for stderr
  ! Local variables
  real, allocatable :: ul(:), ur(:), um(:) ! test values [A]
  type(testing) :: test ! convenience functions
  integer :: k

  call test%set( stdout=stdout ) ! Sets the stdout channel in test
  call test%set( stderr=stderr ) ! Sets the stderr channel in test
  call test%set( verbose=verbose ) ! Sets the verbosity flag in test

  call this%init(3)
  call test%test( this%n /= 3, 'Setting number of levels')
  allocate( um(3), ul(3), ur(3) )

  call this%reconstruct( (/2.,2.,2./), (/1.,3.,5./) )
  call test%real_arr(3, this%u_mean, (/1.,3.,5./), 'Setting cell values')
  do k = 1, 3
    um(k) = this%cell_mean(k)
  enddo
  call test%real_arr(3, um, (/1.,3.,5./), 'Return cell mean')

  do k = 1, 3
    call this%lr_edge(k, ul(k), ur(k))
  enddo
  call test%real_arr(3, ul, (/1.,3.,5./), 'Return left edge')
  call test%real_arr(3, ur, (/1.,3.,5./), 'Return right edge')

  do k = 1, 3
    um(k) = this%average(k, 0.5, 0.75)
  enddo
  call test%real_arr(3, um, (/1.,3.,5./), 'Return interval average')

  do k = 1, 3
    ul(k) = this%inv_f(k, real(2*k-2))
    um(k) = this%inv_f(k, real(2*k-1))
    ur(k) = this%inv_f(k, real(2*k-0))
  enddo
  call test%real_arr(3, ul, (/0.,0.,0./), 'Return position of f<')
  call test%real_arr(3, um, (/0.5,0.5,0.5/), 'Return position of f=')
  call test%real_arr(3, ur, (/1.,1.,1./), 'Return position of f>')

  unit_tests = test%summarize('PCM:unit_tests')

end function unit_tests

!> \namespace recon1d_pcm
!!

end module Recon1d_PCM
