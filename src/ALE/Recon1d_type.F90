!> A generic type for vertical 1D reconstructions
module Recon1d_type

! This file is part of MOM6. See LICENSE.md for the license.

implicit none ; private

public Recon1d
public testing

!> The base class for implementations of 1D reconstructions
type, abstract :: Recon1d

  integer :: n = 0 !< Number of cells in column
  real, allocatable, dimension(:) :: u_mean !< Cell mean [A]
  real :: h_neglect = 0. !< A negligibly small width used in cell reconstructions [same as h, H]

contains

  ! The following functions/subroutines are deferred and must be provided specifically by each scheme

  !> Deferred implementation of initialization
  procedure(i_init), deferred :: init
  !> Deferred implementation of reconstruction function
  procedure(i_reconstruct), deferred :: reconstruct
  !> Deferred implementation of function returning the edge values
  procedure(i_lr_edge), deferred :: lr_edge
  !> Deferred implementation of the average over an interval
  procedure(i_average), deferred :: average
  !> Deferred implementation of the finding position of a value
  procedure(i_inv_f), deferred :: inv_f
  !> Deferred implementation of unit tests for the reconstruction
  procedure(i_unit_tests), deferred :: unit_tests
  !> Deferred implementation of deallocation
  procedure(i_destroy), deferred :: destroy

  ! The following functions/subroutines are shared across all reconstructions and provided by this module

  !> Cell mean of cell k [A]
  procedure :: cell_mean => a_cell_mean

  ! The following functions usually point to the same implementation as above but
  ! for derived secondary children these allow invocation of the parent class function.

  !> Second interface to init(), used to reach the primary class if derived from a primary implementation
  procedure(i_init_parent), deferred :: init_parent
  !> Second interface to reconstruct(), used to reach the primary class if derived from a primary implementation
  procedure(i_reconstruct_parent), deferred :: reconstruct_parent
  !> Second interface to destroy(), used to reach the primary class if derived from a primary implementation
  procedure(i_destroy_parent), deferred :: destroy_parent

end type Recon1d

!> Class to assist in unit tests
type :: testing
  private
  !> True if any fail has been encountered since this instance of "testing" was created
  logical :: state = .false.
  !> Count of tests checked
  integer :: num_tests_checked = 0
  !> Count of tests failed
  integer :: num_tests_failed = 0
  !> If true, be verbose and write results to stdout. Default True.
  logical :: verbose = .true.
  !> Error channel
  integer :: stderr = 0
  !> Standard output channel
  integer :: stdout = 6
  !> If true, stop instantly
  logical :: stop_instantly = .false.
  !> Record instances that fail
  integer :: ifailed(100) = 0.
  !> Record label of first instance that failed
  character(len=:), allocatable :: label_first_fail

  contains
    procedure :: test => test           !< Update the testing state
    procedure :: set => set             !< Set attributes
    procedure :: outcome => outcome     !< Return current outcome
    procedure :: summarize => summarize !< Summarize testing state
!   procedure :: real_scalar => real_scalar !< Compare two reals
    procedure :: real_arr => real_arr   !< Compare array of reals
    procedure :: int_arr => int_arr     !< Compare array of integers
end type

interface

  !> Initialize a 1D reconstruction for n cells
  subroutine i_init(this, n, h_neglect)
    import :: Recon1d
    class(Recon1d), intent(out) :: this !< This reconstruction
    integer,        intent(in)  :: n    !< Number of cells in this column
    real, optional, intent(in)  :: h_neglect !< A negligibly small width used in cell reconstructions [H]
  end subroutine i_init

  !> Calculate a 1D reconstructions based on h(:) and u(:)
  subroutine i_reconstruct(this, h, u)
    import :: Recon1d
    class(Recon1d), intent(inout) :: this !< This reconstruction
    real,           intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
    real,           intent(in)    :: u(*) !< Cell mean values [A]
  end subroutine i_reconstruct

  !> Returns the left/right edge values in cell k of a 1D reconstruction
  subroutine i_lr_edge(this, k, u_left, u_right)
    import :: Recon1d
    class(Recon1d), intent(in)  :: this    !< This reconstruction
    integer,        intent(in)  :: k       !< Cell number
    real,           intent(out) :: u_left  !< Left edge value [A]
    real,           intent(out) :: u_right !< Right edge value [A]
  end subroutine i_lr_edge

  !> Average between xa and xb for cell k of a 1D reconstruction [A]
  !!
  !! It is assumed that 0<=xa<=1, 0<=xb<=1, and xa<=xb
  real function i_average(this, k, xa, xb)
    import :: Recon1d
    class(Recon1d), intent(in) :: this !< This reconstruction
    integer,        intent(in) :: k    !< Cell number
    real,           intent(in) :: xa   !< Start of averaging interval on element (0 to 1)
    real,           intent(in) :: xb   !< End of averaging interval on element (0 to 1)
  end function i_average

  !> Position at which 1D reconstruction has a value f [0..1]
  !!
  !! It is assumed the reconstruction is not decreasing with cell index.
  real function i_inv_f(this, k, f)
    import :: Recon1d
    class(Recon1d), intent(in) :: this !< This reconstruction
    integer,        intent(in) :: k    !< Cell number
    real,           intent(in) :: f    !< Value of reconstruction to solve for [A]
  end function i_inv_f

  !> Runs reconstruction unit tests and returns True for any fails, False otherwise
  !!
  !! Assumes single process/thread context
  logical function i_unit_tests(this, verbose, stdout, stderr)
    import :: Recon1d
    class(Recon1d), intent(inout) :: this    !< This reconstruction
    logical,        intent(in)    :: verbose !< True, if verbose
    integer,        intent(in)    :: stdout  !< I/O channel for stdout
    integer,        intent(in)    :: stderr  !< I/O channel for stderr
  end function i_unit_tests

  !> Deallocate a 1D reconstruction
  subroutine i_destroy(this)
    import :: Recon1d
    class(Recon1d), intent(inout) :: this !< This reconstruction
  end subroutine i_destroy

  !> Second interface to init(), or to parent init()
  subroutine i_init_parent(this, n, h_neglect)
    import :: Recon1d
    class(Recon1d), intent(out) :: this !< This reconstruction
    integer,        intent(in)  :: n    !< Number of cells in this column
    real, optional, intent(in)  :: h_neglect !< A negligibly small width used in cell reconstructions [H]
  end subroutine i_init_parent

  !> Second interface to reconstruct(), or to parent reconstruct()
  subroutine i_reconstruct_parent(this, h, u)
    import :: Recon1d
    class(Recon1d), intent(inout) :: this !< This reconstruction
    real,           intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
    real,           intent(in)    :: u(*) !< Cell mean values [A]
  end subroutine i_reconstruct_parent

  !> Second interface to destroy(), or to parent destroy()
  subroutine i_destroy_parent(this)
    import :: Recon1d
    class(Recon1d), intent(inout) :: this !< This reconstruction
  end subroutine i_destroy_parent

end interface

contains

!> Mean value for cell k [A]
real function a_cell_mean(this, k)
  class(Recon1d), intent(in) :: this !< This reconstruction
  integer,        intent(in) :: k    !< Cell number

  a_cell_mean = this%u_mean(k)

end function a_cell_mean

! =========================================================================================
! The following provide the function for the testing_type helper class
! =========================================================================================

!> Update the state with "test"
subroutine test(this, state, label)
  class(testing),   intent(inout) :: this  !< This testing class
  logical,          intent(in)    :: state !< True to indicate a fail, false otherwise
  character(len=*), intent(in)    :: label !< Message

  this%num_tests_checked = this%num_tests_checked + 1
  if (state) then
    this%state = .true.
    this%num_tests_failed = this%num_tests_failed + 1
    this%ifailed( this%num_tests_failed ) = this%num_tests_checked
    if (this%num_tests_failed == 1) this%label_first_fail = label
  endif
  if (this%stop_instantly .and. this%state) stop 1
end subroutine test

!> Set attributes
subroutine set(this, verbose, stdout, stderr, stop_instantly)
  class(testing), intent(inout) :: this  !< This testing class
  logical, optional, intent(in) :: verbose !< True or false setting to assign to verbosity
  integer, optional, intent(in) :: stdout !< The stdout channel to use
  integer, optional, intent(in) :: stderr !< The stderr channel to use
  logical, optional, intent(in) :: stop_instantly !< If true, stop immediately on error detection

  if (present(verbose)) then
    this%verbose = verbose
  endif
  if (present(stdout)) then
    this%stdout = stdout
  endif
  if (present(stderr)) then
    this%stderr = stderr
  endif
  if (present(stop_instantly)) then
    this%stop_instantly = stop_instantly
  endif
end subroutine set

!> Returns state
logical function outcome(this)
  class(testing), intent(inout) :: this !< This testing class
  outcome = this%state
end function outcome

!> Summarize results
logical function summarize(this, label)
  class(testing),  intent(inout) :: this  !< This testing class
  character(len=*),   intent(in) :: label !< Message
  integer :: i

  if (this%state) then
    write(this%stdout,'(a," : ",a,", ",i4," failed of ",i4," tested")') &
         red('FAIL'), trim(label), this%num_tests_failed, this%num_tests_checked
    write(this%stdout,'(a,100i4)') 'Failed tests:',(this%ifailed(i),i=1,this%num_tests_failed)
    write(this%stdout,'(a,a)') 'First failed test: ',trim(this%label_first_fail)
    write(this%stderr,'(a,100i4)') 'Failed tests:',(this%ifailed(i),i=1,this%num_tests_failed)
    write(this%stderr,'(a,a)') 'First failed test: ',trim(this%label_first_fail)
    write(this%stderr,'(a," : ",a)') trim(label),'FAILED'
  else
    write(this%stdout,'(a," : ",a,", all ",i4," tests passed")') &
         green('Pass'), trim(label), this%num_tests_checked
  endif
  summarize = this%state
end function summarize

!> Pads string with the color codes for red/reset
function red(string) result(newstr)
  character(len=*),  intent(in) :: string !< Input string
  character(len=:), allocatable :: newstr !< The modified output string
  newstr = achar(27)//'[31m'//trim(string)//achar(27)//'[0m'
end function red

!> Pads string with the color codes for green/reset
function green(string) result(newstr)
  character(len=*),  intent(in) :: string !< Input string
  character(len=:), allocatable :: newstr !< The modified output string
  newstr = achar(27)//'[32m'//trim(string)//achar(27)//'[0m'
end function green

!! !> Compare u_test to u_true, report, and return true if a difference larger than tol is measured
!! !!
!! !! If in verbose mode, display results to stdout
!! !! If a difference is measured, display results to stdout and stderr
!! subroutine real_scalar(this, u_test, u_true, label)
!!   class(testing),   intent(inout) :: this   !< This testing class
!!   real,             intent(in)    :: u_test !< Values to test [A]
!!   real,             intent(in)    :: u_true !< Values to test against (correct answer) [A]
!!   character(len=*), intent(in)    :: label  !< Message
!!   ! Local variables
!!   logical :: this_test
!!   real :: err
!!
!!   this_test = .false.
!!
!!   ! Scan for any mismatch between u_test and u_true
!!   err = u_test - u_true
!!   if (abs(err) > 0.) this_test = .true.
!!
!!   ! If either being verbose, or an error was measured then display results
!!   if (this_test .or. this%verbose) then
!!     if (this_test) then
!!       write(this%stdout,'(3(a,1pe24.16,1x),a)') 'Calculated value =',u_test,'correct value =',u_true,'error =',err,label
!!       write(this%stderr,'(3(a,1pe24.16,1x),a)') 'Calculated value =',u_test,'correct value =',u_true,'error =',err,label
!!     else
!!       write(this%stdout,'(2(a,1pe24.16,1x),a)') 'Calculated value =',u_test,'correct value =',u_true,label
!!     endif
!!   endif
!!
!!   call this%test( this_test, label ) ! Updates state and counters in this
!! end subroutine real_scalar

!> Compare u_test to u_true, report, and return true if a difference larger than tol is measured
!!
!! If in verbose mode, display results to stdout
!! If a difference is measured, display results to stdout and stderr
subroutine real_arr(this, n, u_test, u_true, label, tol)
  class(testing),  intent(inout) :: this   !< This testing class
  integer,            intent(in) :: n      !< Number of cells in u
  real, dimension(n), intent(in) :: u_test !< Values to test [A]
  real, dimension(n), intent(in) :: u_true !< Values to test against (correct answer) [A]
  character(len=*),   intent(in) :: label  !< Message
  real,     optional, intent(in) :: tol    !< The tolerance for differences between u and u_true [A]
  ! Local variables
  integer :: k
  logical :: this_test
  real :: tolerance, err

  tolerance = 0.0
  if (present(tol)) tolerance = tol
  this_test = .false.

  ! Scan for any mismatch between u_test and u_true
  do k = 1, n
    if (abs(u_test(k) - u_true(k)) > tolerance) this_test = .true.
  enddo

  ! If either being verbose, or an error was measured then display results
  if (this_test .or. this%verbose) then
    write(this%stdout,'(a4,2a24,1x,a)') 'k','Calculated value','Correct value',label
    if (this_test) write(this%stderr,'(a4,2a24,1x,a)') 'k','Calculated value','Correct value',label
    do k = 1, n
      err = u_test(k) - u_true(k)
      if (abs(err) > tolerance) then
        write(this%stdout,'(i4,1p2e24.16,a,1pe24.16,a)') k, u_test(k), u_true(k), &
                         ' err=', err, red(' <--- WRONG')
        write(this%stderr,'(i4,1p2e24.16,a,1pe24.16,a)') k, u_test(k), u_true(k), &
                         ' err=', err, ' <--- WRONG'
      else
        write(this%stdout,'(i4,1p2e24.16)') k, u_test(k), u_true(k)
      endif
    enddo
  endif

  call this%test( this_test, label ) ! Updates state and counters in this
end subroutine real_arr

!> Compare i_test to i_true and report and return true if a difference is found
!!
!! If in verbose mode, display results to stdout
!! If a difference is measured, display results to stdout and stderr
subroutine int_arr(this, n, i_test, i_true, label)
  class(testing),     intent(inout) :: this   !< This testing class
  integer,               intent(in) :: n      !< Number of cells in u
  integer, dimension(n), intent(in) :: i_test !< Values to test [A]
  integer, dimension(n), intent(in) :: i_true !< Values to test against (correct answer) [A]
  character(len=*),      intent(in) :: label  !< Message
  ! Local variables
  integer :: k
  logical :: this_test

  this_test = .false.

  ! Scan for any mismatch between u_test and u_true
  do k = 1, n
    if (i_test(k) .ne. i_true(k)) this_test = .true.
  enddo

  if (this%verbose) then
     write(this%stdout,'(a12," : calculated =",30i3)') label, i_test
     write(this%stdout,'(12x,"      correct =",30i3)') i_true
     if (this_test) write(this%stdout,'(3x,a,8x,"error =",30i3)') red('FAIL --->'), i_test(:) - i_true(:)
   endif
   if (this_test) then
     write(this%stderr,'(a12," : calculated =",30i3)') label, i_test
     write(this%stderr,'(12x,"      correct =",30i3)') i_true
     write(this%stderr,'("   FAIL --->        error =",30i3)') i_test(:) - i_true(:)
   endif

  call this%test( this_test, label ) ! Updates state and counters in this
end subroutine int_arr

!> \namespace recon1d_type
!!
!! \section section_recon1d_type Generic vertical reconstruction type
!!
!! A class to describe generic reconstruction in 1-D. This module has no implementations
!! but defines the interfaces for members that implement a reconstruction.
!!
!! e.g. a chain of derived reconstructions might look like
!!   Recon1d_type <- Recond1d_XYZ <- Recon1d_XYZ_v2
!! where
!!   Recon1d_type      - defines the interfaces (this module)
!!   Recon1d_XYZ       - extends Recon1d_type, implements the XYZ reconstruction in reconstruct(),
!!                       and reconstruc_parent() -> reconstruct() of the same Recon1d_XYZ module
!!   Recon1d_XYZ_v2    - implements a slight variant of Recon1d_XYZ via reconstruct()
!!                       but reconstruc_parent() is not redefined so that it still is defined by Recon1d_XYZ
!!
!! The module also defines a type "testing" which is used to abbreviate the unit_tests() member.
!! This should no be used outside of these derived reconstructions.
end module Recon1d_type
