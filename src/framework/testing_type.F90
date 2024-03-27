!> Provides convenience functions used in unit tests
!!
!! This module is designed to be very elementary and does not use other parts of MOM6
!! to facilitate rapid compilation that helps during development and testing.
!! The module provides the class "testing" that helps accumulate testing history
!! (such as pass/fail state) in a rudimentary way.
module testing_type

implicit none ; private

! This file is part of MOM6. See LICENSE.md for the license.
! Original module written by Laurent White, 2008.06.09

!> Testing class
type, public :: testing
  private
  !> True if any fail has been encountered since instantiation of "testing"
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
    procedure :: real_arr => real_arr   !< Compare array of reals
    procedure :: int_arr => int_arr     !< Compare array of integers
end type

contains

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

end module testing_type
