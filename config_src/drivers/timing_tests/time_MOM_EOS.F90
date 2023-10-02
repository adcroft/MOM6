program time_MOM_EOS

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_EOS, only : EOS_type
use MOM_EOS, only : EOS_LINEAR, EOS_LINEAR_STRING
use MOM_EOS, only : EOS_UNESCO, EOS_UNESCO_STRING
use MOM_EOS, only : EOS_WRIGHT, EOS_WRIGHT_STRING
use MOM_EOS, only : EOS_WRIGHT_FULL, EOS_WRIGHT_FULL_STRING
use MOM_EOS, only : EOS_WRIGHT_REDUCED, EOS_WRIGHT_RED_STRING
use MOM_EOS, only : EOS_TEOS10, EOS_TEOS10_STRING
use MOM_EOS, only : EOS_ROQUET_RHO, EOS_ROQUET_RHO_STRING
use MOM_EOS, only : EOS_ROQUET_SPV, EOS_ROQUET_SPV_STRING
use MOM_EOS, only : EOS_JACKETT06, EOS_JACKETT06_STRING
use MOM_EOS, only : EOS_manual_init
use MOM_EOS, only : calculate_density
use MOM_EOS, only : calculate_density, calculate_spec_vol

implicit none

integer, parameter :: n_eos = 9
integer :: EOS_list(n_eos)
character(len=20) :: EOS_str(n_eos)
integer, parameter :: n_fns = 4
character(len=40) :: fn_labels(n_fns)

! Testing parameters:
!  nic is number of computations per call
!  halo is a data on either end of the array that should not be used
!  nits is how many times to repeat the call between turning the timer on/off
!  nsamp repeats the timing to collect statistics on the measurement
integer, parameter :: nic=23, halo=4, nits=5000, nsamp=40

real, dimension(n_eos,n_fns) :: timings, tmean, tstd, tmin, tmax
integer :: i, j

EOS_list(1) = EOS_LINEAR
EOS_str(1) = EOS_LINEAR_STRING
EOS_list(2) = EOS_UNESCO
EOS_str(2) = EOS_UNESCO_STRING
EOS_list(3) = EOS_WRIGHT
EOS_str(3) = EOS_WRIGHT_STRING
EOS_list(4) = EOS_WRIGHT_FULL
EOS_str(4) = EOS_WRIGHT_FULL_STRING
EOS_list(5) = EOS_WRIGHT_REDUCED
EOS_str(5) = EOS_WRIGHT_RED_STRING
EOS_list(6) = EOS_TEOS10
EOS_str(6) = EOS_TEOS10_STRING
EOS_list(7) = EOS_ROQUET_RHO
EOS_str(7) = EOS_ROQUET_RHO_STRING
EOS_list(8) = EOS_ROQUET_SPV
EOS_str(8) = EOS_ROQUET_SPV_STRING
EOS_list(9) = EOS_JACKETT06
EOS_str(9) = EOS_JACKETT06_STRING
fn_labels(1) = 'calculate_density_scalar()'
fn_labels(2) = 'calculate_density_array()'
fn_labels(3) = 'calculate_spec_vol_scalar()'
fn_labels(4) = 'calculate_spec_vol_array()'

tmean(:,:) = 0.
tstd(:,:) = 0.
tmin(:,:) = 1.e9
tmax(:,:) = 0.
do i = 1, nsamp
  call run_suite(EOS_list, nic, halo, nits, timings)
  tmean(:,:) = tmean(:,:) + timings(:,:)
  tstd(:,:) = tstd(:,:) + timings(:,:)**2
  tmin(:,:) = min( tmin(:,:), timings(:,:) )
  tmax(:,:) = max( tmax(:,:), timings(:,:) )
enddo
tmean(:,:) = tmean(:,:) / real(nsamp)
tstd(:,:) = tstd(:,:) / real(nsamp) ! Mean squares
tstd(:,:) = tstd(:,:) - tmean(:,:)**2 ! Variance
tstd(:,:) = sqrt( tstd(:,:) * real(nsamp) / real(nsamp-1) ) ! Standard deviation

! Display results in YAML
write(*,'(a,":")') 'timings'
do i = 1, n_eos
  do j = 1, n_fns
    write(*,'(2x,"- name: ",a)') "'MOM_EOS_" // trim(EOS_str(i)) // " " // trim(fn_labels(j))
    write(*,'(4x,"min: ",1pe11.4)') tmin(i,j)
    write(*,'(4x,"mean:",1pe11.4)') tmean(i,j)
    write(*,'(4x,"std: ",1pe11.4)') tstd(i,j)
    write(*,'(4x,"max: ",1pe11.4)') tmax(i,j)
  enddo
enddo

contains

subroutine run_suite(EOS_list, nic, halo, nits, timings)
integer, intent(in)  :: EOS_list(n_eos) !< IDs of EOS forms to loop over
integer, intent(in)  :: nic          !< Width of computational domain
integer, intent(in)  :: halo         !< Width of halo to add on either end
integer, intent(in)  :: nits         !< Number of calls to sample
                                     !! (large enough that the CPU timers can resolve
                                     !! the loop)
real,    intent(out) :: timings(n_eos,n_fns) !< The average time taken for nits calls
                                     !! First index corresponds to EOS
                                     !! Second index: 1 = scalar args,
                                     !! 2 = array args without halo,
                                     !! 3 = array args with halo and "dom".
type(EOS_type) :: EOS
integer :: e, i, dom(2)
real :: start, finish, T, S, P, rho
real, dimension(:), allocatable :: T1, S1, P1, rho1

T = 10.
S = 35.
P = 2000.e4

! Time the scalar interface
do e = 1, n_eos
  call EOS_manual_init(EOS, form_of_EOS=EOS_list(e), &
                       Rho_T0_S0=1030., dRho_dT=0.2, dRho_dS=-0.7)
  call cpu_time(start)
  do i = 1, nits*nic ! Calling nic* to make similar cost to array call
    call calculate_density(T, S, P, rho, EOS)
  enddo
  call cpu_time(finish)
  timings(e,1) = (finish-start)/real(nits)
  call cpu_time(start)
  do i = 1, nits*nic ! Calling nic* to make similar cost to array call
    call calculate_spec_vol(T, S, P, rho, EOS)
  enddo
  call cpu_time(finish)
  timings(e,2) = (finish-start)/real(nits)
enddo

! Time the "dom" interface, 1D array + halos
allocate( T1(nic+2*halo), source=T)
allocate( S1(nic+2*halo), source=S)
allocate( P1(nic+2*halo), source=P)
allocate( rho1(nic+2*halo) )
dom(:) = [1+halo,nic+halo]
do e = 1, n_eos
  call EOS_manual_init(EOS, form_of_EOS=EOS_list(e), &
                       Rho_T0_S0=1030., dRho_dT=0.2, dRho_dS=-0.7)
  call cpu_time(start)
  do i = 1, nits
    call calculate_density(T1, S1, P1, rho1, EOS, dom)
  enddo
  call cpu_time(finish)
  timings(e,3) = (finish-start)/real(nits)
  call cpu_time(start)
  do i = 1, nits
    call calculate_spec_vol(T1, S1, P1, rho1, EOS, dom)
  enddo
  call cpu_time(finish)
  timings(e,4) = (finish-start)/real(nits)
enddo
deallocate( T1, S1, P1, rho1 )

end subroutine run_suite

end program time_MOM_EOS
