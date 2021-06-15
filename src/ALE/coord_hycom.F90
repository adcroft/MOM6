!> Regrid columns for the HyCOM coordinate
module coord_hycom

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL
use MOM_EOS,           only : EOS_type, calculate_density
use regrid_interp,     only : interp_CS_type, build_and_interpolate_grid

implicit none ; private

!> Control structure containing required parameters for the HyCOM coordinate
type, public :: hycom_CS ; private

  !> Number of layers/levels in generated grid
  integer :: nk

  !> Nominal near-surface resolution [Z ~> m]
  real, allocatable, dimension(:) :: coordinateResolution

  !> Nominal density of interfaces [R ~> kg m-3]
  real, allocatable, dimension(:) :: target_density

  !> Use z*-space only in the mixed layer
  logical :: use_z_only_in_ML = .false.

  !> Minimum z*-space region [Z ~> m]
  real :: zspace_min_extent = 0.0

  !> Fractional extra LMB [nondim]
  real :: MLB_extra_frac = 0.0

  !> Maximum extra LMB [Z ~> m]
  real :: MLB_extra_max = 0.0

  !> Resolution to use below mixed-layer base [Z ~> m]
  real :: interior_min_thickness = 0.0

  !> Maximum depths of interfaces [H ~> m or kg m-2]
  real, allocatable, dimension(:) :: max_interface_depths

  !> Maximum thicknesses of layers [H ~> m or kg m-2]
  real, allocatable, dimension(:) :: max_layer_thickness

  !> Interpolation control structure
  type(interp_CS_type) :: interp_CS
end type hycom_CS

public init_coord_hycom, set_hycom_params, build_hycom1_column, end_coord_hycom

contains

!> Initialise a hycom_CS with pointers to parameters
subroutine init_coord_hycom(CS, nk, coordinateResolution, target_density, interp_CS)
  type(hycom_CS),       pointer    :: CS !< Unassociated pointer to hold the control structure
  integer,              intent(in) :: nk !< Number of layers in generated grid
  real, dimension(nk),  intent(in) :: coordinateResolution !< Nominal near-surface resolution [Z ~> m]
  real, dimension(nk+1),intent(in) :: target_density !< Interface target densities [R ~> kg m-3]
  type(interp_CS_type), intent(in) :: interp_CS !< Controls for interpolation

  if (associated(CS)) call MOM_error(FATAL, "init_coord_hycom: CS already associated!")
  allocate(CS)
  allocate(CS%coordinateResolution(nk))
  allocate(CS%target_density(nk+1))

  CS%nk                      = nk
  CS%coordinateResolution(:) = coordinateResolution(:)
  CS%target_density(:)       = target_density(:)
  CS%interp_CS               = interp_CS

end subroutine init_coord_hycom

!> This subroutine deallocates memory in the control structure for the coord_hycom module
subroutine end_coord_hycom(CS)
  type(hycom_CS), pointer :: CS !< Coordinate control structure

  ! nothing to do
  if (.not. associated(CS)) return
  deallocate(CS%coordinateResolution)
  deallocate(CS%target_density)
  if (allocated(CS%max_interface_depths)) deallocate(CS%max_interface_depths)
  if (allocated(CS%max_layer_thickness)) deallocate(CS%max_layer_thickness)
  deallocate(CS)
end subroutine end_coord_hycom

!> This subroutine can be used to set the parameters for the coord_hycom module
subroutine set_hycom_params(CS, max_interface_depths, max_layer_thickness, &
                            use_z_only_in_ML, zspace_min_extent, MLB_extra_frac, MLB_extra_max, &
                            interior_min_thickness, interp_CS)
  type(hycom_CS),                 pointer    :: CS !< Coordinate control structure
  real, dimension(:),   optional, intent(in) :: max_interface_depths !< Maximum depths of interfaces [H ~> m or kg m-2]
  real, dimension(:),   optional, intent(in) :: max_layer_thickness  !< Maximum thicknesses of layers [H ~> m or kg m-2]
  logical,              optional, intent(in) :: use_z_only_in_ML !< If true, only use z-space in the surface mixed layer
  real,                 optional, intent(in) :: zspace_min_extent !< The minimum extent of z*-space region [H ~> m or kg m-2]
  real,                 optional, intent(in) :: MLB_extra_frac !< Extra fraction for MLB estimate [nondim]
  real,                 optional, intent(in) :: MLB_extra_max !< Maximum extra for MLB estimate [H ~> m or kg m-2]
  real,                 optional, intent(in) :: interior_min_thickness !< The minimum thickness to use between interior interfaces [H ~> m or kg m-2]
  type(interp_CS_type), optional, intent(in) :: interp_CS !< Controls for interpolation

  if (.not. associated(CS)) call MOM_error(FATAL, "set_hycom_params: CS not associated")

  if (present(max_interface_depths)) then
    if (size(max_interface_depths) /= CS%nk+1) &
      call MOM_error(FATAL, "set_hycom_params: max_interface_depths inconsistent size")
    allocate(CS%max_interface_depths(CS%nk+1))
    CS%max_interface_depths(:) = max_interface_depths(:)
  endif

  if (present(max_layer_thickness)) then
    if (size(max_layer_thickness) /= CS%nk) &
      call MOM_error(FATAL, "set_hycom_params: max_layer_thickness inconsistent size")
    allocate(CS%max_layer_thickness(CS%nk))
    CS%max_layer_thickness(:) = max_layer_thickness(:)
  endif

  if (present(use_z_only_in_ML)) CS%use_z_only_in_ML = use_z_only_in_ML
  if (present(zspace_min_extent)) CS%zspace_min_extent = zspace_min_extent
  if (present(MLB_extra_max)) CS%MLB_extra_max = MLB_extra_max
  if (present(MLB_extra_frac)) CS%MLB_extra_frac = MLB_extra_frac
  if (present(interior_min_thickness)) CS%interior_min_thickness = interior_min_thickness
  if (present(interp_CS)) CS%interp_CS = interp_CS
end subroutine set_hycom_params

!> Build a HyCOM coordinate column
subroutine build_hycom1_column(CS, eqn_of_state, nz, depth, h, T, S, p_col, &
                               z_col, z_col_new, zScale, h_neglect, h_neglect_edge)
  type(hycom_CS),        intent(in)    :: CS    !< Coordinate control structure
  type(EOS_type),        pointer       :: eqn_of_state !< Equation of state structure
  integer,               intent(in)    :: nz    !< Number of levels
  real,                  intent(in)    :: depth !< Depth of ocean bottom (positive [H ~> m or kg m-2])
  real, dimension(nz),   intent(in)    :: T     !< Temperature of column [degC]
  real, dimension(nz),   intent(in)    :: S     !< Salinity of column [ppt]
  real, dimension(nz),   intent(in)    :: h     !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(nz),   intent(in)    :: p_col !< Layer pressure [R L2 T-2 ~> Pa]
  real, dimension(nz+1), intent(in)    :: z_col !< Interface positions relative to the surface [H ~> m or kg m-2]
  real, dimension(CS%nk+1), intent(inout) :: z_col_new !< Absolute positions of interfaces [H ~> m or kg m-2]
  real, optional,        intent(in)    :: zScale !< Scaling factor from the input coordinate thicknesses in [Z ~> m]
                                                !! to desired units for zInterface, perhaps GV%Z_to_H.
  real, optional,        intent(in)    :: h_neglect !< A negligibly small width for the purpose of
                                                !! cell reconstruction [H ~> m or kg m-2]
  real, optional,        intent(in)    :: h_neglect_edge !< A negligibly small width for the purpose of
                                                !! edge value calculation [H ~> m or kg m-2]

  ! Local variables
  integer :: k
  real, dimension(nz) :: rho_col ! Layer densities in a column [R ~> kg m-3]
  real, dimension(CS%nk) :: h_col_new ! New layer thicknesses
  real :: z_scale    ! A scaling factor from the input thicknesses to the target thicknesses,
                     ! perhaps 1 or a factor in [H Z-1 ~> 1 or kg m-3]
  real :: stretching ! z* stretching, converts z* to z [nondim].
  real :: nominal_z ! Nominal depth of interface when using z* [H ~> m or kg m-2]
  real :: MLB ! Approximate base of the mixed layer [H ~> m or kg m-2]
  real :: NII ! Next interior interface below the mixed layer [H ~> m or kg m-2]
  logical :: maximum_depths_set ! If true, the maximum depths of interface have been set.
  logical :: maximum_h_set      ! If true, the maximum layer thicknesses have been set.
  integer :: KMB, KBZ, KNII ! Indices of specific interfaces: the first interior and the last z-region
  real :: a ! An interpolation weight [nondim]
  real :: dr ! Delta rho [R ~> kg m-3]

  maximum_depths_set = allocated(CS%max_interface_depths)
  maximum_h_set = allocated(CS%max_layer_thickness)

  z_scale = 1.0 ; if (present(zScale)) z_scale = zScale

  ! Work bottom recording potential density
  call calculate_density(T, S, p_col, rho_col, eqn_of_state)
  ! This ensures the potential density profile is monotonic
  ! although not necessarily single valued.
  do k = nz-1, 1, -1
    rho_col(k) = min( rho_col(k), rho_col(k+1) )
  enddo

  ! Interpolates for the target interface position with the rho_col profile
  ! Based on global density profile, interpolate to generate a new grid
  call build_and_interpolate_grid(CS%interp_CS, rho_col, nz, h(:), z_col, &
           CS%target_density, CS%nk, h_col_new, z_col_new, h_neglect, h_neglect_edge)

  ! Sweep down the interfaces bounding to z* space where needed
  stretching = z_col(nz+1) / depth ! Stretches z* to z
  if (CS%use_z_only_in_ML) then
 !  ! Find Mixed layer depth
 !  MLB = z_col(2)
 !  dr = rho_col(1) + 0.05
 !  find_MLD: do K = 2,nz
 !    if ( rho_col(k) >= dr ) then
 !      ! rho_col(k-1) - rho(cool(1) < 0.03
 !      ! rho_col(k) - rho(cool(1) >= 0.03
 !      a = ( dr - rho_col(k-1) ) / ( rho_col(k) - rho_col(k-1) )
 !      MLB = 0.5*( a*(z_col(K)+z_col(K+1)) - (1.0-a)*(z_col(K-1)+z_col(K)) )
 !      exit find_MLD
 !    endif
 !    MLB = z_col(K+1)
 !  enddo find_MLD
 !  ! Given MLB, find a NII and KMB in nominal rho grid
 !  NII = z_col_new(2)
 !  KMB = 1
 !  find_NII: do K = 2, CS%nk+1
 !    NII = z_col_new(K)
 !    KMB = K-1
 !    if (NII > MLB) exit find_NII
 !  enddo find_NII

    ! As a proxy for the ML depth, find the first isopycnal that is not snapped to the surface
    ! (Use the first nominal resolution as a threshold)
    MLB = z_col_new(1)
    NII = z_col_new(2)
    KMB = 1
    find_MLB: do K = 2, CS%nk
      MLB = z_col_new(K)
      NII = z_col_new(K+1)
      KMB = K
      if (MLB - z_col_new(1) > z_scale * CS%coordinateResolution(1) ) exit find_MLB
    enddo find_MLB

    ! At this point MLB is an approximate ML depth. Now we extend and add some constraints
    MLB = MLB + min(MLB * CS%MLB_extra_frac, CS%MLB_extra_max)
    MLB = max(MLB, CS%zspace_min_extent)
    NII = NII + min(NII * CS%MLB_extra_frac, CS%MLB_extra_max)
    NII = max(NII, CS%zspace_min_extent)
    MLB = min(MLB, NII)
    nominal_z = 0.
    KBZ = 1
    do K = 2, CS%nk+1
      ! Here, nominal_z is the depth an interface would have if it is in the z* region
      nominal_z = nominal_z + (z_scale * CS%coordinateResolution(k-1)) * stretching
      if ( z_col_new(K) < NII ) then ! For all interfaces above the NII
        if ( nominal_z <= MLB) then ! z*-space if the nominal z is in the mixed layer
          z_col_new(K) = nominal_z
          KBZ = min(K, KMB) ! Index to the last interface we put in z* region
                            ! (safety limit of KMB might not be needed)
        else
          a = real(K-KBZ) / real(KMB+1-KBZ) ! a=0 @ K=KBZ, a=1 @ K=KMB+1
          z_col_new(K) = (1.0-a)*MLB + a*NII
   !      z_col_new(K) = MLB
          z_col_new(K) = max( z_col_new(K), z_col_new(K-1) + zscale * CS%interior_min_thickness )
        endif
      else
        ! Otherwise we found the target isopycnal in the interior
        z_col_new(K) = max( z_col_new(K), z_col_new(K-1) + zscale * CS%interior_min_thickness )
      endif
   !  nominal_z = min( nominal_z, MLB )
   !  ! If the interpolated position was within the surface mixed-layer and we use z*-space
   !  z_col_new(K) = max( z_col_new(K), nominal_z )
   !  ! Otherwise we found the target isopycnal in the interior
   !  z_col_new(K) = max( z_col_new(K), z_col_new(K-1) + zscale * CS%interior_min_thickness )
      ! Make sure that the interface is not below the bottom topography
      z_col_new(K) = min( z_col_new(K), z_col(nz+1) )
    enddo
  else
    nominal_z = 0.
    do K = 2, CS%nk+1
      ! Make sure that the interface is at least as deep as a nominal target z* grid
      nominal_z = nominal_z + (z_scale * CS%coordinateResolution(k-1)) * stretching
      z_col_new(K) = max( z_col_new(K), nominal_z )
      ! Make sure that the interface is not below the bottom topography
      z_col_new(K) = min( z_col_new(K), z_col(nz+1) )
    enddo
  endif

  ! Optionally limit by a maximum depth or maximum thickness
  if (maximum_depths_set .and. maximum_h_set) then
    do K=2,CS%nk
      ! The loop bounds are 2 & nz so the top and bottom interfaces do not move.
      ! Recall that z_col_new is positive downward.
      z_col_new(K) = min( z_col_new(K), CS%max_interface_depths(K), &
                          z_col_new(K-1) + CS%max_layer_thickness(k-1) )
    enddo
  elseif (maximum_depths_set) then
    do K=2,CS%nk
      z_col_new(K) = min( z_col_new(K), CS%max_interface_depths(K) )
    enddo
  elseif (maximum_h_set) then
    do k=2,CS%nk
      z_col_new(K) = min( z_col_new(K), z_col_new(K-1) + CS%max_layer_thickness(k-1) )
    enddo
  endif

end subroutine build_hycom1_column

end module coord_hycom
