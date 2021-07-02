!> Calculates horizontal viscosity and viscous stresses
module MOM_CNN

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_checksums,             only : hchksum, Bchksum
use MOM_coms,                  only : min_across_PEs
use MOM_diag_mediator,         only : post_data, register_diag_field
use MOM_diag_mediator,         only : diag_ctrl, time_type
use MOM_domains,               only : pass_var, pass_vector, To_ALL, AGRID
use MOM_error_handler,         only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,           only : get_param, log_version, param_file_type
use MOM_grid,                  only : ocean_grid_type

implicit none ; private

#include <MOM_memory.h>

public apply_CNN, CNN_init, CNN_end

!> A single convolution layer
type :: cnv_layer_2d
  integer :: width=0 !< Stencil width of layer (should be odd)
  real, dimension(:,:), allocatable :: uuweights !< Network weights
  real, dimension(:,:), allocatable :: uvweights !< Network weights
  real, dimension(:,:), allocatable :: vuweights !< Network weights
  real, dimension(:,:), allocatable :: vvweights !< Network weights
  real :: ubias !< Offset for layer
  real :: vbias !< Offset for layer
end type cnv_layer_2d

!> Control structure for CNN
type, public :: CNN_CS ; private
  integer :: n_layers=0 !< Number of layers
  type(cnv_layer_2d), allocatable :: cnv_layers(:) !< List of convolution layers
end type CNN_CS

contains

!> Call a CNN in forward mode
subroutine apply_CNN(CS, G, u, v, out_u, out_v)
  type(CNN_CS),                      intent(in)    :: CS    !< A CNN
  type(ocean_grid_type),             intent(inout) :: G     !< The ocean's grid structure
  real, dimension(SZIB_(G),SZJ_(G)), intent(in)    :: u     !< The zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G)), intent(in)    :: v     !< The meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G)), intent(inout) :: out_u !< The output at u-points [unspecified]
  real, dimension(SZI_(G),SZJB_(G)), intent(inout) :: out_v !< The output at v-points [unspecified]
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    u_input, &      ! u at h-points [L-1 T-1 ~> m-1 s-1]
    v_input         ! v at h-points [L-1 T-1 ~> m-1 s-1]

  integer :: i, j, layer, hstencil

  ! Interpolate to cell-centers so u,v are co-located
  do j=G%jsc,G%jec
    do i=G%isc,G%iec
      u_input(i,j) = 0.5 * ( u(I-1,j) + u(I,j) )
      v_input(i,j) = 0.5 * ( v(i,J-1) + v(i,J) )
    enddo
  enddo

  ! Pass data through each layer
  do layer = 1, CS%n_layers
    ! Fill in halos
    call pass_vector(u_input, v_input, G%Domain, To_All, AGRID)

    hstencil = ( CS%cnv_layers(layer)%width - 1 ) / 2
    do j=G%jsc,G%jec
      do i=G%isc,G%iec
        call apply_layer(CS%cnv_layers(layer)%width, CS%cnv_layers(layer), &
                         u_input(i-hstencil:i+hstencil, j-hstencil:j+hstencil), &
                         v_input(i-hstencil:i+hstencil, j-hstencil:j+hstencil), &
                         out_u(i,j), out_v(i,j))
      enddo
    enddo

    ! Copy output of layer into input for next layer
    do j=G%jsc,G%jec
      do i=G%isc,G%iec
      u_input(i,j) = out_u(i,j)
      v_input(i,j) = out_v(i,j)
      enddo
    enddo
    
  enddo

  ! Interpolate back to u,v points
  call pass_vector(out_u, out_v, G%Domain, To_All, AGRID)
  do j=G%jsc,G%jec
    do I=G%isc,G%iec
      out_u(I,j) = 0.5 * ( u(I,j) + u(I,j+1) )
    enddo
  enddo
  do J=G%jsc-1,G%jec
    do i=G%isc,G%iec
      out_v(i,J) = 0.5 * ( v(i,J) + v(i,J+1) )
    enddo
  enddo

end subroutine apply_CNN

!> Apply a layer
subroutine apply_layer(stencil, layer, uin, vin, uout, vout)
  integer,                          intent(in)    :: stencil  !< Stencil of layer
  type(cnv_layer_2d),               intent(in)    :: layer    !< Convolution layer
  real, dimension(stencil,stencil), intent(in)    :: uin      !< Inputs u
  real, dimension(stencil,stencil), intent(in)    :: vin      !< Inputs v
  real,                             intent(out)   :: uout     !< Output u
  real,                             intent(out)   :: vout     !< Output v
  ! Local variables
  integer :: i, j, m

  ! Start with bias
  uout = layer%ubias
  vout = layer%vbias
  ! Tensor multiply weights x inputs
  do j = 1, stencil
    do i = 1, stencil
      uout = uout + layer%uuweights(i,j) * uin(i,j)
      vout = vout + layer%vuweights(i,j) * uin(i,j)
      uout = uout + layer%uvweights(i,j) * vin(i,j)
      vout = vout + layer%vvweights(i,j) * vin(i,j)
    enddo
  enddo
  ! ReLU activation function
  uout = max(0., uout)
  vout = max(0., vout)
end subroutine apply_layer

!> Read parameters and allocate a CNN
subroutine CNN_init(CS, param_file)
  type(CNN_CS),            intent(inout) :: CS   !< CNN control structure
  type(param_file_type),   intent(in)    :: param_file !< Run-time parameters parser
  ! Local variables
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_CNN"  ! module name

  !CS%diag => diag

  ! Read parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "N_LAYERS", CS%n_layers, &
                 "Number of layers in convolution neural network")
end subroutine CNN_init

!> De-allocate any variables allocated in CNN_CS
subroutine CNN_end(CS)
  type(CNN_CS), intent(inout) :: CS !< The control structure returned by a previous call to CNN_init

end subroutine CNN_end

!> \namespace mom_hor_visc
!!
!! About this module
!!
!! ...

end module MOM_CNN
