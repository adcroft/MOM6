! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Creates 1D reconstructions, by name, for a column of a given number of cells
!!
!! This module also provides the base class, Recon1d, and all of the Recon1d classes, so that it is
!! the single point of access to the 1D reconstructions for other code.
module Recon1d_factory

use MOM_error_handler, only : MOM_error, FATAL
use MOM_string_functions, only : uppercase
use Recon1d_type, only : Recon1d
use Recon1d_PCM, only : PCM
use Recon1d_PLM_CW, only : PLM_CW
use Recon1d_PLM_hybgen, only : PLM_hybgen
use Recon1d_PLM_hybgen_poly, only : PLM_hybgen_poly
use Recon1d_PLM_CWK, only : PLM_CWK
use Recon1d_MPLM_CWK, only : MPLM_CWK
use Recon1d_EMPLM_CWK, only : EMPLM_CWK
use Recon1d_MPLM_WA, only : MPLM_WA
use Recon1d_EMPLM_WA, only : EMPLM_WA
use Recon1d_MPLM_WA_poly, only : MPLM_WA_poly
use Recon1d_EMPLM_WA_poly, only : EMPLM_WA_poly
use Recon1d_PPM_CW, only : PPM_CW
use Recon1d_PPM_hybgen, only : PPM_hybgen
use Recon1d_PPM_hybgen_v1, only : PPM_hybgen_v1
use Recon1d_PPM_CWK, only : PPM_CWK
use Recon1d_EPPM_CWK, only : EPPM_CWK
use Recon1d_PPM_H4_2019, only : PPM_H4_2019
use Recon1d_PPM_H4_2018, only : PPM_H4_2018
use Recon1d_PLM_WLS, only : PLM_WLS

implicit none ; private

public new_recon1d
public Recon1d
public PCM, PLM_CW, PLM_hybgen, PLM_hybgen_poly, PLM_CWK, MPLM_CWK, EMPLM_CWK
public MPLM_WA, EMPLM_WA, MPLM_WA_poly, EMPLM_WA_poly, PLM_WLS
public PPM_CW, PPM_hybgen, PPM_hybgen_v1, PPM_CWK, EPPM_CWK, PPM_H4_2019, PPM_H4_2018

contains

!> Allocates and initializes the reconstruction named by scheme for a column of n cells
!!
!! The classes that reproduce OM4-era schemes are replaced, for short columns, by the classes
!! that reproduce the simpler schemes that build_reconstructions_1d() uses in their place.
subroutine new_recon1d(recon, scheme, n, h_neglect)
  class(Recon1d), pointer, intent(out) :: recon  !< The new reconstruction
  character(len=*),        intent(in)  :: scheme !< The name of the reconstruction, e.g. "C_PPM_CW"
  integer,                 intent(in)  :: n      !< Number of cells in the column
  real,          optional, intent(in)  :: h_neglect !< A negligibly small width used in cell
                                                    !! reconstructions [H]
  ! Local variables
  character(len=12) :: n_str ! n, as a label

  select case ( uppercase(trim(scheme)) )
    case ("C_PCM")
      allocate( PCM :: recon )
    case ("C_PLM_CW")
      allocate( PLM_CW :: recon )
    case ("C_PLM_HYBGEN")
      allocate( PLM_hybgen :: recon )
    case ("C_PLM_HYBGEN_POLY")
      if (n <= 3) then ! As for PLM_HYBGEN in build_reconstructions_1d()
        allocate( MPLM_WA_poly :: recon )
      else
        allocate( PLM_hybgen_poly :: recon )
      endif
    case ("C_MPLM_WA")
      allocate( MPLM_WA :: recon )
    case ("C_EMPLM_WA")
      allocate( EMPLM_WA :: recon )
    case ("C_MPLM_WA_POLY")
      allocate( MPLM_WA_poly :: recon )
    case ("C_EMPLM_WA_POLY")
      allocate( EMPLM_WA_poly :: recon )
    case ("C_PLM_CWK")
      allocate( PLM_CWK :: recon )
    case ("C_MPLM_CWK")
      allocate( MPLM_CWK :: recon )
    case ("C_EMPLM_CWK")
      allocate( EMPLM_CWK :: recon )
    case ("C_PPM_CW")
      allocate( PPM_CW :: recon )
    case ("C_PPM_HYBGEN")
      allocate( PPM_hybgen :: recon )
    case ("C_PPM_HYBGEN_V1")
      if (n <= 3) then ! As for PPM_HYBGEN in build_reconstructions_1d()
        allocate( MPLM_WA_poly :: recon )
      elseif (n == 4) then
        allocate( PPM_H4_2019 :: recon )
      else
        allocate( PPM_hybgen_v1 :: recon )
      endif
    case ("C_PPM_CWK")
      allocate( PPM_CWK :: recon )
    case ("C_EPPM_CWK")
      allocate( EPPM_CWK :: recon )
    case ("C_PPM_H4_2019")
      if (n <= 3) then ! As for PPM_H4 in build_reconstructions_1d()
        allocate( MPLM_WA_poly :: recon )
      else
        allocate( PPM_H4_2019 :: recon )
      endif
    case ("C_PPM_H4_2018")
      if (n <= 3) then ! As for PPM_H4 in build_reconstructions_1d()
        allocate( MPLM_WA_poly :: recon )
      else
        allocate( PPM_H4_2018 :: recon )
      endif
    case ("C_PLM_WLS")
      allocate( PLM_WLS :: recon )
    case default
      call MOM_error(FATAL, "Recon1d_factory, new_recon1d: "//&
           "Unrecognized reconstruction scheme ("//trim(scheme)//").")
  end select

  if (n < recon%min_cells()) then
    write(n_str,'(i0)') n
    call MOM_error(FATAL, "Recon1d_factory, new_recon1d: "//trim(scheme)//&
         " can not be used for a column of "//trim(n_str)//" cells.")
  endif

  call recon%init(n, h_neglect=h_neglect)

end subroutine new_recon1d

end module Recon1d_factory
