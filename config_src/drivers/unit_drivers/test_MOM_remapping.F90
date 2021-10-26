!> Tests MOM_remapping module
program test_MOM_remapping

! This file is part of MOM6. See LICENSE.md for the license.
! Original module written by Laurent White, 2008.06.09

use MOM_remapping, only : remapping_unit_tests
use, intrinsic :: iso_fortran_env, only : stdout=>output_unit, stderr=>error_unit

implicit none

if ( remapping_unit_tests(.true.) ) stop 'remapping_unit_tests failed'

end program test_MOM_remapping
