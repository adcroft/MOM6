program test_MOM_regridding

use MOM_regridding, only : regridding_unit_tests

if (regridding_unit_tests(.true.)) stop 1

end program test_MOM_regridding
