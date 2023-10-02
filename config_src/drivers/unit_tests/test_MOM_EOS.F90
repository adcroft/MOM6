program test_MOM_EOS

use MOM_EOS, only : EOS_unit_tests

if ( EOS_unit_tests(.true.) ) error stop "test_MOM_EOS: EOS_unit_tests FAILED"

end program test_MOM_EOS
