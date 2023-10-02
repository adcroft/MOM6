program test_MOM_string_functions

use MOM_string_functions,           only : string_functions_unit_tests

if ( string_functions_unit_tests(.true.) ) error &
   stop "test_MOM_string_functions: string_functions_unit_tests FAILED"

end program test_MOM_string_functions
