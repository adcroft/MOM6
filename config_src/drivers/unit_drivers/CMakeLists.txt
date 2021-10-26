cmake_minimum_required(VERSION 3.10)

# set the project name
project(
  "MOM6-unit-testing"
  LANGUAGES "C" "Fortran"
)

find_package(MPI REQUIRED)
add_definitions(${MPI_Fortran_COMPILE_FLAGS})
include_directories(${MPI_Fortran_INCLUDE_PATH})
link_directories(${MPI_Fortran_LIBRARIES})
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

if("${NETCDF_DIR}" STREQUAL "")
  execute_process(COMMAND nf-config --prefix OUTPUT_VARIABLE NETCDF_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")
find_package(NetCDF REQUIRED COMPONENTS F90 C)
include_directories(${NETCDF_F90_INCLUDE_DIRS})

message("Using Fortran: ${CMAKE_Fortran_COMPILER_ID} and C: ${CMAKE_C_COMPILER_ID}")

if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS
    "${CMAKE_Fortran_FLAGS} -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -i4 -r8 -nowarn -sox -traceback")
  set(CMAKE_Fortran_FLAGS_RELEASE
    "-O2 -debug minimal -fp-model source")
  set(CMAKE_Fortran_FLAGS_DEBUG
    "-g -O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fpe0 -ftrapuv")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS
    "${CMAKE_Fortran_FLAGS} -fcray-pointer -fdefault-real-8 -fdefault-double-8 -Waliasing -ffree-line-length-none -fno-range-check")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fbounds-check")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -W -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow")
endif()

if(CMAKE_C_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -sox -traceback -O2 -debug minimal")
elseif(CMAKE_C_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
endif()

message("Using Fortran flags: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}")
message("NetCDF libs: ${NETCDF_LIBRARIES}")

add_definitions(-Duse_libMPI -Duse_netCDF -Duse_LARGEFILE -DSPMD -D__IFC)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# FMS library
set(FMS_DIR "../../../.testing/deps/fms/src")
file(GLOB_RECURSE FMS_SOURCE LIST_DIRECTORIES false
  ${FMS_DIR}/*.[fF]90 ${FMS_DIR}/*.c)
list(FILTER FMS_SOURCE EXCLUDE REGEX ".*test_.*.F90$")
add_library(FMS ${FMS_SOURCE})
target_include_directories(FMS PRIVATE
  /usr/include ${FMS_DIR}/include ${FMS_DIR}/mosaic ${FMS_DIR}/drifters ${FMS_DIR}/fms ${FMS_DIR}/mpp/include)

# common files for many builds
set(MOM6_DIR "../../../")
file(GLOB_RECURSE MOM6_SOURCE LIST_DIRECTORIES false FOLLOW_SYMLINKS
  ${MOM6_DIR}/src/*.[fF]90 ${MOM6_DIR}/src/*.c ${MOM6_DIR}/config_src/infra/FMS1/*.[fF]90 ${MOM6_DIR}/config_src/external/*.[fF]90 )
set(ALL_INCLUDE_DIRS 
  ${MOM6_DIR}/config_src/memory/dynamic_symmetric ${MOM6_DIR}/src/framework
  /usr/include ${FMS_DIR}/include ${FMS_DIR}/mosaic ${FMS_DIR}/drifters ${FMS_DIR}/fms ${FMS_DIR}/mpp/include)

# MOM6 library
add_library(MOM6 ${MOM6_SOURCE})
target_include_directories(MOM6 PRIVATE ${ALL_INCLUDE_DIRS})
target_link_libraries(MOM6 FMS) 

# MOM6-solo = ocean only executable
file(GLOB_RECURSE MOM6_SOLO_SOURCE LIST_DIRECTORIES false
  ${MOM6_DIR}/config_src/drivers/solo_driver/*.[fF]90)
add_executable(MOM6-solo ${MOM6_SOLO_SOURCE})
target_include_directories(MOM6-solo PRIVATE
  ${ALL_INCLUDE_DIRS})
target_link_libraries(MOM6-solo MOM6 FMS ${NETCDF_LIBRARIES} ${MPI_Fortran_LIBRARIES})

# MOM6_sum_driver
add_executable(MOM_sum_driver MOM_sum_driver.F90)
target_include_directories(MOM_sum_driver PRIVATE ${ALL_INCLUDE_DIRS})
target_link_libraries(MOM_sum_driver MOM6 FMS ${NETCDF_LIBRARIES} ${MPI_Fortran_LIBRARIES})

# test_MOM_remapping
add_executable(test_MOM_remapping test_MOM_remapping.F90)
target_include_directories(test_MOM_remapping PRIVATE ${ALL_INCLUDE_DIRS})
target_link_libraries(test_MOM_remapping MOM6 FMS ${NETCDF_LIBRARIES} ${MPI_Fortran_LIBRARIES})
