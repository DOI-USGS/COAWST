cmake_minimum_required(VERSION 3.18)

project(
  ufs_fire_behavior 
  VERSION 0.0.1
  LANGUAGES Fortran CXX C )

# user defined cmake modules
list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/")

# build configuration (cmake/modules/FireBuildSettings.cmake)
include(FireBuildSettings)
set_fire_defaults()
set_fire_cache()
print_fire_settings()

# preprocessor options
set(CMAKE_CPP_FLAGS "-nostdinc -C -P -w")

# compiler options
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -traceback -O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -traceback -O0 -fno-inline -no-ip -check all -fpe0 -ftrapuv -init=snan,arrays")
  set(CMAKE_Fortran_FLAGS "-ip -fp-model precise -w -ftz -align all -fno-alias -FR -convert big_endian")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -fbacktrace -02")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -fbacktrace -O0 -Wall -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan")
  set(CMAKE_Fortran_FLAGS "-std=f2008 -ffree-form -ffree-line-length-none -fall-intrinsics")
else()
  message(WARNING "${CMAKE_Fortran_COMPILER_ID} Fortran compiler will be used with default options")
endif()

# fortran modules directory
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION ${CMAKE_INSTALL_PREFIX})

# third party libraries
if(${CMAKE_PROJECT_NAME} STREQUAL "WRF")
  find_package(netCDF-Fortran REQUIRED)
else()
  find_package(NetCDF REQUIRED Fortran)
endif()

# turn on DM_PARALLEL preprocessor directive
if(DM_PARALLEL)
  find_package(MPI REQUIRED)
  add_compile_definitions(firelib PUBLIC DM_PARALLEL=1)
endif()

if(OPENMP)
  find_package(OpenMP REQUIRED)
  message(STATUS "OpenMP enabled")
endif()

#------------------------------------------------------------------------------
# Add source files

# State
list(APPEND _state_files state/state_mod.F90
			 state/ignition_line_mod.F90
	                 state/tiles_mod.F90)

list(APPEND _share_files share/proj_lc_mod.F90
                         share/datetime_mod.F90
			 share/emis_mod.F90
			 share/fuel_mod.F90
			 share/ros_mod.F90
			 share/fmc_mod.F90
			 share/interp_mod.F90
			 share/mpi_mod.F90
                         share/constants_mod.F90)

# IO
list(APPEND _io_files io/namelist_mod.F90
                      io/coupling_mod.F90
                      io/geogrid_mod.F90
                      io/wrf_mod.F90
                      io/wrfdata_mod.F90
		      io/netcdf_mod.F90
                      io/stderrout_mod.F90)


# WRF Physics
list(APPEND _wrffire_physics physics/level_set_mod.F90
                             physics/fire_driver_mod.F90
			     physics/fire_model_mod.F90
			     physics/fire_physics_mod.F90
			     physics/fmc_wrffire_mod.F90
			     physics/ros_wrffire_mod.F90
			     physics/fuel_anderson_mod.F90)

# Driver
list(APPEND _driver_files driver/advance_mod.F90
                          driver/initialize_mod.F90)
# note that driver/fire_behavior.F90 is binded below in add_executable(...
# that's because these lists will build the library firelib, which is 
# used by both, standalone and exmx.  

#------------------------------------------------------------------------------

# configure fire behavior library
add_library(firelib STATIC ${_state_files} ${_share_files} ${_io_files} ${_wrffire_physics} ${_driver_files})
if(${CMAKE_PROJECT_NAME} STREQUAL "WRF")
  target_link_libraries(firelib PUBLIC netCDF::netcdff )
else()
  target_include_directories(firelib PUBLIC ${NetCDF_INCLUDE_DIRS})
  target_link_libraries(firelib PUBLIC $<BUILD_INTERFACE:${NetCDF_LIBRARIES}>)
endif()
if (DM_PARALLEL)
  target_link_libraries(firelib PUBLIC $<BUILD_INTERFACE:MPI::MPI_Fortran>)
endif()
if(OPENMP)
  target_link_libraries(firelib PUBLIC $<BUILD_INTERFACE:OpenMP::OpenMP_Fortran>)
endif()
set_target_properties( 
                      firelib
                      PROPERTIES
                        Fortran_FORMAT           FREE
                      )

# configure fire behavior executable
add_executable(fire_behavior.exe driver/fire_behavior.F90)
add_dependencies(fire_behavior.exe firelib)
target_link_libraries(fire_behavior.exe PUBLIC firelib)
if(NOT ${CMAKE_PROJECT_NAME} STREQUAL "WRF")
  target_link_libraries(fire_behavior.exe PUBLIC NetCDF::NetCDF_Fortran)
endif()
install(TARGETS fire_behavior.exe DESTINATION bin)

if(NUOPC OR ESMX)
  add_subdirectory(nuopc)
endif()

add_subdirectory(tests)
