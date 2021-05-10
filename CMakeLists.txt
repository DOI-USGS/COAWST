# svn $Id: CMakeLists.txt 1053 2020-12-29 00:41:48Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2020 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Top-Level CMake File Definitions.

cmake_minimum_required( VERSION 3.12.0 FATAL_ERROR )

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/Compilers;${CMAKE_MODULE_PATH})
Message( STATUS "module path = ${CMAKE_MODULE_PATH}" )

# Figure out if we are using ecbuild.

if( DEFINED ENV{ECBUILD_MODULE_PATH} )
  set( IN_ECBUILD TRUE )
endif()

project(roms VERSION 3.9 LANGUAGES Fortran)

if( IN_ECBUILD )
  set( ECBUILD_DEFAULT_BUILD_TYPE Release )
  set( ENABLE_OS_TESTS           OFF CACHE BOOL "Disable OS tests" FORCE )
  set( ENABLE_LARGE_FILE_SUPPORT OFF CACHE BOOL "Disable testing of large file support" FORCE )
  set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )

  include( ecbuild_system NO_POLICY_SCOPE )

  ecbuild_requires_macro_version( 2.5 )
endif()

# Include ROMS cmake files.

include( roms_functions )
include( roms_compiler_flags )
include( roms_config )

###########################################################################
# Project
###########################################################################

if( IN_ECBUILD )
  ecbuild_declare_project()

  ecbuild_enable_fortran( REQUIRED )
endif()

set( ROMS_LINKER_LANGUAGE Fortran )

if( NOT IN_ECBUILD )

  # Set Fortran modules directory when not using ecbuild. It is set
  # automatically when using ecbuild.

  set( CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules" )
  set( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )
  set( CMAKE_SKIP_BUILD_RPATH FALSE )

  # MacOS requires an extra parameter for rpath to work properly. This is
  # set automatically when using ecbuild.

  if( APPLE )
    set( CMAKE_MACOSX_RPATH ON )
    set( CMAKE_INSTALL_RPATH @loader_path/../lib )
  else()
    set( CMAKE_INSTALL_RPATH \$ORIGIN/../lib )
  endif()

endif()

###########################################################################
# Dependencies
###########################################################################

# NetCDF Library.

if(DEFINED ENV{SINGULARITY_COMMAND} AND IN_ECBUILD )
  find_package( NetCDF REQUIRED COMPONENTS Fortran )
  include_directories( ${NETCDF_INCLUDE_DIRS} )
  list( APPEND netcdf_libs "netcdff" "netcdf" )

else()

  # Find NetCDF outside singularity. "target_link_libraries" is found below
  # because it has to be defined after "add_executable".  The link_netcdf()
  # function can be found in Compilers/roms_functions.cmake.

  link_netcdf()
  include_directories( ${netcdf_idir} )
  link_directories( ${netcdf_ldirs} )
endif()

# MPI Library for use with ecbuild.

if( IN_ECBUILD )
  ecbuild_find_mpi( COMPONENTS Fortran REQUIRED )
  ecbuild_include_mpi()
  link_libraries(${MPI_Fortran_LIBRARIES})
endif()

###########################################################################
# ROMS
###########################################################################

add_subdirectory( "Master" )
add_subdirectory( "ROMS/Drivers" )
add_subdirectory( "ROMS/Functionals" )
add_subdirectory( "ROMS/Modules" )
add_subdirectory( "ROMS/Nonlinear" )
add_subdirectory( "ROMS/Nonlinear/Biology" )
add_subdirectory( "ROMS/Nonlinear/Sediment" )
add_subdirectory( "ROMS/Utility" )

include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  "/usr/local/include"
  "$ENV{NETCDF}/include"
  "Master"
  "ROMS/Drivers"
  "ROMS/Functionals"
  "ROMS/Include"
  "ROMS/Modules"
  "ROMS/Nonlinear"
  "ROMS/Nonlinear/Biology"
  "ROMS/Nonlinear/Sediment"
  "ROMS/Utility"
)

link_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  "/usr/local/lib"
  "$ENV{NETCDF}/lib"
)

list(APPEND srcs
     ${Master_files}
     ${ROMS_Drivers_files}
     ${ROMS_Functionals_files}
     ${ROMS_Modules_files}
     ${ROMS_Nonlinear_files}
     ${ROMS_Nonlinear_Biology_files}
     ${ROMS_Nonlinear_Sediment_files}
     ${ROMS_Utility_files}
)

if( ADJOINT )
  add_subdirectory( "ROMS/Adjoint" )
  add_subdirectory( "ROMS/Adjoint/Biology" )
  include_directories( "ROMS/Adjoint" "ROMS/Adjoint/Biology" )
  list( APPEND srcs
        ${ROMS_Adjoint_files}
        ${ROMS_Adjoint_Biology_files}
  )

  # Free format flags for certain files that may exceed 132 characters per line

  set_property(
    SOURCE ${CMAKE_CURRENT_BINARY_DIR}/f90/ad_biology.f90
    APPEND_STRING PROPERTY COMPILE_FLAGS ${ROMS_FREEFLAGS}
  )
endif()

if( REPRESENTER )
  add_subdirectory( "ROMS/Representer" )
  add_subdirectory( "ROMS/Representer/Biology" )
  include_directories( "ROMS/Representer" "ROMS/Representer/Biology" )
  list( APPEND srcs
        ${ROMS_Representer_files}
        ${ROMS_Representer_Biology_files}
  )

  # Free format flags for certain files that may exceed 132 characters per line

  set_property(
    SOURCE ${CMAKE_CURRENT_BINARY_DIR}/f90/rp_biology.f90
    APPEND_STRING PROPERTY COMPILE_FLAGS ${ROMS_FREEFLAGS}
  )
endif()

if( TANGENT )
  add_subdirectory( "ROMS/Tangent" )
  add_subdirectory( "ROMS/Tangent/Biology" )
  include_directories( "ROMS/Tangent" "ROMS/Tangent/Biology" )
  list( APPEND srcs
        ${ROMS_Tangent_files}
        ${ROMS_Tangent_Biology_files}
  )

  # Free format flags for certain files that may exceed 132 characters per line

  set_property(
    SOURCE ${CMAKE_CURRENT_BINARY_DIR}/f90/tl_biology.f90
    APPEND_STRING PROPERTY COMPILE_FLAGS ${ROMS_FREEFLAGS}
  )
endif()

#--------------------------------------------------------------------------
# ROMS specific rules.
#--------------------------------------------------------------------------

# Special CPP definitions for "mod_strings.F".

set_property(
  SOURCE ROMS/Modules/mod_strings.F
  APPEND_STRING PROPERTY COMPILE_DEFINITIONS
    MY_OS='${my_os}' MY_CPU='${my_cpu}' MY_FORT='${my_fort}'
    MY_FC='${my_fc}'
    MY_FFLAGS='${my_fflags}'
)

# Free format flags for certain files that may exceed 132 characters per line

set_property(
  SOURCE
    ${CMAKE_CURRENT_BINARY_DIR}/f90/mod_ncparam.f90
    ${CMAKE_CURRENT_BINARY_DIR}/f90/mod_strings.f90
    ${CMAKE_CURRENT_BINARY_DIR}/f90/analytical.f90
    ${CMAKE_CURRENT_BINARY_DIR}/f90/biology.f90
  APPEND_STRING PROPERTY COMPILE_FLAGS ${ROMS_FREEFLAGS}
)

###########################################################################
# Clean-up the dependencies
###########################################################################


# Remove the .F files from "DependInfo.cmake" for the Objects target. If
# the "Objects" target changes names this MUST be updated.

add_custom_command(
  OUTPUT fixdep
  COMMAND ${PERL} ${CMAKE_CURRENT_SOURCE_DIR}/ROMS/Bin/FixDependInfo.pl ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Objects.dir/DependInfo.cmake
  VERBATIM
)

# Custom target needed to make the custom command above run

add_custom_target(fix ALL
  DEPENDS fixdep
)

# The "preprocess_fortran" function tells CMake how to create .f90 files.
# It allows CMake crippled Fortran depency tracker to look for the right
# modules using the pre-processed .f90 files instead of the .F files.
#
# The dependency generator is easily confused by CPP if-directives.
# This also allows us to feed .f90 files to CMake as the sources to all
# the ROMS code and eliminate compiler assumptions related to .F files.
#
# The function "preprocess_fortran" is defined in the "roms_functions.cmake"
# file in the "Compilers" directory.

preprocess_fortran( ${srcs} )
set( All_f90s "${f90srcs}" )

# Adding all the sources as an OBJECT library prevents CMake from trying
# to link an executable or collect the object files into a library.
#
# The libXXX.(a|so) libraries will be created below.

add_library( Objects OBJECT ${All_f90s} )

# Add dependency on custom target "fix" so "DependInfo.cmake" is cleaned up
# before the the dependencies are calculated.

add_dependencies( Objects fix )

# This is to set the install prefix to the build directory as long as the
# user has not set their own install prefix.

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "..." FORCE)
endif()

# Shared libraries need PIC. There is perfomance degradation to position
# independent code so we only set this if shared libraries are requested.
# This really has no effect when using ecbuild since Position Indpendent
# Code (PIC) is enabled by default inside ecbuild regardless of the type
# of library requested.

set( deplib )
if( LIBSHARED )
  message( STATUS "Enabling position independent code" )
  set_property(TARGET Objects PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library( ROMS_shared SHARED $<TARGET_OBJECTS:Objects> )
  set_target_properties(ROMS_shared PROPERTIES OUTPUT_NAME ROMS)
  install( TARGETS ROMS_shared LIBRARY )

  # This is needed to add to the depency list of the ROMS executable
  # since it is linking with libROMS.so instead of libROMS_shared.so
  # The executable doesn't see "ROMS" as a target generated by this
  # project so it starts linking before the library is fully assembled.

  list( APPEND deplib ROMS_shared )
endif()

if( LIBSTATIC )
  add_library( ROMS_static STATIC $<TARGET_OBJECTS:Objects> )
  set_target_properties(ROMS_static PROPERTIES OUTPUT_NAME ROMS)
  install( TARGETS ROMS_static ARCHIVE )

  # This is needed to add to the depency list of the ROMS executable
  # since it is linking with libROMS.a instead of libROMS_static.a
  # The executable doesn't see "ROMS" as a target generated by this
  # project so it starts linking before the library is fully assembled.

  list( APPEND deplib ROMS_static )
endif()

# There needs to be at least one source file for the "add_executable" line
# so "Master/master.F" is removed from "Master/CMakeLists.txt" and added
# manually here.
#
# I would have done all the files in "Master" here to simplify things
# but "ocean_control.F", "coupler.F", and "propagator.F" include
# modules needed in other parts of the code outside "Master".

if( ROMS_EXECUTABLE )

  set( CMAKE_BUILD_WITH_INSTALL_RPATH FALSE )

  message( STATUS "CMAKE_INSTALL_RPATH_USE_LINK_PATH = ${CMAKE_INSTALL_RPATH_USE_LINK_PATH}" )
  message( STATUS "CMAKE_SKIP_BUILD_RPATH = ${CMAKE_SKIP_BUILD_RPATH}" )
  message( STATUS "CMAKE_BUILD_WITH_INSTALL_RPATH = ${CMAKE_BUILD_WITH_INSTALL_RPATH}" )
  message( STATUS "CMAKE_INSTALL_RPATH = ${CMAKE_INSTALL_RPATH}" )

  if( APPLE )
    message( STATUS "CMAKE_MACOSX_RPATH = ${CMAKE_MACOSX_RPATH}" )
  endif()

  # If PARPACK/ARPACK is needed add directory for the linker. The
  # "link_directories" and "target_link_libraries" for PARPACK/ARPACK
  # are separated because "link_libraries" has to be set before
  # "add_executable" and "target_link_libraries" has to be set after.

  if( ARPACK )
    link_directories( ${PARPACK_LIBDIR} ${ARPACK_LIBDIR} )
  endif()

  preprocess_fortran( Master/master.F )
  set( master_f90 "${f90srcs}" )
  add_executable( "${BIN}" ${master_f90} )
  add_dependencies( "${BIN}" Objects ${deplib} )

  if( LIBSHARED AND LIBSTATIC )
    Message( STATUS "Both LIBSHARED and LIBSTATIC are enabled, the linker will probably default to the shared ROMS library." )
    target_link_libraries( "${BIN}" ROMS ${netcdf_libs} )
  elseif( LIBSHARED )
    Message( STATUS "Linking ROMS using shared ROMS library." )
    target_link_libraries( "${BIN}" ROMS ${netcdf_libs} )
  elseif( LIBSTATIC )
    Message( STATUS "Linking ROMS using static ROMS library." )
    target_link_libraries( "${BIN}" ROMS ${netcdf_libs} )
  else()
    Message( FATAL_ERROR "Unknown configuration!" )
  endif()

  # If ARPACK is needed, add to linking library list. The path to the library is set above.

  if( ARPACK )
    target_link_libraries( "${BIN}" parpack arpack )
  endif()

  install( TARGETS "${BIN}" RUNTIME )

endif()
