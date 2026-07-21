# a CMake script to build SWAN
#
# Copyright (C) 2026  Delft University of Technology

# use of CMake version 3.20 or greater
cmake_minimum_required( VERSION 3.20 )

# by default, the SWAN package will be installed in user's local or home folder
if( WIN32 )
  cmake_path( SET LOCALSWAN $ENV{LOCALAPPDATA}/Programs/wavemodels/swan )
  set( CMAKE_INSTALL_PREFIX "${LOCALSWAN}" CACHE STRING "" )
else()
  set( CMAKE_INSTALL_PREFIX "$ENV{HOME}/wavemodels/swan" CACHE STRING "" )
endif()

# Fortran project for SWAN
project( SWAN Fortran )

# set the SWAN version
set( VERSION 41.51 )

# enable common build
set( CMAKE_BUILD_TYPE None CACHE STRING "" FORCE )

# set the executable names
if( WIN32 )
   set( EXE swan )
   set( HEXE hcat )
else()
   set( EXE swan.exe )
   set( HEXE hcat.exe )
endif()

# set the Fortran77 extension
if( WIN32 )
   set( EXTF for )
else()
   set( EXTF f   )
endif()

# extra build options for the user
option( OPENMP "" OFF )
option( MPI "" OFF )
option( JAC "" OFF )
option( METIS "" OFF )
option( NETCDF "" OFF )

# define bin/lib//mod/src directories
set( BIN ${CMAKE_BINARY_DIR}/bin )
set( LIB ${CMAKE_BINARY_DIR}/lib )
set( MOD ${CMAKE_BINARY_DIR}/mod )
set( SRC ${CMAKE_SOURCE_DIR}/src_coawst )

# have the bin and lib folders placed in the build folder
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN} )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LIB} )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB} )

# have the modules placed in the mod folder
set( CMAKE_Fortran_MODULE_DIRECTORY ${MOD} )

# add find modules to the module path
list( APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" )

# Perl is required for setting switches in the source code
find_package( Perl REQUIRED )

# find Metis
if( METIS )
  find_package( Metis REQUIRED )
endif()

# find netCDF
if( NETCDF )
  # netCDF Fortran library is required
  find_package( NetCDF REQUIRED COMPONENTS Fortran )
endif()

# add the source to the library and place it in the lib folder
add_subdirectory( ${SRC} ${LIB} )

# create the SWAN executable and place it in the bin folder
add_subdirectory( ${SRC} ${BIN} )

# if desired, create the HCAT executable
add_subdirectory( ${SRC}/hcat )
