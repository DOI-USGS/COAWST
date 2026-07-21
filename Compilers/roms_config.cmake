# git $Id$
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2026 The ROMS Group                                :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# CMake configuration for any ROMS application.

string( TOLOWER "${ROMS_APP}.h" ROMS_APP_HEADER )
add_compile_definitions( ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}" ${ROMS_APP} HEADER="${ROMS_APP_HEADER}" )

# Do you want a shared or static library?

if( NOT DEFINED LIBTYPE )
  option( LIBSHARED "Build ROMS as a Shared Library" OFF )
  option( LIBSTATIC "Build ROMS as a Static Library" ON )
  Message( STATUS "No LIBTYPE chosen, defaulting to static. Valid types are SHARED, STATIC, and BOTH" )
else()
  string( TOLOWER "${LIBTYPE}" libtype )
  if( ${libtype} STREQUAL "shared" )
    option( LIBSHARED "Build ROMS as a Shared Library" ON )
    option( LIBSTATIC "Build ROMS as a Static Library" OFF )
    Message( STATUS "Building with shared ROMS library." )
  elseif( ${libtype} STREQUAL "static" )
    option( LIBSHARED "Build ROMS as a Shared Library" OFF )
    option( LIBSTATIC "Build ROMS as a Static Library" ON )
    Message( STATUS "Building with static ROMS library." )
  elseif( ${libtype} STREQUAL "both" )
    option( LIBSHARED "Build ROMS as a Shared Library" ON )
    option( LIBSTATIC "Build ROMS as a Static Library" ON )
    Message( STATUS "Building both static and shared ROMS libraries." )
  else()
    Message( FATAL_ERROR "Invalid value for LIBTYPE. Valid types are SHARED, STATIC, and BOTH" )
  endif()
endif()

# Do you want a ROMS executable?

if( NOT DEFINED ROMS_EXECUTABLE )
  option( ROMS_EXECUTABLE "Turn on/off building ROMS executable" ON )
  Message( STATUS "ROMS_EXECUTABLE not set, defaulting to ON. The ROMS executable will be built." )
elseif( ROMS_EXECUTABLE )
  option( ROMS_EXECUTABLE"Turn on/off building ROMS executable" ON )
  Message( STATUS "The ROMS executable will be built." )
elseif( NOT ROMS_EXECUTABLE )
  option( ROMS_EXECUTABLE "Turn on/off building ROMS executable" OFF )
  Message( STATUS "The ROMS executable will NOT be built." )
else()
  Message( FATAL_ERROR "Invalid value for ROMS_EXECUTABLE. Valid values are ON and OFF" )
endif()

# Where is your <application>.h file?

set( HEADER_DIR  ${MY_HEADER_DIR} )

# Any custom analytical files?

if( NOT DEFINED ENV{MY_ANALYTICAL_DIR} )
  set( ANALYTICAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ROMS/Functionals" )
  add_compile_definitions( ANALYTICAL_DIR="${ANALYTICAL_DIR}" )
  # no need for "include_directories" since it is already included
else()
  include_directories( $ENV{MY_ANALYTICAL_DIR} )
endif()

# Location(s) of ARPACK/PARPACK libraries. If both libarpack.a and libparpack.a
# both live in the same directory, leave ARPACK_LIBDIR blank.
#
# The decision about whether to use them in linking is computed below.
# This CMake setup will NOT build ARPACK/PARPACK for you.

if( DEFINED PARPACK_LIBDIR )
  set( PARPACK_LIBDIR "${PARPACK_LIBDIR}" )
else()
  set( PARPACK_LIBDIR "" )
endif()
if( DEFINED ARPACK_LIBDIR )
  set( ARPACK_LIBDIR "${ARPACK_LIBDIR}" )
else()
  set( ARPACK_LIBDIR "" )
endif()

Message( STATUS "PARPACK_LIBDIR = ${PARPACK_LIBDIR}" )
Message( STATUS "ARPACK_LIBDIR  = ${ARPACK_LIBDIR}" )

set( ROMS_HEADER ${HEADER_DIR}/${ROMS_APP_HEADER} )

add_compile_definitions( ROMS_HEADER="${ROMS_HEADER}" )

# Set ROMS Executable Name.

if( ${CMAKE_BUILD_TYPE} MATCHES "Debug" )
  set( BIN "romsG" )
  if( MPI )
    add_compile_definitions( MPI )
  endif()
elseif( MPI )
  set( BIN "romsM" )
  add_compile_definitions( MPI )
else()
  set( BIN "romsS" )
endif()

if( MY_CPP_FLAGS )
  foreach( flag ${MY_CPP_FLAGS} )
    add_compile_definitions( ${flag} )
  endforeach()
  # get_options is found in roms_functions.cmake and populates the $defs
  # variable used below. 
  get_options( ${ROMS_HEADER} ${MY_CPP_FLAGS} )
else()
  get_options( ${ROMS_HEADER} )
endif()

# If use_4dvar returns a line containing "CPPDEFS" then somehow ROMS_HEADER
# was not defined and the build script configuation or CMake command should
# be checked.

if( "${defs}" MATCHES "CPPDEFS" )
  message( FATAL_ERROR "The C-preprocessor was unable to find your application's header file please check you configuation" )
endif()

if( "${defs}" MATCHES "ADJOINT" )
  option( ADJOINT "Turn on/off Adjoint Model" ON )
  Message( STATUS "ROMS Adjoint Model ENABLED" )
endif()

if( "${defs}" MATCHES "REPRESENTER" )
  option( REPRESENTER "Turn on/off Representer Model" ON )
  Message( STATUS "ROMS Representer Model ENABLED" )
endif()

if( "${defs}" MATCHES "TANGENT" )
  option( TANGENT "Turn on/off Tangent Linear Model" ON )
  message( STATUS "ROMS Tangent Linear Model ENABLED" )
endif()

if( "${defs}" MATCHES "ARPACK" )
  option( ARPACK "ARPACK/PARPACK Library" ON )
  message( STATUS "ROMS Link With ARPACK/PARPACK ENABLED" )
endif()

# Locations of PIO and PNetCDF libraries.
#
# The decision about whether to use them in linking is computed below.
# This CMake setup will NOT build PIO or PNetCDF for you.

if( "${defs}" MATCHES "PIO" )
  option( PIO "Link with PIO Libraries" ON )
  Message( STATUS "ROMS link With PIO ENABLED" )
  if ( DEFINED PIO_LIBDIR AND DEFINED PIO_INCDIR )
    set( PIO_LIBDIR "${PIO_LIBDIR}" )
    set( PIO_INCDIR "${PIO_INCDIR}" )
    Message( STATUS "    PIO_LIBDIR = ${PIO_LIBDIR}" )
    Message( STATUS "    PIO_INCDIR = ${PIO_INCDIR}" )

    if( DEFINED PNETCDF_LIBDIR AND DEFINED PNETCDF_INCDIR )
      set( PNETCDF_LIBDIR "${PNETCDF_LIBDIR}" )
      set( PNETCDF_INCDIR "${PNETCDF_INCDIR}" )
      Message( STATUS "PNETCDF_LIBDIR = ${PNETCDF_LIBDIR}" )
      Message( STATUS "PNETCDF_INCDIR = ${PNETCDF_INCDIR}" )
    else()
      set( PNETCDF_LIBDIR "" )
      set( PNETCDF_INCDIR "" )
      Message( STATUS "Using System Library and Include Locations for PNetCDF" )
    endif()
  else()
    set( PIO_LIBDIR "" )
    set( PIO_INCDIR "" )
    Message( STATUS "Using System Library and Include Locations for PIO" )
  endif()
endif()

