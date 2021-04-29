# svn $Id: roms_config.cmake 1054 2021-03-06 19:47:12Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2021 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# CMake configuration for any ROMS application.

string( TOLOWER "${APP}.h" HEADER )
string( TOLOWER "${APP}" application )
add_definitions( -DROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}" -D${APP} -DHEADER="${HEADER}" )

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

if( NOT DEFINED MY_ANALYTICAL_DIR )
  set( ANALYTICAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ROMS/Functionals" )
  add_definitions( -DANALYTICAL_DIR="${ANALYTICAL_DIR}" )
  # no need for "include_directories" since it is already included
else()
  add_definitions( -DANALYTICAL_DIR="${MY_ANALYTICAL_DIR}" )
  include_directories( ${MY_ANALYTICAL_DIR} )
endif()

# Location(s) of ARPACK/PARPACK libraries. If both libarpack.a and libparpack.a
# both live in the same directory, leave ARPACK_LIBDIR blank.
#
# The decision about whether to use them in linking is computed below.
# This CMake setup will NOT build ARPACK/PARPACK for you.

if(DEFINED PARPACK_LIBDIR)
  set( PARPACK_LIBDIR "${PARPACK_LIBDIR}" )
else()
  set( PARPACK_LIBDIR "" )
endif()
if(DEFINED ARPACK_LIBDIR)
  set( ARPACK_LIBDIR "${ARPACK_LIBDIR}" )
else()
  set( ARPACK_LIBDIR "" )
endif()

Message( STATUS "PARPACK_LIBDIR = ${PARPACK_LIBDIR}" )
Message( STATUS "ARPACK_LIBDIR  = ${ARPACK_LIBDIR}" )

# Set ROMS SVN repository information.

set( SVN_URL "${MY_SVN_URL}" )
set( SVN_REV "${MY_SVN_REV}" )

set( ROMS_HEADER ${HEADER_DIR}/${HEADER} )

add_definitions(
  -DROMS_HEADER="${ROMS_HEADER}"
  -DSVN_URL="${SVN_URL}"
  -DSVN_REV="${SVN_REV}"
)

# Set ROMS Executable Name.

if( ${CMAKE_BUILD_TYPE} MATCHES "Debug" )
  set( BIN "romsG" )
  if( MPI )
    add_definitions( -DMPI )
  endif()
elseif( MPI )
  set( BIN "romsM" )
  add_definitions( -DMPI )
else()
  set( BIN "romsS" )
endif()

if( MY_CPP_FLAGS )
  foreach( flag ${MY_CPP_FLAGS} )
    add_definitions( -D${flag} )
  endforeach()
  use_4dvar( ${ROMS_HEADER} ${MY_CPP_FLAGS} )
else()
  use_4dvar( ${ROMS_HEADER} )
endif()

# If use_4dvar returns a line containing "CPPDEFS" then somehow ROMS_HEADER
# was not defined and the build script configuation or CMake command should
# be checked.

if( "${defs}" MATCHES "CPPDEFS" )
  message( FATAL_ERROR "The C-preprocessor was unable to find your application's header file please check you configuation" )
endif()

if( "${defs}" MATCHES "ADJOINT" )
  option( ADJOINT "Turn on/off Adjoint Model" ON )
  Message( STATUS "Adjoint Model ENABLED" )
endif()

if( "${defs}" MATCHES "REPRESENTER" )
  option( REPRESENTER "Turn on/off Representer Model" ON )
  Message( STATUS "Representer Model ENABLED" )
endif()

if( "${defs}" MATCHES "TANGENT" )
  option( TANGENT "Turn on/off Tangent Linear Model" ON )
  message( STATUS "Tangent Linear Model ENABLED" )
endif()

if( "${defs}" MATCHES "ARPACK" )
  option( ARPACK "ARPACK/PARPACK Library" ON )
  message( STATUS "ROMS Link With ARPACK/PARPACK ENABLED" )
endif()


