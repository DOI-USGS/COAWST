# git $Id$
# svn $Id: compiler_flags_Intel_Fortran.cmake 1151 2023-02-09 03:08:53Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# CMake Flags for the Intel Fortran Compiler.

###########################################################################
# Set the name of the Fortran compiler based on FORT and USE_MPI
###########################################################################

if( MPI )
  if( ${COMM} MATCHES "intel")
    set( CMAKE_Fortran_COMPILER mpiifort )
  else()
    set( CMAKE_Fortran_COMPILER mpif90 )
  endif()
else()
  set( CMAKE_Fortran_COMPILER ifort )
endif()

###########################################################################
# FLAGS COMMON TO ALL BUILD TYPES
###########################################################################

set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fp-model precise" )

###########################################################################
# RELEASE FLAGS
###########################################################################

set( CMAKE_Fortran_FLAGS_RELEASE "-ip -O3 -traceback -check uninit" )

###########################################################################
# DEBUG FLAGS
###########################################################################

set( CMAKE_Fortran_FLAGS_DEBUG   "-g -check all -check bounds -traceback -warn interfaces,nouncalled -gen-interfaces" )

###########################################################################
# BIT REPRODUCIBLE FLAGS
###########################################################################

set( CMAKE_Fortran_FLAGS_BIT     "-O2" )

###########################################################################
# LINK FLAGS
###########################################################################

if( APPLE )
  message( STATUS "MacOS so setting -Wl,-undefined dynamic_lookup linker flag" )
  set( CMAKE_SHARED_LINKER_FLAGS    "-Wl,-undefined,dynamic_lookup" )
endif()

set( CMAKE_Fortran_LINK_FLAGS    "" )

###########################################################################
# ROMS Definitions
###########################################################################

# Special flags for files that may contain very long lines.

set( ROMS_FREEFLAGS "-free" )

# Special flag for if compiler is confused by "dimension(*)"

set( ROMS_NOBOUNDSFLAG "" )

# Special defines for "mod_strings.F"

set( my_fort   "ifort" )
set( my_fc     "${CMAKE_Fortran_COMPILER}" )

if( ${CMAKE_BUILD_TYPE} MATCHES "Debug" )
  set( my_fflags "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}" )
else()
  set( my_fflags "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE}" )
endif()

# Flags for the C-preprocessor

set( CPPFLAGS  "-P" "--traditional-cpp" "-w" )

###########################################################################
