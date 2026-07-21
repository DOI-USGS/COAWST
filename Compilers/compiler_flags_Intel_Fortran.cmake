# git $Id$
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2026 The ROMS Group                                :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# CMake Flags for the Intel Fortran Compiler.

###########################################################################
# Set the name of the Fortran compiler based on FORT and USE_MPI
###########################################################################

if( MPI )
  if( ${COMM} MATCHES "intel")
    execute_process( COMMAND which mpiifort
                     OUTPUT_VARIABLE CMAKE_Fortran_COMPILER
                     OUTPUT_STRIP_TRAILING_WHITESPACE )
  else()
    execute_process( COMMAND which mpif90
                     OUTPUT_VARIABLE CMAKE_Fortran_COMPILER
                     OUTPUT_STRIP_TRAILING_WHITESPACE )
  endif()
else()
  execute_process( COMMAND which ifort
                   OUTPUT_VARIABLE CMAKE_Fortran_COMPILER
                   OUTPUT_STRIP_TRAILING_WHITESPACE )
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
# RELEASE WITH DEBUG INFORMATION FLAGS
###########################################################################

set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-ip -O3 -g -traceback -check all -check bounds" )

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

# Flags for the C-preprocessor

set( CPPFLAGS  "-P" "--traditional-cpp" "-w" )

###########################################################################
