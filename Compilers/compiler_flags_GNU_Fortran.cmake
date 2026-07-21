# git $Id$
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2026 The ROMS Group                                :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# CMake Flags for the GNU Fortran Compiler.

###########################################################################
# Set the name of the Fortran compiler based on FORT and USE_MPI
###########################################################################

if( MPI )
  execute_process( COMMAND which mpif90
                   OUTPUT_VARIABLE CMAKE_Fortran_COMPILER
                   OUTPUT_STRIP_TRAILING_WHITESPACE )
else()
  execute_process( COMMAND which gfortran
                   OUTPUT_VARIABLE CMAKE_Fortran_COMPILER
                   OUTPUT_STRIP_TRAILING_WHITESPACE )
endif()

###########################################################################
# FLAGS COMMON TO ALL BUILD TYPES
###########################################################################

set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -frepack-arrays -fallow-argument-mismatch" )

###########################################################################
# RELEASE FLAGS
###########################################################################

set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -ffast-math" )

###########################################################################
# RELEASE WITH DEBUG INFORMATION FLAGS
###########################################################################

set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O3 -g -ffast-math -fbounds-check -fbacktrace -fcheck=all" )

###########################################################################
# DEBUG FLAGS
###########################################################################

set( CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fbounds-check -fbacktrace -fcheck=all -finit-real=nan -ffpe-trap=invalid,zero,overflow" )

###########################################################################
# BIT REPRODUCIBLE FLAGS
###########################################################################

set( CMAKE_Fortran_FLAGS_BIT     "-O2 -fbounds-check" )

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

set( ROMS_FREEFLAGS "-ffree-form -ffree-line-length-none" )

# Special flag for def_var since gfortran is confused by "dimension(*)"

set( ROMS_NOBOUNDSFLAG "-fno-bounds-check" )

# Special defines for "mod_strings.F"

set( my_fort   "gfortran" )
set( my_fc     "${CMAKE_Fortran_COMPILER}" )

# Flags for the C-preprocessor

set( CPPFLAGS  "-P" "--traditional-cpp" "-w" )

###########################################################################
