# svn $Id: roms_compiler_flags.cmake 1054 2021-03-06 19:47:12Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2021 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# CMake Fortran Compiler Flags.

set( my_os     "${CMAKE_SYSTEM_NAME}" )
set( my_cpu    "${CMAKE_SYSTEM_PROCESSOR}" )

###########################################################################
# Fortran
###########################################################################

if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU" )
  include( compiler_flags_GNU_Fortran )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
  include( compiler_flags_Intel_Fortran )
else()
  message( STATUS "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()

###########################################################################
# C-Preprocessor Definitions
###########################################################################

string( TOUPPER ${my_os} OS )
string( TOUPPER ${my_cpu} CPU )
string( TOUPPER ${my_fort} FORT )
add_definitions ( -D${OS} -D${CPU} -D${FORT} )

