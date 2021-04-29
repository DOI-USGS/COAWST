#!/bin/csh -f
#
# svn $Id: build_roms.csh 1054 2021-03-06 19:47:12Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2021 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# ROMS Compiling CSH Script                                             :::
#                                                                       :::
# Script to compile an user application where the application-specific  :::
# files are kept separate from the ROMS source code.                    :::
#                                                                       :::
# Q: How/why does this script work?                                     :::
#                                                                       :::
# A: The ROMS makefile configures user-defined options with a set of    :::
#    flags such as ROMS_APPLICATION. Browse the makefile to see these.  :::
#    If an option in the makefile uses the syntax ?= in setting the     :::
#    default, this means that make will check whether an environment    :::
#    variable by that name is set in the shell that calls make. If so   :::
#    the environment variable value overrides the default (and the      :::
#    user need not maintain separate makefiles, or frequently edit      :::
#    the makefile, to run separate applications).                       :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    ./build_roms.csh [options]                                         :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    -j [N]      Compile in parallel using N CPUs                       :::
#                  omit argument for all available CPUs                 :::
#                                                                       :::
#    -p macro    Prints any Makefile macro value. For example,          :::
#                                                                       :::
#                  build_roms.csh -p FFLAGS                             :::
#                                                                       :::
#    -noclean    Do not clean already compiled objects                  :::
#                                                                       :::
# Notice that sometimes the parallel compilation fail to find MPI       :::
# include file "mpif.h".                                                :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

setenv which_MPI openmpi                      #  default, overwritten below

set parallel = 0
set clean = 1
set dprint = 0

setenv MY_CPP_FLAGS ''

while ( ($#argv) > 0 )
  switch ($1)
    case "-noclean"
      shift
      set clean = 0
    breaksw

    case "-p"
      shift
      set clean = 0
      set dprint = 1
      set debug = "print-$1"
      shift
    breaksw

    case "-j"
      shift
      set parallel = 1
      if (`echo $1 | grep '^[0-9]\+$'` != "" ) then
        set NCPUS = "-j $1"
        shift
      else
        set NCPUS = "-j"
      endif
    breaksw

    case "-*":
      echo ""
      echo "$0 : Unknown option [ $1 ]"
      echo ""
      echo "Available Options:"
      echo ""
      echo "-j [N]      Compile in parallel using N CPUs"
      echo "              omit argument for all avaliable CPUs"
      echo ""
      echo "-p macro    Prints any Makefile macro value"
      echo "              For example:  build_roms.csh -p FFLAGS"
      echo ""
      echo "-noclean    Do not clean already compiled objects"
      echo ""
      exit 1
    breaksw

  endsw
end

# Set the CPP option defining the particular application. This will
# determine the name of the ".h" header file with the application
# CPP definitions.

setenv ROMS_APPLICATION      UPWELLING

# Set a local environmental variable to define the path to the directories
# where all this project's files are kept.

setenv MY_ROOT_DIR           ${HOME}/ocean/repository
setenv MY_PROJECT_DIR        ${PWD}

# The path to the user's local current ROMS source code.
#
# If using svn locally, this would be the user's Working Copy Path (WCPATH).
# Note that one advantage of maintaining your source code locally with svn
# is that when working simultaneously on multiple machines (e.g. a local
# workstation, a local cluster and a remote supercomputer) you can checkout
# the latest release and always get an up-to-date customized source on each
# machine. This script is designed to more easily allow for differing paths
# to the code and inputs on differing machines.

 setenv MY_ROMS_SRC          ${MY_ROOT_DIR}/trunk

# Set path of the directory containing makefile configuration (*.mk) files.
# The user has the option to specify a customized version of these files
# in a different directory than the one distributed with the source code,
# ${MY_ROMS_SRC}/Compilers. If this is the case, you need to keep these
# configurations files up-to-date.

 setenv COMPILERS            ${MY_ROMS_SRC}/Compilers
#setenv COMPILERS            ${HOME}/Compilers/ROMS

#--------------------------------------------------------------------------
# Set tunable CPP options.
#--------------------------------------------------------------------------
#
# Sometimes it is desirable to activate one or more CPP options to run
# different variants of the same application without modifying its header
# file. If this is the case, specify each options here using the -D syntax.
# Notice also that you need to use shell's quoting syntax to enclose the
# definition.  Both single or double quotes work. For example,
#
#    setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DAVERAGES"
#    setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DDEBUGGING"
#
# can be used to write time-averaged fields. Notice that you can have as
# many definitions as you want by appending values.

#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -D"

#--------------------------------------------------------------------------
# Compiler options.
#--------------------------------------------------------------------------
#
# Other user defined environmental variables. See the ROMS makefile for
# details on other options the user might want to set here. Be sure to
# leave the switches meant to be off set to an empty string or commented
# out. Any string value (including off) will evaluate to TRUE in
# conditional if-statements.

 setenv USE_MPI             on          # distributed-memory parallelism
 setenv USE_MPIF90          on          # compile with mpif90 script
#setenv which_MPI           mpich       # compile with MPICH library
#setenv which_MPI           mpich2      # compile with MPICH2 library
#setenv which_MPI           mvapich2    # compile with MVAPICH2 library
 setenv which_MPI           openmpi     # compile with OpenMPI library

#setenv USE_OpenMP          on          # shared-memory parallelism

 setenv FORT                ifort
#setenv FORT                gfortran
#setenv FORT                pgi

#setenv USE_DEBUG           on          # use Fortran debugging flags
 setenv USE_LARGE           on          # activate 64-bit compilation
#setenv USE_NETCDF4         on          # compile with NetCDF-4 library
#setenv USE_PARALLEL_IO     on          # Parallel I/O with NetCDF-4/HDF5

#--------------------------------------------------------------------------
# If coupling Earth Systems Models (ESM), set the location of the ESM
# component libraries and modules. The strategy is to compile and link
# each ESM component separately first, and then ROMS since it is driving
# the coupled system. Only the ESM components activated are considered
# and the rest are ignored.  Some components like WRF cannot be built
# in a directory specified by the user but in it is the root directory,
# and cannot be moved when debugging with tools like TotalView.
#--------------------------------------------------------------------------

setenv WRF_SRC_DIR         ${HOME}/ocean/repository/WRF

if ($?USE_DEBUG) then
  setenv CICE_LIB_DIR      ${MY_PROJECT_DIR}/Build_ciceG
  setenv COAMPS_LIB_DIR    ${MY_PROJECT_DIR}/Build_coampsG
  setenv REGCM_LIB_DIR     ${MY_PROJECT_DIR}/Build_regcmG
  setenv WAM_LIB_DIR       ${MY_PROJECT_DIR}/Build_wamG
# setenv WRF_LIB_DIR       ${MY_PROJECT_DIR}/Build_wrfG
  setenv WRF_LIB_DIR       ${WRF_SRC_DIR}
else
  setenv CICE_LIB_DIR      ${MY_PROJECT_DIR}/Build_cice
  setenv COAMPS_LIB_DIR    ${MY_PROJECT_DIR}/Build_coamps
  setenv REGCM_LIB_DIR     ${MY_PROJECT_DIR}/Build_regcm
  setenv WAM_LIB_DIR       ${MY_PROJECT_DIR}/Build_wam
  setenv WRF_LIB_DIR       ${MY_PROJECT_DIR}/Build_wrf
# setenv WRF_LIB_DIR       ${WRF_SRC_DIR}
endif

#--------------------------------------------------------------------------
# If applicable, use my specified library paths.
#--------------------------------------------------------------------------

 setenv USE_MY_LIBS no           # use system default library paths
#setenv USE_MY_LIBS yes          # use my customized library paths

set MY_PATHS = ${COMPILERS}/my_build_paths.csh

if ($USE_MY_LIBS == 'yes') then
  source ${MY_PATHS} ${MY_PATHS}
endif

#--------------------------------------------------------------------------
# The rest of this script sets the path to the users header file and
# analytical source files, if any. See the templates in User/Functionals.
#--------------------------------------------------------------------------
#
# If applicable, use the MY_ANALYTICAL_DIR directory to place your
# customized biology model header file (like fennel.h, nemuro.h, ecosim.h,
# etc).

 setenv MY_HEADER_DIR       ${MY_PROJECT_DIR}

 setenv MY_ANALYTICAL_DIR   ${MY_PROJECT_DIR}

# Put the binary to execute in the following directory.

 setenv BINDIR              ${MY_PROJECT_DIR}

# Put the f90 files in a project specific Build directory to avoid conflict
# with other projects.

if ($?USE_DEBUG) then
  setenv SCRATCH_DIR        ${MY_PROJECT_DIR}/Build_romsG
else
  setenv SCRATCH_DIR        ${MY_PROJECT_DIR}/Build_roms
endif

# Go to the users source directory to compile. The options set above will
# pick up the application-specific code from the appropriate place.

 cd ${MY_ROMS_SRC}

# Stop if activating both MPI and OpenMP at the same time.

if ( ${?USE_MPI} & ${?USE_OpenMP} ) then
  echo "You cannot activate USE_MPI and USE_OpenMP at the same time!"
  exit 1
endif

#--------------------------------------------------------------------------
# Compile.
#--------------------------------------------------------------------------

# Remove build directory.

if ( $clean == 1 ) then
  make clean
endif

# Compile (the binary will go to BINDIR set above).

if ( $dprint == 1 ) then
  make $debug
else
  if ( $parallel == 1 ) then
    make $NCPUS
  else
    make
  endif
endif
