#!/bin/csh -ef
#
# svn $Id: cbuild_roms.csh 1054 2021-03-06 19:47:12Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2021 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# ROMS ecbuild (CMake) Compiling CSH Script                             :::
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
#    ./cbuild_roms.csh [options]                                        :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    -j [N]      Compile in parallel using N CPUs                       :::
#                  omit argument for all available CPUs                 :::
#                                                                       :::
#    -p macro    Prints any Makefile macro value. For example,          :::
#                                                                       :::
#                  cbuild_roms.csh -p MY_CPP_FLAGS                      :::
#                                                                       :::
#    -noclean    Do not clean already compiled objects                  :::
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
      set debug = "$1"
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
      echo "              For example:  cbuild_roms.csh -p FFLAGS"
      echo ""
      echo "-noclean    Do not clean already compiled objects"
      echo ""
      exit 1
    breaksw

  endsw
end

# Set the CPP option defining the particular application. This will
# determine the name of the ".h" header file with the application
# CPP definitions. REQUIRED

 setenv ROMS_APPLICATION     UPWELLING

# Set a local environmental variable to define the path to the directories
# where all this project's files are kept.

 setenv MY_ROOT_DIR          ${HOME}/ocean/repository
 setenv MY_PROJECT_DIR       ${PWD}

# The path to the user's local current ROMS source code.

 setenv MY_ROMS_SRC          ${MY_ROOT_DIR}/trunk

# Which type(s) of libraries would you like?
# NOTE: If you choose both and also choose to build the ROMS executable,
#       it will be linked to the static version of the library.
#
# Valid options are SHARED, STATIC, and BOTH.

 setenv LIBTYPE              STATIC

# Do you want to build the ROMS executable?
# Valid values are: ON (build the executable) and OFF (do NOT build the
# executable). If you comment this out the executable WILL be built.

 setenv ROMS_EXECUTABLE      ON

# Set path of the directory containing my_build_paths.csh
# The user has the option to specify a customized version of this file
# in a different directory than the one distributed with the source code,
# ${MY_ROMS_SRC}/Compilers. If this is the case, you need to keep this
# configurations files up-to-date.

 setenv COMPILERS            ${MY_ROMS_SRC}/Compilers
#setenv COMPILERS            ${HOME}/Compilers/ROMS

#--------------------------------------------------------------------------
# Set tunable CPP options.
#--------------------------------------------------------------------------
#
# Sometimes it is desirable to activate one or more CPP options to run
# different variants of the same application without modifying its header
# file. If this is the case, specify each option here.
#
# Notice also that you need to use shell's quoting syntax to enclose the
# definition.
#
#    setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DAVERAGES"
#    setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DDEBUGGING"

#--------------------------------------------------------------------------
# Compilation options and paths.
#--------------------------------------------------------------------------

#setenv USE_ECBUILD          on              # use "ecbuild" wrapper

 setenv USE_MPI              on              # distributed-memory
 setenv USE_MPIF90           on              # compile with mpif90 script
#setenv which_MPI            mpich           # compile with MPICH library
#setenv which_MPI            mpich2          # compile with MPICH2 library
#setenv which_MPI            mvapich2        # compile with MVAPICH2 library
 setenv which_MPI            openmpi         # compile with OpenMPI library

 setenv FORT                 ifort
#setenv FORT                 gfortran
#setenv FORT                 pgi

#setenv USE_DEBUG            on              # use Fortran debugging flags
#setenv USE_NETCDF4          on              # use NetCDF4

#--------------------------------------------------------------------------
# If applicable, use my specified library paths.
#--------------------------------------------------------------------------

 setenv USE_MY_LIBS no           # use system default library paths
#setenv USE_MY_LIBS yes          # use my customized library paths

set MY_PATHS = ${COMPILERS}/my_build_paths.csh

if ($USE_MY_LIBS == 'yes') then
  source ${MY_PATHS} ${MY_PATHS}
endif

# Set location of the application header file.

 setenv MY_HEADER_DIR        ${MY_PROJECT_DIR}

# If you have custom analytical functions to include, enter the path here.

 setenv MY_ANALYTICAL_DIR    ${MY_PROJECT_DIR}

# Put the CMake files in a project specific Build directory to avoid conflict
# with other projects.

if ( $?USE_DEBUG ) then
  if ( "${USE_DEBUG}" = "on" ) then
    setenv SCRATCH_DIR       ${MY_PROJECT_DIR}/CBuild_romsG
  else
    setenv SCRATCH_DIR       ${MY_PROJECT_DIR}/CBuild_roms
  endif
else
  setenv SCRATCH_DIR         ${MY_PROJECT_DIR}/CBuild_roms
endif

# Create the build directory specified above and change into it.

if ( -d ${SCRATCH_DIR} ) then
  if ( $clean == 1 ) then
    rm -rf ${SCRATCH_DIR}
    mkdir ${SCRATCH_DIR}
    cd ${SCRATCH_DIR}
  else
    cd ${SCRATCH_DIR}
  endif
else
  if ( $clean == 1 ) then
    mkdir ${SCRATCH_DIR}
    cd ${SCRATCH_DIR}
  else
    echo "-noclean option activated when the build directory didn't exist"
    echo "creating the directory and disabling -noclean"
    clean=1
    mkdir ${SCRATCH_DIR}
    cd ${SCRATCH_DIR}
  endif
endif

#--------------------------------------------------------------------------
# Configure.
#--------------------------------------------------------------------------

# Construct the ecbuild command.

if ( $?LIBTYPE ) then
  set ltype="-DLIBTYPE=${LIBTYPE}"
else
  set ltype=""
endif

if ( $?MY_CPP_FLAGS ) then
  set tmp=`echo ${MY_CPP_FLAGS} | sed 's/^ *-D//' | sed 's/ *-D/;/g'`
  set extra_flags="-DMY_CPP_FLAGS=${tmp}"
else
  set extra_flags=""
endif

if ( $?PARPACK_LIBDIR ) then
  set parpack_ldir="-DPARPACK_LIBDIR=${PARPACK_LIBDIR}"
else
  set parpack_ldir=""
endif

if ( $?ARPACK_LIBDIR ) then
  set arpack_ldir="-DARPACK_LIBDIR=${ARPACK_LIBDIR}"
else
  set arpack_ldir=""
endif

# The nested ifs are required to avoid breaking the script, as tcsh
# apparently does not short-circuit if when the first truth is found

if ( $?USE_MPI ) then
  if ( "${USE_MPI}" == "on" ) then
    set mpi="-DMPI=ON"
  else
    set mpi=""
  endif
else
  set mpi=""
endif

if ( $?ROMS_EXECUTABLE ) then
  if ( "${ROMS_EXECUTABLE}" == "ON" ) then
    set roms_exec="-DROMS_EXECUTABLE=ON"
  else
    set roms_exec="-DROMS_EXECUTABLE=OFF"
  endif
else
  set roms_exec=""
endif

if ( $?USE_DEBUG ) then
  if ( "${USE_DEBUG}" == "on" ) then
    set dbg="-DCMAKE_BUILD_TYPE=Debug"
  else
    set dbg="-DCMAKE_BUILD_TYPE=Release"
  endif
else
  set dbg="-DCMAKE_BUILD_TYPE=Release"
endif

#--------------------------------------------------------------------------
# Run the chosen build command.
#--------------------------------------------------------------------------

set my_hdir="-DMY_HEADER_DIR=${MY_HEADER_DIR}"

if ( $dprint == 0 ) then
  if ( $?USE_ECBUILD ) then
    if ( "${USE_ECBUILD}" == "on" ) then
      set conf_com = "ecbuild"
    else if ( "${USE_ECBUILD}" == "off" ) then
      set conf_com = "cmake"
    else
      set conf_com = "Unknown"
    endif
  else
    set conf_com = "cmake"
  endif

  if ( "${conf_com}" == "cmake" ) then
    cmake -DAPP=${ROMS_APPLICATION} \
                ${my_hdir} \
                ${ltype} \
                ${extra_flags} \
                ${parpack_ldir} \
                ${arpack_ldir} \
                ${mpi} \
                ${roms_exec} \
                ${dbg} \
                ${MY_ROMS_SRC}
  else if ( "${conf_com}" == "ecbuild" ) then
    ecbuild -DAPP=${ROMS_APPLICATION} \
                  ${my_hdir} \
                  ${ltype} \
                  ${extra_flags} \
                  ${parpack_ldir} \
                  ${arpack_ldir} \
                  ${mpi} \
                  ${roms_exec} \
                  ${dbg} \
                  ${MY_ROMS_SRC}
  else
    echo "Unrecognized value, '${USE_ECBUILD}' set for USE_ECBUILD"
    exit 1
  endif
endif

#--------------------------------------------------------------------------
# Compile.
#--------------------------------------------------------------------------

if ( $dprint == 1 ) then
  set val = `eval echo \$${debug}`
  echo "${debug}:$val"
else
  if ( $parallel == 1 ) then
    make $NCPUS
  else
    make
  endif
  make install
endif

cd ${MY_PROJECT_DIR}

# Create symlink to executable. This should work even if ROMS was
# linked to the shared library (libROMS.{so|dylib}) because
# CMAKE_BUILD_WITH_INSTALL_RPATH is set to FALSE so that
# RPATH/RUNPATH are set correctly for both the build tree and
# installed locations of the ROMS executable.

if ( $dprint == 0 ) then
  if ( $?USE_DEBUG ) then
    if ( "${USE_DEBUG}" == "on" ) then
      ln -sfv ${SCRATCH_DIR}/romsG
    endif
  else if ( $?USE_MPI ) then
    if ( "${USE_MPI}" == "on" ) then
      ln -sfv ${SCRATCH_DIR}/romsM
    endif
  else
    ln -sfv ${SCRATCH_DIR}/romsS
  endif
endif
