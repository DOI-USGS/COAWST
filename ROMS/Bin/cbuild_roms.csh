#!/bin/csh -ef
#
# git $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2026 The ROMS Group                                :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# ROMS CMake Compiling CSH Script                                       :::
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
#    -b          Compile a specific ROMS GitHub branch                  :::
#                                                                       :::
#                  cbuild_roms.csh -j 5 -b feature/kernel               :::
#                                                                       :::
#    -g          Compile with debug flag (slower code)                  :::
#                                                                       :::
#                  cbuild_roms.csh -g -j 10                             :::
#                                                                       :::
#    -j [N]      Compile in parallel using N CPUs                       :::
#                  omit argument for all available CPUs                 :::
#                                                                       :::
#    -pio        Compile with PIO (Parallel I/O) NetCDF library         :::
#                  Otherwise, it used standard NetCDF library (slower)  :::
#                                                                       :::
#                  cbuild_roms.csh -pio -j 10                           :::
#                                                                       :::
#    -p macro    Prints any Makefile macro value. For example,          :::
#                                                                       :::
#                  cbuild_roms.csh -p MY_CPP_FLAGS                      :::
#                                                                       :::
#    -noclean    Do not clean already compiled objects                  :::
#                                                                       :::
#    -v             Compile in verbose mode (VERBOSE=1)                 :::
#                                                                       :::
# The branch option -b is only possible for ROMS source code from       :::
# https://github.com/myroms. Such versions are under development        :::
# and targeted to advanced users, superusers, and beta testers.         :::
# Regular and novice users must use the default 'develop' branch.       :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

setenv which_MPI openmpi                      #  default, overwritten below

set g_flags = 0
set parallel = 0
set pio_lib = 0
set clean = 1
set dprint = 0
set Verbose = 0
set branch = 0

set command = "cbuild_roms.csh $argv[*]"

set separator = `perl -e "print '<>' x 50;"`

setenv MY_CPP_FLAGS ''

while ( ($#argv) > 0 )
  switch ($1)
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

    case "-g"
      shift
      set g_flags = 1
    breaksw

    case "-pio"
      shift
      set pio_lib = 1
    breaksw

    case "-p"
      shift
      set clean = 0
      set dprint = 1
      set debug = "$1"
      shift
    breaksw

    case "-v"
      shift
      set Verbose = 1
    breaksw

    case "-noclean"
      shift
      set clean = 0
    breaksw

    case "-b"
      shift
      set branch = 1
      set branch_name = `echo $1 | grep -v '^-'`
      if ( "$branch_name" == "" ) then
        echo "Please enter a branch name."
        exit 1
      endif
      shift
    breaksw

    case "-*":
      echo ""
      echo "${separator}"
      echo "$0 : Unknown option [ $1 ]"
      echo ""
      echo "Available Options:"
      echo ""
      echo "-b branch_name  Compile specific ROMS GitHub branch name"
      echo "                  For example:  cbuild_roms.csh -b feature/kernel"
      echo ""
      echo "-g              Compile with debugging flags, slower code"
      echo ""
      echo "-j [N]          Compile in parallel using N CPUs"
      echo "                  omit argument for all avaliable CPUs"
      echo ""
      echo "-pio            Compile with the PIO NetCDF Library"
      echo ""
      echo "-p macro        Prints any Makefile macro value"
      echo "                  For example:  cbuild_roms.csh -p FFLAGS"
      echo ""
      echo "-noclean        Do not clean already compiled objects"
      echo ""
      echo "-v              Compile in verbose mode"
      echo "${separator}"
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
# where the ROMS source code is located (MY_ROOT_DIR), and this project's
# configuration and files are kept (MY_PROJECT_DIR). Notice that if the
# User sets the ROMS_ROOT_DIR environment variable in their computer logging
# script describing the location from where the ROMS source code was cloned
# or downloaded, it uses that value.

if ( $?ROMS_ROOT_DIR ) then
  if ( "${ROMS_ROOT_DIR}" != "" ) then
    setenv MY_ROOT_DIR       ${ROMS_ROOT_DIR}
  else
    setenv MY_ROOT_DIR       ${HOME}/ocean/repository/git
  endif
else
  setenv MY_ROOT_DIR         ${HOME}/ocean/repository/git
endif

setenv MY_PROJECT_DIR        ${PWD}

# The path to the user's local current ROMS source code.
#
# If downloading ROMS locally, this would be the user's Working Copy Path.
# One advantage of maintaining your source code copy is that when working
# simultaneously on multiple machines (e.g., a local workstation, a local
# cluster, and a remote supercomputer), you can update with the latest ROMS
# release and always get an up-to-date customized source on each machine.
# This script allows for differing paths to the code and inputs on other
# computers.

 setenv MY_ROMS_SRC          ${MY_ROOT_DIR}/roms

# Which type(s) of libraries would you like?
#
# NOTE: If you choose both and also choose to build the ROMS executable,
#       it will be linked to the static version of the library.
#
# Valid options are SHARED, STATIC, and BOTH.

 setenv LIBTYPE              STATIC

# Do you want to build the ROMS executable?
#
# Valid values are: ON (build the executable) and OFF (do NOT build the
# executable). If you comment this out the executable WILL be built.

 setenv ROMS_EXECUTABLE      ON

# Set path of the directory containing "my_build_paths.csh".
#
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
#
# can be used to write time-averaged fields. Notice that you can have as
# many definitions as you want by appending values.

if ( $pio_lib == 1 ) then
  setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DPIO_LIB"
endif

#--------------------------------------------------------------------------
# Compilation options.
#--------------------------------------------------------------------------

 setenv USE_MPI              on              # distributed-memory
 setenv USE_MPIF90           on              # compile with mpif90 script
#setenv which_MPI            intel           # compile with mpiifort library
#setenv which_MPI            mpich           # compile with MPICH library
#setenv which_MPI            mpich2          # compile with MPICH2 library
#setenv which_MPI            mvapich2        # compile with MVAPICH2 library
 setenv which_MPI            openmpi         # compile with OpenMPI library

#setenv FORT                 ifx
 setenv FORT                 ifort
#setenv FORT                 gfortran
#setenv FORT                 pgi

if ( $g_flags == 1 ) then
 setenv USE_DEBUG            on              # use Fortran debugging flags
endif

# ROMS I/O choices and combinations. A more complete description of the
# available options can be found in the wiki (https://myroms.org/wiki/IO).
# Most users will want to enable at least USE_NETCDF4 because that will
# instruct the ROMS build system to use nf-config to determine the
# necessary libraries and paths to link into the ROMS executable.

 setenv USE_NETCDF4          on              # compile with NetCDF4 library
#setenv USE_PARALLEL_IO      on              # Parallel I/O with NetCDF-4/HDF5

if ( $pio_lib == 1 ) then
  setenv USE_PIO             on              # Parallel I/O with PIO library
endif

# If any of the coupling component use the HDF5 Fortran API for primary
# I/O, we need to compile the main driver with the HDF5 library.

#setenv USE_HDF5             on              # compile with HDF5 library

#--------------------------------------------------------------------------
# If coupling Earth Systems Models (ESM), set the location of the ESM
# component libraries and modules.
#--------------------------------------------------------------------------

source ${MY_ROMS_SRC}/ESM/esm_libs.csh ${MY_ROMS_SRC}/ESM/esm_libs.csh

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

 echo ""
 echo "${separator}"

# Put the CMake files in a project specific Build directory to avoid conflict
# with other projects.

if ( $?USE_DEBUG ) then
  if ( "${USE_DEBUG}" == "on" ) then
    setenv BUILD_DIR         ${MY_PROJECT_DIR}/CBuild_romsG
  else if ( $?USE_MPI ) then
    if ( "${USE_MPI}" == "on" ) then
      setenv BUILD_DIR       ${MY_PROJECT_DIR}/CBuild_romsM
    else
      setenv BUILD_DIR       ${MY_PROJECT_DIR}/CBuild_roms
    endif
  else
    setenv BUILD_DIR         ${MY_PROJECT_DIR}/CBuild_roms
  endif
else if ($?USE_MPI) then
  if ( "${USE_MPI}" == "on" ) then
    setenv BUILD_DIR         ${MY_PROJECT_DIR}/CBuild_romsM
  else
    setenv BUILD_DIR         ${MY_PROJECT_DIR}/CBuild_roms
  endif
else
  setenv BUILD_DIR           ${MY_PROJECT_DIR}/CBuild_roms
endif

# For backward compatibility, set deprecated SCRATCH_DIR to compile
# older released versions of ROMS.

setenv SCRATCH_DIR ${BUILD_DIR}

# If necessary, create ROMS build directory.

if ( $dprint == 0 ) then
  if ( -d ${BUILD_DIR} ) then
    if ( $clean == 1 ) then
      echo ""
      echo "Removing ROMS build directory: ${BUILD_DIR}"
      echo ""
      rm -rf ${BUILD_DIR}
      echo ""
      echo "Creating ROMS build directory: ${BUILD_DIR}"
      echo ""
      mkdir ${BUILD_DIR}
    endif
  else
    if ( $clean == 1 ) then
      mkdir ${BUILD_DIR}
      cd ${BUILD_DIR}
    else
      echo ""
      echo "Option -noclean activated when the ROMS build directory didn't exist"
      echo "Creating ROMS build directory and disabling -noclean"
      echo ""
      set clean = 1
      mkdir ${BUILD_DIR}
      cd ${BUILD_DIR}
    endif
  endif

  # If requested, check out requested branch from ROMS GitHub

  if ( $branch == 1 ) then
    if ( ! -d ${MY_PROJECT_DIR}/src ) then
      echo ""
      echo "Downloading ROMS source code from GitHub: https://github.com/myroms"
      echo ""
      git clone https://github.com/myroms/roms.git src
    endif
    echo ""
    echo "Checking out ROMS GitHub branch: $branch_name"
    echo ""
    cd src
    git checkout $branch_name

    # If we are using the COMPILERS from the ROMS source code
    # overide the value set above

    if ( ${COMPILERS} =~ ${MY_ROMS_SRC}* ) then
      setenv COMPILERS ${MY_PROJECT_DIR}/src/Compilers
    endif
    setenv MY_ROMS_SRC ${MY_PROJECT_DIR}/src

  else
    echo ""
    echo "Using ROMS source code from: ${MY_ROMS_SRC}"
    echo ""
    cd ${MY_ROMS_SRC}
  endif
endif

#--------------------------------------------------------------------------
# Add environmental variables constructed in 'makefile' to MY_CPP_FLAGS
# so can be passed to ROMS.
#--------------------------------------------------------------------------

set ANALYTICAL_DIR = "ANALYTICAL_DIR='${MY_ANALYTICAL_DIR}'"
set HEADER = `echo ${ROMS_APPLICATION} | tr '[:upper:]' '[:lower:]'`.h
set HEADER_DIR = "HEADER_DIR='${MY_HEADER_DIR}'"
set ROOT_DIR = "ROOT_DIR='${MY_ROMS_SRC}'"

set mycppflags = "${MY_CPP_FLAGS}"

setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -D${ANALYTICAL_DIR}"
setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -D${HEADER_DIR}"
setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -D${ROOT_DIR}"

cd ${BUILD_DIR}

#--------------------------------------------------------------------------
# Configure.
#--------------------------------------------------------------------------

# Construct the cmake command.

if ( $?LIBTYPE ) then
  set ltype="-DLIBTYPE=${LIBTYPE}"
else
  set ltype=""
endif

if ( $?FORT ) then
  if ( "${FORT}" == "ifort" ) then
    set compiler="-DCMAKE_Fortran_COMPILER=ifort"
  elif ( "${FORT}" == "ifx" ) then
    set compiler="-DCMAKE_Fortran_COMPILER=ifx"
  else if ( "${FORT}" == "gfortran" ) then
    set compiler="-DCMAKE_Fortran_COMPILER=gfortran"
  else
    set compiler=""
  endif
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

if ( $?PIO_LIBDIR && $?PIO_INCDIR ) then
  set pio_ldir="-DPIO_LIBDIR=${PIO_LIBDIR}"
  set pio_idir="-DPIO_INCDIR=${PIO_INCDIR}"
  if ( $?PNETCDF_LIBDIR && $?PNETCDF_INCDIR ) then
    set pnetcdf_ldir="-DPNETCDF_LIBDIR=${PNETCDF_LIBDIR}"
    set pnetcdf_idir="-DPNETCDF_INCDIR=${PNETCDF_INCDIR}"
  else
    set pnetcdf_ldir=""
    set pnetcdf_idir=""
  endif
else
  set pio_ldir=""
  set pio_idir=""
  set pnetcdf_ldir=""
  set pnetcdf_idir=""
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

if ( $?USE_MPIF90 ) then
  if ( "${USE_MPIF90}" == "on" ) then
    set comm="-DCOMM=${which_MPI}"
  else
    set comm=""
  endif
else
  set comm=""
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
  if ( $clean == 1 ) then
    echo ""
    echo "Configuring CMake for ROMS application:"
    echo ""
    cmake -DROMS_APP=${ROMS_APPLICATION} \
                     ${my_hdir} \
                     ${ltype} \
                     ${compiler} \
                     ${extra_flags} \
                     ${parpack_ldir} \
                     ${arpack_ldir} \
                     ${pio_ldir} \
                     ${pio_idir} \
                     ${pnetcdf_ldir} \
                     ${pnetcdf_idir} \
                     ${mpi} \
                     ${comm} \
                     ${roms_exec} \
                     ${dbg} \
                     ${MY_ROMS_SRC}
  endif
endif

#--------------------------------------------------------------------------
# Compile.
#--------------------------------------------------------------------------

if ( $dprint == 1 ) then
  set val = `eval echo \$${debug}`
  echo "${debug}:$val"
else

  echo ""
  echo "Compiling ROMS source code:"
  echo ""

  if ( $parallel == 1 ) then
    if ( $Verbose == 1 ) then
      make VERBOSE=1 $NCPUS
    else
      make $NCPUS
    endif
  else
    if ( $Verbose == 1 ) then
      make VERBOSE=1
    else
      make
    endif
  endif
  make install

  set HEADER = `echo ${ROMS_APPLICATION} | tr '[:upper:]' '[:lower:]'`.h

  echo ""
  echo "${separator}"
  echo "CMake Build script command:    ${command}"
  echo "ROMS source directory:         ${MY_ROMS_SRC}"
  echo "ROMS header file:              ${MY_HEADER_DIR}/${HEADER}"
  echo "ROMS build  directory:         ${BUILD_DIR}"
  if ( $branch == 1 ) then
    echo "ROMS downloaded from:          https://github.com/myroms/roms.git"
    echo "ROMS compiled branch:          $branch_name"
  endif
  echo "ROMS Application:              ${ROMS_APPLICATION}"
  set FFLAGS = `cat fortran_flags`
  echo "Fortran compiler:              ${FORT}"
  echo "Fortran flags:                 ${FFLAGS}"
  if ($?mycppflags) then
    echo "Added CPP Options:            ${mycppflags}"
  endif
  echo "${separator}"
  echo ""
endif

cd ${MY_PROJECT_DIR}

# If ROMS_EXECUTABLE is set to OFF remove the symlink from
# previous build if present.

if ( $?ROMS_EXECUTABLE ) then
  if ( "${ROMS_EXECUTABLE}" == "OFF" ) then
    if ( $?USE_DEBUG ) then
      if ( "${USE_DEBUG}" == "on" ) then
        if { test -L romsG } then
          rm -f romsG
        endif
      endif
    else if ( $?USE_MPI ) then
      if ( "${USE_MPI}" == "on" ) then
        if { test  -L romsM } then
          rm -f romsM
        endif
      endif
    else
      if { test -L romsS } then
        rm -f romsS
      endif
    endif
  endif
endif

# Copy executable to project directory. This should work even
# if ROMS was linked to the shared library (libROMS.{so|dylib})
# because CMAKE_BUILD_WITH_INSTALL_RPATH is set to FALSE so that
# RPATH/RUNPATH are set correctly for both the build tree and
# installed locations of the ROMS executable.

if ( $dprint == 0 ) then
  if ( ! $?ROMS_EXECUTABLE ) then
    if ( $?USE_DEBUG ) then
      if ( "${USE_DEBUG}" == "on" ) then
        cp -pfv ${BUILD_DIR}/romsG .
      endif
    else if ( $?USE_MPI ) then
      if ( "${USE_MPI}" == "on" ) then
        cp -pfv ${BUILD_DIR}/romsM .
      endif
    else
      cp -pfv ${BUILD_DIR}/romsS .
    endif
  else
    if ( "${ROMS_EXECUTABLE}" == "ON" ) then
      if ( $?USE_DEBUG ) then
        if ( "${USE_DEBUG}" == "on" ) then
          cp -pfv ${BUILD_DIR}/romsG .
        endif
      else if ( $?USE_MPI ) then
        if ( "${USE_MPI}" == "on" ) then
          cp -pfv ${BUILD_DIR}/romsM .
        endif
      else
        cp -pfv ${BUILD_DIR}/romsS .
      endif
    endif
  endif
endif
