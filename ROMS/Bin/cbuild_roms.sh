#!/bin/bash
#
# git $Id$
# svn $Id: cbuild_roms.sh 1151 2023-02-09 03:08:53Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
#                                                                       :::
# ROMS CMake Compiling BASH Script                                      :::
#                                                                       :::
# Script to configure and compile a user application where the          :::
# application-specific files are kept separate from the ROMS            :::
# source code.                                                          :::
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
#    ./cbuild_roms.sh [options]                                         :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    -j [N]      Compile in parallel using N CPUs                       :::
#                  omit argument for all available CPUs                 :::
#                                                                       :::
#    -p macro    Prints any Makefile macro value. For example,          :::
#                                                                       :::
#                  cbuild_roms.sh -p MY_CPP_FLAGS                       :::
#                                                                       :::
#    -noclean    Do not clean already compiled objects                  :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

export which_MPI=openmpi                       # default, overwritten below

parallel=0
clean=1
dprint=0

export MY_CPP_FLAGS=

while [ $# -gt 0 ]
do
  case "$1" in
    -j )
      shift
      parallel=1
      test=`echo $1 | grep '^[0-9]\+$'`
      if [ "$test" != "" ]; then
        NCPUS="-j $1"
        shift
      else
        NCPUS="-j"
      fi
      ;;

    -p )
      shift
      clean=0
      dprint=1
      debug="$1"
      shift
      ;;

    -noclean )
      shift
      clean=0
      ;;

    * )
      echo ""
      echo "$0 : Unknown option [ $1 ]"
      echo ""
      echo "Available Options:"
      echo ""
      echo "-j [N]      Compile in parallel using N CPUs"
      echo "              omit argument for all avaliable CPUs"
      echo ""
      echo "-p macro    Prints any Makefile macro value"
      echo "              For example:  cbuild_roms.sh -p FFLAGS"
      echo ""
      echo "-noclean    Do not clean already compiled objects"
      echo ""
      exit 1
      ;;
  esac
done

# Set the CPP option defining the particular application. This will
# determine the name of the ".h" header file with the application
# CPP definitions. REQUIRED

export   ROMS_APPLICATION=UPWELLING

# Set a local environmental variable to define the path to the directories
# where all this project's files are kept.

export        MY_ROOT_DIR=${HOME}/ocean/repository
export     MY_PROJECT_DIR=${PWD}

# The path to the user's local current ROMS source code.

 export       MY_ROMS_SRC=${MY_ROOT_DIR}/trunk

# Which type(s) of libraries would you like?
# NOTE: If you choose both and also choose to build the ROMS executable,
#       it will be linked to the shared version of the library.
#
# Valid options are SHARED, STATIC, and BOTH.

 export           LIBTYPE=STATIC

# Do you want to build the ROMS executable?
# Valid values are: ON (build the executable) and OFF (do NOT build the
# executable). If you comment this out the executable WILL be built.

 export   ROMS_EXECUTABLE=ON

# Set path of the directory containing my_build_paths.sh.
# The user has the option to specify a customized version of this file
# in a different directory than the one distributed with the source code,
# ${MY_ROMS_SRC}/Compilers. If this is the case, you need to keep these
# configurations files up-to-date.

 export         COMPILERS=${MY_ROMS_SRC}/Compilers
#export         COMPILERS=${HOME}/Compilers/ROMS

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
#export      MY_CPP_FLAGS="${MY_CPP_FLAGS} -DAVERAGES"
#export      MY_CPP_FLAGS="${MY_CPP_FLAGS} -DDEBUGGING"

#--------------------------------------------------------------------------
# Compilation options.
#--------------------------------------------------------------------------

# Set this option to "on" if you wish to use the "ecbuild" CMake wrapper.
# Setting this to "off" or commenting it out will use cmake directly.

 export       USE_ECBUILD=off              # don't use "ecbuild" wrapper
#export       USE_ECBUILD=on               # use "ecbuild" wrapper

 export           USE_MPI=on               # distributed-memory parallelism
 export        USE_MPIF90=on               # compile with mpif90 script
#export         which_MPI=intel            # compile with mpiifort library
#export         which_MPI=mpich            # compile with MPICH library
#export         which_MPI=mpich2           # compile with MPICH2 library
#export         which_MPI=mvapich2         # compile with MVAPICH2 library
 export         which_MPI=openmpi          # compile with OpenMPI library

 export              FORT=ifort
#export              FORT=gfortran
#export              FORT=pgi

#export         USE_DEBUG=on               # use Fortran debugging flags

# ROMS I/O choices and combinations. A more complete description of the
# available options can be found in the wiki (https://myroms.org/wiki/IO).
# Most users will want to enable at least USE_NETCDF4 because that will
# instruct the ROMS build system to use nf-config to determine the
# necessary libraries and paths to link into the ROMS executable.

 export       USE_NETCDF4=on               # compile with NetCDF-4 library
#export   USE_PARALLEL_IO=on               # Parallel I/O with NetCDF-4/HDF5
#export           USE_PIO=on               # Parallel I/O with PIO library
#export       USE_SCORPIO=on               # Parallel I/O with SCORPIO library

# If any of the coupling component use the HDF5 Fortran API for primary
# I/O, we need to compile the main driver with the HDF5 library.

#export          USE_HDF5=on               # compile with HDF5 library

#--------------------------------------------------------------------------
# If applicable, use my specified library paths.
#--------------------------------------------------------------------------

 export USE_MY_LIBS=no            # use system default library paths
#export USE_MY_LIBS=yes           # use my customized library paths

MY_PATHS=${COMPILERS}/my_build_paths.sh

if [ "${USE_MY_LIBS}" = "yes" ]; then
  source ${MY_PATHS} ${MY_PATHS}
fi

# Set location of the application header file.

 export     MY_HEADER_DIR=${MY_PROJECT_DIR}

# If you have custom analytical functions to include, enter the path here.

 export MY_ANALYTICAL_DIR=${MY_PROJECT_DIR}

# Put the CMake files in a project specific Build directory to avoid conflict
# with other projects.

if [ -n "${USE_DEBUG:+1}" ]; then
  export      SCRATCH_DIR=${MY_PROJECT_DIR}/CBuild_romsG
else
  export      SCRATCH_DIR=${MY_PROJECT_DIR}/CBuild_roms
fi

# Create the build directory specified above and change into it.

if [ -d ${SCRATCH_DIR} ]; then
  if [ $clean -eq 1 ]; then
    rm -rf ${SCRATCH_DIR}
    mkdir ${SCRATCH_DIR}
    cd ${SCRATCH_DIR}
  else
    cd ${SCRATCH_DIR}
  fi
else
  if [ $clean -eq 1 ]; then
    mkdir ${SCRATCH_DIR}
    cd ${SCRATCH_DIR}
  else
    echo "-noclean option activated when the build directory didn't exist"
    echo "creating the directory and disabling -noclean"
    clean=1
    mkdir ${SCRATCH_DIR}
    cd ${SCRATCH_DIR}
  fi
fi

#--------------------------------------------------------------------------
# Add enviromental variables constructed in 'makefile' to MY_CPP_FLAGS
# so can be passed to ROMS.
#--------------------------------------------------------------------------

ANALYTICAL_DIR="ANALYTICAL_DIR='${MY_ANALYTICAL_DIR}'"
HEADER=`echo ${ROMS_APPLICATION} | tr '[:upper:]' '[:lower:]'`.h
HEADER_DIR="HEADER_DIR='${MY_HEADER_DIR}'"
ROOT_DIR="ROOT_DIR='${MY_ROMS_SRC}'"

export       MY_CPP_FLAGS="${MY_CPP_FLAGS} -D${ANALYTICAL_DIR}"
export       MY_CPP_FLAGS="${MY_CPP_FLAGS} -D${HEADER_DIR}"
export       MY_CPP_FLAGS="${MY_CPP_FLAGS} -D${ROOT_DIR}"

if [[ -d "${MY_ROMS_SRC}/.git" ]]; then
  cd ${MY_ROMS_SRC}
  GITURL=$(git config --get remote.origin.url)
  GITREV=$(git rev-parse --verify HEAD)
  GIT_URL="GIT_URL='${GITURL}'"
  GIT_REV="GIT_REV='${GITREV}'"
  SVN_URL="SVN_URL='https://www.myroms.org/svn/src'"

  export     MY_CPP_FLAGS="${MY_CPP_FLAGS} -D${GIT_URL}"
  export     MY_CPP_FLAGS="${MY_CPP_FLAGS} -D${GIT_REV}"
  export     MY_CPP_FLAGS="${MY_CPP_FLAGS} -D${SVN_URL}"
  cd ${SCRATCH_DIR}
else
  cd ${MY_ROMS_SRC}
  SVNURL=$(svn info | grep '^URL:' | sed 's/URL: //')
  SVNREV=$(svn info | grep '^Revision:' | sed 's/Revision: //')
  SVN_URL="SVN_URL='${SVNURL}'"
  SVN_REV="SVN_REV='${SVNREV}'"

  export     MY_CPP_FLAGS="${MY_CPP_FLAGS} -D${SVN_URL}"
  export     MY_CPP_FLAGS="${MY_CPP_FLAGS} -D${SVN_REV}"
  cd ${SCRATCH_DIR}
fi

#--------------------------------------------------------------------------
# Configure.
#--------------------------------------------------------------------------

# Construct the "ecbuild" command.

if [ ! -z "${LIBTYPE}" ]; then
  ltype="-DLIBTYPE=${LIBTYPE}"
else
  ltype=""
fi

if [ ! -z "${MY_CPP_FLAGS}" ]; then
  tmp=`echo ${MY_CPP_FLAGS} | sed 's/^ *-D//' | sed 's/ *-D/;/g'`
  extra_flags="-DMY_CPP_FLAGS=${tmp}"
else
  extra_flags=""
fi

if [ ! -z "${PARPACK_LIBDIR}" ]; then
  parpack_ldir="-DPARPACK_LIBDIR=${PARPACK_LIBDIR}"
else
  parpack_ldir=""
fi

if [ ! -z "${ARPACK_LIBDIR}" ]; then
  arpack_ldir="-DARPACK_LIBDIR=${ARPACK_LIBDIR}"
else
  arpack_ldir=""
fi

if [ ! -z "${USE_SCORPIO}" ]; then
  if [[ ! -z "${PIO_LIBDIR}" && ! -z "${PIO_INCDIR}" ]]; then
    pio_ldir="-DPIO_LIBDIR=${PIO_LIBDIR}"
    pio_idir="-DPIO_INCDIR=${PIO_INCDIR}"
    if [[ ! -z "${PNETCDF_LIBDIR}" && ! -z "${PNETCDF_INCDIR}" ]]; then
      pnetcdf_ldir="-DPNETCDF_LIBDIR=${PNETCDF_LIBDIR}"
      pnetcdf_idir="-DPNETCDF_INCDIR=${PNETCDF_INCDIR}"
    else
      pnetcdf_ldir=""
      pnetcdf_idir=""
    fi
  else
    pio_ldir=""
    pio_idir=""
    pnetcdf_ldir=""
    pnetcdf_idir=""
  fi
fi

if [[ ! -z "${USE_MPI}" && "${USE_MPI}" == "on" ]]; then
  mpi="-DMPI=ON"
else
  mpi=""
fi

if [[ ! -z "${USE_MPIF90}" && "${USE_MPIF90}" == "on" ]]; then
  comm="-DCOMM=${which_MPI}"
else
  comm=""
fi

if [ ! -z "${ROMS_EXECUTABLE}" ]; then
  if [[ "${ROMS_EXECUTABLE}" == "ON" ]]; then
    roms_exec="-DROMS_EXECUTABLE=ON"
  else
    roms_exec="-DROMS_EXECUTABLE=OFF"
  fi
else
  roms_exec=""
fi

if [[ ! -z "${USE_DEBUG}" && "${USE_DEBUG}" == "on" ]]; then
  dbg="-DCMAKE_BUILD_TYPE=Debug"
else
  dbg="-DCMAKE_BUILD_TYPE=Release"
fi

#--------------------------------------------------------------------------
# Run the chosen build command.
#--------------------------------------------------------------------------

my_hdir="-DMY_HEADER_DIR=${MY_HEADER_DIR}"

if [[ $dprint -eq 0 && $clean -eq 1 ]]; then
  if [[ -z ${USE_ECBUILD+x} || "${USE_ECBUILD}" == "off" ]]; then
    conf_com="cmake"
    cmake -DAPP=${ROMS_APPLICATION} \
                ${my_hdir} \
                ${ltype} \
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
  elif [[ "${USE_ECBUILD}" == "on" ]]; then
    conf_com="ecbuild"
    ecbuild -DAPP=${ROMS_APPLICATION} \
                  ${my_hdir} \
                  ${ltype} \
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
  else
    echo "Unrecognized value, '${USE_ECBUILD}' set for USE_ECBUILD"
    exit 1
  fi
fi

if [ $? -ne 0 ]; then
  echo "$conf_com did not complete successfully"
  exit 1
fi

#--------------------------------------------------------------------------
# Compile.
#--------------------------------------------------------------------------

if [ $dprint -eq 1 ]; then
  echo $debug:"${!debug}"
else
  if [ $parallel -eq 1 ]; then
    make $NCPUS
  else
    make
  fi
  make install
fi

cd ${MY_PROJECT_DIR}

# Create symlink to executable. This should work even if ROMS was
# linked to the shared library (libROMS.{so|dylib}) because
# CMAKE_BUILD_WITH_INSTALL_RPATH is set to FALSE so that
# RPATH/RUNPATH are set correctly for both the build tree and
# installed locations of the ROMS executable.

if [ $dprint -eq 0 ]; then
  if [[ ! -z "${USE_DEBUG}" && "${USE_DEBUG}" == "on" ]]; then
    ln -sfv ${SCRATCH_DIR}/romsG
  elif [[ ! -z "${USE_MPI}" && "${USE_MPI}" == "on" ]]; then
    ln -sfv ${SCRATCH_DIR}/romsM
  else
    ln -sfv ${SCRATCH_DIR}/romsS
  fi
fi
