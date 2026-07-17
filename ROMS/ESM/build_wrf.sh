#!/bin/bash
#
# git $Id$
# svn $Id: build_wrf.sh 1151 2023-02-09 03:08:53Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# WRF Compiling BASH Script: WRF Versions 4.1 and up                    :::
#                                                                       :::
# Script to compile WRF where source code files are kept separate       :::
# from the application configuration and build objects.                 :::
#                                                                       :::
# Q: How/why does this script work?                                     :::
#                                                                       :::
# A:                                                                    :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    ./build_wrf.sh [options]                                           :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    -j [N]      Compile in parallel using N CPUs                       :::
#                  omit argument for all available CPUs                 :::
#                                                                       :::
#    -move       Move compiled objects to build directory               :::
#                                                                       :::
#    -noclean    Do not run clean -a script                             :::
#                                                                       :::
#    -noconfig   Do not run configure compilation script                :::
#                                                                       :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

export which_MPI=openmpi                      # default, overwritten below

#Initialize.

separator=`perl -e "print ':' x 100;"`

parallel=0
clean=1
config=1
move=0

export CPLFLAG=''
export MY_CPP_FLAGS=''

while [ $# -gt 0 ]
do
  case "$1" in
    -j )
      shift
      parallel=1
      test=`echo $1 | grep '^[0-9]\+$'`
      if [ "$test" != "" ]; then
        NCPUS="$1"
        shift
      else
        NCPUS="2"
      fi
      ;;

    -move )
      shift
      move=1
      ;;

    -noclean )
      shift
      clean=0
      ;;

    -noconfig )
      shift
      config=0
      ;;

    * )
      echo ""
      echo "$0 : Unknown option [ $1 ]"
      echo ""
      echo "Available Options:"
      echo ""
      echo "-j [N]        Compile in parallel using N CPUs"
      echo ""
      echo "-move         Move compiled objects to build directory"
      echo ""
      echo "-noclean      Do not run clean script"
      echo ""
      echo "-nocconfig    Do not run configure compilation script"
      echo ""
      exit 1
      ;;
  esac
done

#--------------------------------------------------------------------------
# Set a local environmental variable to define the root path to the
# directories where ROMS and WRF source files are located.  The ROMS
# source directory is needed for replacing several WRF files for
# ESMF/NUOPC Coupling (see below).
#--------------------------------------------------------------------------

 export ROMS_SRC_DIR=${HOME}/ocean/repository/trunk

#export WRF_ROOT_DIR=${HOME}/ocean/repository/WRF.4.1.2
#export WRF_ROOT_DIR=${HOME}/ocean/repository/WRF.4.1.3
#export WRF_ROOT_DIR=${HOME}/ocean/repository/WRF.4.2.2
#export WRF_ROOT_DIR=${HOME}/ocean/repository/WRF.4.3
 export WRF_ROOT_DIR=${HOME}/ocean/repository/WRF

 export WRF_SRC_DIR=${WRF_ROOT_DIR}

# Decode WRF version from its README file to decide the appropriate data
# file links needed.

 wrf_ver=`grep 'WRF Model Version ' $WRF_ROOT_DIR/README | sed -e 's/[^0-9]*\([0-9]\.[0-9]\).*/\1/'`
 export WRF_VERSION=${wrf_ver}

#--------------------------------------------------------------------------
# Set a local environmental variable to define the path of the working
# application directory where all this project's files are kept.
#--------------------------------------------------------------------------

 export MY_PROJECT_DIR=${PWD}
 export MY_PROJECT_DATA=`dirname ${PWD}`/Data

#--------------------------------------------------------------------------
# COAMPS configuration CPP options.
#--------------------------------------------------------------------------

# Sometimes it is desirable to activate one or more CPP options to
# configure a particular application. If it is the case, specify each
# option here using the -D syntax. Notice also that you need to use
# shell's quoting syntax to enclose the definition. Both single or
# double quotes work. For example,
#
# export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DHDF5"
#
# Notice that you can have as many definitions as you want by appending
# values.

# export MY_CPP_FLAGS="${MY_CPP_FLAGS} -D"

#--------------------------------------------------------------------------
# Set Fortran compiler and MPI library to use.
#--------------------------------------------------------------------------

# Set path of the directory containing makefile configuration (*.mk) files.
# The user has the option to specify a customized version of these files
# in a different directory than the one distributed with the source code,
# ${COAMPS_ROOT_DIR}/compilers. If this is the case, the you need to keep
# these configurations files up-to-date.

 export USE_MPI=on             # distributed-memory parallelism
 export USE_MPIF90=on          # compile with mpif90 script
#export which_MPI=intel        # compile with mpiifort library
#export which_MPI=mpich        # compile with MPICH library
#export which_MPI=mpich2       # compile with MPICH2 library
#export which_MPI=mvapich2     # compile with MVAPICH2 library
 export which_MPI=openmpi      # compile with OpenMPI library

 export FORT=ifort
#export FORT=gfortran
#export FORT=pgi

#export USE_REAL_DOUBLE=on          # use real double precision (-r8)
#export USE_DEBUG=on                # use Fortran debugging flags
 export USE_NETCDF=on               # compile with NetCDF
 export USE_NETCDF4=on              # compile with NetCDF-4 library
                                    # (Must also set USE_NETCDF)

 export -n PNETCDF                  # disable compiling with pNetCDF

#--------------------------------------------------------------------------
# Use my specified library paths. It is not needed but it is added for
# for checking in the future.
#--------------------------------------------------------------------------

 export USE_MY_LIBS=no           # use system default library paths
#export USE_MY_LIBS=yes          # use my customized library paths

 MY_PATHS=${ROMS_SRC_DIR}/Compilers/my_build_paths.sh
#MY_PATHS=${HOME}/Compilers/ROMS/my_build_paths.sh

if [ "${USE_MY_LIBS}" == 'yes' ]; then
  source ${MY_PATHS} ${MY_PATHS}
fi

#--------------------------------------------------------------------------
# WRF build and executable directory.
#--------------------------------------------------------------------------

# Put the *.a, .f and .f90 files in a project specific Build directory to
# avoid conflict with other projects.

if [ -n "${USE_DEBUG:+1}" ]; then
  export  WRF_BUILD_DIR=${MY_PROJECT_DIR}/Build_wrfG
else
  export  WRF_BUILD_DIR=${MY_PROJECT_DIR}/Build_wrf
fi

# Put WRF executables in the following directory.

export  WRF_BIN_DIR=${WRF_BUILD_DIR}/Bin

# Go to the users source directory to compile. The options set above will
# pick up the application-specific code from the appropriate place.

 cd ${WRF_ROOT_DIR}

#--------------------------------------------------------------------------
# Configure. It creates configure.wrf script used for compilation.
#
#   -d    build with debugging information and no optimization
#   -D    build with -d AND floating traps, traceback, uninitialized variables
#   -r8   build with 8 byte reals
#
# During configuration the WRF/arch/Config.pl perl script is executed and
# we need to interactively select the combination of compiler and parallel
# comunications option. For example, for Darwin operating system we get:
#
# Please select from among the following Darwin ARCH options:
#
#  1. (serial)   2. (smpar)   3. (dmpar)   4. (dm+sm)   PGI (pgf90/pgcc)
#  5. (serial)   6. (smpar)   7. (dmpar)   8. (dm+sm)   INTEL (ifort/icc)
#  9. (serial)  10. (smpar)  11. (dmpar)  12. (dm+sm)   INTEL (ifort/clang)
# 13. (serial)               14. (dmpar)                GNU (g95/gcc)
# 15. (serial)  16. (smpar)  17. (dmpar)  18. (dm+sm)   GNU (gfortran/gcc)
# 19. (serial)  20. (smpar)  21. (dmpar)  22. (dm+sm)   GNU (gfortran/clang)
# 23. (serial)               24. (dmpar)                IBM (xlf90_r/cc)
# 25. (serial)  26. (smpar)  27. (dmpar)  28. (dm+sm)   PGI (pgf90/pgcc): -f90=pgf90
# 29. (serial)  30. (smpar)  31. (dmpar)  32. (dm+sm)   INTEL (ifort/icc): Open MPI
# 33. (serial)  34. (smpar)  35. (dmpar)  36. (dm+sm)   INTEL/GNU (ifort/gcc): Open MPI
# 37. (serial)  38. (smpar)  39. (dmpar)  40. (dm+sm)   GNU (gfortran/gcc): Open MPI
#
# Enter selection [1-40] : ??
#
# For coupling with ESMF/NUOPC, we need select an option from the (dmpar)
# column for distributed-memory configuration.
#--------------------------------------------------------------------------

# Clean source code.

if [ "$clean" -eq "1" ]; then
  echo ""
  echo "${separator}"
  echo "Cleaning WRF source code:  ${WRF_ROOT_DIR}/clean -a"
  echo "${separator}"
  echo ""
  ${WRF_ROOT_DIR}/clean -a
fi

if [ -n "${USE_DEBUG:+1}" ]; then
#  DEBUG_FLAG="-d"
   DEBUG_FLAG="-D"
fi

export CONFIG_FLAGS=''

if [ "$config" -eq "1" ]; then
  if [ -n "${USE_DEBUG:+1}" -a -n "${USE_REAL_DOUBLE:+1}" ]; then
    export CONFIG_FLAGS="${DEBUG_FLAG} -r8"
  elif [ -n "${USE_DEBUG:+1}" ]; then
    export CONFIG_FLAGS="${DEBUG_FLAG}"
  elif [ "${USE_REAL_DOUBLE:+1}" ]; then
    export CONFIG_FLAGS="-r8"
  fi

# Check if WRF needs to be patched and do so if necessary.

  ${ROMS_SRC_DIR}/ESM/wrf_patch.sh

  echo ""
  echo "${separator}"
  echo "Configuring WRF code:  ${WRF_ROOT_DIR}/configure ${CONFIG_FLAGS}"
  echo "${separator}"
  echo ""

  ${WRF_ROOT_DIR}/configure ${CONFIG_FLAGS}

# If which_MPI is "intel" then we need to replace DM_FC and DM_CC in configure.wrf

  if [ "${which_MPI}" == "intel" ]; then
    perl -i -pe 's/^DM_FC(\s*)=(\s*)mpif90/DM_FC$1=$2mpiifort/' ${WRF_SRC_DIR}/configure.wrf
    perl -i -pe 's/^DM_CC(\s*)=(\s*)mpicc/DM_CC$1=$2mpiicc/' ${WRF_SRC_DIR}/configure.wrf
  fi
fi

#--------------------------------------------------------------------------
# WRF Compile script:
#
# Usage:
#
#     compile [-j n] wrf   compile wrf in run dir (NOTE: no real.exe,
#                          ndown.exe, or ideal.exe generated)
#
#   or choose a test case (see README_test_cases for details) :
#
#     compile [-j n] em_b_wave
#     compile [-j n] em_convrad
#     compile [-j n] em_esmf_exp
#     compile [-j n] em_fire
#     compile [-j n] em_grav2d_x
#     compile [-j n] em_heldsuarez
#     compile [-j n] em_hill2d_x
#     compile [-j n] em_les
#     compile [-j n] em_quarter_ss
#     compile [-j n] em_real
#     compile [-j n] em_scm_xy
#     compile [-j n] em_seabreeze2d_x
#     compile [-j n] em_squall2d_x
#     compile [-j n] em_squall2d_y
#     compile [-j n] em_tropical_cyclone
#     compile [-j n] nmm_real
#     compile [-j n] nmm_tropical_cyclone
#
#     compile -j n         parallel make using n tasks if supported
#                          (default 2)
#
#     compile -h           help message
#--------------------------------------------------------------------------

export WRF_DA_CORE=0             # no Data Assimilation core
export WRF_EM_CORE=1             # Eurelian Mass-coordinate core
export WRF_NMM_CORE=0            # Nonhydrostatic Mesoscale Model core

# Remove existing build directory.

if [[ "$move" -eq "1" && "$clean" -eq "1" ]]; then
  echo ""
  echo "${separator}"
  echo "Removing WRF build directory:  ${WRF_BUILD_DIR}"
  echo "${separator}"
  echo ""
  /bin/rm -rf ${WRF_BUILD_DIR}
fi

# Compile (if -move is set, the binaries will go to WRF_BIN_DIR set above).

#export WRF_CASE=wrf
 export WRF_CASE=em_real

if [ "$parallel" -eq "1" ]; then
  export J="-j ${NCPUS}"
else
  export J="-j 2"
fi

echo ""
echo "${separator}"
echo "Compiling WRF using  ${MY_PROJECT_DIR}/${0}:"
echo ""
echo "   ${WRF_ROOT_DIR}/compile ${WRF_CASE}"
echo "        WRF_DA_CORE = ${WRF_DA_CORE},    Data Assimilation core"
echo "        WRF_EM_CORE = ${WRF_EM_CORE},    Eurelian Mass-coordinate core"
echo "        WRF_NMM_CORE = ${WRF_NMM_CORE}, Nonhydrostatic Mesoscale Model core"
echo "        J = ${J},          number of compiling CPUs"
echo "${separator}"
echo ""

${WRF_ROOT_DIR}/compile ${WRF_CASE}

#--------------------------------------------------------------------------
# Move WRF objects and executables to working project directory.
#--------------------------------------------------------------------------

if [ "$move" -eq "1" ]; then
  ${ROMS_SRC_DIR}/ESM/wrf_move.sh
fi

#--------------------------------------------------------------------------
# Create WRF data links in the working project directory.
#--------------------------------------------------------------------------

if [[ "$move" -eq "1" && "$WRF_CASE" == "em_real" ]]; then
  ${ROMS_SRC_DIR}/ESM/wrf_links.sh
fi
