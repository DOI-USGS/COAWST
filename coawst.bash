#!/bin/bash
#
# svn $Id: build.bash 429 2009-12-20 17:30:26Z jcwarner $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2018 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: John C. Warner   :::
#                                                                       :::
# ROMS/TOMS Compiling Script                                            :::
# Modified to configure the COAWST Modeling System                      :::
#                                                                       :::
# Script to identify locations of application-specific files for        :::
# compiling the modeling sytem.                                         :::
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
#    ./build.bash [options]                                             :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    -j [N]       Compile in parallel using N CPUs                      :::
#                  omit argument for all available CPUs                 :::
#    -p macro     Prints any Makefile macro value. For example,          :::
#                  build.bash -p FFLAGS                                 :::
#    -noclean     Do not clean already compiled roms objects            :::
#    -nocleanwrf  Do not clean already compiled wrf objects             :::
#    -nocleanww3  Do not clean already compiled ww3 objects             :::
#                                                                       :::
# Notice that sometimes the parallel compilation fail to find MPI       :::
# include file "mpif.h".                                                :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
############################################################################
# Top area here is to set flags from calling this routine. Do not change.
#
parallel=0
clean=1
dprint=0
cleanwrf=1
cleanwrfhydro=1
cleanww3=1

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
      debug="print-$1"
      shift
      ;;

    -noclean )
      shift
      clean=0
      ;;

    -nocleanwrf )
      shift
      cleanwrf=0
      ;;

    -nocleanwrfhydro )
      shift
      cleanwrfhydro=0
      ;;

    -nocleanww3 )
      shift
      cleanww3=0
      ;;

    * )
      echo ""
      echo "$0 : Unknown option [ $1 ]"
      echo ""
      echo "Available Options:"
      echo ""
      echo "-j [N]      Compile in parallel using N CPUs"
      echo "              omit argument for all avaliable CPUs"
      echo "-p macro    Prints any Makefile macro value"
      echo "              For example:  build.bash -p FFLAGS"
      echo "-noclean       Do not clean already compiled objects"
      echo "-nocleanwrf    Do not clean already compiled wrf objects"
      echo "-nocleanwrfhydro Do not clean already compiled wrf objects"
      echo "-nocleanww3    Do not clean already compiled ww3 objects"
      echo ""
      exit 1
      ;;
  esac
done
############################################################################
# Start of USER definitions area:
#
# Set the CPP option defining the particular application. This will
# determine the name of the ".h" header file with the application
# CPP definitions. Also this will activate the switch file for WW3.
export   COAWST_APPLICATION=INLET_TEST

# Set the ROMS_APPLICATION to be the same as the COAWST_APP.
# Do not change this. We use the COAWST APP for other checks.
export   ROMS_APPLICATION=${COAWST_APPLICATION}

# Set a local environmental variable to define the path to the directories
# where all this project's files are kept.
export   MY_ROOT_DIR=/cygdrive/e/data/models/COAWST
export   MY_PROJECT_DIR=${MY_ROOT_DIR}

# The path to the user's local current ROMS source code.
# If using svn locally, this would be the user's Working Copy Path (WCPATH).
# Note that one advantage of maintaining your source code locally with svn
# is that when working simultaneously on multiple machines (e.g. a local
# workstation, a local cluster and a remote supercomputer) you can checkout
# the latest release and always get an up-to-date customized source on each
# machine. This script is designed to more easily allow for differing paths
# to the code and inputs on differing machines.
export   MY_ROMS_SRC=${MY_ROOT_DIR}/

############################################################################
# WRF : Needs to have the env variable NETCDF set.
#export  NETCDF=${NETCDF_INCDIR}/../
#
############################################################################
# Wave Watch 3: Here we provide 5 environment variables for WW3.
#
# 1) COAWST_WW3_DIR is a pointer to root WW3 code, do not change.
export   COAWST_WW3_DIR=${MY_ROOT_DIR}/WW3/model
#
# 2) WWATCH3_NETCDF can be NC3 or NC4. We need NC4 for COAWST. do not change.
export   WWATCH3_NETCDF=NC4
#
# 3) WWATCH_ENV points to WW3 environment listing. do not change.
export   WWATCH_ENV=${COAWST_WW3_DIR}/wwatch.env
#
# 4) NETCDF_CONFIG is needed by WW3. You need to set this:
#export   NETCDF_CONFIG=${NETCDF_LIBDIR}/../bin/nc-config
#    This may require nf-config, depending on your system.
export   NETCDF_CONFIG=/usr/bin/nf-config
#export   NETCDF_CONFIG=/vortexfs1/apps/impistack-1.0/bin/nf-config
#
# 5) WW3_SWITCH_FILE is like cpp options for WW3. You need to create it and
#    list the name here.  You need to have COAWST listed in the switch file.
 export  WW3_SWITCH_FILE=switch_sandy

############################################################################
# Compiler selections.
#
# Set path of the directory containing makefile configuration (*.mk) files.
# The user has the option to specify a customized version of these files
# in a different directory than the one distributed with the source code,
# ${MY_ROMS_SRC}/Compilers. If this is the case, the you need to keep
# these configurations files up-to-date.

#export         COMPILERS=${MY_ROMS_SRC}/Compilers

# Set tunable CPP options.
#
# Sometimes it is desirable to activate one or more CPP options to run
# different variants of the same application without modifying its header
# file. If this is the case, specify each options here using the -D syntax.
# Notice also that you need to use shell's quoting syntax to enclose the
# definition.  Both single or double quotes work. For example,
#
#export      MY_CPP_FLAGS="-DAVERAGES"
#export      MY_CPP_FLAGS="${MY_CPP_FLAGS} -DDEBUGGING"
#
# can be used to write time-averaged fields. Notice that you can have as
# many definitions as you want by appending values.

#export      MY_CPP_FLAGS="-D"

# Other user defined environmental variables. See the ROMS makefile for
# details on other options the user might want to set here. Be sure to
# leave the switches meant to be off set to an empty string or commented
# out. Any string value (including off) will evaluate to TRUE in
# conditional if-statements.

 export           USE_MPI=on            # distributed-memory parallelism
 export        USE_MPIF90=on            # compile with mpif90 script
#export         which_MPI=mpich         # compile with MPICH library
#export         which_MPI=mpich2        # compile with MPICH2 library
 export         which_MPI=openmpi       # compile with OpenMPI library

#export        USE_OpenMP=on            # shared-memory parallelism

 export              FORT=ifort
#export              FORT=gfortran
#export              FORT=pgi

 export         USE_DEBUG=              # use Fortran debugging flags
 export         USE_LARGE=on            # activate 64-bit compilation
 export       USE_NETCDF4=on            # compile with NetCDF-4 library
#export   USE_PARALLEL_IO=on            # Parallel I/O with Netcdf-4/HDF5

#export       USE_MY_LIBS=on            # use my library paths below

# There are several MPI libraries available. Here, we set the desired
# "mpif90" script to use during compilation. This only works if the make
# configuration file (say, Linux-pgi.mk) in the "Compilers" directory
# has the following definition for FC (Fortran Compiler) in the USE_MPI
# section:
#
#              FC := mpif90
#
# that is, "mpif90" defined without any path. Notice that the path
# where the MPI library is installed is computer dependent. Recall
# that you still need to use the appropriate "mpirun" to execute.

if [ -n "${USE_MPIF90:+1}" ]; then
  case "$FORT" in
    ifort )
      if [ "${which_MPI}" = "mpich" ]; then
        export PATH=/opt/intelsoft/mpich/bin:$PATH
      elif [ "${which_MPI}" = "mpich2" ]; then
        export PATH=/opt/intelsoft/mpich2/bin:$PATH
      elif [ "${which_MPI}" = "openmpi" ]; then
        export PATH=/opt/intelsoft/openmpi/bin:$PATH
      fi
      ;;

    pgi )
      if [ "${which_MPI}" = "mpich" ]; then
        export PATH=/opt/pgisoft/mpich/bin:$PATH
      elif [ "${which_MPI}" = "mpich2" ]; then
        export PATH=/opt/pgisoft/mpich2/bin:$PATH
      elif [ "${which_MPI}" = "openmpi" ]; then
        export PATH=/opt/pgisoft/openmpi/bin:$PATH
      fi
      ;;

    gfortran )
      if [ "${which_MPI}" = "mpich2" ]; then
        export PATH=/opt/gfortransoft/mpich2/bin:$PATH
      elif [ "${which_MPI}" = "openmpi" ]; then
        export PATH=/opt/gfortransoft/openmpi/bin:$PATH
      fi
      ;;

  esac
fi
############################################################################
# Additional libraries selections.
#
# If the USE_MY_LIBS is activated above, the path of the libraries
# required by ROMS can be set here using environmental variables
# which take precedence to the values specified in the make macro
# definitions file (Compilers/*.mk). For most applications, only
# the location of the NetCDF library is needed during compilation.
#
# Notice that when the USE_NETCDF4 macro is activated, we need the
# serial or parallel version of the NetCDF-4/HDF5 library. The
# configuration script NF_CONFIG (available since NetCDF 4.0.1)
# is used to set up all the required libraries according to the
# installed options (openDAP, netCDF4/HDF5 file format). The
# parallel library uses the MPI-I/O layer (usually available
# in MPICH2 and OpenMPI) requiring compiling with the selected
# MPI library.
#
# In ROMS distributed-memory applications, you may use either the
# serial or parallel version of the NetCDF-4/HDF5 library. The
# parallel version is required when parallel I/O is activated
# (ROMS cpp option PARALLEL_IO and HDF5).
#
# However, in serial or shared-memory ROMS applications, we need
# to use the serial version of the NetCDF-4/HDF5 to avoid conflicts
# with the compiler. We cannot activate MPI constructs in serial
# or shared-memory ROMS code. Hybrid parallelism is not possible.
#
# Recall also that the MPI library comes in several flavors:
# MPICH, MPICH2, OpenMPI, etc.

if [ -n "${USE_MY_LIBS:+1}" ]; then
  case "$FORT" in
    ifort )
      export             ESMF_OS=Linux
      export       ESMF_COMPILER=ifort
      export           ESMF_BOPT=O
      export            ESMF_ABI=64
      export           ESMF_COMM=mpich
      export           ESMF_SITE=default

      export       ARPACK_LIBDIR=/opt/intelsoft/serial/ARPACK
      if [ -n "${USE_MPI:+1}" ]; then
        if [ "${which_MPI}" = "mpich" ]; then
          export        ESMF_DIR=/opt/intelsoft/mpich/esmf
          export      MCT_INCDIR=/opt/intelsoft/mpich/mct/include
          export      MCT_LIBDIR=/opt/intelsoft/mpich/mct/lib
          export  PARPACK_LIBDIR=/opt/intelsoft/mpich/PARPACK
        elif [ "${which_MPI}" = "mpich2" ]; then
          export        ESMF_DIR=/opt/intelsoft/mpich2/esmf
          export      MCT_INCDIR=/opt/intelsoft/mpich2/mct/include
          export      MCT_LIBDIR=/opt/intelsoft/mpich2/mct/lib
          export  PARPACK_LIBDIR=/opt/intelsoft/mpich2/PARPACK
        elif [ "${which_MPI}" = "openmpi" ]; then
          export        ESMF_DIR=/opt/intelsoft/openmpi/esmf
          export      MCT_INCDIR=/opt/intelsoft/openmpi/mct/include
          export      MCT_LIBDIR=/opt/intelsoft/openmpi/mct/lib
          export  PARPACK_LIBDIR=/opt/intelsoft/openmpi/PARPACK
        fi
      fi

      if [ -n "${USE_NETCDF4:+1}" ]; then
        if [ -n "${USE_PARALLEL_IO:+1}" ] && [ -n "${USE_MPI:+1}" ]; then
          if [ "${which_MPI}" = "mpich" ]; then
            export     NF_CONFIG=/opt/intelsoft/mpich/netcdf4/bin/nf-config
            export NETCDF_INCDIR=/opt/intelsoft/mpich/netcdf4/include
          elif [ "${which_MPI}" = "mpich2" ]; then
            export     NF_CONFIG=/opt/intelsoft/mpich2/netcdf4/bin/nf-config
            export NETCDF_INCDIR=/opt/intelsoft/mpich2/netcdf4/include
          elif [ "${which_MPI}" = "openmpi" ]; then
            export     NF_CONFIG=/opt/intelsoft/openmpi/netcdf4/bin/nf-config
            export NETCDF_INCDIR=/opt/intelsoft/openmpi/netcdf4/include
          fi
        else
          export       NF_CONFIG=/opt/intelsoft/serial/netcdf4/bin/nf-config
          export   NETCDF_INCDIR=/opt/intelsoft/serial/netcdf4/include
        fi
      else
        export     NETCDF_INCDIR=/opt/intelsoft/serial/netcdf3/include
        export     NETCDF_LIBDIR=/opt/intelsoft/serial/netcdf3/lib
      fi
      ;;

    pgi )
      export             ESMF_OS=Linux
      export       ESMF_COMPILER=pgi
      export           ESMF_BOPT=O
      export            ESMF_ABI=64
      export           ESMF_COMM=mpich
      export           ESMF_SITE=default

      export       ARPACK_LIBDIR=/opt/pgisoft/serial/ARPACK
      if [ -n "${USE_MPI:+1}" ]; then
        if [ "${which_MPI}" = "mpich" ]; then
          export        ESMF_DIR=/opt/pgisoft/mpich/esmf
          export      MCT_INCDIR=/opt/pgisoft/mpich/mct/include
          export      MCT_LIBDIR=/opt/pgisoft/mpich/mct/lib
          export  PARPACK_LIBDIR=/opt/pgisoft/mpich/PARPACK
        elif [ "${which_MPI}" = "mpich2" ]; then
          export        ESMF_DIR=/opt/pgisoft/mpich2/esmf
          export      MCT_INCDIR=/opt/pgisoft/mpich2/mct/include
          export      MCT_LIBDIR=/opt/pgisoft/mpich2/mct/lib
          export  PARPACK_LIBDIR=/opt/pgisoft/mpich2/PARPACK
        elif [ "${which_MPI}" = "openmpi" ]; then
          export        ESMF_DIR=/opt/pgisoft/openmpi/esmf
          export      MCT_INCDIR=/opt/pgisoft/openmpi/mct/include
          export      MCT_LIBDIR=/opt/pgisoft/openmpi/mct/lib
          export  PARPACK_LIBDIR=/opt/pgisoft/openmpi/PARPACK
        fi
      fi

      if [ -n "${USE_NETCDF4:+1}" ]; then
        if [ -n "${USE_PARALLEL_IO:+1}" ] && [ -n "${USE_MPI:+1}" ]; then
          if [ "${which_MPI}" = "mpich" ]; then
            export     NF_CONFIG=/opt/pgisoft/mpich/netcdf4/bin/nf-config
            export NETCDF_INCDIR=/opt/pgisoft/mpich/netcdf4/include
          elif [ "${which_MPI}" = "mpich2" ]; then
            export     NF_CONFIG=/opt/pgisoft/mpich2/netcdf4/bin/nf-config
            export NETCDF_INCDIR=/opt/pgisoft/mpich2/netcdf4/include
          elif [ "${which_MPI}" = "openmpi" ]; then
            export     NF_CONFIG=/opt/pgisoft/openmpi/netcdf4/bin/nf-config
            export NETCDF_INCDIR=/opt/pgisoft/openmpi/netcdf4/include
          fi
        else
          export       NF_CONFIG=/opt/pgisoft/serial/netcdf4/bin/nf-config
          export   NETCDF_INCDIR=/opt/pgisoft/serial/netcdf4/include
        fi
      else
        export     NETCDF_INCDIR=/opt/pgisoft/serial/netcdf3/include
        export     NETCDF_LIBDIR=/opt/pgisoft/serial/netcdf3/lib
      fi
      ;;

    gfortran )
      export             ESMF_OS=Linux
      export       ESMF_COMPILER=gfortran
      export           ESMF_BOPT=O
      export            ESMF_ABI=64
      export           ESMF_COMM=mpich
      export           ESMF_SITE=default

      export       ARPACK_LIBDIR=/opt/gfortransoft/serial/ARPACK
      if [ -n "${USE_MPI:+1}" ]; then
        if [ "${which_MPI}" = "mpich2" ]; then
          export        ESMF_DIR=/opt/gfortransoft/mpich2/esmf
          export      MCT_INCDIR=/opt/gfortransoft/mpich2/mct/include
          export      MCT_LIBDIR=/opt/gfortransoft/mpich2/mct/lib
          export  PARPACK_LIBDIR=/opt/gfortransoft/mpich2/PARPACK
        elif [ "${which_MPI}" = "openmpi" ]; then
          export        ESMF_DIR=/opt/gfortransoft/openmpi/esmf
          export      MCT_INCDIR=/opt/gfortransoft/openmpi/mct/include
          export      MCT_LIBDIR=/opt/gfortransoft/openmpi/mct/lib
          export  PARPACK_LIBDIR=/opt/gfortransoft/openmpi/PARPACK
        fi
      fi

      if [ -n "${USE_NETCDF4:+1}" ]; then
        if [ -n "${USE_PARALLEL_IO:+1}" ] && [ -n "${USE_MPI:+1}" ]; then
          if [ "${which_MPI}" = "mpich2" ]; then
            export     NF_CONFIG=/opt/gfortransoft/mpich2/netcdf4/bin/nf-config
            export NETCDF_INCDIR=/opt/gfortransoft/mpich2/netcdf4/include
          elif [ "${which_MPI}" = "openmpi" ]; then
            export     NF_CONFIG=/opt/gfortransoft/openmpi/netcdf4/bin/nf-config
            export NETCDF_INCDIR=/opt/gfortransoft/openmpi/netcdf4/include
          fi
        else
          export       NF_CONFIG=/opt/gfortransoft/serial/netcdf4/bin/nf-config
          export   NETCDF_INCDIR=/opt/gfortransoft/serial/netcdf4/include
        fi
      else
        export     NETCDF_INCDIR=/opt/gfortransoft/serial/netcdf3/include
        export     NETCDF_LIBDIR=/opt/gfortransoft/serial/netcdf3/lib
      fi
      ;;

  esac
fi

############################################################################
# Header and other source directories selections.
#
# The rest of this script sets the path to the users header file and
# analytical source files, if any. See the templates in User/Functionals.
#
# If applicable, use the MY_ANALYTICAL_DIR directory to place your
# customized biology model header file (like fennel.h, nemuro.h, ecosim.h,
# etc).

#  export     MY_HEADER_DIR=${MY_PROJECT_DIR}/ROMS/Include
#  export MY_ANALYTICAL_DIR=${MY_PROJECT_DIR}/ROMS/Functionals
   export     MY_HEADER_DIR=${MY_PROJECT_DIR}/Projects/Inlet_test/Coupled
   export MY_ANALYTICAL_DIR=${MY_PROJECT_DIR}/Projects/Inlet_test/Coupled

# Put the binary to execute in the following directory.

# export            BINDIR=${MY_PROJECT_DIR}
  export            BINDIR=./

# Put the f90 files in a project specific Build directory to avoid conflict
# with other projects.

# export       SCRATCH_DIR=${MY_PROJECT_DIR}/Build
  export       SCRATCH_DIR=./Build

# Go to the users source directory to compile. The options set above will
# pick up the application-specific code from the appropriate place.

############################################################################
# End of USER definitions area.  You really should not change anything 
# down here.
#
 cd ${MY_ROMS_SRC}

#--------------------------------------------------------------------------
# Compile.
#--------------------------------------------------------------------------
# Remove build directory.

if [ $clean -eq 1 ]; then
  make clean
fi

# Compile (the binary will go to BINDIR set above).

# first go to make some coupler files
if [ $clean -eq 1 ] && [ $cleanwrf -eq 1 ] && [ $cleanwrfhydro -eq 1 ] && [ $cleanww3 -eq 1 ] ; then
  make mct_params
fi
  cd ${SCRATCH_DIR}
  export MCT_PARAMS_DIR=${PWD}
  cd ${MY_ROMS_SRC}

  export WRF_DIR=${MY_ROMS_SRC}/WRF
if [ $cleanwrf -eq 1 ]; then
  make wrfclean
  cd ${MY_ROMS_SRC}
fi
  make wrf

  export WRFHYDRO_DIR=${MY_ROMS_SRC}/WRF/hydro_v5.0
if [ $cleanwrfhydro -eq 1 ]; then
  make wrfhydroclean
  cd ${MY_ROMS_SRC}
fi
  make wrfhydro

if [ $cleanww3 -eq 1 ]; then
  make ww3clean
  cd ${MY_ROMS_SRC}
fi
  make ww3

if [ $dprint -eq 1 ]; then
  make $debug
else
  if [ $parallel -eq 1 ]; then
    make $NCPUS
  else
    make
  fi
fi
