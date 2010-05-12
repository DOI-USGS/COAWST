#!/bin/bash
#
# svn $Id: build.bash 429 2009-12-20 17:30:26Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2010 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# ROMS/TOMS Compiling Script                                            :::
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
#    ./build.bash [options]                                             :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    -j [N]      Compile in parallel using N CPUs                       :::
#                  omit argument for all available CPUs                 :::
#    -noclean    Do not clean already compiled objects                  :::
#                                                                       :::
# Notice that sometimes the parallel compilation fail to find MPI       :::
# include file "mpif.h".                                                :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

parallel=0
clean=1

while [ $# -gt 0 ]
do
  case "$1" in
    -j )
      shift
      parallel=1
      test=`echo $1 | grep -P '^\d+$'`
      if [ "$test" != "" ]; then
        NCPUS="-j $1"
        shift
      else
        NCPUS="-j"
      fi
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
      echo "-noclean    Do not clean already compiled objects"
      echo ""
      exit 1
      ;;
  esac
done

# Set the CPP option defining the particular application. This will
# determine the name of the ".h" header file with the application
# CPP definitions.

export   ROMS_APPLICATION=UPWELLING

# Set number of nested/composed/mosaic grids.  Currently, only one grid
# is supported.

export        NestedGrids=1

# Set a local environmental variable to define the path to the directories
# where all this project's files are kept.

export        MY_ROOT_DIR=/home/arango/ocean/toms/repository
export     MY_PROJECT_DIR=${PWD}

# The path to the user's local current ROMS source code.
#
# If using svn locally, this would be the user's Working Copy Path (WCPATH).
# Note that one advantage of maintaining your source code locally with svn
# is that when working simultaneously on multiple machines (e.g. a local
# workstation, a local cluster and a remote supercomputer) you can checkout
# the latest release and always get an up-to-date customized source on each
# machine. This script is designed to more easily allow for differing paths
# to the code and inputs on differing machines.

export        MY_ROMS_SRC=${MY_ROOT_DIR}/branches/arango

# Set path of the directory containing makefile configuration (*.mk) files.
# The user has the option to specify a customized version of these files
# in a different directory than the one distributed with the source code,
# ${MY_ROMS_SCR}/Compilers. If this is the case, the you need to keep
# these configurations files up-to-date.

#export         COMPILERS=${MY_ROMS_SRC}/Compilers

# Set tunable CPP options.
#
# Sometimes it is desirable to activate one or more CPP options to run
# different variants of the same application without modifying its header
# file. If this is the case, specify each options here using the -D syntax.
# Notice also that you need to use shell's quoting syntax to enclose the
# definition.  Both single or double quotes works. For example, to write
# time-averaged fields set:
#
#export      MY_CPP_FLAGS="-DAVERAGES"

# Other user defined environmental variables. See the ROMS makefile for
# details on other options the user might want to set here. Be sure to
# leave the switched meant to be off set to an empty string or commented
# out. Any string value (including off) will evaluate to TRUE in
# conditional if-stamentents.

 export           USE_MPI=on
 export        USE_MPIF90=on
 export              FORT=ifort

#export        USE_OpenMP=on

#export         USE_DEBUG=on
 export         USE_LARGE=on
#export       USE_NETCDF4=on

# There are several MPI libraries out there. The user can select here the
# appropriate "mpif90" script to compile, provided that the makefile
# macro file (say, Linux-pgi.mk) in the Compilers directory has:
#
#              FC := mpif90
#
# "mpif90" defined without any path. Recall that you still need to use the
# appropriate "mpirun" to execute. Also notice that the path where the
# MPI library is installed is computer dependent.

if [ -n "${USE_MPIF90:+1}" ]; then
  case "$FORT" in
    ifort )
#     export PATH=/opt/intelsoft/mpich/bin:$PATH
      export PATH=/opt/intelsoft/mpich2/bin:$PATH
#     export PATH=/opt/intelsoft/openmpi/bin:$PATH
      ;;

    pgi )
#     export PATH=/opt/pgisoft/mpich/bin:$PATH
      export PATH=/opt/pgisoft/mpich2/bin:$PATH
#     export PATH=/opt/pgisoft/openmpi/bin:$PATH
      ;;

    gfortran )
#     export PATH=/opt/gfortransoft/mpich2/bin:$PATH
      export PATH=/opt/gfortansoft/openmpi/bin:$PATH
      ;;

    g95 )
#     export PATH=/opt/g95soft/mpich2/bin:$PATH
      export PATH=/opt/g95soft/openmpi/bin:$PATH
      ;;

  esac
fi

# The path of the libraries required by ROMS can be set here using
# environmental variables which take precedence to the values
# specified in the makefile macro definitions file (Compilers/*.mk).
# If so desired, uncomment the local USE_MY_LIBS definition below
# and edit the paths to your values. For most applications, only
# the location of the NetCDF library (NETCDF_LIBDIR) and include
# directorry (NETCDF_INCDIR) are needed!
#
# Notice that when the USE_NETCDF4 macro is activated, we need a
# serial and parallel version of the NetCDF-4/HDF5 library. The
# parallel library uses parallel I/O through MPI-I/O so we need
# compile also with the MPI library. This is fine in ROMS
# distributed-memory applications.  However, in serial or
# shared-memory ROMS applications we need to use the serial
# version of the NetCDF-4/HDF5 to avoid conflicts with the
# compiler. Recall also that the MPI library comes in several
# flavors: MPICH, MPICH2, OpenMPI.

#export           USE_MY_LIBS=on

if [ -n "${USE_MY_LIBS:+1}" ]; then
  case "$FORT" in
    ifort )
      export           ESMF_DIR=/opt/intelsoft/esmf-3.1.0
      export            ESMF_OS=Linux
      export      ESMF_COMPILER=ifort
      export          ESMF_BOPT=O
      export           ESMF_ABI=64
      export          ESMF_COMM=mpich
      export          ESMF_SITE=default
      export         MCT_INCDIR=/opt/intelsoft/mct/include
      export         MCT_LIBDIR=/opt/intelsoft/mct/lib

      export      ARPACK_LIBDIR=/opt/intelsoft/serial/ARPACK
      if [ -n "${USE_MPI:+1}" ]; then
#       export   PARPACK_LIBDIR=/opt/intelsoft/mpich/PARPACK
        export   PARPACK_LIBDIR=/opt/intelsoft/mpich2/PARPACK
#       export   PARPACK_LIBDIR=/opt/intelsoft/openmpi/PARPACK
      fi

      if [ -n "${USE_NETCDF4:+1}" ]; then
        if [ -n "${USE_MPI:+1}" ]; then
#         export  NETCDF_INCDIR=/opt/intelsoft/mpich/netcdf4/include
#         export  NETCDF_LIBDIR=/opt/intelsoft/mpich/netcdf4/lib
#         export    HDF5_LIBDIR=/opt/intelsoft/mpich/hdf5/lib

          export  NETCDF_INCDIR=/opt/intelsoft/mpich2/netcdf4/include
          export  NETCDF_LIBDIR=/opt/intelsoft/mpich2/netcdf4/lib
          export    HDF5_LIBDIR=/opt/intelsoft/mpich2/hdf5/lib

#         export  NETCDF_INCDIR=/opt/intelsoft/openmpi/netcdf4/include
#         export  NETCDF_LIBDIR=/opt/intelsoft/openmpi/netcdf4/lib
#         export    HDF5_LIBDIR=/opt/intelsoft/openmpi/hdf5/lib
        else
          export  NETCDF_INCDIR=/opt/intelsoft/serial/netcdf4/include
          export  NETCDF_LIBDIR=/opt/intelsoft/serial/netcdf4/lib
          export    HDF5_LIBDIR=/opt/intelsoft/serial/hdf5/lib
        fi
      else
        export    NETCDF_INCDIR=/opt/intelsoft/serial/netcdf3/include
        export    NETCDF_LIBDIR=/opt/intelsoft/serial/netcdf3/lib
      fi
      ;;

    pgi )
      export           ESMF_DIR=/opt/pgisoft/esmf-3.1.0
      export            ESMF_OS=Linux
      export      ESMF_COMPILER=pgi
      export          ESMF_BOPT=O
      export           ESMF_ABI=64
      export          ESMF_COMM=mpich
      export          ESMF_SITE=default
      export         MCT_INCDIR=/opt/pgisoft/mct/include
      export         MCT_LIBDIR=/opt/pgisoft/mct/lib

      export      ARPACK_LIBDIR=/opt/pgisoft/serial/ARPACK
      if [ -n "${USE_MPI:+1}" ]; then
        export   PARPACK_LIBDIR=/opt/pgisoft/mpich/PARPACK
#       export   PARPACK_LIBDIR=/opt/pgisoft/mpich2/PARPACK
#       export   PARPACK_LIBDIR=/opt/pgisoft/openmpi/PARPACK
      fi

      if [ -n "${USE_NETCDF4:+1}" ]; then
        if [ -n "${USE_MPI:+1}" ]; then
#         export  NETCDF_INCDIR=/opt/pgisoft/mpich/netcdf4/include
#         export  NETCDF_LIBDIR=/opt/pgisoft/mpich/netcdf4/lib
#         export    HDF5_LIBDIR=/opt/pgisoft/mpich/hdf5/lib

          export  NETCDF_INCDIR=/opt/pgisoft/mpich2/netcdf4/include
          export  NETCDF_LIBDIR=/opt/pgisoft/mpich2/netcdf4/lib
          export    HDF5_LIBDIR=/opt/pgisoft/mpich2/hdf5/lib

#         export  NETCDF_INCDIR=/opt/pgisoft/openmpi/netcdf4/include
#         export  NETCDF_LIBDIR=/opt/pgisoft/openmpi/netcdf4/lib
#         export    HDF5_LIBDIR=/opt/pgisoft/openmpi/hdf5/lib
        else
          export  NETCDF_INCDIR=/opt/pgisoft/serial/netcdf4/include
          export  NETCDF_LIBDIR=/opt/pgisoft/serial/netcdf4/lib
          export    HDF5_LIBDIR=/opt/pgisoft/serial/hdf5/lib
        fi
      else
        export    NETCDF_INCDIR=/opt/pgisoft/serial/netcdf3/include
        export    NETCDF_LIBDIR=/opt/pgisoft/serial/netcdf3/lib
      fi
      ;;

    gfortran )
      export         MCT_INCDIR=/opt/gfortransoft/mct/include
      export         MCT_LIBDIR=/opt/gfortransoft/mct/lib

      export      ARPACK_LIBDIR=/opt/gfortransoft/serial/ARPACK
      if [ -n "${USE_MPI:+1}" ]; then
#       export   PARPACK_LIBDIR=/opt/gfortransoft/mpich/PARPACK
        export   PARPACK_LIBDIR=/opt/gfortransoft/mpich2/PARPACK
#       export   PARPACK_LIBDIR=/opt/gfortransoft/openmpi/PARPACK
      fi

      if [ -n "${USE_NETCDF4:+1}" ]; then
        if [ -n "${USE_MPI:+1}" ]; then
#         export  NETCDF_INCDIR=/opt/gfortransoft/mpich/netcdf4/include
#         export  NETCDF_LIBDIR=/opt/gfortransoft/mpich/netcdf4/lib
#         export    HDF5_LIBDIR=/opt/gfortransoft/mpich/hdf5/lib

          export  NETCDF_INCDIR=/opt/gfortransoft/mpich2/netcdf4/include
          export  NETCDF_LIBDIR=/opt/gfortransoft/mpich2/netcdf4/lib
          export    HDF5_LIBDIR=/opt/gfortransoft/mpich2/hdf5/lib

#         export  NETCDF_INCDIR=/opt/gfortransoft/openmpi/netcdf4/include
#         export  NETCDF_LIBDIR=/opt/gfortransoft/openmpi/netcdf4/lib
#         export    HDF5_LIBDIR=/opt/gfortransoft/openmpi/hdf5/lib
        else
          export  NETCDF_INCDIR=/opt/gfortransoft/serial/netcdf4/include
          export  NETCDF_LIBDIR=/opt/gfortransoft/serial/netcdf4/lib
          export    HDF5_LIBDIR=/opt/gfortransoft/serial/hdf5/lib
        fi
      else
          export  NETCDF_INCDIR=/opt/gfortransoft/serial/netcdf3/include
          export  NETCDF_LIBDIR=/opt/gfortransoft/serial/netcdf3/lib
      fi
      ;;

    g95 )
      export         MCT_INCDIR=/opt/g95soft/mct/include
      export         MCT_LIBDIR=/opt/g95soft/mct/lib

      export      ARPACK_LIBDIR=/opt/g95soft/serial/ARPACK
      if [ -n "${USE_MPI:+1}" ]; then
#       export   PARPACK_LIBDIR=/opt/g95soft/mpich/PARPACK
        export   PARPACK_LIBDIR=/opt/g95soft/mpich2/PARPACK
        export   PARPACK_LIBDIR=/opt/g95soft/openmpi/PARPACK
      fi

      if [ -n "${USE_NETCDF4:+1}" ]; then
        if [ -n "${USE_MPI:+1}" ]; then
          export  NETCDF_INCDIR=/opt/g95soft/mpich2/netcdf4/include
          export  NETCDF_LIBDIR=/opt/g95soft/mpich2/netcdf4/lib
          export    HDF5_LIBDIR=/opt/g95soft/mpich2/hdf5/lib

#         export  NETCDF_INCDIR=/opt/g95soft/openmpi/netcdf4/include
#         export  NETCDF_LIBDIR=/opt/g95soft/openmpi/netcdf4/lib
#         export    HDF5_LIBDIR=/opt/g95soft/openmpi/hdf5/lib
        else
          export  NETCDF_INCDIR=/opt/g95soft/serial/netcdf4/include
          export  NETCDF_LIBDIR=/opt/g95soft/serial/netcdf4/lib
          export    HDF5_LIBDIR=/opt/g95soft/serial/hdf5/lib
        fi
      else
        export    NETCDF_INCDIR=/opt/g95soft/serial/netcdf3/include
        export    NETCDF_LIBDIR=/opt/g95soft/serial/netcdf3/lib
      fi
      ;;

  esac
fi

# The rest of this script sets the path to the users header file and
# analytical source files, if any. See the templates in User/Functionals.
#
# If applicable, use the MY_ANALYTICAL_DIR directory to place your
# customized biology model header file (like fennel.h, nemuro.h, ecosim.h,
# etc).

 export     MY_HEADER_DIR=${MY_PROJECT_DIR}

 export MY_ANALYTICAL_DIR=${MY_PROJECT_DIR}

# Put the binary to execute in the following directory.

 export            BINDIR=${MY_PROJECT_DIR}

# Put the f90 files in a project specific Build directory to avoid conflict
# with other projects.

 export       SCRATCH_DIR=${MY_PROJECT_DIR}/Build

# Go to the users source directory to compile. The options set above will
# pick up the application-specific code from the appropriate place.

 cd ${MY_ROMS_SRC}

# Remove build directory.

if [ $clean -eq 1 ]; then
  make clean
fi

# Compile (the binary will go to BINDIR set above).

if [ $parallel -eq 1 ]; then
  make $NCPUS
else
  make
fi
