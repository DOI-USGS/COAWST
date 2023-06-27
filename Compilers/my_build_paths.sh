#!/bin/bash
#
# git $Id$
# svn $Id: my_build_paths.sh 1151 2023-02-09 03:08:53Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# ROMS Customized Compiling Libraries BASH Script:                      :::
#                                                                       :::
# This bash script sets the customized library paths needed by the      :::
# build script when the enviromental variable USE_MY_LIBS has a 'yes'   :::
# value.                                                                :::
#                                                                       :::
# For example, in build_roms.sh we have:                                :::
#                                                                       :::
#       if [ "${USE_MY_LIBS}" = "yes" ]; then                           :::
#         source ${COMPILERS}/my_build_paths.sh                         :::
#       fi                                                              :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

separator=`perl -e "print ':' x 100;"`

echo ""
echo "${separator}"
echo "Using customized library paths from:  $1"
echo "${separator}"
echo ""

#--------------------------------------------------------------------------
# Add MPI library to compile.
#--------------------------------------------------------------------------
#
# Recall also that the MPI library comes in several flavors:
# MPICH, MPICH2, OpenMPI, etc.

# There are several MPI libraries available (MPICH, MPICH2, OpenMPI,
# etc.). Here, we set the desired  "mpif90" script to use during
# compilation. This only works if the make configuration file (say,
# Linux-ifort.mk) in the "Compilers" directory has the following
# definition for FC (Fortran Compiler) in the USE_MPI section:
#
#              FC := mpif90
#
# that is, "mpif90" defined without any path. Notice that the path
# where the MPI library is installed is computer dependent. Recall
# that you still need to use the appropriate "mpirun" to execute.

export               MPI_ROOT=""

if [ -n "${USE_MPIF90:+1}" ]; then
  case "$FORT" in
    ifort )
      if [ "${which_MPI}" = "mpich" ]; then
        export       MPI_ROOT=/opt/intelsoft/mpich
      elif [ "${which_MPI}" = "mpich2" ]; then
        export       MPI_ROOT=/opt/intelsoft/mpich2
      elif [ "${which_MPI}" = "openmpi" ]; then
        export       MPI_ROOT=/opt/intelsoft/openmpi
      elif [ "${which_MPI}" = "mvapich2" ]; then
        export       MPI_ROOT=/opt/intelsoft/mvapich2
      fi
      export             PATH=${MPI_ROOT}/bin:$PATH
      export       MPI_INCDIR=${MPI_ROOT}/include
      ;;

    pgi )
      if [ "${which_MPI}" = "mpich" ]; then
        export       MPI_ROOT=/opt/pgisoft/mpich
      elif [ "${which_MPI}" = "mpich2" ]; then
        export       MPI_ROOT=/opt/pgisoft/mpich2
      elif [ "${which_MPI}" = "openmpi" ]; then
        export       MPI_ROOT=/opt/pgisoft/openmpi
      elif [ "${which_MPI}" = "mvapich2" ]; then
        export       MPI_ROOT=/opt/pgisoft/mvapich2
      fi
      export             PATH=${MPI_ROOT}/bin:$PATH
      export       MPI_INCDIR=${MPI_ROOT}/include
      ;;

    gfortran )
      if [ "${which_MPI}" = "mpich2" ]; then
        export       MPI_ROOT=/opt/gfortransoft/mpich2
      elif [ "${which_MPI}" = "openmpi" ]; then
        export       MPI_ROOT=/opt/gfortransoft/openmpi
      elif [ "${which_MPI}" = "mvapich2" ]; then
        export       MPI_ROOT=/opt/gfortransoft/mvapich2
      fi
      export             PATH=${MPI_ROOT}/bin:$PATH
      export       MPI_INCDIR=${MPI_ROOT}/include
      ;;

  esac
fi

#--------------------------------------------------------------------------
# Set libraries to compile and link.
#--------------------------------------------------------------------------
#
# The path of the libraries required by ROMS can be set here using
# environmental variables which take precedence to the values
# specified in the make macro definitions file (Compilers/*.mk).
# For most applications, only the location of the NetCDF library
# is needed during compilation.
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

export               MPI_SOFT=""

case "$FORT" in
  ifort )
    export      ESMF_COMPILER=intelgcc
#   export      ESMF_COMPILER=intel
    if [ -n "${USE_DEBUG:+1}" ]; then
      export        ESMF_BOPT=g
    else
      export        ESMF_BOPT=O
    fi
    export           ESMF_ABI=64
    export          ESMF_COMM=${which_MPI}
    export          ESMF_SITE=default

    export      ARPACK_LIBDIR=/opt/intelsoft/serial/ARPACK
    if [ -n "${USE_MPI:+1}" ]; then
      if [ "${which_MPI}" = "mpich" ]; then
        export       MPI_SOFT=/opt/intelsoft/mpich
      elif [ "${which_MPI}" = "mpich2" ]; then
        export       MPI_SOFT=/opt/intelsoft/mpich2
      elif [ "${which_MPI}" = "openmpi" ]; then
        export       MPI_SOFT=/opt/intelsoft/openmpi
      elif [ "${which_MPI}" = "mvapich2" ]; then
        export       MPI_SOFT=/opt/intelsoft/mvapich2
      fi
      export       MCT_INCDIR=${MPI_SOFT}/mct/include
      export       MCT_LIBDIR=${MPI_SOFT}/mct/lib
      export   PARPACK_LIBDIR=${MPI_SOFT}/PARPACK
    fi

    if [ ! -n "${SINGULARITY_COMMAND:+1}" ]; then
      if [ -n "${USE_NETCDF4:+1}" ]; then
        if [ -n "${USE_PARALLEL_IO:+1}" ] && [ -n "${USE_MPI:+1}" ]; then
          export       ESMF_DIR=${MPI_SOFT}/esmf_nc4
          export         NETCDF=${MPI_SOFT}/netcdf4
          export      NF_CONFIG=${NETCDF}/bin/nf-config
          export  NETCDF_INCDIR=${NETCDF}/include
          export        NETCDF4=1
        else
          export       ESMF_DIR=${MPI_SOFT}/esmf_nc4
          export         NETCDF=/opt/intelsoft/serial/netcdf4
          export      NF_CONFIG=${NETCDF}/bin/nf-config
          export  NETCDF_INCDIR=${NETCDF}/include
          export        NETCDF4=1
        fi
      else
        export         ESMF_DIR=${MPI_SOFT}/esmf_nc3
        export           NETCDF=/opt/intelsoft/serial/netcdf3
        export    NETCDF_INCDIR=${NETCDF}/include
        export    NETCDF_LIBDIR=${NETCDF}/lib
        export   NETCDF_classic=1
      fi
    fi

    if [ -n "${USE_PIO:+1}" ]; then
      export              PIO=${MPI_SOFT}/pio
      export       PIO_LIBDIR=${PIO}/lib
      export       PIO_INCDIR=${PIO}/include
      export          PNETCDF=${MPI_SOFT}/pnetcdf
      export   PNETCDF_LIBDIR=${PNETCDF}/lib
      export   PNETCDF_INCDIR=${PNETCDF}/include
    fi

    if [ -n "${USE_SCORPIO:+1}" ]; then
      export          SCORPIO=${MPI_SOFT}/scorpio
      export       PIO_LIBDIR=${SCORPIO}/lib
      export       PIO_INCDIR=${SCORPIO}/include
      export          PNETCDF=${MPI_SOFT}/pnetcdf
      export   PNETCDF_LIBDIR=${PNETCDF}/lib
      export   PNETCDF_INCDIR=${PNETCDF}/include
    fi

    if [ -n "${USE_HDF5:+1}" ]; then
      if [ -n "${USE_PARALLEL_IO:+1}" ] && [ -n "${USE_MPI:+1}" ]; then
        export           HDF5=${MPI_SOFT}/hdf5
        export    HDF5_LIBDIR=${HDF5}/lib
        export    HDF5_INCDIR=${HDF5}/include
      else
        export           HDF5=/opt/intelsoft/serial/hdf5
        export    HDF5_LIBDIR=${HDF5}/lib
        export    HDF5_INCDIR=${HDF5}/include
      fi
    fi
    ;;

  pgi )
    export      ESMF_COMPILER=pgi
    if [ -n "${USE_DEBUG:+1}" ]; then
      export        ESMF_BOPT=g
    else
      export        ESMF_BOPT=O
    fi
    export           ESMF_ABI=64
    export          ESMF_COMM=${which_MPI}
    export          ESMF_SITE=default

    export      ARPACK_LIBDIR=/opt/pgisoft/serial/ARPACK
    if [ -n "${USE_MPI:+1}" ]; then
      if [ "${which_MPI}" = "mpich" ]; then
        export       MPI_SOFT=/opt/pgisoft/mpich
      elif [ "${which_MPI}" = "mpich2" ]; then
        export       MPI_SOFT=/opt/pgisoft/mpich2
      elif [ "${which_MPI}" = "openmpi" ]; then
        export       MPI_SOFT=/opt/pgisoft/openmpi
      elif [ "${which_MPI}" = "mvapich2" ]; then
        export       MPI_SOFT=/opt/pgisoft/mvapich2
      fi
      export       MCT_INCDIR=${MPI_SOFT}/mct/include
      export       MCT_LIBDIR=${MPI_SOFT}/mct/lib
      export   PARPACK_LIBDIR=${MPI_SOFT}/PARPACK
    fi

    if [ ! -n "${SINGULARITY_COMMAND:+1}" ]; then
      if [ -n "${USE_NETCDF4:+1}" ]; then
        if [ -n "${USE_PARALLEL_IO:+1}" ] && [ -n "${USE_MPI:+1}" ]; then
          export       ESMF_DIR=${MPI_SOFT}/esmf_nc4
          export         NETCDF=${MPI_SOFT}/netcdf4
          export      NF_CONFIG=${NETCDF}/bin/nf-config
          export  NETCDF_INCDIR=${NETCDF}/include
          export        NETCDF4=1
        else
          export       ESMF_DIR=${MPI_SOFT}/esmf_nc4
          export         NETCDF=/opt/pgisoft/serial/netcdf4
          export      NF_CONFIG=${NETCDF}/bin/nf-config
          export  NETCDF_INCDIR=${NETCDF}/include
          export        NETCDF4=1
        fi
      else
        export         ESMF_DIR=${MPI_SOFT}/esmf_nc3
        export           NETCDF=/opt/pgisoft/serial/netcdf3
        export    NETCDF_INCDIR=${NETCDF}/include
        export    NETCDF_LIBDIR=${NETCDF}/lib
        export   NETCDF_classic=1
      fi
    fi

    if [ -n "${USE_PIO:+1}" ]; then
      export              PIO=${MPI_SOFT}/pio
      export       PIO_LIBDIR=${PIO}/lib
      export       PIO_INCDIR=${PIO}/include
      export          PNETCDF=${MPI_SOFT}/pnetcdf
      export   PNETCDF_LIBDIR=${PNETCDF}/lib
      export   PNETCDF_INCDIR=${PNETCDF}/include
    fi

    if [ -n "${USE_SCORPIO:+1}" ]; then
      export          SCORPIO=${MPI_SOFT}/scorpio
      export       PIO_LIBDIR=${SCORPIO}/lib
      export       PIO_INCDIR=${SCORPIO}/include
      export          PNETCDF=${MPI_SOFT}/pnetcdf
      export   PNETCDF_LIBDIR=${PNETCDF}/lib
      export   PNETCDF_INCDIR=${PNETCDF}/include
    fi

    if [ -n "${USE_HDF5:+1}" ]; then
      if [ -n "${USE_PARALLEL_IO:+1}" ] && [ -n "${USE_MPI:+1}" ]; then
        export           HDF5=${MPI_SOFT}/hdf5
        export    HDF5_LIBDIR=${HDF5}/lib
        export    HDF5_INCDIR=${HDF5}/include
      else
        export           HDF5=/opt/pgisoft/serial/hdf5
        export    HDF5_LIBDIR=${HDF5}/lib
        export    HDF5_INCDIR=${HDF5}/include
      fi
    fi
    ;;

  gfortran )
    export      ESMF_COMPILER=gfortran
    if [ -n "${USE_DEBUG:+1}" ]; then
      export        ESMF_BOPT=g
    else
      export        ESMF_BOPT=O
    fi
    export           ESMF_ABI=64
    export          ESMF_COMM=${which_MPI}
    export          ESMF_SITE=default

    export      ARPACK_LIBDIR=/opt/gfortransoft/serial/ARPACK
    if [ -n "${USE_MPI:+1}" ]; then
      if [ "${which_MPI}" = "mpich2" ]; then
        export       MPI_SOFT=/opt/gfortransoft/mpich2
      elif [ "${which_MPI}" = "openmpi" ]; then
        export       MPI_SOFT=/opt/gfortransoft/openmpi
      elif [ "${which_MPI}" = "mvapich2" ]; then
        export       MPI_SOFT=/opt/gfortransoft/mvapich2
      fi
      export       MCT_INCDIR=${MPI_SOFT}/mct/include
      export       MCT_LIBDIR=${MPI_SOFT}/mct/lib
      export   PARPACK_LIBDIR=${MPI_SOFT}/PARPACK
    fi

    if [ ! -n "${SINGULARITY_COMMAND:+1}" ]; then
      if [ -n "${USE_NETCDF4:+1}" ]; then
        if [ -n "${USE_PARALLEL_IO:+1}" ] && [ -n "${USE_MPI:+1}" ]; then
          export       ESMF_DIR=${MPI_SOFT}/esmf_nc4
          export         NETCDF=${MPI_SOFT}/netcdf4
          export      NF_CONFIG=${NETCDF}/bin/nf-config
          export  NETCDF_INCDIR=${NETCDF}/include
          export        NETCDF4=1
        else
          export       ESMF_DIR=${MPI_SOFT}/esmf_nc4
          export         NETCDF=/opt/gfortransoft/serial/netcdf4
          export      NF_CONFIG=${NETCDF}/bin/nf-config
          export  NETCDF_INCDIR=${NETCDF}/include
          export        NETCDF4=1
        fi
      else
        export         ESMF_DIR=${MPI_SOFT}/esmf_nc3
        export           NETCDF=/opt/gfortransoft/serial/netcdf3
        export    NETCDF_INCDIR=${NETCDF}/include
        export    NETCDF_LIBDIR=${NETCDF}/lib
        export   NETCDF_classic=1
      fi
    fi

    if [ -n "${USE_PIO:+1}" ]; then
      export              PIO=${MPI_SOFT}/pio
      export       PIO_LIBDIR=${PIO}/lib
      export       PIO_INCDIR=${PIO}/include
      export          PNETCDF=${MPI_SOFT}/pnetcdf
      export   PNETCDF_LIBDIR=${PNETCDF}/lib
      export   PNETCDF_INCDIR=${PNETCDF}/include
    fi

    if [ -n "${USE_SCORPIO:+1}" ]; then
      export          SCORPIO=${MPI_SOFT}/scorpio
      export       PIO_LIBDIR=${SCORPIO}/lib
      export       PIO_INCDIR=${SCORPIO}/include
      export          PNETCDF=${MPI_SOFT}/pnetcdf
      export   PNETCDF_LIBDIR=${PNETCDF}/lib
      export   PNETCDF_INCDIR=${PNETCDF}/include
    fi

    if [ -n "${USE_HDF5:+1}" ]; then
      if [ -n "${USE_PARALLEL_IO:+1}" ] && [ -n "${USE_MPI:+1}" ]; then
        export           HDF5=${MPI_SOFT}/hdf5
        export    HDF5_LIBDIR=${HDF5}/lib
        export    HDF5_INCDIR=${HDF5}/include
      else
        export           HDF5=/opt/gfortransoft/serial/hdf5
        export    HDF5_LIBDIR=${HDF5}/lib
        export    HDF5_INCDIR=${HDF5}/include
      fi
    fi
    ;;

esac
