#!/bin/csh -f
#
# git $Id$
# svn $Id: my_build_paths.csh 1151 2023-02-09 03:08:53Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# ROMS Customized Compiling Libraries CSH Script                        :::
#                                                                       :::
# This C-shell script sets the customized library paths needed by the   :::
# build script when the enviromental variable USE_MY_LIBS has a 'yes'   :::
# value.                                                                :::
#                                                                       :::
# For example, in build_roms.csh we have:                               :::
#                                                                       :::
#       if ($USE_MY_LIBS == 'yes') then                                 :::
#          source ${COMPILERS}/my_build_paths.csh                       :::
#       endif                                                           :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set separator = `perl -e "print ':' x 100;"`

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
#
# Some models like COAMPS need to know where the include file
# "mpif.h" is located.

setenv MPI_ROOT ""

if ($?USE_MPIF90) then
  switch ($FORT)

    case "ifort"
      if ($which_MPI == "mpich" ) then
        setenv MPI_ROOT /opt/intelsoft/mpich
      else if ($which_MPI == "mpich2" ) then
        setenv MPI_ROOT /opt/intelsoft/mpich2
      else if ($which_MPI == "openmpi" ) then
        setenv MPI_ROOT /opt/intelsoft/openmpi
      else if ($which_MPI == "mvapich2" ) then
        setenv MPI_ROOT /opt/intelsoft/mvapich2
      endif
      setenv PATH       ${MPI_ROOT}/bin:$PATH
      setenv MPI_INCDIR ${MPI_ROOT}/include
    breaksw

    case "pgi"
      if ($which_MPI == "mpich" ) then
        setenv MPI_ROOT /opt/pgisoft/mpich
      else if ($which_MPI == "mpich2" ) then
        setenv MPI_ROOT /opt/pgisoft/mpich2
      else if ($which_MPI == "openmpi" ) then
        setenv MPI_ROOT /opt/pgisoft/openmpi
      else if ($which_MPI == "mvapich2" ) then
        setenv MPI_ROOT /opt/pgisoft/mvapich2
      endif
      setenv PATH       ${MPI_ROOT}/bin:$PATH
      setenv MPI_INCDIR ${MPI_ROOT}/include
    breaksw

    case "gfortran"
      if ($which_MPI == "mpich2" ) then
        setenv MPI_ROOT /opt/gfortransoft/mpich2
      else if ($which_MPI == "openmpi" ) then
        setenv MPI_ROOT /opt/gfortransoft/openmpi
      else if ($which_MPI == "mvapich2" ) then
        setenv MPI_ROOT /opt/gfortransoft/mvapich2
      endif
      setenv PATH       ${MPI_ROOT}/bin:$PATH
      setenv MPI_INCDIR ${MPI_ROOT}/include
    breaksw

  endsw
endif

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

setenv MPI_SOFT ""

switch ($FORT)

# Intel Compiler:

  case "ifort"
    setenv ESMF_COMPILER        intelgcc
#   setenv ESMF_COMPILER        intel
    if ($?USE_DEBUG) then
      setenv ESMF_BOPT          g
    else
      setenv ESMF_BOPT          O
    endif
    setenv ESMF_ABI             64
    setenv ESMF_COMM            ${which_MPI}
    setenv ESMF_SITE            default

    setenv ARPACK_LIBDIR        /opt/intelsoft/serial/ARPACK
    if ($?USE_MPI) then
      if ($which_MPI == "mpich" ) then
        setenv MPI_SOFT         /opt/intelsoft/mpich
      else if ($which_MPI == "mpich2" ) then
        setenv MPI_SOFT         /opt/intelsoft/mpich2
      else if ($which_MPI == "openmpi" ) then
        setenv MPI_SOFT         /opt/intelsoft/openmpi
      else if ($which_MPI == "mvapich2" ) then
        setenv MPI_SOFT         /opt/intelsoft/mvapich2
      endif
      setenv MCT_INCDIR         ${MPI_SOFT}/mct/include
      setenv MCT_LIBDIR         ${MPI_SOFT}/mct/lib
      setenv PARPACK_LIBDIR     ${MPI_SOFT}/PARPACK
    endif

    if (! $?SINGULARITY_COMMAND) then
        if ($?USE_NETCDF4) then
          if ($?USE_PARALLEL_IO && $?USE_MPI) then
              setenv ESMF_DIR       ${MPI_SOFT}/esmf_nc4
              setenv NETCDF         ${MPI_SOFT}/netcdf4
              setenv NF_CONFIG      ${NETCDF}/bin/nf-config
              setenv NETCDF_INCDIR  ${NETCDF}/include
              setenv NETCDF4        1
          else
            setenv ESMF_DIR         ${MPI_SOFT}/esmf_nc4
            setenv NETCDF           /opt/intelsoft/serial/netcdf4
            setenv NF_CONFIG        ${NETCDF}/bin/nf-config
            setenv NETCDF_INCDIR    ${NETCDF}/include
            setenv NETCDF4          1
          endif
        else
          setenv ESMF_DIR           ${MPI_SOFT}/esmf_nc3
          setenv NETCDF             /opt/intelsoft/serial/netcdf3
          setenv NETCDF_INCDIR      ${NETCDF}/include
          setenv NETCDF_LIBDIR      ${NETCDF}/lib
          setenv NETCDF_classic     1
        endif
      endif
    endif

    if ($?USE_PIO) then
      setenv PIO                ${MPI_SOFT}/pio
      setenv PIO_LIBDIR         ${PIO}/lib
      setenv PIO_INCDIR         ${PIO}/include
      setenv PNETCDF            ${MPI_SOFT}/pnetcdf
      setenv PNETCDF_LIBDIR     ${PNETCDF}/lib
      setenv PNETCDF_INCDIR     ${PNETCDF}/include
    endif

    if ($?USE_SCORPIO) then
      setenv SCORPIO            ${MPI_SOFT}/scorpio
      setenv PIO_LIBDIR         ${SCORPIO}/lib
      setenv PIO_INCDIR         ${SCORPIO}/include
      setenv PNETCDF            ${MPI_SOFT}/pnetcdf
      setenv PNETCDF_LIBDIR     ${PNETCDF}/lib
      setenv PNETCDF_INCDIR     ${PNETCDF}/include
    endif

    if ($?USE_HDF5) then
      if ($?USE_PARALLEL_IO && $?USE_MPI) then
        setenv HDF5             ${MPI_SOFT}/hdf5
        setenv HDF5_LIBDIR      ${HDF5}/lib
        setenv HDF5_INCDIR      ${HDF5}/include
      else
        setenv HDF5             /opt/intelsoft/serial/hdf5
        setenv HDF5_LIBDIR      ${HDF5}/lib
        setenv HDF5_INCDIR      ${HDF5}/include
      endif
    endif

  breaksw

# PGI Compiler:

  case "pgi"
    setenv ESMF_COMPILER        pgi
    if ($?USE_DEBUG) then
      setenv ESMF_BOPT          g
    else
      setenv ESMF_BOPT          O
    endif
    setenv ESMF_ABI             64
    setenv ESMF_COMM            ${which_MPI}
    setenv ESMF_SITE            default

    setenv ARPACK_LIBDIR        /opt/pgisoft/serial/ARPACK
    if ($?USE_MPI) then
      if ($which_MPI == "mpich" ) then
        setenv MPI_SOFT         /opt/pgisoft/mpich
      else if ($which_MPI == "mpich2" ) then
        setenv MPI_SOFT         /opt/pgisoft/mpich2
      else if ($which_MPI == "openmpi" ) then
        setenv MPI_SOFT         /opt/pgisoft/openmpi
      else if ($which_MPI == "mvapich2" ) then
        setenv MPI_SOFT         /opt/pgisoft/mvapich2
      endif
      setenv MCT_INCDIR         ${MPI_SOFT}/mct/include
      setenv MCT_LIBDIR         ${MPI_SOFT}/mct/lib
      setenv PARPACK_LIBDIR     ${MPI_SOFT}/PARPACK
    endif

    if (! $?SINGULARITY_COMMAND) then
      if ($?USE_NETCDF4) then
        if ($?USE_PARALLEL_IO && $?USE_MPI) then
            setenv ESMF_DIR       ${MPI_SOFT}/esmf_nc4
            setenv NETCDF         ${MPI_SOFT}/netcdf4
            setenv NF_CONFIG      ${NETCDF}/bin/nf-config
            setenv NETCDF_INCDIR  ${NETCDF}/include
            setenv NETCDF4        1
        else
          setenv ESMF_DIR         ${MPI_SOFT}/esmf_nc4
          setenv NETCDF           /opt/pgisoft/serial/netcdf4
          setenv NF_CONFIG        ${NETCDF}/bin/nf-config
          setenv NETCDF_INCDIR    ${NETCDF}/include
          setenv NETCDF4          1
        endif
      else
        setenv ESMF_DIR           ${MPI_SOFT}/esmf_nc3
        setenv NETCDF             /opt/pgisoft/serial/netcdf3
        setenv NETCDF_INCDIR      ${NETCDF}/include
        setenv NETCDF_LIBDIR      ${NETCDF}/lib
        setenv NETCDF_classic     1
      endif
    endif

    if ($?USE_PIO) then
      setenv PIO                ${MPI_SOFT}/pio
      setenv PIO_LIBDIR         ${PIO}/lib
      setenv PIO_INCDIR         ${PIO}/include
      setenv PNETCDF            ${MPI_SOFT}/pnetcdf
      setenv PNETCDF_LIBDIR     ${PNETCDF}/lib
      setenv PNETCDF_INCDIR     ${PNETCDF}/include
    endif

    if ($?USE_SCORPIO) then
      setenv SCORPIO            ${MPI_SOFT}/scorpio
      setenv PIO_LIBDIR         ${SCORPIO}/lib
      setenv PIO_INCDIR         ${SCORPIO}/include
      setenv PNETCDF            ${MPI_SOFT}/pnetcdf
      setenv PNETCDF_LIBDIR     ${PNETCDF}/lib
      setenv PNETCDF_INCDIR     ${PNETCDF}/include
    endif

    if ($?USE_HDF5) then
      if ($?USE_PARALLEL_IO && $?USE_MPI) then
        setenv HDF5             ${MPI_SOFT}/hdf5
        setenv HDF5_LIBDIR      ${HDF5}/lib
        setenv HDF5_INCDIR      ${HDF5}/include
      else
        setenv HDF5             /opt/pgisoft/serial/hdf5
        setenv HDF5_LIBDIR      ${HDF5}/lib
        setenv HDF5_INCDIR      ${HDF5}/include
      endif
    endif

  breaksw

# GNU Compiler:

  case "gfortran"
    setenv ESMF_COMPILER        gfortran
    if ($?USE_DEBUG) then
      setenv ESMF_BOPT          g
    else
      setenv ESMF_BOPT          O
    endif
    setenv ESMF_ABI             64
    setenv ESMF_COMM            ${which_MPI}
    setenv ESMF_SITE            default

    setenv ARPACK_LIBDIR        /opt/gfortransoft/serial/ARPACK
    if ($?USE_MPI) then
      if ($which_MPI == "mpich2" ) then
        setenv MPI_SOFT         /opt/gfortransoft/mpich2
      else if ($which_MPI == "openmpi" ) then
        setenv MPI_SOFT         /opt/gfortransoft/openmpi
      else if ($which_MPI == "mvapich2" ) then
        setenv MPI_SOFT         /opt/gfortransoft/mvapich2
      endif
      setenv MCT_INCDIR         ${MPI_SOFT}/mct/include
      setenv MCT_LIBDIR         ${MPI_SOFT}/mct/lib
      setenv PARPACK_LIBDIR     ${MPI_SOFT}/PARPACK
    endif

    if (! $?SINGULARITY_COMMAND) then
      if ($?USE_NETCDF4) then
        if ($?USE_PARALLEL_IO && $?USE_MPI) then
            setenv ESMF_DIR       ${MPI_SOFT}/esmf_nc4
            setenv NETCDF         ${MPI_SOFT}/netcdf4
            setenv NF_CONFIG      ${NETCDF}/bin/nf-config
            setenv NETCDF_INCDIR  ${NETCDF}/include
            setenv NETCDF4        1
        else
          setenv ESMF_DIR         ${MPI_SOFT}/esmf_nc4
          setenv NETCDF           /opt/gfortransoft/serial/netcdf4
          setenv NF_CONFIG        ${NETCDF}/bin/nf-config
          setenv NETCDF_INCDIR    ${NETCDF}/include
          setenv NETCDF4          1
        endif
      else
        setenv ESMF_DIR           ${MPI_SOFT}/esmf_nc3
        setenv NETCDF             /opt/gfortransoft/serial/netcdf3
        setenv NETCDF_INCDIR      ${NETCDF}/include
        setenv NETCDF_LIBDIR      ${NETCDF}/lib
        setenv NETCDF_classic     1
      endif
    endif

    if ($?USE_PIO) then
      setenv PIO                ${MPI_SOFT}/pio
      setenv PIO_LIBDIR         ${PIO}/lib
      setenv PIO_INCDIR         ${PIO}/include
      setenv PNETCDF            ${MPI_SOFT}/pnetcdf
      setenv PNETCDF_LIBDIR     ${PNETCDF}/lib
      setenv PNETCDF_INCDIR     ${PNETCDF}/include
    endif

    if ($?USE_SCORPIO) then
      setenv SCORPIO            ${MPI_SOFT}/scorpio
      setenv PIO_LIBDIR         ${SCORPIO}/lib
      setenv PIO_INCDIR         ${SCORPIO}/include
      setenv PNETCDF            ${MPI_SOFT}/pnetcdf
      setenv PNETCDF_LIBDIR     ${PNETCDF}/lib
      setenv PNETCDF_INCDIR     ${PNETCDF}/include
    endif

    if ($?USE_HDF5) then
      if ($?USE_PARALLEL_IO && $?USE_MPI) then
        setenv HDF5             ${MPI_SOFT}/hdf5
        setenv HDF5_LIBDIR      ${HDF5}/lib
        setenv HDF5_INCDIR      ${HDF5}/include
      else
        setenv HDF5             /opt/gfortransoft/serial/hdf5
        setenv HDF5_LIBDIR      ${HDF5}/lib
        setenv HDF5_INCDIR      ${HDF5}/include
      endif
    endif

  breaksw

endsw
