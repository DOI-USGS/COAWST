#!/bin/csh -f
#
# git $Id$
# svn $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# WPS Compiling CSH Script                                              :::
#                                                                       :::
# Script to compile WPS where source code files are kept separate       :::
# from the application configuration and build objects.                 :::
#                                                                       :::
# If WRF I/O libraries are needed to compile WPS, it is expected to     :::
# have the WRF model compiled one directory level up named 'WRF',       :::
# 'WRF-4.1', ..., 'WRF-4.3', and so on. Alternatevely, specify the      :::
# full path of the compiled WRF model with the environmental variable   :::
# WRF_DIR. Notice that the WRF I/O libraries is needed to create        :::
# "metgrid.exe" and "real.exe".                                         :::
#                                                                       :::
# WARNING:                                                              :::
#                                                                       :::
# If WRF is built with parallel enable NetCDF4/HDF5, you must choose    :::
# the appropriate (dmpar) option for your compiler during configure.    :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    ./build_wps.csh [options]                                          :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    -move       Move compiled objects to build directory               :::
#                                                                       :::
#    -noclean    Do not run clean -a script                             :::
#                                                                       :::
#    -noconfig   Do not run configure compilation script                :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

setenv which_MPI openmpi                        # default, overwriten below

# Initialize.

set separator = `perl -e "print ':' x 100;"`

set parallel = 0
set clean = 1
set config = 1
set move = 0

setenv CPLFLAG ''
setenv MY_CPP_FLAGS ''

while ( ($#argv) > 0 )
  switch ($1)
    case "-noclean"
      shift
      set clean = 0
    breaksw

    case "-noconfig"
      shift
      set config = 0
    breaksw

    case "-move"
      shift
      set move = 1
    breaksw

    case "-*":
      echo ""
      echo "$0 : Unknown option [ $1 ]"
      echo ""
      echo "Available Options:"
      echo ""
      echo "-move         Move compiled objects to build directory"
      echo ""
      echo "-noclean      Do not run clean script"
      echo ""
      echo "-noconfig     Do not run configure compilation script"
      echo ""
     exit 1
    breaksw

  endsw
end

#--------------------------------------------------------------------------
# Set a local environmental variable to define the path of the working
# application directory where all this project's files are kept.
#--------------------------------------------------------------------------

 setenv MY_PROJECT_DIR       ${PWD}

#--------------------------------------------------------------------------
# Set a local environmental variable to define the root path to the
# directories where WPS source and WRF library and objects files are
# located. WRF library and objects files are needed for building and
# linking WPS executables.
#--------------------------------------------------------------------------

 setenv WPS_ROOT_DIR         ${HOME}/ocean/repository/WPS
 setenv WPS_SRC_DIR          ${WPS_ROOT_DIR}

 setenv WRF_DIR              ${HOME}/ocean/repository/WRF
#setenv WRF_DIR              ${MY_PROJECT_DIR}/Build_wrf

#--------------------------------------------------------------------------
# Set Fortran compiler and MPI library to use.
#--------------------------------------------------------------------------

 setenv USE_MPI              on          # distributed-memory parallelism
 setenv USE_MPIF90           on          # compile with mpif90 script
#setenv which_MPI            intel       # compile with mpiifort library
#setenv which_MPI            mpich       # compile with MPICH library
#setenv which_MPI            mpich2      # compile with MPICH2 library
#setenv which_MPI            mvapich2    # compile with MVAPICH2 library
 setenv which_MPI            openmpi     # compile with OpenMPI library

 setenv FORT                 ifort
#setenv FORT                 gfortran
#setenv FORT                 pgi

#setenv USE_REAL_DOUBLE      on          # use real double precision (-r8)
#setenv USE_DEBUG            on          # use Fortran debugging flags
 setenv USE_NETCDF           on          # compile with NetCDF
 setenv USE_NETCDF4          on          # compile with NetCDF-4 library
                                         # (Must also set USE_NETCDF)

#--------------------------------------------------------------------------
# Use my specified library paths.
#--------------------------------------------------------------------------

 setenv USE_MY_LIBS no           # use system default library paths
#setenv USE_MY_LIBS yes          # use my customized library paths

#set MY_PATHS = ${ROMS_SRC_DIR}/Compilers/my_build_paths.csh
 set MY_PATHS = ${HOME}/Compilers/ROMS/my_build_paths.csh

if ($USE_MY_LIBS == 'yes') then
  source ${MY_PATHS} ${MY_PATHS}
endif

#--------------------------------------------------------------------------
# WPS build and executable directory.
#--------------------------------------------------------------------------

# Put the *.a, .f and .f90 files in a project specific Build directory to
# avoid conflict with other projects.

if ($?USE_DEBUG) then
  setenv  WPS_BUILD_DIR   ${MY_PROJECT_DIR}/Build_wpsG
else
  setenv  WPS_BUILD_DIR   ${MY_PROJECT_DIR}/Build_wps
endif

# Put WPS executables in the following directory.

setenv  WPS_BIN_DIR       ${WPS_BUILD_DIR}/Bin

# Go to the users source directory to compile. The options set above will
# pick up the application-specific code from the appropriate place.

cd ${WPS_ROOT_DIR}

# Patch util/src/Makefile so int2nc.exe using parallel enabled NetCDF/HDF5
# can link without errors.

echo ""
echo "${separator}"
echo "Patching ${WPS_SRC_DIR}/util/src/Makefile"
echo "${separator}"
echo ""

# Replace $(SFC) with $(FC) in recipe for int2nc.exe

perl -i -pe 's/^(\s*)\$\(SFC\) (.*\$\(WRF_LIB\)\s*)/$1\$\(FC\) $2/' ${WPS_SRC_DIR}/util/src/Makefile

#--------------------------------------------------------------------------
# Configure. It creates configure.wps script used for compilation.
#
#   -d    build with debugging information and no optimization
#   -D    build with -d AND floating traps, traceback, uninitialized variables
#   -r8   build with 8 byte reals
#
# During configuration the WPS/arch/Config.pl perl script is executed and
# we need to interactively select the combination of compiler and parallel
# comunications option. For example, for Darwin operating system we get:
#
#Please select from among the following supported platforms.
#
#   1.  Darwin Intel PGI compiler   (serial)
#   2.  Darwin Intel PGI compiler   (serial_NO_GRIB2)
#   3.  Darwin Intel PGI compiler   (dmpar)
#   4.  Darwin Intel PGI compiler   (dmpar_NO_GRIB2)
#   5.  Darwin Intel PGI compiler; optional DM -f90=pgf90   (serial)
#   6.  Darwin Intel PGI compiler; optional DM -f90=pgf90   (serial_NO_GRIB2)
#   7.  Darwin Intel PGI compiler; optional DM -f90=pgf90   (dmpar)
#   8.  Darwin Intel PGI compiler; optional DM -f90=pgf90   (dmpar_NO_GRIB2)
#   9.  Darwin Intel Intel compiler   (serial)
#  10.  Darwin Intel Intel compiler   (serial_NO_GRIB2)
#  11.  Darwin Intel Intel compiler   (dmpar)
#  12.  Darwin Intel Intel compiler   (dmpar_NO_GRIB2)
#  13.  Darwin Intel g95 compiler     (serial)
#  14.  Darwin Intel g95 compiler     (serial_NO_GRIB2)
#  15.  Darwin Intel g95 compiler     (dmpar)
#  16.  Darwin Intel g95 compiler     (dmpar_NO_GRIB2)
#  17.  Darwin Intel gfortran/gcc    (serial)
#  18.  Darwin Intel gfortran/gcc    (serial_NO_GRIB2)
#  19.  Darwin Intel gfortran/gcc    (dmpar)
#  20.  Darwin Intel gfortran/gcc    (dmpar_NO_GRIB2)
#  21.  Darwin Intel gfortran/clang    (serial)
#  22.  Darwin Intel gfortran/clang    (serial_NO_GRIB2)
#  23.  Darwin Intel gfortran/clang    (dmpar)
#  24.  Darwin Intel gfortran/clang    (dmpar_NO_GRIB2)
#  25.  Darwin PPC xlf    (serial)
#  26.  Darwin PPC xlf    (serial_NO_GRIB2)
#  27.  Darwin PPC xlf gcc3.3 SystemStubs   (serial)
#  28.  Darwin PPC xlf gcc3.3 SystemStubs   (serial_NO_GRIB2)
#  29.  Darwin PPC g95    (serial)
#  30.  Darwin PPC g95    (serial_NO_GRIB2)
#  31.  Darwin PPC g95    (dmpar)
#  32.  Darwin PPC g95    (dmpar_NO_GRIB2)
#
#Enter selection [1-32] :
#
#--------------------------------------------------------------------------

# Clean source code.

if ( $clean == 1 ) then
  echo ""
  echo "${separator}"
  echo "Cleaning WPS source code:  ${WPS_ROOT_DIR}/clean -a"
  echo "${separator}"
  echo ""
  ${WPS_ROOT_DIR}/clean -a            # clean source code
endif

if ($?USE_DEBUG) then
   set DEBUG_FLAG = -d
#  set DEBUG_FLAG = -D
endif

setenv CONFIG_FLAGS ''

if ( $config == 1 ) then
  if ($?USE_DEBUG && $?USE_REAL_DOUBLE) then
    setenv CONFIG_FLAGS "${CONFIG_FLAGS} ${DEBUG_FLAG} -r8"
  else if ($?USE_DEBUG) then
    setenv CONFIG_FLAGS "${CONFIG_FLAGS} ${DEBUG_FLAG}"
  else if ($?USE_REAL_DOUBLE) then
    setenv CONFIG_FLAGS "${CONFIG_FLAGS} -r8"
  endif

  echo ""
  echo "${separator}"
  echo "Configuring WPS code:  ${WPS_ROOT_DIR}/configure ${CONFIG_FLAGS}"
  echo "${separator}"
  echo ""

  ${WPS_ROOT_DIR}/configure ${CONFIG_FLAGS}

# Add the NetCDF-4 dependecies to WRF_LIB.

  if ( -f ${NETCDF}/bin/nf-config ) then
    set NF_LIBS=`${NETCDF}/bin/nf-config --flibs`
    setenv NF_LIBS "${NF_LIBS}"

    echo ""
    echo "${separator}"
    echo "Appending WRF_LIB += ${NF_LIBS}"
    echo "to the end of configure.wps"
    echo ""
    echo "WRF_LIB += ${NF_LIBS}" >> ${WPS_SRC_DIR}/configure.wps
  else
    echo "${NETCDF}/bin/nf-config not found. Please check your configuration."
    exit 1
  endif

# If which_MPI is "intel" then we need to replace DM_FC and DM_CC in configure.wps

  if ( ${which_MPI} == "intel" ) then
    perl -i -pe 's/^DM_FC(\s*)=(\s*)mpif90/DM_FC$1=$2mpiifort/' ${WPS_SRC_DIR}/configure.wps
    perl -i -pe 's/^DM_CC(\s*)=(\s*)mpicc/DM_CC$1=$2mpiicc/' ${WPS_SRC_DIR}/configure.wps
  endif
endif

#--------------------------------------------------------------------------
# WPS Compile script:
#
# Usage:
#
#     compile [target]
#
#   where target is one of
#
#      wps
#      util
#      geogrid
#      ungrib
#      metgrid
#      g1print
#      g2print
#      plotfmt
#      rd_intermediate
#      plotgrids
#      mod_levs
#      avg_tsfc
#      calc_ecmwf_p
#      height_ukmo
#      int2nc
#
#    or just run compile with no target to build everything.
#--------------------------------------------------------------------------

# Remove existing build directory.

if ( $move == 1 ) then
  if ( $clean == 1 ) then
    echo ""
    echo "${separator}"
    echo "Removing WPS build directory: ${WPS_BUILD_DIR}"
    echo "${separator}"
    echo ""
    /bin/rm -rf ${WPS_BUILD_DIR}
  endif
endif

# Compile (if -move is set, the binares will go to WPS_BIN_DIR set above).

echo ""
echo "${separator}"
echo "Compiling WPS using  ${MY_PROJECT_DIR}/${0}:"
echo ""
echo "   ${WPS_ROOT_DIR}/compile"
echo "${separator}"
echo ""

${WPS_ROOT_DIR}/compile

#--------------------------------------------------------------------------
# Move WPS objects and executables.
#--------------------------------------------------------------------------

if ( $move == 1 ) then

  echo ""
  echo "${separator}"
  echo "Moving WPS configuration and executables to Build directory  ${WPS_BUILD_DIR}:"
  echo "${separator}"
  echo ""

  if ( ! -d ${WPS_BUILD_DIR} ) then
    /bin/mkdir -pv ${WPS_BUILD_DIR}
    /bin/mkdir -pv ${WPS_BIN_DIR}
    echo ""
  endif

  /bin/cp -pfv configure.wps ${WPS_BUILD_DIR}

  echo "WPS_ROOT_DIR = ${WPS_ROOT_DIR}"
  find -H ${WPS_ROOT_DIR} -type f -name "*.exe" -exec /bin/mv -fv {} ${WPS_BIN_DIR} \;

  echo ""
  echo "${separator}"
  echo "Removing symbolic links to executables in  ${WPS_ROOT_DIR}:"
  echo "${separator}"
  echo ""

  find -H ${WPS_ROOT_DIR} -type l -name "*.exe" -exec /bin/rm -fv {} \;

endif
