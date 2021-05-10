#!/bin/csh -f
#
# svn $Id: build_cice.csh 1054 2021-03-06 19:47:12Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2021 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# CICE Compiling CSH Script                                             :::
#                                                                       :::
# Script to compile an user application where the application-specific  :::
# files are kept separate from the CICE source code.                    :::
#                                                                       :::
# Q: How/why does this script work?                                     :::
#                                                                       :::
# A: The CICE Makefile configures user-defined options with a set of    :::
#    flags such as CPPDEFS. Browse the CICE Makefile to see these.      :::
#                                                                       :::
#    If an option in the Makefile uses the syntax ?= in setting the     :::
#    default, this means that make will check whether an environment    :::
#    variable by that name is set in the shell that calls make. If so   :::
#    the environment variable value overrides the default (and the      :::
#    user need not maintain separate Makefiles, or frequently edit      :::
#    the Makefile, to run separate applications).                       :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    ./build_cice.csh [options]                                         :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    -j [N]      Compile in parallel using N CPUs                       :::
#                  omit argument for all available CPUs                 :::
#                                                                       :::
#    -p macro    Prints any Makefile macro value. For example,          :::
#                                                                       :::
#                  build_cice.csh -p CPPDEFS                            :::
#                                                                       :::
#    -noclean    Do not clean already compiled objects                  :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

setenv which_MPI  openmpi                     #  default, overwritten below

set parallel = 0
set clean = 1
set dprint = 0

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
      echo "              For example:  build_cice.csh -p CPPDEFS"
      echo ""
      echo "-noclean    Do not clean already compiled objects"
      echo ""
      exit 1
    breaksw

  endsw
end

#--------------------------------------------------------------------------
# CICE Application configuration.
#--------------------------------------------------------------------------

# Set the CPP option defining the particular application. This will
# determine the name of the ".h" header file with the application
# CPP definitions.

 setenv CICE_APPLICATION      BERING

# Set application grid size and decode values in X- and Y-directions.

 setenv GRID                  237x277
 setenv NXGLOB               `echo $GRID | sed s/x.\*//`
 setenv NYGLOB               `echo $GRID | sed s/.\*x//`

# Set number of number of processors (NTASK) and number of points
# per domain partition (blocks/tiles) in each X-direction (BLCKX)
# and Y-direction (BLCKY). Here, NTASK is equal to "nprocs" in
# namelist file "ice_in". If BLCKX (BLCKY) does not divide NXGLOB
# (NYGLOB) evenly, padding will be used on the right (top) of the
# grid. The grid partitions BLCKX and BLCKY does not include the
# ghost points.

 setenv PARTITION '1x4'

 if ($PARTITION == '1x1') then
   setenv NTASK 1
   setenv BLCKX $NXGLOB
   setenv BLCKY $NYGLOB
 else if ($PARTITION == '1x2') then
   setenv NTASK 2
   setenv BLCKX 237
   setenv BLCKY 139
 else if ($PARTITION == '1x4') then
   setenv NTASK 4
   setenv BLCKX 237
   setenv BLCKY 70
 else if ($PARTITION == '2x2') then
   setenv NTASK 4
   setenv BLCKX 120
   setenv BLCKY 139
 else if ($PARTITION == '4x2') then
   setenv NTASK 8
   setenv BLCKX 120
   setenv BLCKY 70
 else if ($PARTITION == '2x6') then
   setenv NTASK 12
   setenv BLCKX 120
   setenv BLCKY 47
 else if ($PARTITION == '2x8') then
   setenv NTASK 16
   setenv BLCKX 120
   setenv BLCKY 35
 else if ($PARTITION == '4x16') then
   setenv NTASK 64
   setenv BLCKX 60
   setenv BLCKY 18
 endif

# May need to increase MXBLCKS with rake distribution or padding.

 @ a = $NXGLOB * $NYGLOB
 @ b = $BLCKX * $BLCKY * $NTASK
 @ m = $a / $b

 setenv MXBLCKS $m
 if ($MXBLCKS == 0) then
   setenv MXBLCKS 1
 endif

# Set number of vertical layers and thickness categories.

 setenv NICELYR    7       # number of vertical layers in the ice
 setenv NSNWLYR    1       # number of vertical layers in the snow
 setenv NICECAT    5       # number of ice thickness categories

# Set tracers parameters. To conserve memory match values with namelist
# "tracer_nml" in file "ice_in".

 setenv TRAGE   1          # set to 1 for ice age tracer
 setenv TRFY    1          # set to 1 for first-year ice area tracer
 setenv TRLVL   1          # set to 1 for level and deformed ice tracers
 setenv TRPND   1          # set to 1 for melt pond tracers
 setenv NTRAERO 0          # number of aerosol tracers
                           # (up to max_aero in ice_domain_size.F90)
                           # CESM uses 3 aerosol tracers
 setenv TRBRI   0          # set to 1 for brine height tracer
 setenv NBGCLYR 7          # number of zbgc layers
 setenv TRBGCS  0          # number of skeletal layer bgc tracers
                           # TRBGCS=0 or 2<=TRBGCS<=9)

# Set specialty code.

 setenv CAM_ICE  no        # set to yes for CAM runs (single column)
 setenv SHRDIR   csm_share # location of CCSM shared code
 setenv IO_TYPE  netcdf    # set to none if netcdf library is unavailable
                           # set to pio for parallel netcdf
 setenv DITTO    no        # reproducible diagnostics
 setenv BARRIERS no        # prevent MPI buffer overflow during gather/scatter

# Set share-memory application.

 setenv THRD     no        # set to yes for OpenMP threading

 if ( $THRD == 'yes') then
   setenv OMP_NUM_THREADS $NTASK
 endif

# Set file unit numbers.

 setenv NUMIN 11           # minimum file unit number
 setenv NUMAX 99           # maximum file unit number

#--------------------------------------------------------------------------
# Set source directory path and compilation file.
#--------------------------------------------------------------------------

# Set a local environmental variable to define the path to the directories
# where all this project's files are kept.

 setenv CICE_ROOT_DIR         ${HOME}/ocean/repository
 setenv ROMS_ROOT_DIR         ${HOME}/ocean/repository

# The path to the user's CICE and ROMS source code.
#
# If using svn locally, this would be the user's Working Copy Path (WCPATH).
# Note that one advantage of maintaining your source code locally with svn
# is that when working simultaneously on multiple machines (e.g. a local
# workstation, a local cluster and a remote supercomputer) you can checkout
# the latest release and always get an up-to-date customized source on each
# machine. This script is designed to more easily allow for differing paths
# to the code and inputs on differing machines.

 setenv MY_CICE_SRC          ${CICE_ROOT_DIR}/cice
 setenv MY_ROMS_SRC          ${ROMS_ROOT_DIR}/coupling

# Set CICE sub-directories to compile.

 setenv DRVDIR cice

 if ($IO_TYPE == 'netcdf') then
   setenv IODIR io_netcdf
 else if ($IO_TYPE == 'pio') then
   setenv IODIR io_pio
 else
   setenv IODIR io_binary
 endif

 if ($NTASK == 1) then
   setenv COMMDIR serial
 else
   setenv COMMDIR mpi
 endif

 setenv SRCDIR  ${MY_CICE_SRC}
 setenv EXEDIR  ${PWD}
 setenv OBJDIR  ${PWD}/Build_cice

 setenv VPATH  "${MY_CICE_SRC}/drivers/${DRVDIR}"
 setenv VPATH  "${VPATH} ${MY_CICE_SRC}/source"
 setenv VPATH  "${VPATH} ${MY_CICE_SRC}/$COMMDIR"
 setenv VPATH  "${VPATH} ${MY_CICE_SRC}/$IODIR"
 setenv VPATH  "${VPATH} ${MY_CICE_SRC}/$SHRDIR"

 setenv SCRATCH_DIR $OBJDIR           # ROMS make macros compatibility

#--------------------------------------------------------------------------
# Set CPP options passed to Makefile.
#--------------------------------------------------------------------------

 setenv MY_CPP_FLAGS "-D${CICE_APPLICATION}"
 setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DNXGLOB=${NXGLOB} -DNYGLOB=${NYGLOB}"
 setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DBLCKX=${BLCKX} -DBLCKY=${BLCKY}"
 setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DMXBLCKS=${MXBLCKS}"
 setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DNICELYR=${NICELYR}"
 setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DNSNWLYR=${NSNWLYR}"
 setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DNICECAT=${NICECAT}"
 setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DTRAGE=${TRAGE}"
 setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DTRFY=${TRFY}"
 setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DTRLVL=${TRLVL}"
 setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DTRPND=${TRPND}"
 setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DTRBRI=${TRBRI}"
 setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DNTRAERO=${NTRAERO}"
 setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DNBGCLYR=${NBGCLYR}"
 setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DTRBGCS=${TRBGCS}"
 setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DNUMIN=${NUMIN}"
 setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DNUMAX=${NUMAX}"

 if ($DITTO == 'yes') then
   setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DREPRODUCIBLE"
 endif

 if ($IO_TYPE == 'netcdf') then
   setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -Dncdf"
 endif

#--------------------------------------------------------------------------
# Compiler options.
#--------------------------------------------------------------------------
#
# Set path of the directory containing makefile configuration (*.mk) files.
# The user has the option to specify a customized version of these files
# in a different directory than the one distributed with the source code,
# ${MY_ROMS_SRC}/Compilers. If this is the case, the you need to keep
# these configurations files up-to-date.

 setenv COMPILERS            ${MY_ROMS_SRC}/Compilers
#setenv COMPILERS            ${HOME}/Compilers

# Set make dependencies perl script.

 setenv SFMAKEDEPEND         ${MY_ROMS_SRC}/ROMS/Bin/sfmakedepend

# Other user defined environmental variables. See the ROMS makefile for
# details on other options the user might want to set here. Be sure to
# leave the switches meant to be off set to an empty string or commented
# out. Any string value (including off) will evaluate to TRUE in
# conditional if-statements.

 setenv USE_MPI             on          # distributed-memory parallelism
 setenv USE_MPIF90          on          # compile with mpif90 script
#setenv which_MPI           mpich       # compile with MPICH library
#setenv which_MPI           mpich2      # compile with MPICH2 library
 setenv which_MPI           openmpi     # compile with OpenMPI library

#setenv USE_OpenMP          on          # shared-memory parallelism

#setenv FORT                ifort
#setenv FORT                gfortran
 setenv FORT                pgi

#setenv USE_DEBUG           on          # use Fortran debugging flags
 setenv USE_LARGE           on          # activate 64-bit compilation
#setenv USE_NETCDF4         on          # compile with NetCDF-4 library
#setenv USE_PARALLEL_IO     on          # Parallel I/O with NetCDF-4/HDF5

# Stop if activating both MPI and OpenMP at the same time.

 if ( ${?USE_MPI} & ${?USE_OpenMP} ) then
   echo "You cannot activate USE_MPI and USE_OpenMP at the same time!"
   exit 1
 endif

#--------------------------------------------------------------------------
# Use my specified library paths.
#--------------------------------------------------------------------------

 setenv USE_MY_LIBS no           # use system default library paths
#setenv USE_MY_LIBS yes          # use my customized library paths

if ($USE_MY_LIBS == 'yes') then
  source ${COMPILERS}/my_build_paths.csh
endif

#--------------------------------------------------------------------------
# Compile.
#--------------------------------------------------------------------------

# Set CICE Makefile. Use the one in ROMS Repository.

setenv CICE_MAKEFILE $MY_ROMS_SRC/Compilers/Makefile.cice

# Execute Makefile.

if ( $clean == 1 ) then
  make clean -f $CICE_MAKEFILE
endif

if ( $dprint == 1 ) then
  make $debug -f $CICE_MAKEFILE
else
 if ( $parallel == 1 ) then
   make $NCPUS -f $CICE_MAKEFILE
 else
   make -f $CICE_MAKEFILE
 endif

 echo ""
 echo "Application: $CICE_APPLICATION, Grid: $GRID, Partition: $PARTITION"
 echo ""
 echo "Number of global X-points,              NXGLOB  = $NXGLOB"
 echo "Number of global Y-points,              NYGLOB  = $NYGLOB"
 echo "Number of processors,                   NTASK   = $NTASK"
 echo "Number of X-points per tile,            BLCKX   = $BLCKX"
 echo "Number of Y-points per tile,            BLCKY   = $BLCKY"
 echo "Maximum number of blocks per processor, MXBLCKS = $MXBLCKS"
 echo "Number of vertical layers in the ice,   NICELYR = $NICELYR"
 echo "Number of vertical layers in the snow,  NSNWLYR = $NSNWLYR"
 echo "Number of ice thickness categories,     NICECAT = $NICECAT"
 echo "Ice age tracer switch,                  TRAGE   = $TRAGE"
 echo "First-year ice tracer switch,           TRFY    = $TRFY"
 echo "Level ice tracers switch,               TRLVL   = $TRLVL"
 echo "Melt pond tracers switch,               TRPND   = $TRPND"
 echo "Number of aerosol tracers,              NTRAERO = $NTRAERO"
 echo "Brine height tracer switch,             TRBRI   = $TRBRI"
 echo "Number of biogeochemical grid layers,   NBGCLYR = $NBGCLYR"
 echo "Number of biogeochemical tracers,       TRBGCS  = $TRBGCS"
 echo "Minimum file unit number,               NUMIN   = $NUMIN"
 echo "Maximum file unit number,               NUMAX   = $NUMAX"
 echo ""
 echo "MY_BUILD = $PWD/build_cice.csh"
 echo "MAKEFILE = $CICE_MAKEFILE"
 echo "FORT     = $FORT"
 echo "MDEPEND  = $SFMAKEDEPEND"
 echo "MPI      = $which_MPI"
 echo "CMACROS  = $COMPILERS"
 echo "EXEDIR   = $EXEDIR"
 echo "OBJDIR   = $OBJDIR"
 echo "SRCDIR   = $SRCDIR"
 echo ""
 echo "VPATH    = $VPATH"
 echo ""
 echo "CPPDEFS  = $MY_CPP_FLAGS"
 echo ""
endif
