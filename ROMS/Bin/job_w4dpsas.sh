#!/bin/csh -f
#
# svn $Id: job_w4dpsas.sh 503 2008-01-10 00:11:51Z arango $
#######################################################################
# Copyright (c) 2002-2008 The ROMS/TOMS Group                         #
#   Licensed under a MIT/X style license                              #
#   See License_ROMS.txt                                              #
#######################################################################
#                                                                     #
# Weak constraint, PSAS job script:                                   #
#                                                                     #
# This script NEEDS to be run before any run:                         #
#                                                                     #
#   (1) It copies a new clean nonlinear model initial conditions      #
#       file. The nonlinear model is initialized from the             #
#       background or reference state.                                #
#   (2) It copies a new clean tangent linear model initial conditions #
#       file. The tangent linear model is started from rest. That is  #
#       all the state variables are zero.                             #
#   (3) Set model error covariance standard deviation and             #
#       normalization coefficients NetCDF files.                      #
#   (4) Copy a clean copy of the observations NetCDF file.            #
#   (5) Create 4DVAR input script from a template and specify the     #
#       the background-error standard deviation, normalization, and   #
#       observation files to be used.                                 #
#                                                                     #
#######################################################################

# Set working directory root.

 set MYROOT=${MYHOME}/ocean/toms/adjoint/Test/Atoy

# Set application prefix.

 set PREFIX="atoy"

# Set string manipulations perl script.

 set SUBSTITUTE=${MYHOME}/ocean/toms/adjoint/src/ROMS/Bin/substitute

# Set ROMS data assimilation standard input scripts.

 set DA_TEMPLATE="da_ocean.tmp"

 set DA_STDINP="w4dpsas.in"

# Copy nonlinear model initial conditions file, use background or
# first guess state.

 cp -p ${MYROOT}/Data/${PREFIX}_bck.nc ${PREFIX}_ini.nc

# Copy tangent linear model initial conditions file.

 cp -p ${MYROOT}/Data/${PREFIX}_ini_zero.nc ${PREFIX}_itl.nc

# Set model-error covariance standard deviations file.

 set STDname=${MYROOT}/Data/${PREFIX}_std.nc

# Set background-error covariance normalization factor file

 set NRMname=${MYROOT}/Data/${PREFIX}_nrm.nc

# Set observations file.

 set OBSname=${PREFIX}_obs.nc

# Get a clean copy of the observation file.  This is really
# important since this file is modified to compute the
# fractional vertical position of the observations when
# they are specified as depth in meter (negative values).

 cp -p ${MYROOT}/OBS/$OBSname .

# Build data assimilation standard input script, specify above files.

 if (-e $DA_STDINP) then
   /bin/rm $DA_STDINP
 endif
 cp $DA_TEMPLATE $DA_STDINP

 $SUBSTITUTE $DA_STDINP ocean_std.nc $STDname
 $SUBSTITUTE $DA_STDINP ocean_nrm.nc $NRMname
 $SUBSTITUTE $DA_STDINP ocean_obs.nc $OBSname
 $SUBSTITUTE $DA_STDINP ocean_mod.nc ${PREFIX}_mod.nc
