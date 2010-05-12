#!/bin/csh -f
#
# svn $Id: job_is4dvar.sh 429 2009-12-20 17:30:26Z arango $
#######################################################################
# Copyright (c) 2002-2010 The ROMS/TOMS Group                         #
#   Licensed under a MIT/X style license                              #
#   See License_ROMS.txt                                              #
#######################################################################
#                                                                     #
# Incremental strong constraint, IS4DVAR job script:                  #
#                                                                     #
# This script NEEDS to be run before any run:                         #
#                                                                     #
#   (1) It copies a new clean nonlinear model initial conditions      #
#       file. The nonlinear model is initialized from the             #
#       background or reference state.                                #
#   (2) Specify initial conditions, boundary conditions, and surface  #
#       forcing error convariance input standard deviations files.    #
#   (3) Specify initial conditions, boundary conditions, and surface  #
#       forcing error convariance input/output normalization factors  #
#       files.                                                        #
#   (4) Copy a clean copy of the observations NetCDF file.            #
#   (5) Create 4DVAR input script from a template and specify the     #
#       the error covariance standard deviation, error covariance     #
#       normalization factors, and observation files to be used.      #
#                                                                     #
#######################################################################

# Set working directory root.

 set MYROOT=${MYHOME}/ocean/toms/adjoint/Test/gyre3d

# Set application prefix.

 set PREFIX="gyre3d"

# Set string manipulations perl script.

 set SUBSTITUTE=${MYHOME}/ocean/toms/adjoint/src/ROMS/Bin/substitute

# Set ROMS data assimilation standard input scripts.

 set DA_TEMPLATE="s4dvar.in"

 set DA_STDINP="is4dpsas.in"

# Copy nonlinear model initial conditions file, use background or
# first guess state.

 cp -p ${MYROOT}/Data/${PREFIX}_bck.nc ${PREFIX}_ini.nc

# Set initial conditions, boundary conditions, and surface forcing
# error covariance input standard deviations files.

 set STDnameI=${MYROOT}/Data/${PREFIX}_std_i.nc
 set STDnameB=${MYROOT}/Data/${PREFIX}_std_b.nc
 set STDnameF=${MYROOT}/Data/${PREFIX}_std_f.nc

# Set initial conditions, boundary conditions, and surface forcing
# error covariance input or output normalization factors files.

 set NRMnameI=${MYROOT}/Data/${PREFIX}_nrm_i.nc
 set NRMnameB=${MYROOT}/Data/${PREFIX}_nrm_b.nc
 set NRMnameF=${MYROOT}/Data/${PREFIX}_nrm_f.nc

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

 $SUBSTITUTE $DA_STDINP ocean_std_i.nc $STDnameI
 $SUBSTITUTE $DA_STDINP ocean_std_b.nc $STDnameB
 $SUBSTITUTE $DA_STDINP ocean_std_f.nc $STDnameF
 $SUBSTITUTE $DA_STDINP ocean_nrm_i.nc $NRMnameI
 $SUBSTITUTE $DA_STDINP ocean_nrm_b.nc $NRMnameB
 $SUBSTITUTE $DA_STDINP ocean_nrm_f.nc $NRMnameF
 $SUBSTITUTE $DA_STDINP ocean_obs.nc $OBSname
 $SUBSTITUTE $DA_STDINP ocean_hss.nc ${PREFIX}_hss.nc
 $SUBSTITUTE $DA_STDINP ocean_lcz.nc ${PREFIX}_lcz.nc
 $SUBSTITUTE $DA_STDINP ocean_mod.nc ${PREFIX}_mod.nc
