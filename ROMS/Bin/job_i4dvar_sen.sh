#!/bin/csh -f
#
# svn $Id$
#######################################################################
# Copyright (c) 2002-2016 The ROMS/TOMS Group                         #
#   Licensed under a MIT/X style license                              #
#   See License_ROMS.txt                                              #
#######################################################################
#                                                                     #
# Incremental I4D-Var observations impact job script:                 #
#                                                                     #
# This script NEEDS to be run before any run:                         #
#                                                                     #
#   (1) It copies a new clean nonlinear model initial conditions      #
#       file. The nonlinear model is initialized from the             #
#       background or reference state.                                #
#   (2) It copies the Lanczos vectors from a previos I4D-Var run.     #
#       They are stored in the output adjoint NetCDF file.            #
#   (3) It copies the adjoint sensitivy functional file for the       #
#       observation impact.                                           #
#   (4) Specify initial conditions, boundary conditions, and surface  #
#       forcing error convariance input standard deviations files.    #
#   (5) Specify initial conditions, boundary conditions, and surface  #
#       forcing error convariance input/output normalization factors  #
#       files.                                                        #
#   (6) Copy a clean copy of the observations NetCDF file.            #
#   (7) Create 4D-Var input script "i4dvar.in" from template and      #
#       specify the error covariance standard deviation, error        #
#       covariance normalization factors, and observation files to    #
#       be used.                                                      #
#                                                                     #
#######################################################################

# Set path definition to one directory up in the tree.

 set Dir=`dirname ${PWD}`

# Set string manipulations perl script.

 set SUBSTITUTE=${ROMS_ROOT}/ROMS/Bin/substitute

# Copy nonlinear model initial conditions file.

 cp -p ${Dir}/Data/wc13_ini.nc wc13_ini.nc

# Copy Lanczos vectors (adjoint Netcdf file) from a previous I4D-Var
# run.

 cp -p ${Dir}/I4DVAR/wc13_adj_001.nc wc13_lcz.nc

# Copy adjoint sensitivity functional.

 cp -p ${Dir}/Data/wc13_ads.nc wc13_ads.nc

# Set initial conditions, boundary conditions and surface forcing
# error covariance standard deviations files.

 set STDnameI=${Dir}/Data/wc13_std_i.nc
 set STDnameB=${Dir}/Data/wc13_std_b.nc
 set STDnameF=${Dir}/Data/wc13_std_f.nc

# Set initial conditions, boundary conditions and surface forcing
# error covariance normalization factors files.

 set NRMnameI=${Dir}/Data/wc13_nrm_i.nc
 set NRMnameB=${Dir}/Data/wc13_nrm_b.nc
 set NRMnameF=${Dir}/Data/wc13_nrm_f.nc

# Set observations file.

 set OBSname=wc13_obs.nc

# Get a clean copy of the observation file.  This is really
# important since this file is modified.

 cp -p ${Dir}/Data/${OBSname} .

# Modify 4D-Var template input script and specify above files.

 set I4DVAR=i4dvar.in
 if (-e $I4DVAR) then
   /bin/rm $I4DVAR
 endif
 cp s4dvar.in $I4DVAR

 $SUBSTITUTE $I4DVAR ocean_std_i.nc $STDnameI
 $SUBSTITUTE $I4DVAR ocean_std_b.nc $STDnameB
 $SUBSTITUTE $I4DVAR ocean_std_f.nc $STDnameF
 $SUBSTITUTE $I4DVAR ocean_nrm_i.nc $NRMnameI
 $SUBSTITUTE $I4DVAR ocean_nrm_b.nc $NRMnameB
 $SUBSTITUTE $I4DVAR ocean_nrm_f.nc $NRMnameF
 $SUBSTITUTE $I4DVAR ocean_obs.nc $OBSname
 $SUBSTITUTE $I4DVAR ocean_hss.nc wc13_hss.nc
 $SUBSTITUTE $I4DVAR ocean_lcz.nc wc13_lcz.nc
 $SUBSTITUTE $I4DVAR ocean_mod.nc wc13_mod.nc
 $SUBSTITUTE $I4DVAR ocean_err.nc wc13_err.nc
