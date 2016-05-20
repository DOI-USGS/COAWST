#!/bin/csh -f
#
# svn $Id$
#######################################################################
# Copyright (c) 2002-2016 The ROMS/TOMS Group                         #
#   Licensed under a MIT/X style license                              #
#   See License_ROMS.txt                                              #
#######################################################################
#                                                                     #
# Strong/Weak constraint 4D-PSAS job script:                          #
#                                                                     #
# This script NEEDS to be run before any run:                         #
#                                                                     #
#   (1) It copies a new clean nonlinear model initial conditions      #
#       file. The nonlinear model is initialized from the             #
#       background or reference state.                                #
#   (2) Specify model, initial conditions, boundary conditions, and   #
#       surface forcing error convariance input standard deviations   #
#       files.                                                        #
#   (3) Specify model, initial conditions, boundary conditions, and   #
#       surface forcing error convariance input/output normalization  #
#       factors files.                                                #
#   (4) Copy a clean copy of the observations NetCDF file.            #
#   (5) Create 4D-Var input script "psas.in" from template and        #
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

# Set model, initial conditions, boundary conditions and surface
# forcing error covariance standard deviations files.

 set STDnameM=${Dir}/Data/wc13_std_m.nc
 set STDnameI=${Dir}/Data/wc13_std_i.nc
 set STDnameB=${Dir}/Data/wc13_std_b.nc
 set STDnameF=${Dir}/Data/wc13_std_f.nc

# Set model, initial conditions, boundary conditions and surface
# forcing error covariance normalization factors files.

 set NRMnameM=${Dir}/Data/wc13_nrm_m.nc
 set NRMnameI=${Dir}/Data/wc13_nrm_i.nc
 set NRMnameB=${Dir}/Data/wc13_nrm_b.nc
 set NRMnameF=${Dir}/Data/wc13_nrm_f.nc

# Set observations file.

 set OBSname=wc13_obs.nc

# Get a clean copy of the observation file.  This is really
# important since this file is modified.

 cp -p ${Dir}/Data/${OBSname} .

# Modify 4D-Var template input script and specify above files.

 set PSAS=psas.in
 if (-e $PSAS) then
   /bin/rm $PSAS
 endif
 cp s4dvar.in $PSAS

 $SUBSTITUTE $PSAS ocean_std_m.nc $STDnameM
 $SUBSTITUTE $PSAS ocean_std_i.nc $STDnameI
 $SUBSTITUTE $PSAS ocean_std_b.nc $STDnameB
 $SUBSTITUTE $PSAS ocean_std_f.nc $STDnameF
 $SUBSTITUTE $PSAS ocean_nrm_m.nc $NRMnameM
 $SUBSTITUTE $PSAS ocean_nrm_i.nc $NRMnameI
 $SUBSTITUTE $PSAS ocean_nrm_b.nc $NRMnameB
 $SUBSTITUTE $PSAS ocean_nrm_f.nc $NRMnameF
 $SUBSTITUTE $PSAS ocean_obs.nc $OBSname
 $SUBSTITUTE $PSAS ocean_hss.nc wc13_hss.nc
 $SUBSTITUTE $PSAS ocean_lcz.nc wc13_lcz.nc
 $SUBSTITUTE $PSAS ocean_mod.nc wc13_mod.nc
 $SUBSTITUTE $PSAS ocean_err.nc wc13_err.nc
