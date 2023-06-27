#!/bin/bash
#
# git $Id$
# svn $Id: job_rbl4dvar.sh 1151 2023-02-09 03:08:53Z arango $
#######################################################################
# Copyright (c) 2002-2023 The ROMS/TOMS Group                         #
#   Licensed under a MIT/X style license                              #
#   See License_ROMS.txt                                              #
#######################################################################
#                                                                     #
# Strong/Weak constraint RBL4D-Var job BASH script:                   #
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
#   (5) Create 4D-Var input script "rbl4dvar.in" from template and    #
#       specify the error covariance standard deviation, error        #
#       covariance normalization factors, and observation files to    #
#       be used.                                                      #
#                                                                     #
#######################################################################

# Set path definition to one directory up in the tree.

 Dir=`dirname ${PWD}`

# Set string manipulations perl script.

 SUBSTITUTE=${ROMS_ROOT}/ROMS/Bin/substitute

# Copy nonlinear model initial conditions file.

 cp -p ${Dir}/Data/wc13_ini.nc wc13_ini.nc

# Set model, initial conditions, boundary conditions and surface
# forcing error covariance standard deviations files.

 STDnameM=${Dir}/Data/wc13_std_m.nc
 STDnameI=${Dir}/Data/wc13_std_i.nc
 STDnameB=${Dir}/Data/wc13_std_b.nc
 STDnameF=${Dir}/Data/wc13_std_f.nc

# Set model, initial conditions, boundary conditions and surface
# forcing error covariance normalization factors files.

 NRMnameM=${Dir}/Data/wc13_nrm_m.nc
 NRMnameI=${Dir}/Data/wc13_nrm_i.nc
 NRMnameB=${Dir}/Data/wc13_nrm_b.nc
 NRMnameF=${Dir}/Data/wc13_nrm_f.nc

# Set observations file.

 OBSname=wc13_obs.nc

# Get a clean copy of the observation file.  This is really
# important since this file is modified.

 cp -p ${Dir}/Data/${OBSname} .

# Modify 4D-Var template input script and specify above files.

 RBL4DVAR=rbl4dvar.in
 if [ -f $RBL4DVAR ]; then
   /bin/rm $RBL4DVAR
 fi
 cp s4dvar.in $RBL4DVAR

 $SUBSTITUTE $RBL4DVAR roms_std_m.nc $STDnameM
 $SUBSTITUTE $RBL4DVAR roms_std_i.nc $STDnameI
 $SUBSTITUTE $RBL4DVAR roms_std_b.nc $STDnameB
 $SUBSTITUTE $RBL4DVAR roms_std_f.nc $STDnameF
 $SUBSTITUTE $RBL4DVAR roms_nrm_m.nc $NRMnameM
 $SUBSTITUTE $RBL4DVAR roms_nrm_i.nc $NRMnameI
 $SUBSTITUTE $RBL4DVAR roms_nrm_b.nc $NRMnameB
 $SUBSTITUTE $RBL4DVAR roms_nrm_f.nc $NRMnameF
 $SUBSTITUTE $RBL4DVAR roms_obs.nc $OBSname
 $SUBSTITUTE $RBL4DVAR roms_hss.nc wc13_hss.nc
 $SUBSTITUTE $RBL4DVAR roms_lcz.nc wc13_lcz.nc
 $SUBSTITUTE $RBL4DVAR roms_mod.nc wc13_mod.nc
 $SUBSTITUTE $RBL4DVAR roms_err.nc wc13_err.nc
