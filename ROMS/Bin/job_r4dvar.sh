#!/bin/csh -f
#
# svn $Id: job_r4dvar.sh 995 2020-01-10 04:01:28Z arango $
#######################################################################
# Copyright (c) 2002-2020 The ROMS/TOMS Group                         #
#   Licensed under a MIT/X style license                              #
#   See License_ROMS.txt                                              #
#######################################################################
#                                                                     #
# Strong/Weak constraint R4D-Var job script:                          #
#                                                                     #
# This script NEEDS to be run before any run:                         #
#                                                                     #
#   (1) It copies a new clean nonlinear model initial conditions      #
#       file. The nonlinear model is initialized from the             #
#       background or reference state.                                #
#   (2) It copies representer model initial condition, same as        #
#       nonlinear model.                                              #
#   (3) Specify model, initial conditions, boundary conditions, and   #
#       surface forcing error convariance input standard deviations   #
#       files.                                                        #
#   (4) Specify model, initial conditions, boundary conditions, and   #
#       surface forcing error convariance input/output normalization  #
#       factors files.                                                #
#   (5) Copy a clean copy of the observations NetCDF file.            #
#   (6) Create 4D-Var input script "r4dvar.in" from template and      #
#       specify the error covariance standard deviation, error        #
#       covariance normalization factors, and observation files to    #
#       be used.                                                      #
#                                                                     #
#######################################################################

# Set path definition to one directory up in the tree.

 set Dir=`dirname ${PWD}`

# Set string manipulations perl script.

 set SUBSTITUTE=${ROMS_ROOT}/ROMS/Bin/substitute

# Copy nonlinear model initial conditions file, use background or
# first guess state.

 cp -p ${Dir}/Data/wc13_ini.nc wc13_ini.nc

# Copy representer model initial conditions file, same as nonlinear
# model.

 cp -p ${Dir}/Data/wc13_ini.nc wc13_irp.nc

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

 set R4DVAR=r4dvar.in
 if (-e $R4DVAR) then
   /bin/rm $R4DVAR
 endif
 cp s4dvar.in $R4DVAR

 $SUBSTITUTE $R4DVAR roms_std_m.nc $STDnameM
 $SUBSTITUTE $R4DVAR roms_std_i.nc $STDnameI
 $SUBSTITUTE $R4DVAR roms_std_b.nc $STDnameB
 $SUBSTITUTE $R4DVAR roms_std_f.nc $STDnameF
 $SUBSTITUTE $R4DVAR roms_nrm_m.nc $NRMnameM
 $SUBSTITUTE $R4DVAR roms_nrm_i.nc $NRMnameI
 $SUBSTITUTE $R4DVAR roms_nrm_b.nc $NRMnameB
 $SUBSTITUTE $R4DVAR roms_nrm_f.nc $NRMnameF
 $SUBSTITUTE $R4DVAR roms_obs.nc $OBSname
 $SUBSTITUTE $R4DVAR roms_hss.nc wc13_hss.nc
 $SUBSTITUTE $R4DVAR roms_lcz.nc wc13_lcz.nc
 $SUBSTITUTE $R4DVAR roms_mod.nc wc13_mod.nc
 $SUBSTITUTE $R4DVAR roms_err.nc wc13_err.nc
