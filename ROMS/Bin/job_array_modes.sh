#!/bin/csh -f
#
# svn $Id$
#######################################################################
# Copyright (c) 2002-2016 The ROMS/TOMS Group                         #
#   Licensed under a MIT/X style license                              #
#   See License_ROMS.txt                                              #
#######################################################################
#                                                                     #
# Stabilized representer matrix array modes job script:               #
#                                                                     #
# This script NEEDS to be run before any run:                         #
#                                                                     #
#   (1) It copies a new clean nonlinear model initial conditions      #
#       file. The nonlinear model is initialized from the             #
#       background or reference state.                                #
#   (2) It copies representer model initial condition, same as        #
#       nonlinear model.                                              #
#   (3) It copies Lanczos vectors from previous R4D-Var run. They     #
#       are stored in 4D-Var data assimilation file.                  #
#   (4) Specify model, initial conditions, boundary conditions, and   #
#       surface forcing error convariance input standard deviations   #
#       files.                                                        #
#   (5) Specify model, initial conditions, boundary conditions, and   #
#       surface forcing error convariance input/output normalization  #
#       factors files.                                                #
#   (6) Copy a clean copy of the observations NetCDF file.            #
#   (7) Create 4D-Var input script "r4dvar.in" from template and      #
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

# Copies Lanczos vectors from previous R4D-Var run. They are stored
# in 4D-Var data assimilation file.

 cp -p ${Dir}/R4DVAR/wc13_mod.nc wc13_lcz.nc

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

 $SUBSTITUTE $R4DVAR ocean_std_m.nc $STDnameM
 $SUBSTITUTE $R4DVAR ocean_std_i.nc $STDnameI
 $SUBSTITUTE $R4DVAR ocean_std_b.nc $STDnameB
 $SUBSTITUTE $R4DVAR ocean_std_f.nc $STDnameF
 $SUBSTITUTE $R4DVAR ocean_nrm_m.nc $NRMnameM
 $SUBSTITUTE $R4DVAR ocean_nrm_i.nc $NRMnameI
 $SUBSTITUTE $R4DVAR ocean_nrm_b.nc $NRMnameB
 $SUBSTITUTE $R4DVAR ocean_nrm_f.nc $NRMnameF
 $SUBSTITUTE $R4DVAR ocean_obs.nc $OBSname
 $SUBSTITUTE $R4DVAR ocean_hss.nc wc13_hss.nc
 $SUBSTITUTE $R4DVAR ocean_lcz.nc wc13_lcz.nc
 $SUBSTITUTE $R4DVAR ocean_mod.nc wc13_mod.nc
 $SUBSTITUTE $R4DVAR ocean_err.nc wc13_err.nc
