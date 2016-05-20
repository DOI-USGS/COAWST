#!/bin/csh -f
#
# svn $Id$
#######################################################################
# Copyright (c) 2002-2016 The ROMS/TOMS Group                         #
#   Licensed under a MIT/X style license                              #
#   See License_ROMS.txt                                              #
#######################################################################
#                                                                     #
# 4D-Var error covariance normalization coefficients job script:      #
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
#       factors filenames.                                            #
#   (4) Create 4D-Var input script "c4dvar.in" from a template and    #
#       specify the error covariance standard deviation, and error    #
#       covariance normalization factors files to be used.            #
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
# forcing error covariance normalization factors filenames.

 set NRMnameM=wc13_nrm_m.nc
 set NRMnameI=wc13_nrm_i.nc
 set NRMnameB=wc13_nrm_b.nc
 set NRMnameF=wc13_nrm_f.nc

# Modify 4D-Var template input script and specify above files.

 set NORM=c4dvar.in
 if (-e $NORM) then
   /bin/rm $NORM
 endif
 cp s4dvar.in $NORM

 $SUBSTITUTE $NORM ocean_std_m.nc $STDnameM
 $SUBSTITUTE $NORM ocean_std_i.nc $STDnameI
 $SUBSTITUTE $NORM ocean_std_b.nc $STDnameB
 $SUBSTITUTE $NORM ocean_std_f.nc $STDnameF
 $SUBSTITUTE $NORM ocean_nrm_m.nc $NRMnameM
 $SUBSTITUTE $NORM ocean_nrm_i.nc $NRMnameI
 $SUBSTITUTE $NORM ocean_nrm_b.nc $NRMnameB
 $SUBSTITUTE $NORM ocean_nrm_f.nc $NRMnameF
 $SUBSTITUTE $NORM ocean_obs.nc wc13_obs.nc
 $SUBSTITUTE $NORM ocean_hss.nc wc13_hss.nc
 $SUBSTITUTE $NORM ocean_lcz.nc wc13_lcz.nc
 $SUBSTITUTE $NORM ocean_mod.nc wc13_mod.nc
 $SUBSTITUTE $NORM ocean_err.nc wc13_err.nc
