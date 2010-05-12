#!/bin/csh -f
#
# svn $Id: job_op.sh 429 2009-12-20 17:30:26Z arango $
#######################################################################
# Copyright (c) 2002-2010 The ROMS/TOMS Group                         #
#   Licensed under a MIT/X style license                              #
#   See License_ROMS.txt                                              #
#######################################################################
#                                                                     #
#  Generalized Stability Theory: Optimal Perturbations                #
#                                                                     #
#  This script is used to run the ROMS/TOMS optimal perturbations     #
#  algorithm.                                                         #
#                                                                     #
#######################################################################

# Set ROOT of the directory to run Optimal Perturbations.

set MYROOT="/home/arango/ocean/toms/adjoint/Test/EAC4"

# Set application prefix.

set PREFIX="eac4"

# Set basic state trajectory, forward file:

set HISname=${MYROOT}/Forward/${PREFIX}_his.nc

set FWDname=${PREFIX}_fwd.nc

if (-e $FWDname) then
  /bin/rm $FWDname
endif
ln -s $HISname $FWDname

# Set zero fields initial condition file

set ZEROname=${MYROOT}/Data/${PREFIX}_ini_zero.nc

# Set tangent linear model initial conditions file: zero fields.

set ITLname=${PREFIX}_itl.nc

if (-e $ITLname) then
  /bin/rm $ITLname
endif
ln -s $ZEROname $ITLname

# Set adjoint model initial conditions file: zero fields.

set IADname=${PREFIX}_iad.nc

if (-e $IADname) then
  /bin/rm $IADname
endif
ln -s $ZEROname $IADname
