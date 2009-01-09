#!/bin/csh -f
#
# svn $Id: job_afte.sh 503 2008-01-10 00:11:51Z arango $
#######################################################################
# Copyright (c) 2002-2008 The ROMS/TOMS Group                         #
#   Licensed under a MIT/X style license                              #
#   See License_ROMS.txt                                              #
#######################################################################
#                                                                     #
#  Generalized Stability Theory: Adjoint Finite Time Eigenmodes       #
#                                                                     #
#  This script is used to run the ROMS/TOMS Adjoint Finite Time       #
#  Eigenmodes algorithm.                                              #
#                                                                     #
#######################################################################

# Set ROOT of the directory to run Optimal Perturbations.

set MYROOT="/home/arango/Work/EAC4"

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
