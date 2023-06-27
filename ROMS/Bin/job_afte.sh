#!/bin/bash
#
# git $Id$
# svn $Id: job_afte.sh 1151 2023-02-09 03:08:53Z arango $
#######################################################################
# Copyright (c) 2002-2023 The ROMS/TOMS Group                         #
#   Licensed under a MIT/X style license                              #
#   See License_ROMS.txt                                              #
#######################################################################
#                                                                     #
#  Generalized Stability Theory: Adjoint Finite Time Eigenmodes       #
#                                                                     #
#  This script is used to set-up ROMS Adjoint Finite Time Eigenmodes  #
#  algorithm.                                                         #
#                                                                     #
#######################################################################

# Set ROOT of the directory to run application.  The following
# "dirname" command returns a path by removing any suffix from
# the last slash ('/').  It returns a path above current diretory.

Dir=`dirname ${PWD}`

# Set basic state trajectory, forward file:

#HISname=${Dir}/Forward/gyre3d_his_00.nc
 HISname=${Dir}/Forward/gyre3d_his_01.nc

FWDname=gyre3d_fwd.nc

if [ -f $FWDname ]; then
  /bin/rm $FWDname
fi
ln -s -v $HISname $FWDname

# Set adjoint and tangent linear model initial conditions file:
# zero fields.

IADname=gyre3d_iad.nc

if [ -f $IADname ]; then
  /bin/rm $IADname
fi

ln -s -v ${Dir}/Data/gyre3d_ini_zero.nc $IADname
