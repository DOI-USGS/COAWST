#!/bin/csh -f
#
# svn $Id: job_afte.sh 429 2009-12-20 17:30:26Z arango $
#######################################################################
# Copyright (c) 2002-2016 The ROMS/TOMS Group                         #
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

# Set ROOT of the directory to run application.  The following
# "dirname" command returns a path by removing any suffix from
# the last slash ('/').  It returns a path above current diretory.

set Dir=`dirname ${PWD}`

# Set basic state trajectory, forward file:

#set HISname=${Dir}/Forward/gyre3d_his_00.nc
 set HISname=${Dir}/Forward/gyre3d_his_01.nc

set FWDname=gyre3d_fwd.nc

if (-e $FWDname) then
  /bin/rm $FWDname
endif
ln -s -v $HISname $FWDname

# Set adjoint and tangent linear model initial conditions file:
# zero fields.

set IADname=gyre3d_iad.nc

if (-e $IADname) then
  /bin/rm $IADname
endif

ln -s -v ${Dir}/Data/gyre3d_ini_zero.nc $IADname
