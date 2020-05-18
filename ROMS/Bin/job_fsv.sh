#!/bin/csh -f
#
# svn $Id: job_fsv.sh 995 2020-01-10 04:01:28Z arango $
#######################################################################
# Copyright (c) 2002-2020 The ROMS/TOMS Group                         #
#   Licensed under a MIT/X style license                              #
#   See License_ROMS.txt                                              #
#######################################################################
#                                                                     #
#  Generalized Stability Theory: Forcing Singular Vectors             #
#                                                                     #
#  This script is used to run the ROMS/TOMS Forcing Singular Vectors  #
#  algorithm.                                                         #
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

# Set tangent linear model initial conditions file: zero fields.

set ITLname=gyre3d_itl.nc
if (-e $ITLname) then
  /bin/rm $ITLname
endif
ln -s -v ${Dir}/Data/gyre3d_ini_zero.nc $ITLname
