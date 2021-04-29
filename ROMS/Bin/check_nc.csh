#!/bin/csh -f
#
# svn $Id: check_nc.csh 1054 2021-03-06 19:47:12Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2021 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
#                                                                       :::
# ROMS/TOMS NetCDF file checker between simulations:                    :::
#                                                                       :::
# This script compares the binary differences of all ROMS/TOMS output   :::
# NetCDF files between two simulations of the same application.  The    :::
# DEBUGGING and POSITIVE_ZERO options need to be activated to avoid     :::
# time marks in output NetCDF files and other header information.       :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    check_nc.csh <dir1> <dir2>                                         :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set dir1=$1
set dir2=$2
set diffs=0

# diff the file pairs one by one and increment $diffs if necessary.
#
# $? stores the exit code of the previous command, in this case diff.
# If $? is non-zero it means that differences were found so we should
# increment $diffs.

foreach file ($dir1/*.nc)
  set f = `basename $file`
  diff -q ${dir1}/${f} ${dir2}/${f} >& /dev/null
  if( $? != 0 ) then
    echo "${dir1}/${f} and ${dir2}/${f} differ!"
    @ diffs +=1
  endif
end

# Exit and set exit code to $diffs so we can sum total differences
# (hopefully 0) in the calling script.

exit ${diffs}

