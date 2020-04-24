#!/bin/csh -f
#
# svn $Id: copyright.sh 996 2020-01-10 04:28:56Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2020 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
#                                                                       :::
# ROMS/TOMS Copyright Update Script                                     :::
#                                                                       :::
# Script to update the copyright information on 'matlab' source files.  :::
# This script replaces the copyright string in the source files and     :::
# updates the copyright svn property. This script must be executed      :::
# from the top level 'matlab' source code.                              :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    ./bin/copyright.sh [options]                                       :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    -nosvn    skip updating of the svn properties. Meant for systems   :::
#                that don't have the comandline svn tools.              :::
#                                                                       :::
#    -verbose  list files that are modified.                            :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set search = "2002-2019 The ROMS/TOMS"
set replace = "2002-2020 The ROMS/TOMS"

# Directories to search for replacements.

set c_dirs = "4dvar"
set c_dirs = "$c_dirs bathymetry"
set c_dirs = "$c_dirs bin"
set c_dirs = "$c_dirs boundary"
set c_dirs = "$c_dirs coastlines"
set c_dirs = "$c_dirs colormaps"
set c_dirs = "$c_dirs coupling"
set c_dirs = "$c_dirs forcing"
set c_dirs = "$c_dirs grid"
set c_dirs = "$c_dirs initial"
set c_dirs = "$c_dirs landmask"
set c_dirs = "$c_dirs netcdf"
set c_dirs = "$c_dirs utility"

set setsvn = 1

# verbose is a csh command to print all lines of the script so I changed
# this variable to "verb".

set verb = 0

while ( ($#argv) > 0 )
  switch ($1)
    case "-nosvn":
      shift
      set setsvn = 0
    breaksw

    case "-verbose":
      shift
      set verb = 1
    breaksw

    case "-*":
      echo ""
      echo "$0 : Unknown option [ $1 ]"
      echo ""
      echo "Available Options:"
      echo ""
      echo "-nosvn    skip updating of the svn properties. Meant for systems"
      echo "            that don't have the comandline svn tools."
      echo ""
      echo "-verbose  list files that are modified"
      echo ""
      exit 1
    breaksw

  endsw
end

echo ""
echo "Replacing Copyright String in Files ..."
echo ""

# The "! -path '*/.svn/*'" is there to keep it from messing with
# files in the .svn directories. The "! -name 'copyright.*'" is to
# keep it from messing with the file that's making the reaplacements.
# There is no way to redirect only stderr with csh.

foreach FILE ( `find ${c_dirs} ! -path '*/.svn/*' ! -name 'copyright.*' -type f -print` )

# Double check that we're not changing a file in a .svn folder.

  if ( `echo $FILE | grep -vc '.svn/'` ) then
    if ( $verb == 1 ) then
      grep -l "${search}" $FILE && sed -i -e "s|${search}|${replace}|g" $FILE
    else
      grep -l "${search}" $FILE > /dev/null && sed -i -e "s|${search}|${replace}|g" $FILE
    endif
  else
    echo "There is a .svn in the path: $FILE skipped"
  endif

end

echo ""
echo "Done."
echo ""

if ( $setsvn == 1 ) then
  svn propset -R copyright '(c) 2002-2020 The ROMS/TOMS Group' 4dvar
  svn propset -R copyright '(c) 2002-2020 The ROMS/TOMS Group' bathymetry
  svn propset -R copyright '(c) 2002-2020 The ROMS/TOMS Group' bin
  svn propset -R copyright '(c) 2002-2020 The ROMS/TOMS Group' boundary
  svn propset -R copyright '(c) 2002-2020 The ROMS/TOMS Group' coastlines
  svn propset -R copyright '(c) 2002-2020 The ROMS/TOMS Group' colormaps
  svn propset -R copyright '(c) 2002-2020 The ROMS/TOMS Group' coupling
  svn propset -R copyright '(c) 2002-2020 The ROMS/TOMS Group' forcing
  svn propset -R copyright '(c) 2002-2020 The ROMS/TOMS Group' grid
  svn propset -R copyright '(c) 2002-2020 The ROMS/TOMS Group' initial
  svn propset -R copyright '(c) 2002-2020 The ROMS/TOMS Group' landmask
  svn propset -R copyright '(c) 2002-2020 The ROMS/TOMS Group' netcdf
  svn propset -R copyright '(c) 2002-2020 The ROMS/TOMS Group' utility

  svn propset copyright '(c) 2002-2020 The ROMS/TOMS Group' mex
  svn propset copyright '(c) 2002-2020 The ROMS/TOMS Group' mex/Contents.m
  svn propset copyright '(c) 2002-2020 The ROMS/TOMS Group' mex/mexinside
  svn propset copyright '(c) 2002-2020 The ROMS/TOMS Group' mex/mexrect
  svn propset copyright '(c) 2002-2020 The ROMS/TOMS Group' mex/mexsepeli
  svn propset copyright '(c) 2002-2020 The ROMS/TOMS Group' seagrid
  svn propset copyright '(c) 2002-2020 The ROMS/TOMS Group' seagrid/presto
  svn propset copyright '(c) 2002-2020 The ROMS/TOMS Group' seagrid/presto/@presto
  svn propset copyright '(c) 2002-2020 The ROMS/TOMS Group' seagrid/presto/@ps
  svn propset copyright '(c) 2002-2020 The ROMS/TOMS Group' seagrid/presto/@pst
  svn propset copyright '(c) 2002-2020 The ROMS/TOMS Group' seagrid/@seagrid
  svn propset copyright '(c) 2002-2020 The ROMS/TOMS Group' seagrid/test_data
  svn propset copyright '(c) 2002-2020 The ROMS/TOMS Group' seawater
  svn propset copyright '(c) 2002-2020 The ROMS/TOMS Group' t_tide
  svn propset copyright '(c) 2002-2020 The ROMS/TOMS Group' . startup.m
else
  echo ""
  echo "Not updating svn properties."
  echo ""
endif

