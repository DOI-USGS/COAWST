#!/bin/bash
#
# svn $Id: copyright.bash 996 2020-01-10 04:28:56Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2019 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
#                                                                       :::
# ROMS/TOMS Copyright Update Script                                     :::
#                                                                       :::
# Script to update the copyright information on 'matlab' source files.  :::
# This script replaces the copyright string in the source files and     :::
# updates the copyright svn property. This script must be executed      :::
# from top level of the 'matlab' source code.                           :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    ./bin/copyright.bash [options]                                     :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    -nosvn    skip updating of the svn properties. Meant for systems   :::
#                that don't have the comandline svn tools.              :::
#                                                                       :::
#    -verbose  list files that are modified.                            :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

search="2002-2019 The ROMS/TOMS"
replace="2002-2020 The ROMS/TOMS"

# Directories to search for replacements.

c_dirs="4dvar"
c_dirs="$c_dirs bathymetry"
c_dirs="$c_dirs bin"
c_dirs="$c_dirs boundary"
c_dirs="$c_dirs coastlines"
c_dirs="$c_dirs colormaps"
c_dirs="$c_dirs coupling"
c_dirs="$c_dirs forcing"
c_dirs="$c_dirs grid"
c_dirs="$c_dirs initial"
c_dirs="$c_dirs landmask"
c_dirs="$c_dirs netcdf"
c_dirs="$c_dirs utility"

setsvn=1
verbose=0

while [ $# -gt 0 ]
do
  case "$1" in
    -nosvn )
      shift
      setsvn=0
      ;;

    -verbose )
      shift
      verbose=1
      ;;

    * )
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
      ;;
  esac
done

echo -e "\nReplacing Copyright String in Files ...\n"

# The "! -path '*/.svn/*'" is there to keep it from messing with
# files in the .svn directories. The "! -name 'copyright.*'" is to
# keep it from messing with the file that's making the reaplacements.
# The "2>" redirects stderr so errors don't get put in FILE.

for FILE in `find ${c_dirs} ! -path '*/.svn/*' ! -name 'copyright.*' -type f -print 2> /dev/null`
do

# Double check that we're not changing a file in a .svn folder.

  if [ `echo $FILE | grep -vc '.svn/'` -gt 0 ]; then
    if [ $verbose -eq 1 ]; then
      grep -l "${search}" $FILE && sed -i -e "s|${search}|${replace}|g" $FILE
    else
      grep -l "${search}" $FILE > /dev/null && sed -i -e "s|${search}|${replace}|g" $FILE
    fi
  else
    echo "There is a .svn in the path: $FILE skipped"
  fi
done

echo -e "\nDone.\n"

if [ $setsvn -eq 1 ]; then
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
  echo -e "Not updating svn properties.\n"
fi

