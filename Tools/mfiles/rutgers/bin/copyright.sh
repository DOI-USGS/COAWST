#!/bin/bash
#
# git $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2025 The ROMS Group                                :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
#                                                                       :::
# ROMS Copyright Update Script                                          :::
#                                                                       :::
# Script to update the copyright information on 'matlab' source files.  :::
# This script replaces the copyright string in the source files. It     :::
# must be executed from top level of the 'matlab' source code.          :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    ./bin/copyright.sh [options]                                       :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    -verbose  list files that are modified.                            :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

search="2002-2024 The ROMS"
replace="2002-2025 The ROMS"

# Directories to search for replacements.

c_dirs="./"
c_dirs="$c_dirs 4dvar"
c_dirs="$c_dirs bathymetry"
c_dirs="$c_dirs bin"
c_dirs="$c_dirs boundary"
c_dirs="$c_dirs coastlines"
c_dirs="$c_dirs colormaps"
c_dirs="$c_dirs coupling"
c_dirs="$c_dirs forcing"
c_dirs="$c_dirs grid"
c_dirs="$c_dirs initial"
c_dirs="$c_dirs ioda"
c_dirs="$c_dirs landmask"
c_dirs="$c_dirs mex"
c_dirs="$c_dirs m_map"
c_dirs="$c_dirs netcdf"
c_dirs="$c_dirs seagrid"
c_dirs="$c_dirs tidal_ellipse"
c_dirs="$c_dirs t_tide"
c_dirs="$c_dirs utility"

verbose=0

while [ $# -gt 0 ]
do
  case "$1" in
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
      echo "-verbose  list files that are modified"
      echo ""
      exit 1
      ;;
  esac
done

echo -e "\nReplacing Copyright String in Files ...\n"

# The "! -path '*/.git/*'" is there to keep it from messing with
# files in the .svn directories. The "! -name 'copyright.*'" is to
# keep it from messing with the file that's making the reaplacements.
# The "2>" redirects stderr so errors don't get put in FILE.

for FILE in `find ${c_dirs} ! -path '*/.git/*' ! -name 'copyright.*' -type f -print 2> /dev/null`
do

# Double check that we're not changing a file in a .svn folder.

  if [ `echo $FILE | grep -vc '.git/'` -gt 0 ]; then
    if [ $verbose -eq 1 ]; then
      grep -l "${search}" $FILE && sed -i -e "s|${search}|${replace}|g" $FILE
    else
      grep -l "${search}" $FILE > /dev/null && sed -i -e "s|${search}|${replace}|g" $FILE
    fi
  else
    echo "There is a .git in the path: $FILE skipped"
  fi
done

echo -e "\nDone.\n"


