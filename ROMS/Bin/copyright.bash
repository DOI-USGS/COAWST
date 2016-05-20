#!/bin/bash
#
# svn $Id: copyright.bash 795 2016-05-11 01:42:43Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2016 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
#                                                                       :::
# ROMS/TOMS Copyright Update Script                                     :::
#                                                                       :::
# Script to update the copyright information on ROMS source files.      :::
# This script replaces the copyright string in the source files and     :::
# updates the copyright svn property. This script must be executed      :::
# from top level of the ROMS source code.                               :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    ./ROMS/Bin/copyright.bash [options]                                :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    -nosvn    skip updating of the svn properties. Meant for systems   :::
#                that don't have the comandline svn tools.              :::
#                                                                       :::
#    -verbose  list files that are modified.                            :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

search="2002-2015 The ROMS/TOMS"
replace="2002-2016 The ROMS/TOMS"

# Directories to search for replacements.

c_dirs="Apps Compilers Master ROMS User SeaIce"

# Specific files not in the "c_dirs".

special_files="makefile Waves/SWAN/Src/Module.mk Waves/SWAN/Src/waves_coupler.F Waves/SWAN/Src/swancpp.h"

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

# Replace the string in the "special_files" separately.

for FILE in $special_files
do
  if [ $verbose -eq 1 ]; then
    grep -l "${search}" $FILE && sed -i -e "s|${search}|${replace}|g" $FILE
  else
    grep -l "${search}" $FILE > /dev/null && sed -i -e "s|${search}|${replace}|g" $FILE
  fi
done

echo -e "\nDone.\n"

if [ $setsvn -eq 1 ]; then
  svn propset -R copyright '(c) 2002-2016 The ROMS/TOMS Group' Compilers
  svn propset -R copyright '(c) 2002-2016 The ROMS/TOMS Group' Data
  svn propset -R copyright '(c) 2002-2016 The ROMS/TOMS Group' Master
  svn propset -R copyright '(c) 2002-2016 The ROMS/TOMS Group' ROMS
  svn propset -R copyright '(c) 2002-2016 The ROMS/TOMS Group' User
  svn propset copyright '(c) 2002-2016 The ROMS/TOMS Group' . makefile
  svn propset copyright '(c) 2002-2016 The ROMS/TOMS Group' Waves/SWAN/Src/Module.mk
  svn propset copyright '(c) 2002-2016 The ROMS/TOMS Group' Waves/SWAN/Src/waves_coupler.F
else
  echo -e "Not updating svn properties.\n"
fi

