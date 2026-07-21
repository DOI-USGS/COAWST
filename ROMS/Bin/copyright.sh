#!/bin/bash
#
# git $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2026 The ROMS Group                                :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
#                                                                       :::
# ROMS Copyright Update BASH Script                                     :::
#                                                                       :::
# Script to update the copyright information on ROMS source files.      :::
# This script replaces the copyright string in the source files and     :::
# updates the copyright svn property. This script must be executed      :::
# from top level of the ROMS source code.                               :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    ./ROMS/Bin/copyright.sh [options]                                  :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    -nosvn    skip updating of the svn properties. Meant for systems   :::
#                that don't have the comandline svn tools.              :::
#                                                                       :::
#    -verbose  list files that are modified.                            :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

search="2002-2025 The ROMS"
replace="2002-2026 The ROMS"

# Special case for License_ROMS.md

searchL="2002-2025 Regional"
replaceL="2002-2026 Regional"


# Directories to search for replacements.

c_dirs="Compilers Data docs ESM Master ROMS User"

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

sed -i -e "s|${search}|${replace}|g" CMakeLists.txt
sed -i -e "s|${searchL}|${replaceL}|g" License_ROMS.md
sed -i -e "s|${search}|${replace}|g" makefile

echo -e "\nDone.\n"

if [ $setsvn -eq 1 ]; then
  svn propset -R copyright '(c) 2002-2025 The ROMS Group' Compilers
  svn propset -R copyright '(c) 2002-2025 The ROMS Group' Data
  svn propset -R copyright '(c) 2002-2025 The ROMS Group' ESM
  svn propset -R copyright '(c) 2002-2025 The ROMS Group' Master
  svn propset -R copyright '(c) 2002-2025 The ROMS Group' ROMS
  svn propset -R copyright '(c) 2002-2025 The ROMS Group' User
  svn propset copyright '(c) 2002-2025 The ROMS Group' . makefile CMakeLists.txt
else
  echo -e "Not updating svn properties.\n"
fi

