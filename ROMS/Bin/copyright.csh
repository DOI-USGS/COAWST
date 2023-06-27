#!/bin/csh -f
#
# git $Id$
# svn $Id: copyright.csh 1152 2023-02-09 03:12:48Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
#                                                                       :::
# ROMS Copyright Update CSH Script                                      :::
#                                                                       :::
# Script to update the copyright information on ROMS source files.      :::
# This script replaces the copyright string in the source files and     :::
# updates the copyright svn property. This script must be executed      :::
# from top level of the ROMS source code.                               :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    ./ROMS/Bin/copyright.csh [options]                                 :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    -nosvn    skip updating of the svn properties. Meant for systems   :::
#                that don't have the comandline svn tools.              :::
#                                                                       :::
#    -verbose  list files that are modified.                            :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set search = "2002-2022 The ROMS/TOMS"
set replace = "2002-2023 The ROMS/TOMS"

# Directories to search for replacements.

set c_dirs = "Compilers ESM Master ROMS User"

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
  svn propset -R copyright '(c) 2002-2023 The ROMS/TOMS Group' Compilers
  svn propset -R copyright '(c) 2002-2023 The ROMS/TOMS Group' Data
  svn propset -R copyright '(c) 2002-2023 The ROMS/TOMS Group' ESM
  svn propset -R copyright '(c) 2002-2023 The ROMS/TOMS Group' Master
  svn propset -R copyright '(c) 2002-2023 The ROMS/TOMS Group' ROMS
  svn propset -R copyright '(c) 2002-2023 The ROMS/TOMS Group' User
  svn propset copyright '(c) 2002-2023 The ROMS/TOMS Group' . makefile CMakeLists.txt
else
  echo ""
  echo "Not updating svn properties."
  echo ""
endif
