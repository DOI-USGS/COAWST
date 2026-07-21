#!/bin/csh -f
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
# must be executed from the top level 'matlab' source code.             :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    ./bin/copyright.csh [options]                                      :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    -verbose  list files that are modified.                            :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set search = "2002-2024 The ROMS"
set replace = "2002-2025 The ROMS"

# Directories to search for replacements.

set c_dirs = "./"
set c_dirs = "$c_dirs 4dvar"
set c_dirs = "$c_dirs bathymetry"
set c_dirs = "$c_dirs bin"
set c_dirs = "$c_dirs boundary"
set c_dirs = "$c_dirs coastlines"
set c_dirs = "$c_dirs colormaps"
set c_dirs = "$c_dirs coupling"
set c_dirs = "$c_dirs forcing"
set c_dirs = "$c_dirs grid"
set c_dirs = "$c_dirs initial"
set c_dirs = "$c_dirs ioda"
set c_dirs = "$c_dirs landmask"
set c_dirs = "$c_dirs mex"
set c_dirs = "$c_dirs m_map"
set c_dirs = "$c_dirs netcdf"
set c_dirs = "$c_dirs segrid"
set c_dirs = "$c_dirs tidal_ellipse"
set c_dirs = "$c_dirs t_tide"
set c_dirs = "$c_dirs utility"

# verbose is a csh command to print all lines of the script so I changed
# this variable to "verb".

set verb = 0

while ( ($#argv) > 0 )
  switch ($1)
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
      echo "-verbose  list files that are modified"
      echo ""
      exit 1
    breaksw

  endsw
end

echo ""
echo "Replacing Copyright String in Files ..."
echo ""

# The "! -path '*/.git/*'" is there to keep it from messing with
# files in the .svn directories. The "! -name 'copyright.*'" is to
# keep it from messing with the file that's making the reaplacements.
# There is no way to redirect only stderr with csh.

foreach FILE ( `find ${c_dirs} ! -path '*/.git/*' ! -name 'copyright.*' -type f -print` )

# Double check that we're not changing a file in a .git folder.

  if ( `echo $FILE | grep -vc '.git/'` ) then
    if ( $verb == 1 ) then
      grep -l "${search}" $FILE && sed -i -e "s|${search}|${replace}|g" $FILE
    else
      grep -l "${search}" $FILE > /dev/null && sed -i -e "s|${search}|${replace}|g" $FILE
    endif
  else
    echo "There is a .git in the path: $FILE skipped"
  endif

end

echo ""
echo "Done."
echo ""

