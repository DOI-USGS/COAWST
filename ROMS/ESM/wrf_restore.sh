#!/bin/bash
#
# git $Id$
# svn $Id: wrf_restore.sh 1151 2023-02-09 03:08:53Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# WRF Restore BASH Script: WRF Versions 4.1 and up                      :::
#                                                                       :::
# This script restores WRF directory to its original checkout state.    :::
# It undo the changes made by "wrf_patch.sh".                           :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    ./build_restore.sh wrf_root_dir                                    :::
#                                                                       :::
# where "wrf_root_dir" is the local path for the installed WRF code.    :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Initialize.

separator=`perl -e "print ':' x 100;"`

if [ $# -eq 1 ]; then
  if [ -d $1 ]; then
    WRF_ROOT_DIR=$1
  else
    echo "Provided directory does not exist."
    exit 1
  fi
else
  echo "You must provide only the WRF source code root directory."
  exit 1
fi

# Clean source code.

if [ -f ${WRF_ROOT_DIR}/clean ]; then
  echo ""
  echo "${separator}"
  echo "Cleaning WRF source code:  ${WRF_ROOT_DIR}/clean -a"
  echo "${separator}"
  echo ""
  ${WRF_ROOT_DIR}/clean -a
else
  echo "${WRF_ROOT_DIR} does not appear to contain WRF source code"
  exit 1
fi

# Undo ROMS patches

if [ -f ${WRF_ROOT_DIR}/configure.orig ]; then
  mv -v ${WRF_ROOT_DIR}/configure.orig ${WRF_ROOT_DIR}/configure
else
  echo "   No need to restore: ${WRF_ROOT_DIR}/configure"
fi

if [ -f ${WRF_ROOT_DIR}/Makefile.orig ]; then
  mv -v  ${WRF_ROOT_DIR}/Makefile.orig ${WRF_ROOT_DIR}/Makefile
else
  echo "   No need to restore: ${WRF_ROOT_DIR}/Makefile"
fi

if [ -f ${WRF_ROOT_DIR}/arch/postamble.orig ]; then
  mv -v  ${WRF_ROOT_DIR}/arch/postamble.orig ${WRF_ROOT_DIR}/arch/postamble
else
  echo "   No need to restore: ${WRF_ROOT_DIR}/arch/postamble"
fi

if [ -f ${WRF_ROOT_DIR}/arch/configure.defaults.orig ]; then
  mv -v  ${WRF_ROOT_DIR}/arch/configure.defaults.orig ${WRF_ROOT_DIR}/arch/configure.defaults
else
  echo "   No need to restore: ${WRF_ROOT_DIR}/arch/configure.defaults"
fi

if [ -f ${WRF_ROOT_DIR}/external/esmf_time_f90/Makefile.orig ]; then
  mv -v  ${WRF_ROOT_DIR}/external/esmf_time_f90/Makefile.orig ${WRF_ROOT_DIR}/external/esmf_time_f90/Makefile
else
  echo "   No need to restore: ${WRF_ROOT_DIR}/external/esmf_time_f90/Makefile"
fi

if [ -f ${WRF_ROOT_DIR}/external/esmf_time_f90/Test1.F90.orig ]; then
  mv -v  ${WRF_ROOT_DIR}/external/esmf_time_f90/Test1.F90.orig ${WRF_ROOT_DIR}/external/esmf_time_f90/Test1.F90
else
  echo "   No need to restore: ${WRF_ROOT_DIR}/external/esmf_time_f90/Test1.F90"
fi

