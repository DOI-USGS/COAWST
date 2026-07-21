#!/bin/bash
#
# git $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2026 The ROMS Group                                :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
#                                                                       :::
# ROMS-UFS Git sparse checkout BASH script:                             :::
#                                                                       :::
# Script to checkout just the directories needed for configuring and    :::
# running ROMS test applications (currently only IRENE) in the UFS      :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    ./roms_test_ufs.sh [options]                                       :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    -b branch      Checkout a specific roms_test GitHub branch         :::
#                                                                       :::
#                     roms_test_ufs.sh -b feature/ufs                   :::
#                                                                       :::
#    -o outdir      Clone into the specified directory. If omitted,     :::
#                     roms_ufs will be used                             :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

outdir="roms_ufs"
branch="main"

command="roms_test_ufs.sh $@"

separator=`perl -e "print '<>' x 50;"`

while [ $# -gt 0 ]
do
  case "$1" in
    -b )
      shift
      branch=`echo $1 | grep -v '^-'`
      if [ "$branch" == "" ]; then
        echo "Please enter a roms_test GitHub branch name or omit the -b option for main."
        exit 1
      fi
      shift
      ;;

    -o )
      shift
      outdir=`echo $1 | grep -v '^-'`
      if [ "outdir" == "" ]; then
        echo "Please enter an output directory or omit the -o option for roms_ufs."
        exit 1
      fi
      shift
      ;;

    * )
      echo ""
      echo "${separator}"
      echo "$0 : Unknown option [ $1 ]"
      echo ""
      echo "Available Options:"
      echo ""
      echo "-b branch       Checkout a specific roms_test GitHub branch"
      echo "                  For example:  roms_test_ufs.sh -b feature/ufs"
      echo ""
      echo "-o outdir       Clone into the specified directory"
      echo "                  For example:  roms_test_ufs.sh -o roms_ufs_tests"
      echo ""
      exit 1
      ;;
  esac
done

# Make sure the requested output directory does not exist. This would indicate either a
# mistake or a canceled or failed previous run of this script.

if [ -d ${outdir} ]; then
  echo "Directory '${outdir}' already exist; move, remove, or choose a different directory."
  exit 1
fi

# Clone just the repository info but do not download any code or data

echo "${separator}"
echo "Cloning repository structure into $outdir"
echo "${separator}"
git clone --filter=blob:none --no-checkout https://github.com/myroms/roms_test $outdir

cd $outdir

# Initialize the sparse-checkout. The --cone option improves performance by excluding
# checking for patterns below directories we are not interested in.

echo "${separator}"
echo "Initializing sparse-checkout"
echo ""
git sparse-checkout init --cone

# Set the list of directories that we want included in our sparse-checkout cone

echo "Setting directories for download"
echo ""
git sparse-checkout set External IRENE/Data IRENE/Coupling/roms_data_cmeps IRENE/Coupling/roms_data_cdeps

# This is the command that actually filters out the extra directories you did not request.
# This typically takes about a minute but future pulls and branch switches should be much quicker

echo "Downloading and filtering data from requested branch."
echo "${separator}"
git checkout $branch

echo "${separator}"
echo "Checkout complete"
echo ""
echo "Script command: ${command}"
echo "${separator}"
echo ""
