#!/bin/bash
#
# svn $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2016 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
#                                                                       :::
# ROMS/TOMS applications parallel partitions checker:                   :::
#                                                                       :::
# Script to compile an user application with different tile partitions  :::
# (1x1, 2x2, 3x3) in serial, shared-memory, and distributed-memory to   :::
# find parallel bugs. It can compare two versions of the code.          :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    verify.bash [options]                                              :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    [-nobuild]   Do not compile the code (optional)                    :::
#    [-j [N]]     Compile in parallel using N CPUs (optional)           :::
#                   omit argument for all available CPUs                :::
#    [-cleanup]   Remove all traces that comparison tests were run      :::
#    -in          ROMS/TOMS input script (ocean.in)                     :::
#    <src_dir1>   ROMS/TOMS source directory 1                          :::
#    <src_dir2>   ROMS/TOMS source directory 2, used when comparing     :::
#                   two different versions of the code (optional)       :::
#                                                                       :::
#    The -j and -noclean options are passed directly to a template      :::
#    build script "build_tpl.sh", which is created from the application :::
#    "build.sh" script.  The default ROMS build command is:             :::
#                                                                       :::
#    ./build_tpl.sh -j 4                                                :::
#                                                                       :::
#    If "src_dir1" and "src_dir2" does not start with a / then it is    :::
#    appended to the string "${MY_ROOT_DIR}/" for insertion in the      :::
#    build script. If it starts with a / then the full path is used     :::
#    ("${MY_ROOT_DIR}/" is not pre-pended).                             :::
#                                                                       :::
# Examples:                                                             :::
#                                                                       :::
# (1) verify.bash -in ocean_upwelling.in branches/arango nesting        :::
#                                                                       :::
#                 It will compile, run and compare all configurations   :::
#                   of branches/arango and nesting                      :::
#                                                                       :::
# (2) verify.bash -nobuild -in ocean_upwelling.in                       :::
#                 /Users/arango/ocean/repository/trunk                  :::
#                                                                       :::
#                 This run and compare all configuration of trunk       :::
#                   without recompiling the executables. However, if    :::
#                   an executable is missing it will rebuilt it.        :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# Variables that need to be passed to perl for build script replacements

declare -x src
declare -x config
declare -x cpps

# Set variables that get passed to the build script that will not change
# for the duration of the tests

fort="ifort"
which_mpi="openmpi"

# Serial build configuration

serial_config="#####STARTCONFIG#####
 setenv FORT          $fort
 setenv USE_LARGE     on
#####ENDCONFIG#####"

# MPI build configuation

mpi_config="#####STARTCONFIG#####
 setenv USE_MPI       on
 setenv USE_MPIF90    on
 setenv which_MPI     $which_mpi
 setenv FORT          $fort
 setenv USE_LARGE     on
#####ENDCONFIG#####"

# OpenMP build configuration

openmp_config="#####STARTCONFIG#####
 setenv USE_OpenMP    on
 setenv FORT          $fort
 setenv USE_LARGE     on
#####ENDCONFIG#####"

# Set unchanging variables local to this script.

log="log.01"
infile=""


# Explaination of perl lines:
#
# perl -p0777 -i -e
#
#    -p tells perl to print results of the -e statements to stdout.
#
#    0777 tells perl to treat the entire file as if it were one long string.
#    This facilitates replacing blocks of text.
#
#    -i tells perl to operate on the file in place
#
#    -e tells per to execute the following command
#
# 's|^\s*[^#]\s*(setenv MY_CPP_FLAGS)(.*)$| $1$2\n$ENV{cpps}\n|m' build_tpl.sh
#
#    ^ means start of line, \s* means 0 or more white space, and [^#] means
#    not # (not commented)
#
#    $1 and $2 hold "setenv MY_CPP_FLAGS" and the first CPP flag respectively
#
#    $ENV{cpps} holds needed additional MY_CPP_FLAGS.
#
#    The "m" modifier allows ^ and $ to correspond to the begining and end of
#    lines instead of the entire string (the whole file).
#
#    build_tpl.sh is the file you are modifying in place
#
# 's/^(#| |)setenv USE_MPI(.*|\n)*setenv USE_MY_LIBS(.*)$/#####STARTCONFIG#####\n
# #####ENDCONFIG#####/m' build_tpl.sh
#
#

# 's/^[^#](setenv\s+MY_ROMS_SRC\s+).*$| $1ENV{src}/m;s/^(#| |)setenv USE_MPI(.*|\n)*
#  setenv USE_MY_LIBS(.*)$/$ENV{config}/m' build.sh > build_tpl.sh
#
# Notice that the semicolon at the end of the first -e statement is
# needed to allow the second statement to be executed.
#
# -p tells perl to print results of the -e statement to stdout.
#
# 0777 tells perl to treat the entire file as if it were one long string.
# This facilitates replacing blocks of text.
#
# The "m" modifier allows ^ and $ to correspond to the begining and end of
# lines instead of the entire string (the whole file).
#
# ^ means start of line, [^#] means not # (not commented), $ENV{src}
# holds the source path.
#
# ^(#| |) matches lines that start with #, space or nothing (no space).
# This will match the likely forms of the setenv USE_MPI line.
#
# (.*|\n)* matches anything including new lines (\n), the second * allows
# repeat matches to match multiple lines.
#
# $ENV{config} holds the block of compile options that choose whether to
# enable serial, MPI, or OpenMP and choses the compiler.


# Processes command line options.

parallel=1
build=1
cleanup=0
cnt=0
ncpus=""
dir1=""
dir2=""

while [ $# -gt 0 ]
do
  case "$1" in
    -nobuild )
      shift
      build=0
      ;;

    -cleanup )
      shift
      cleanup=1
      ;;

    -j )
      shift
      parallel=1
      test=`echo $1 | grep '^[0-9]\+$'`
      if [ "$test" != "" ]; then
        ncpus="-j $1"
        shift
      else
        ncpus="-j"
      fi
      ;;

    -in )
      shift
      infile=$1
      shift
      ;;

    -* )
      echo ""
      echo "$0 : Unknown option [ $1 ]"
      echo ""
      echo "Available Options:"
      echo ""
      echo "-nobuild    Do not build roms if executable already"
      echo "              exist, only run and compare"
      echo "-j [N]      Compile in parallel using N CPUs"
      echo "              set N to 1 for serial compile"
      echo "              omit argument for all avaliable CPUs"
      echo "-cleanup    Remove all traces that comparison tests"
      echo "              were run"
      echo "-in infile  Set the ocean.in file to infile"
      echo ""
      exit 1
      ;;

    * )
      cnt=`expr ${cnt} + 1`
      if [[ $dir1 == "" ]]; then
        # Strip trailing slashes if present
        dir1=`echo $1 | sed 's|/*$||'`
        # if path is not absolute prepend '${MY_ROOT_DIR}/'
        if [[ $dir1 != /* ]]; then
          dir1='${MY_ROOT_DIR}/'$dir1
        fi
        shift
      else
        # Strip trailing slashes if present
        dir2=`echo $1 | sed 's|/*$||'`
        # if path is not absolute prepend '${MY_ROOT_DIR}/'
        if [[ $dir2 != /* ]]; then
          dir2='${MY_ROOT_DIR}/'$dir2
        fi
        shift
      fi

  esac
done

if [[ $cnt == 0 || $cnt > 2 ]]; then
  echo "You must specify 1 or 2 source directories"
  exit 1
fi

if [[ $infile == "" ]]; then
  echo "You must specify the ocean.in file with the -in flag"
  exit 1
fi

# Determine build script options.

if [[ $ncpus != "" ]]; then
  bld_opts=$ncpus
else
  bld_opts="-j 4"
fi

#************************************************************************
#***********************   Build script template   **********************
#************************************************************************

# Copy build script for manipulation

cp -p build.sh build_tpl.sh

# These lines determine what MY_CPP_FLAGS are already set in the build script

my_cpp=`grep -c '^\s*[^#]\s*setenv MY_CPP_FLAGS' build_tpl.sh`
dbg=`grep -c '^\s*[^#].*-DDEBUGGING' build_tpl.sh`
doub=`grep -c '^\s*[^#].*-DOUT_DOUBLE' build_tpl.sh`
p0=`grep -c '^\s*[^#].*-DPOSITIVE_ZERO' build_tpl.sh`

cpps=''

# If there are no uncommented setenv MY_CPP_FLAGS lines use export to
# set them outside of the build script. These values will persist for
# all model comparison builds.

if [[ $my_cpp == 0 ]]; then
  export MY_CPP_FLAGS="-DDEBUGGING -DOUT_DOUBLE -DPOSITIVE_ZERO"

# If there are active MY_CPP_FLAGS lines assemble debugging/comparison flags

else

  # if -DDEBUGGING is not set add to $cpps
  if [[ $dbg == 0 ]]; then
    cpps='
 setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DDEBUGGING"'
  fi

  # if -DOUT_DOUBLE is not set add to $cpps
  if [[ $doub == 0 ]]; then
    cpps=${cpps}'
 setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DOUT_DOUBLE"'
  fi

  # if -DPOSITIVE_ZERO is not set add to $cpps
  if [[ $p0 == 0 ]]; then
    cpps=${cpps}'
 setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DPOSITIVE_ZERO"'
  fi

  # Set proper CPP flags and setup config section
  perl -p0777 -i -e 's|^\s*[^#]\s*(setenv MY_CPP_FLAGS)(.*)$| $1$2\n$ENV{cpps}\n|m' build_tpl.sh
  perl -p0777 -i -e 's/^(#| |)setenv USE_MPI(.*|\n)*setenv USE_MY_LIBS(.*)$/#####STARTCONFIG#####\n#####ENDCONFIG#####/m' build_tpl.sh
fi

#************************************************************************
#****************************   $dir1 RUNS   ****************************
#************************************************************************

#********  Serial runs with $dir1  ********

# determine first letter of $run_pfx and save for tests at the end.

pre1=`echo ${dir1} | perl -pe 's|.*/(?!.*/)([a-zA-Z]).*|$1|'`

# Set directory and executable identifier.

run_pfx="${pre1}s"

# Configure build script

src=$dir1
config=$serial_config

# Set ROMS source code location and configuration

perl -p0777 -i -e 's|^[^#](setenv\s+MY_ROMS_SRC\s+).*$| $1$ENV{src}|m;s/^#####STARTCONFIG#####(.*|\n)*#####ENDCONFIG#####$/$ENV{config}/m' build_tpl.sh

# Build $dir1 in serial then rename the resulting executable.

if [ $build -eq 0 ]; then
  if [ -f "oceanS_${run_pfx}" ]; then
    echo "oceanS_${run_pfx} will not be rebuilt"
  else
    echo "oceanS_${run_pfx} does not exist so build is forced"
    echo -n "Building ${dir1} in serial . . . "
    ./build_tpl.sh ${bld_opts} &> build_${run_pfx}.log
    mv oceanS oceanS_${run_pfx}
    echo "Done."
  fi
else
  echo -n "Building ${dir1} in serial . . . "
  ./build_tpl.sh ${bld_opts} &> build_${run_pfx}.log
  mv oceanS oceanS_${run_pfx}
  echo "Done."
fi

# Set tiles and the directory the output files will be moved to.

odir="${run_pfx}1"
i=1
j=1
sed -i '' -e "s|NtileI == [0-9+]|NtileI == $i|g" ${infile}
sed -i '' -e "s|NtileJ == [0-9+]|NtileJ == $j|g" ${infile}

# Run application.

echo -n "   Running ${i}x${j} . . . "
./oceanS_${run_pfx} < ${infile} &> ${log}
echo "Done."

# Move model output to the correct directory, creating directory if necessary.

echo -n "   Moving files to ${odir} . . . "
if [ ! -d ${odir} ]; then
  mkdir ${odir}
fi
mv -f *.nc ${log} ${odir}
echo "Done."
echo

# Set tiles and the directory the output files will be moved to.

odir="${run_pfx}2"
i=2
j=2
sed -i '' -e "s|NtileI == [0-9+]|NtileI == $i|g" ${infile}
sed -i '' -e "s|NtileJ == [0-9+]|NtileJ == $j|g" ${infile}

# Run test case

echo -n "   Running ${i}x${j} . . . "
./oceanS_${run_pfx} < ${infile} &> ${log}
echo "Done."

# Move model output to the correct directory, creating directory if necessary.

echo -n "   Moving files to ${odir} . . . "
if [ ! -d ${odir} ]; then
  mkdir ${odir}
fi
mv -f *.nc ${log} ${odir}
echo "Done."
echo

# Set tiles and the directory the output files will be moved to.

odir="${run_pfx}3"
i=3
j=3
sed -i '' -e "s|NtileI == [0-9+]|NtileI == $i|g" ${infile}
sed -i '' -e "s|NtileJ == [0-9+]|NtileJ == $j|g" ${infile}

# Run application.

echo -n "   Running ${i}x${j} . . . "
./oceanS_${run_pfx} < ${infile} &> ${log}
echo "Done."

# Move model output to the correct directory, creating directory if necessary.

echo -n "   Moving files to ${odir} . . . "
if [ ! -d ${odir} ]; then
  mkdir ${odir}
fi
mv -f *.nc ${log} ${odir}
echo "Done."
echo


#********  MPI runs with $dir1  ********

# Set directory and executable identifier.

run_pfx="${pre1}m"

# Configure build script.

config=$mpi_config

perl -p0777 -i -e 's/^#####STARTCONFIG#####(.*|\n)*#####ENDCONFIG#####$/$ENV{config}/m' build_tpl.sh

# Build $dir1 with MPI then rename the resulting executable.

if [ $build -eq 0 ]; then
  if [ -f "oceanM_${run_pfx}" ]; then
    echo "oceanM_${run_pfx} will not be rebuilt"
  else
    echo "oceanM_${run_pfx} does not exist so build is forced"
    echo -n "Building ${dir1} with MPI . . . "
    ./build_tpl.sh ${bld_opts} &> build_${run_pfx}.log
    mv oceanM oceanM_${run_pfx}
    echo "Done."
  fi
else
  echo -n "Building ${dir1} with MPI . . . "
  ./build_tpl.sh ${bld_opts} &> build_${run_pfx}.log
  mv oceanM oceanM_${run_pfx}
  echo "Done."
fi

# Set tiles and the directory the output files will be moved to.

odir="${run_pfx}1"
i=1
j=1
proc=$[i*j]
sed -i '' -e "s|NtileI == [0-9+]|NtileI == $i|g" ${infile}
sed -i '' -e "s|NtileJ == [0-9+]|NtileJ == $j|g" ${infile}

# Run application.

echo -n "   Running ${i}x${j} . . . "
mpirunO -np $proc oceanM_${run_pfx} ${infile} &> ${log}
echo "Done."

# Move model output to the correct directory, creating directory if necessary.

echo -n "   Moving files to ${odir} . . . "
if [ ! -d ${odir} ]; then
  mkdir ${odir}
fi
mv -f *.nc ${log} ${odir}
echo "Done."
echo

# Set tiles and the directory the output files will be moved to.

odir="${run_pfx}2"
i=2
j=2
proc=$[i*j]
sed -i '' -e "s|NtileI == [0-9+]|NtileI == $i|g" ${infile}
sed -i '' -e "s|NtileJ == [0-9+]|NtileJ == $j|g" ${infile}

# Run application.

echo -n "   Running ${i}x${j} . . . "
mpirunO -np $proc oceanM_${run_pfx} ${infile} &> ${log}
echo "Done."

# Move model output to the correct directory, creating directory if necessary.

echo -n "   Moving files to ${odir} . . . "
if [ ! -d ${odir} ]; then
  mkdir ${odir}
fi
mv -f *.nc ${log} ${odir}
echo "Done."
echo

# Set tiles and the directory the output files will be moved to.

odir="${run_pfx}3"
i=3
j=3
proc=$[i*j]
sed -i '' -e "s|NtileI == [0-9+]|NtileI == $i|g" ${infile}
sed -i '' -e "s|NtileJ == [0-9+]|NtileJ == $j|g" ${infile}

# Run application.

echo -n "   Running ${i}x${j} . . . "
mpirunO -np $proc oceanM_${run_pfx} ${infile} &> ${log}
echo "Done."

# Move model output to the correct directory, creating directory if necessary.

echo -n "   Moving files to ${odir} . . . "
if [ ! -d ${odir} ]; then
  mkdir ${odir}
fi
mv -f *.nc ${log} ${odir}
echo "Done."
echo

#********  OpenMP runs with $dir1  ********

# Set directory and executable identifier.

run_pfx="${pre1}o"

# Configure build script.

config=$openmp_config

perl -p0777 -i -e 's/^#####STARTCONFIG#####(.*|\n)*#####ENDCONFIG#####$/$ENV{config}/m' build_tpl.sh

# Build $dir1 with OpenMP then rename the resulting executable.

if [ $build -eq 0 ]; then
  if [ -f "oceanO_${run_pfx}" ]; then
    echo "oceanO_${run_pfx} will not be rebuilt"
  else
    echo "oceanO_${run_pfx} does not exist so build is forced"
    echo -n "Building ${dir1} with OpenMP . . . "
    ./build_tpl.sh ${bld_opts} &> build_${run_pfx}.log
    mv oceanO oceanO_${run_pfx}
    echo "Done."
  fi
else
  echo -n "Building ${dir1} with OpenMP . . . "
  ./build_tpl.sh ${bld_opts} &> build_${run_pfx}.log
  mv oceanO oceanO_${run_pfx}
  echo "Done."
fi

# Set tiles and the directory the output files will be moved to
# and set number of OpenMP Threads.

odir="${run_pfx}1"
i=1
j=1
proc=$[i*j]
export OMP_NUM_THREADS=${proc}
sed -i '' -e "s|NtileI == [0-9+]|NtileI == $i|g" ${infile}
sed -i '' -e "s|NtileJ == [0-9+]|NtileJ == $j|g" ${infile}

# Run application.

echo -n "   Running ${i}x${j} . . . "
./oceanO_${run_pfx} < ${infile} &> ${log}
echo "Done."

# Move model output to the correct directory, creating directory if necessary.

echo -n "      Moving files to ${odir} . . . "
if [ ! -d ${odir} ]; then
  mkdir ${odir}
fi
mv -f *.nc ${log} ${odir}
echo "Done."
echo

# Set tiles and the directory the output files will be moved to
# and set number of OpenMP Threads.

odir="${run_pfx}2"
i=2
j=2
proc=$[i*j]
export OMP_NUM_THREADS=${proc}
sed -i '' -e "s|NtileI == [0-9+]|NtileI == $i|g" ${infile}
sed -i '' -e "s|NtileJ == [0-9+]|NtileJ == $j|g" ${infile}

# Run application.

echo -n "   Running ${i}x${j} . . . "
./oceanO_${run_pfx} < ${infile} &> ${log}
echo "Done."

# Move model output to the correct directory, creating directory if necessary.

echo -n "      Moving files to ${odir} . . . "
if [ ! -d ${odir} ]; then
  mkdir ${odir}
fi
mv -f *.nc ${log} ${odir}
echo "Done."
echo

# Set tiles and the directory the output files will be moved to
# and set number of OpenMP Threads.

odir="${run_pfx}3"
i=3
j=3
proc=$[i*j]
export OMP_NUM_THREADS=${proc}
sed -i '' -e "s|NtileI == [0-9+]|NtileI == $i|g" ${infile}
sed -i '' -e "s|NtileJ == [0-9+]|NtileJ == $j|g" ${infile}

# Run application.

echo -n "   Running ${i}x${j} . . . "
./oceanO_${run_pfx} < ${infile} &> ${log}
echo "Done."

# Move model output to the correct directory, creating directory if necessary.

echo -n "      Moving files to ${odir} . . . "
if [ ! -d ${odir} ]; then
  mkdir ${odir}
fi
mv -f *.nc ${log} ${odir}
echo "Done."
echo

#************************************************************************
#****************************   $dir2 RUNS   ****************************
#************************************************************************

if [[ $dir2 != "" ]]; then

  #********  Serial runs with $dir2  ********

  # determine first letter of $run_pfx and save for tests at the end.

  pre2=`echo ${dir2} | perl -pe 's|.*/(?!.*/)([a-zA-Z]).*|$1|'`

  # Set directory and executable identifier.

  run_pfx="${pre2}s"

  # Configure build script.

  src=$dir2
  config=$serial_config

  # Set ROMS source code location and configuration
  perl -p0777 -i -e 's|^[^#](setenv\s+MY_ROMS_SRC\s+).*$| $1$ENV{src}|m;s/^#####STARTCONFIG#####(.*|\n)*#####ENDCONFIG#####$/$ENV{config}/m' build_tpl.sh

  # Build $dir2 in serial then rename the resulting executable.

  if [ $build -eq 0 ]; then
    if [ -f "oceanS_${run_pfx}" ]; then
      echo "oceanS_${run_pfx} will not be rebuilt"
    else
      echo "oceanS_${run_pfx} does not exist so build is forced"
      echo -n "Building ${dir2} in serial . . . "
      ./build_tpl.sh ${bld_opts} &> build_${run_pfx}.log
      mv oceanS oceanS_${run_pfx}
      echo "Done."
    fi
  else
    echo -n "Building ${dir2} in serial . . . "
    ./build_tpl.sh ${bld_opts} &> build_${run_pfx}.log
    mv oceanS oceanS_${run_pfx}
    echo "Done."
  fi

  # Set tiles and the directory the output files will be moved to.

  odir="${run_pfx}1"
  i=1
  j=1
  sed -i '' -e "s|NtileI == [0-9+]|NtileI == $i|g" ${infile}
  sed -i '' -e "s|NtileJ == [0-9+]|NtileJ == $j|g" ${infile}

  # Run application.

  echo -n "   Running ${i}x${j} . . . "
  ./oceanS_${run_pfx} < ${infile} &> ${log}
  echo "Done."

  # Move model output to the correct directory, creating directory if necessary.

  echo -n "   Moving files to ${odir} . . . "
  if [ ! -d ${odir} ]; then
    mkdir ${odir}
  fi
  mv -f *.nc ${log} ${odir}
  echo "Done."
  echo

  # Set tiles and the directory the output files will be moved to.

  odir="${run_pfx}2"
  i=2
  j=2
  sed -i '' -e "s|NtileI == [0-9+]|NtileI == $i|g" ${infile}
  sed -i '' -e "s|NtileJ == [0-9+]|NtileJ == $j|g" ${infile}

  # Run application.

  echo -n "   Running ${i}x${j} . . . "
  ./oceanS_${run_pfx} < ${infile} &> ${log}
  echo "Done."

  # Move model output to the correct directory, creating directory if necessary.

  echo -n "   Moving files to ${odir} . . . "
  if [ ! -d ${odir} ]; then
    mkdir ${odir}
  fi
  mv -f *.nc ${log} ${odir}
  echo "Done."
  echo

  # Set tiles and the directory the output files will be moved to.

  odir="${run_pfx}3"
  i=3
  j=3
  sed -i '' -e "s|NtileI == [0-9+]|NtileI == $i|g" ${infile}
  sed -i '' -e "s|NtileJ == [0-9+]|NtileJ == $j|g" ${infile}

  # Run application.

  echo -n "   Running ${i}x${j} . . . "
  ./oceanS_${run_pfx} < ${infile} &> ${log}
  echo "Done."

  # Move model output to the correct directory, creating directory if necessary.

  echo -n "   Moving files to ${odir} . . . "
  if [ ! -d ${odir} ]; then
    mkdir ${odir}
  fi
  mv -f *.nc ${log} ${odir}
  echo "Done."
  echo

  #********  MPI runs with $dir2  ********

  # Set directory and executable identifier.

  run_pfx="${pre2}m"

  # Configure build script.

  config=$mpi_config

  perl -p0777 -i -e 's/^#####STARTCONFIG#####(.*|\n)*#####ENDCONFIG#####$/$ENV{config}/m' build_tpl.sh

  # Build $dir2 with MPI then rename the resulting executable.

  if [ $build -eq 0 ]; then
    if [ -f "oceanM_${run_pfx}" ]; then
      echo "oceanM_${run_pfx} will not be rebuilt"
    else
      echo "oceanM_${run_pfx} does not exist so build is forced"
      echo -n "Building ${dir2} with MPI . . . "
      ./build_tpl.sh ${bld_opts} &> build_${run_pfx}.log
      mv oceanM oceanM_${run_pfx}
      echo "Done."
    fi
  else
    echo -n "Building ${dir2} with MPI . . . "
    ./build_tpl.sh ${bld_opts} &> build_${run_pfx}.log
    mv oceanM oceanM_${run_pfx}
    echo "Done."
  fi

  # Set tiles and the directory the output files will be moved to.

  odir="${run_pfx}1"
  i=1
  j=1
  proc=$[i*j]
  sed -i '' -e "s|NtileI == [0-9+]|NtileI == $i|g" ${infile}
  sed -i '' -e "s|NtileJ == [0-9+]|NtileJ == $j|g" ${infile}

  # Run application.

  echo -n "   Running ${i}x${j} . . . "
  mpirunO -np $proc oceanM_${run_pfx} ${infile} &> ${log}
  echo "Done."

  # Move model output to the correct directory, creating directory if necessary.

  echo -n "   Moving files to ${odir} . . . "
  if [ ! -d ${odir} ]; then
    mkdir ${odir}
  fi
  mv -f *.nc ${log} ${odir}
  echo "Done."
  echo

  # Set tiles and the directory the output files will be moved to.

  odir="${run_pfx}2"
  i=2
  j=2
  proc=$[i*j]
  sed -i '' -e "s|NtileI == [0-9+]|NtileI == $i|g" ${infile}
  sed -i '' -e "s|NtileJ == [0-9+]|NtileJ == $j|g" ${infile}

  # Run application.

  echo -n "   Running ${i}x${j} . . . "
  mpirunO -np $proc oceanM_${run_pfx} ${infile} &> ${log}
  echo "Done."

  # Move model output to the correct directory, creating directory if necessary.

  echo -n "   Moving files to ${odir} . . . "
  if [ ! -d ${odir} ]; then
    mkdir ${odir}
  fi
  mv -f *.nc ${log} ${odir}
  echo "Done."
  echo

  # Set tiles and the directory the output files will be moved to

  odir="${run_pfx}3"
  i=3
  j=3
  proc=$[i*j]
  sed -i '' -e "s|NtileI == [0-9+]|NtileI == $i|g" ${infile}
  sed -i '' -e "s|NtileJ == [0-9+]|NtileJ == $j|g" ${infile}

  # Run application.

  echo -n "   Running ${i}x${j} . . . "
  mpirunO -np $proc oceanM_${run_pfx} ${infile} &> ${log}
  echo "Done."

  # Move model output to the correct directory, creating directory if necessary.

  echo -n "   Moving files to ${odir} . . . "
  if [ ! -d ${odir} ]; then
    mkdir ${odir}
  fi
  mv -f *.nc ${log} ${odir}
  echo "Done."
  echo

  #********  OpenMP runs with $dir2  ********

  # Set directory and executable identifier.

  run_pfx="${pre2}o"

  # Configure build script.

  config=$openmp_config

  perl -p0777 -i -e 's/^#####STARTCONFIG#####(.*|\n)*#####ENDCONFIG#####$/$ENV{config}/m' build_tpl.sh

  # Build $dir2 with OpenMP then rename the resulting executable.

  if [ $build -eq 0 ]; then
    if [ -f "oceanO_${run_pfx}" ]; then
      echo "oceanO_${run_pfx} will not be rebuilt"
    else
      echo "oceanO_${run_pfx} does not exist so build is forced"
      echo -n "Building ${dir2} with OpenMP . . . "
      ./build_tpl.sh ${bld_opts} &> build_${run_pfx}.log
      mv oceanO oceanO_${run_pfx}
      echo "Done."
    fi
  else
    echo -n "Building ${dir2} with OpenMP . . . "
    ./build_tpl.sh ${bld_opts} &> build_${run_pfx}.log
    mv oceanO oceanO_${run_pfx}
    echo "Done."
  fi

  # Set tiles and the directory the output files will be moved to
  # and set number of OpenMP Threads.

  odir="${run_pfx}1"
  i=1
  j=1
  proc=$[i*j]
  export OMP_NUM_THREADS=${proc}
  sed -i '' -e "s|NtileI == [0-9+]|NtileI == $i|g" ${infile}
  sed -i '' -e "s|NtileJ == [0-9+]|NtileJ == $j|g" ${infile}

  # Run application.

  echo -n "   Running ${i}x${j} . . . "
  ./oceanO_${run_pfx} < ${infile} &> ${log}
  echo "Done."

  # Move model output to the correct directory, creating directory if necessary.

  echo -n "   i   Moving files to ${odir} . . . "
  if [ ! -d ${odir} ]; then
    mkdir ${odir}
  fi
  mv -f *.nc ${log} ${odir}
  echo "Done."
  echo

  # Set tiles and the directory the output files will be moved to
  # and set number of OpenMP Threads.

  odir="${run_pfx}2"
  i=2
  j=2
  proc=$[i*j]
  export OMP_NUM_THREADS=${proc}
  sed -i '' -e "s|NtileI == [0-9+]|NtileI == $i|g" ${infile}
  sed -i '' -e "s|NtileJ == [0-9+]|NtileJ == $j|g" ${infile}

  # Run application.

  echo -n "   Running ${i}x${j} . . . "
  ./oceanO_${run_pfx} < ${infile} &> ${log}
  echo "Done."

  # Move model output to the correct directory, creating directory if necessary.

  echo -n "      Moving files to ${odir} . . . "
  if [ ! -d ${odir} ]; then
    mkdir ${odir}
  fi
  mv -f *.nc ${log} ${odir}
  echo "Done."
  echo

  # Set tiles and the directory the output files will be moved to
  # and set number of OpenMP Threads.

  odir="${run_pfx}3"
  i=3
  j=3
  proc=$[i*j]
  export OMP_NUM_THREADS=${proc}
  sed -i '' -e "s|NtileI == [0-9+]|NtileI == $i|g" ${infile}
  sed -i '' -e "s|NtileJ == [0-9+]|NtileJ == $j|g" ${infile}

  # Run application.

  echo -n "   Running ${i}x${j} . . . "
  ./oceanO_${run_pfx} < ${infile} &> ${log}
  echo "Done."

  # Move model output to the correct directory, creating directory if necessary.

  echo -n "      Moving files to ${odir} . . . "
  if [ ! -d ${odir} ]; then
    mkdir ${odir}
  fi
  mv -f *.nc ${log} ${odir}
  echo "Done."
  echo
fi # if $dir2 is specified run this block

#************************************************************************
#************************************************************************
#*******************                               **********************
#*****************        Result Comparisons         ********************
#*******************                               **********************
#************************************************************************
#************************************************************************

# Set the difference counter to zero

diffs=0

# Check each configuration against the <src_dir1>, serial, single
# tile configuration and count the sum diferences

# NOTE 1: The check_nc.sh script will report any differing files to
#         standard out as it finds them $diffs is used to sum the total
# number of differing files.

# NOTE 2: $? holds the exit code for the previous command (check_nc.sh
#         arg1 arg2 in this case) I am manually setting the exit code in
# check_nc.sh to be the total number of differing files between dir1 and
# dir2.  See check_nc.sh for details.

check_nc.sh ${pre1}s1 ${pre1}s2
diffs=`expr ${diffs} + $?`
check_nc.sh ${pre1}s1 ${pre1}s3
diffs=`expr ${diffs} + $?`
check_nc.sh ${pre1}s1 ${pre1}m1
diffs=`expr ${diffs} + $?`
check_nc.sh ${pre1}s1 ${pre1}m2
diffs=`expr ${diffs} + $?`
check_nc.sh ${pre1}s1 ${pre1}m3
diffs=`expr ${diffs} + $?`
check_nc.sh ${pre1}s1 ${pre1}o1
diffs=`expr ${diffs} + $?`
check_nc.sh ${pre1}s1 ${pre1}o2
diffs=`expr ${diffs} + $?`
check_nc.sh ${pre1}s1 ${pre1}o3
diffs=`expr ${diffs} + $?`

# If directory 2 was specified compare them as well.

if [[ $dir2 != "" ]]; then
  check_nc.sh ${pre1}s1 ${pre2}s1
  diffs=`expr ${diffs} + $?`
  check_nc.sh ${pre1}s1 ${pre2}s2
  diffs=`expr ${diffs} + $?`
  check_nc.sh ${pre1}s1 ${pre2}s3
  diffs=`expr ${diffs} + $?`
  check_nc.sh ${pre1}s1 ${pre2}m1
  diffs=`expr ${diffs} + $?`
  check_nc.sh ${pre1}s1 ${pre2}m2
  diffs=`expr ${diffs} + $?`
  check_nc.sh ${pre1}s1 ${pre2}m3
  diffs=`expr ${diffs} + $?`
  check_nc.sh ${pre1}s1 ${pre2}o1
  diffs=`expr ${diffs} + $?`
  check_nc.sh ${pre1}s1 ${pre2}o2
  diffs=`expr ${diffs} + $?`
  check_nc.sh ${pre1}s1 ${pre2}o3
  diffs=`expr ${diffs} + $?`
fi

# Output the total number of differing files. The actual files that
# differ, if any, will be displayed in standard out as they are found

echo "There were ${diffs} differing files found"

# clean temporary files.

/bin/rm build_tpl.sh

if [[ $diffs == 0 ]]; then
  /bin/rm build_*.log

  # if cleanup flag set remove all traces that test were run
  if [[ $cleanup == 1 ]]; then
    echo -n "  Cleaning up . . . "
    /bin/rm -rf Build ${pre1}s1 ${pre1}s2 ${pre1}s3 ${pre1}m1 ${pre1}m2 ${pre1}m3 ${pre1}o1 ${pre1}o2 ${pre1}o3
    /bin/rm oceanS_${pre1}s oceanM_${pre1}m oceanO_${pre1}o

    if [[ $dir2 != "" ]]; then
      /bin/rm -rf ${pre2}s1 ${pre2}s2 ${pre2}s3 ${pre2}m1 ${pre2}m2 ${pre2}m3 ${pre2}o1 ${pre2}o2 ${pre2}o3
      /bin/rm oceanS_${pre2}s oceanM_${pre2}m oceanO_${pre2}o
    fi
    echo "Done."
    echo
  fi
fi
