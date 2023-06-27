#!/bin/bash
#
# git $Id$
# svn $Id: wrf_patch.sh 1151 2023-02-09 03:08:53Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# WRF Patching BASH Script: WRF Versions 4.1 and up                     :::
#                                                                       :::
# Script called from "build_wrf.sh" to check whether WRF source code    :::
# has been patched for NetCDF4 library depencies, added configure       :::
# options, creating clean .f90 files for debugging, rename modules to   :::
# WRF_ESMF_*,correct optional argment from defaultCalendar to           :::
# defaultCalKind in ESMF_Initialize call, and remove limits on the      :::
# moisture and latent heat flux.                                        :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    This script is called from "build_wrf.sh"                          :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Initialize.

# ${WRF_ROOT_DIR} and ${ROMS_SRC_DIR} are environment variables set in
# "build_wrf.sh".

separator=`perl -e "print ':' x 100;"`
CHECK_STRING='WRF-ROMS ESMF-NUOPC Coupling'

echo ""
echo "${separator}"
echo "If applicable, replacing several WRF files for WRF ESMF/NUOPC Coupling"
echo "${separator}"
echo ""

# Reworking linking NetCDF4 library dependencies

if [ `grep -c "${CHECK_STRING}" ${WRF_ROOT_DIR}/configure` -eq "0" ]; then
  mv -v ${WRF_ROOT_DIR}/configure ${WRF_ROOT_DIR}/configure.orig
  cp -fv ${ROMS_SRC_DIR}/ESM/wrf_configure ${WRF_ROOT_DIR}/configure
else
  echo "   No need to replace: ${WRF_ROOT_DIR}/configure"
fi

if [ `grep -c "${CHECK_STRING}" ${WRF_ROOT_DIR}/Makefile` -eq "0" ]; then
  mv -v  ${WRF_ROOT_DIR}/Makefile ${WRF_ROOT_DIR}/Makefile.orig
  cp -fv ${ROMS_SRC_DIR}/ESM/wrf_Makefile  ${WRF_ROOT_DIR}/Makefile
else
  echo "   No need to replace: ${WRF_ROOT_DIR}/Makefile"
fi

if [ `grep -c "${CHECK_STRING}" ${WRF_ROOT_DIR}/arch/postamble` -eq "0" ]; then
  mv -v  ${WRF_ROOT_DIR}/arch/postamble ${WRF_ROOT_DIR}/arch/postamble.orig
  cp -fv ${ROMS_SRC_DIR}/ESM/wrf_postamble ${WRF_ROOT_DIR}/arch/postamble
else
  echo "   No need to replace: ${WRF_ROOT_DIR}/arch/postamble"  else
fi

# Adding MacOS Intel/GNU with OpenMPI

if [ `grep -c "${CHECK_STRING}" ${WRF_ROOT_DIR}/arch/configure.defaults` -eq "0" ]; then
  mv -v  ${WRF_ROOT_DIR}/arch/configure.defaults ${WRF_ROOT_DIR}/arch/configure.defaults.orig
  cp -fv ${ROMS_SRC_DIR}/ESM/wrf_configure.defaults ${WRF_ROOT_DIR}/arch/configure.defaults
else
  echo "   No need to replace: ${WRF_ROOT_DIR}/arch/configure.defaults"
fi

# Create clean .f90 files for debugging and rename modules to WRF_ESMF_*

if [ `grep -c "${CHECK_STRING}" ${WRF_ROOT_DIR}/external/esmf_time_f90/Makefile` -eq "0" ]; then
  mv -v  ${WRF_ROOT_DIR}/external/esmf_time_f90/Makefile ${WRF_ROOT_DIR}/external/esmf_time_f90/Makefile.orig
  cp -fv ${ROMS_SRC_DIR}/ESM/wrf_Makefile.esmf ${WRF_ROOT_DIR}/external/esmf_time_f90/Makefile
else
  echo "   No need to replace: ${WRF_ROOT_DIR}/external/esmf_time_f90/Makefile"
fi

# Correcting optional argument from defaultCalendar to defaultCalKind in
# ESMF_Initialize call

if [ `grep -c "${CHECK_STRING}" ${WRF_ROOT_DIR}/external/esmf_time_f90/Test1.F90` -eq "0" ]; then
  mv -v  ${WRF_ROOT_DIR}/external/esmf_time_f90/Test1.F90 ${WRF_ROOT_DIR}/external/esmf_time_f90/Test1.F90.orig
  cp -fv ${ROMS_SRC_DIR}/ESM/wrf_Test1.F90 ${WRF_ROOT_DIR}/external/esmf_time_f90/Test1.F90
else
  echo "   No need to replace: ${WRF_ROOT_DIR}/external/esmf_time_f90/Test1.F90"
fi

# We are not sure why moisture flux and sensible heat flux are limited in
# phys/module_sf_sfclay.F and phys/module_sf_sfclayrev.F, but it causes
# incorrect results in our hurricane, upwelling and sea-breeze simiulations.
#
# The -i.orig option saves the original file before performing the replacement.

if [ ! -f ${WRF_ROOT_DIR}/phys/module_sf_sfclay.F.orig ]; then
  echo "Patching ${WRF_ROOT_DIR}/phys/module_sf_sfclay.F"
  perl -i.orig -pe 's/^(\s*)\s([QH]FX\(I\)=AMAX1)/\!$1$2/' ${WRF_ROOT_DIR}/phys/module_sf_sfclay.F
else
  echo "   No need to patch:   ${WRF_ROOT_DIR}/phys/module_sf_sfclay.F"
fi
if [ ! -f ${WRF_ROOT_DIR}/phys/module_sf_sfclayrev.F.orig ]; then
  echo "Patching ${WRF_ROOT_DIR}/phys/module_sf_sfclayrev.F"
  perl -i.orig -pe 's/^(\s*)\s([QH]FX\(I\)=AMAX1)/\!$1$2/' ${WRF_ROOT_DIR}/phys/module_sf_sfclayrev.F
else
  echo "   No need to patch:   ${WRF_ROOT_DIR}/phys/module_sf_sfclayrev.F"
fi

