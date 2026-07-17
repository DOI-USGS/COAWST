#!/bin/csh -f
#
# git $Id$
# svn $Id: wrf_links.csh 1151 2023-02-09 03:08:53Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# WRF em_real symlinks CSH Script: WRF Versions 4.1 and up              :::
#                                                                       :::
# Script called from "build_wrf.csh" if ${WRF_CASE} is set to           :::
# "em_real". This script creates the data links for running the         :::
# 'em_real' executable.                                                 :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    This script is called from "build_wrf.csh".                        :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Initialize.

# ${MY_PROJECT_DIR}, ${WRF_ROOT_DIR}, ${WRF_VERSION}, ${USE_REAL_DOUBLE},
# ${WRF_BUILD_DIR} are environment variables set in "build_wrf.csh".

set separator = `perl -e "print ':' x 100;"`

echo ""
echo "${separator}"
echo "Creating WRF data links:  Case em_real"
echo "${separator}"
echo ""

cd ${MY_PROJECT_DIR}

find ./ -type l -exec /bin/rm -f {} \;

ln -sfv ${WRF_ROOT_DIR}/run/ETAMPNEW_DATA .
ln -sfv ${WRF_ROOT_DIR}/run/ETAMPNEW_DATA.expanded_rain .
ln -sfv ${WRF_ROOT_DIR}/run/RRTM_DATA .
ln -sfv ${WRF_ROOT_DIR}/run/RRTMG_LW_DATA .
ln -sfv ${WRF_ROOT_DIR}/run/RRTMG_SW_DATA .
ln -sfv ${WRF_ROOT_DIR}/run/CAM_ABS_DATA .
ln -sfv ${WRF_ROOT_DIR}/run/CAM_AEROPT_DATA .
ln -sfv ${WRF_ROOT_DIR}/run/CAMtr_volume_mixing_ratio.RCP4.5 .
ln -sfv ${WRF_ROOT_DIR}/run/CAMtr_volume_mixing_ratio.RCP6 .
ln -sfv ${WRF_ROOT_DIR}/run/CAMtr_volume_mixing_ratio.RCP8.5 CAMtr_volume_mixing_ratio
ln -sfv ${WRF_ROOT_DIR}/run/CAMtr_volume_mixing_ratio.A1B .
ln -sfv ${WRF_ROOT_DIR}/run/CAMtr_volume_mixing_ratio.A2 .
ln -sfv ${WRF_ROOT_DIR}/run/CLM_ALB_ICE_DFS_DATA .
ln -sfv ${WRF_ROOT_DIR}/run/CLM_ALB_ICE_DRC_DATA .
ln -sfv ${WRF_ROOT_DIR}/run/CLM_ASM_ICE_DFS_DATA .
ln -sfv ${WRF_ROOT_DIR}/run/CLM_ASM_ICE_DRC_DATA .
ln -sfv ${WRF_ROOT_DIR}/run/CLM_DRDSDT0_DATA .
ln -sfv ${WRF_ROOT_DIR}/run/CLM_EXT_ICE_DFS_DATA .
ln -sfv ${WRF_ROOT_DIR}/run/CLM_EXT_ICE_DRC_DATA .
ln -sfv ${WRF_ROOT_DIR}/run/CLM_KAPPA_DATA .
ln -sfv ${WRF_ROOT_DIR}/run/CLM_TAU_DATA .
ln -sfv ${WRF_ROOT_DIR}/run/ozone.formatted .
ln -sfv ${WRF_ROOT_DIR}/run/ozone_lat.formatted .
ln -sfv ${WRF_ROOT_DIR}/run/ozone_plev.formatted .
ln -sfv ${WRF_ROOT_DIR}/run/aerosol.formatted .
ln -sfv ${WRF_ROOT_DIR}/run/aerosol_lat.formatted .
ln -sfv ${WRF_ROOT_DIR}/run/aerosol_lon.formatted .
ln -sfv ${WRF_ROOT_DIR}/run/aerosol_plev.formatted .
ln -sfv ${WRF_ROOT_DIR}/run/capacity.asc .
ln -sfv ${WRF_ROOT_DIR}/run/coeff_p.asc .
ln -sfv ${WRF_ROOT_DIR}/run/coeff_q.asc .
ln -sfv ${WRF_ROOT_DIR}/run/constants.asc .
ln -sfv ${WRF_ROOT_DIR}/run/masses.asc .
ln -sfv ${WRF_ROOT_DIR}/run/termvels.asc .
ln -sfv ${WRF_ROOT_DIR}/run/kernels.asc_s_0_03_0_9 .
ln -sfv ${WRF_ROOT_DIR}/run/kernels_z.asc .
ln -sfv ${WRF_ROOT_DIR}/run/bulkdens.asc_s_0_03_0_9 .
ln -sfv ${WRF_ROOT_DIR}/run/bulkradii.asc_s_0_03_0_9 .
ln -sfv ${WRF_ROOT_DIR}/run/CCN_ACTIVATE.BIN .
if ( $WRF_VERSION == "4.3" ) then
  ln -sfv ${WRF_ROOT_DIR}/run/p3_lookupTable_1.dat-2momI_v5.1.6_oldDimax .
  ln -sfv ${WRF_ROOT_DIR}/run/p3_lookupTable_1.dat-3momI_v5.1.6 .
  ln -sfv ${WRF_ROOT_DIR}/run/p3_lookupTable_2.dat-4.1 .
else
  ln -sfv ${WRF_ROOT_DIR}/run/p3_lookup_table_1.dat-v4.1 .
  ln -sfv ${WRF_ROOT_DIR}/run/p3_lookup_table_2.dat-v4.1 .
endif

if ( $?USE_REAL_DOUBLE ) then
  ln -sfv ${WRF_ROOT_DIR}/run/ETAMPNEW_DATA_DBL ETAMPNEW_DATA
  ln -sfv ${WRF_ROOT_DIR}/run/ETAMPNEW_DATA.expanded_rain_DBL ETAMPNEW_DATA.expanded_rain
  ln -sfv ${WRF_ROOT_DIR}/run/RRTM_DATA_DBL RRTM_DATA
  ln -sfv ${WRF_ROOT_DIR}/run/RRTMG_LW_DATA_DBL RRTMG_LW_DATA
  ln -sfv ${WRF_ROOT_DIR}/run/RRTMG_SW_DATA_DBL RRTMG_SW_DATA
else
  ln -sfv ${WRF_ROOT_DIR}/run/ETAMPNEW_DATA ETAMPNEW_DATA
  ln -sfv ${WRF_ROOT_DIR}/run/ETAMPNEW_DATA.expanded_rain ETAMPNEW_DATA.expanded_rain
  ln -sfv ${WRF_ROOT_DIR}/run/RRTM_DATA RRTM_DATA
  ln -sfv ${WRF_ROOT_DIR}/run/RRTMG_LW_DATA RRTMG_LW_DATA
  ln -sfv ${WRF_ROOT_DIR}/run/RRTMG_SW_DATA RRTMG_SW_DATA
endif

ln -sfv ${WRF_ROOT_DIR}/run/GENPARM.TBL .
ln -sfv ${WRF_ROOT_DIR}/run/LANDUSE.TBL .
ln -sfv ${WRF_ROOT_DIR}/run/SOILPARM.TBL .
ln -sfv ${WRF_ROOT_DIR}/run/URBPARM.TBL .
ln -sfv ${WRF_ROOT_DIR}/run/VEGPARM.TBL .
ln -sfv ${WRF_ROOT_DIR}/run/MPTABLE.TBL .
ln -sfv ${WRF_ROOT_DIR}/run/tr49t67 .
ln -sfv ${WRF_ROOT_DIR}/run/tr49t85 .
ln -sfv ${WRF_ROOT_DIR}/run/tr67t85 .
ln -sfv ${WRF_ROOT_DIR}/run/gribmap.txt .
ln -sfv ${WRF_ROOT_DIR}/run/grib2map.tbl .

# Set directory structure needed to compile WPS.

echo ""
echo "${separator}"
echo "Creating WPS directory layout."
echo "${separator}"
echo ""

/bin/mkdir -vp ${WRF_BUILD_DIR}/external/io_grib1
/bin/mkdir -vp ${WRF_BUILD_DIR}/external/io_grib_share
/bin/mkdir -vp ${WRF_BUILD_DIR}/external/io_int
/bin/mkdir -vp ${WRF_BUILD_DIR}/external/io_netcdf
/bin/mkdir -vp ${WRF_BUILD_DIR}/frame

# We do not copy files out of the inc directory so we create a symlink
# do that directory.

echo "cd ${WRF_BUILD_DIR}"
cd ${WRF_BUILD_DIR}
ln -sfv ${WRF_ROOT_DIR}/inc .

# io_grib1

echo "cd ${WRF_BUILD_DIR}/external/io_grib1"
cd ${WRF_BUILD_DIR}/external/io_grib1
ln -sfv ../../libio_grib1.a .

# io_grib_share

echo "cd ${WRF_BUILD_DIR}/external/io_grib_share"
cd ${WRF_BUILD_DIR}/external/io_grib_share
ln -sfv ../../libio_grib_share.a .
ln -sfv ../../wrf_io_flags.h     .
ln -sfv ../../wrf_status_codes.h .

# io_int

echo "cd ${WRF_BUILD_DIR}/external/io_int"
cd ${WRF_BUILD_DIR}/external/io_int
ln -sfv ../../libwrfio_int.a                  .
ln -sfv ../../module_internal_header_util.mod .

# io_netcdf

echo "cd ${WRF_BUILD_DIR}/external/io_netcdf"
cd ${WRF_BUILD_DIR}/external/io_netcdf
ln -sfv ../../libwrfio_nf.a .

# frame

echo "cd ${WRF_BUILD_DIR}/frame"
cd ${WRF_BUILD_DIR}/frame
ln -sfv ../pack_utils.o                  .
ln -sfv ../module_machine.o              .
ln -sfv ../module_internal_header_util.o .
ln -sfv ../module_driver_constants.o     .

# Remove links in WRF/test/em_real sub-directory

echo ""
echo "${separator}"
echo "Removing WRF data links from  ${WRF_ROOT_DIR}/test/em_real:"
echo "${separator}"
echo ""

# The -H flag will only dereference paths on the commandline
# (${WRF_ROOT_DIR}/test/em_real) but will allow -type l to work
# properly for all "found" files.

find -H ${WRF_ROOT_DIR}/test/em_real -type l -exec /bin/rm -fv {} \;

