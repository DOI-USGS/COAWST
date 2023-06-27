#!/bin/bash
#
# git $Id$
# svn $Id: wrf_move.sh 1151 2023-02-09 03:08:53Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# WRF objects move BASH Script: WRF Versions 4.1 and up                 :::
#                                                                       :::
# Script called from "build_wrf.sh" if -move flag is set. This script   :::
# moves the WRF objects and executables needed to run WRF in the        :::
# coupled ESMF/NUOPC system.                                            :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    This script is called from "build_wrf.sh".                         :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Initialize.

# ${WRF_ROOT_DIR}, ${WRF_BIN_DIR}, and ${WRF_BUILD_DIR} are environment
# variables set in "build_wrf.sh".

separator=`perl -e "print ':' x 100;"`

echo ""
echo "${separator}"
echo "Moving WRF objects to Build directory  ${WRF_BUILD_DIR}:"
echo "${separator}"
echo ""

if [ ! -d ${WRF_BUILD_DIR} ]; then
  /bin/mkdir -pv ${WRF_BUILD_DIR}
  /bin/mkdir -pv ${WRF_BUILD_DIR}/Bin
  echo ""
fi

/bin/cp -pfv configure.wrf ${WRF_BUILD_DIR}
/bin/cp -pfv Registry/Registry ${WRF_BUILD_DIR}
/bin/cp -pfv Registry/io_boilerplate_temporary.inc ${WRF_BUILD_DIR}

/bin/mv -fv run/namelist.input ${WRF_BUILD_DIR}
/bin/mv -fv test/em_real/README.namelist ${WRF_BUILD_DIR}

/bin/mv -fv frame/md_calls.inc ${WRF_BUILD_DIR}
/bin/mv -fv frame/module_dm.F ${WRF_BUILD_DIR}
/bin/mv -fv frame/module_state_description.F ${WRF_BUILD_DIR}
/bin/mv -fv external/io_grib1/io_grib1.f90 ${WRF_BUILD_DIR}
/bin/mv -fv external/io_grib_share/io_grib_share.f90 ${WRF_BUILD_DIR}

/bin/mv -fv tools/gen_comms.c ${WRF_BUILD_DIR}

/bin/mv -fv inc/dm_comm_cpp_flags ${WRF_BUILD_DIR}
/bin/mv -fv inc/wrf_io_flags.h ${WRF_BUILD_DIR}
/bin/mv -fv inc/wrf_status_codes.h ${WRF_BUILD_DIR}

/bin/mv -fv external/io_int/io_int_idx_tags.h ${WRF_BUILD_DIR}

# The "arch", "chem",  "run" and "var" directories have source code files
# (*.inc, *.F90, *.f90) that are part of repository. We cannot use the
# compact "find" function for these file extensions. If WRF_ROOT_DIR is
# NOT a real path but a link to the disk location, we need to dereference
# just that link (-H option) but leave the "found" links alone for "find"
# to work correctly.

for DIR in frame chem share dyn_em dyn_exp dyn_nmm phys main tools wrftladj
do
  if [ -d $DIR ]; then
    /bin/mv -fv ${DIR}/*.f90 ${WRF_BUILD_DIR}
    /bin/mv -fv ${DIR}/*.F90 ${WRF_BUILD_DIR}
    /bin/mv -fv ${DIR}/*.inc ${WRF_BUILD_DIR}
  fi
done

find -H ${WRF_ROOT_DIR} -type f -name "*.mod" -exec /bin/mv -fv {} ${WRF_BUILD_DIR} \;
find -H ${WRF_ROOT_DIR} -type f -name "*.o"   -exec /bin/mv -fv {} ${WRF_BUILD_DIR} \;
find -H ${WRF_ROOT_DIR} -type f -name "*.obj" -exec /bin/mv -fv {} ${WRF_BUILD_DIR} \;
find -H ${WRF_ROOT_DIR} -type f -name "*.a"   -exec /bin/mv -fv {} ${WRF_BUILD_DIR} \;

find -H ${WRF_ROOT_DIR} -type f -name "*.exe" -exec /bin/mv -fv {} ${WRF_BIN_DIR} \;
find -H ${WRF_ROOT_DIR} -type l -name "*.exe" -exec /bin/rm -fv {} \;

/bin/mv -fv external/esmf_time_f90/*.f ${WRF_BUILD_DIR}
/bin/cp -pv external/esmf_time_f90/ESMF*.inc ${WRF_BUILD_DIR}

/bin/mv -fv external/io_int/diffwrf ${WRF_BIN_DIR}/diffwrf_int
/bin/mv -fv external/io_int/test_io_idx ${WRF_BIN_DIR}
/bin/mv -fv external/io_netcdf/diffwrf ${WRF_BIN_DIR}/diffwrf_nc

/bin/mv -fv tools/fseeko_test ${WRF_BIN_DIR}
/bin/mv -fv tools/registry ${WRF_BIN_DIR}
