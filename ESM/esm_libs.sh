#!/bin/bash
#
# git $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2026 The ROMS Group                                :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# ESM Component Configuration BASH Script: sourced from build_roms.sh   :::
#                                                                       :::
# The strategy is to compile and link each ESM component separately     :::
# first, and then ROMS since it is driving the coupled system. Only     :::
# the ESM components activated are considered and the rest are ignored. :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

echo ""
echo "${separator}"
echo "Using ESM library paths from:  $1"
echo "${separator}"
echo ""

if [ -n "${USE_DEBUG:+1}" ]; then
  export     CICE_LIB_DIR=${MY_PROJECT_DIR}/Build_ciceG
  export   COAMPS_LIB_DIR=${MY_PROJECT_DIR}/Build_coampsG
  export    REGCM_LIB_DIR=${MY_PROJECT_DIR}/Build_regcmG
  export      WAM_LIB_DIR=${MY_PROJECT_DIR}/Build_wamG
  export      WRF_LIB_DIR=${MY_PROJECT_DIR}/Build_wrfG
else
  export     CICE_LIB_DIR=${MY_PROJECT_DIR}/Build_cice
  export   COAMPS_LIB_DIR=${MY_PROJECT_DIR}/Build_coamps
  export    REGCM_LIB_DIR=${MY_PROJECT_DIR}/Build_regcm
  export      WAM_LIB_DIR=${MY_PROJECT_DIR}/Build_wam
  export      WRF_LIB_DIR=${MY_PROJECT_DIR}/Build_wrf
fi
