#!/bin/csh -f
#
# git $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2026 The ROMS Group                                :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# ESM Component Configuration BASH Script: sourced from build_roms.csh  :::
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

if ($?USE_DEBUG) then
  setenv CICE_LIB_DIR      ${MY_PROJECT_DIR}/Build_ciceG
  setenv COAMPS_LIB_DIR    ${MY_PROJECT_DIR}/Build_coampsG
  setenv REGCM_LIB_DIR     ${MY_PROJECT_DIR}/Build_regcmG
  setenv WAM_LIB_DIR       ${MY_PROJECT_DIR}/Build_wamG
  setenv WRF_LIB_DIR       ${MY_PROJECT_DIR}/Build_wrfG
else
  setenv CICE_LIB_DIR      ${MY_PROJECT_DIR}/Build_cice
  setenv COAMPS_LIB_DIR    ${MY_PROJECT_DIR}/Build_coamps
  setenv REGCM_LIB_DIR     ${MY_PROJECT_DIR}/Build_regcm
  setenv WAM_LIB_DIR       ${MY_PROJECT_DIR}/Build_wam
  setenv WRF_LIB_DIR       ${MY_PROJECT_DIR}/Build_wrf
endif
