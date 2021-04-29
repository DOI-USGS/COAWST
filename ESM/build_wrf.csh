#!/bin/csh -f
#
# svn $Id: build_wrf.csh 1054 2021-03-06 19:47:12Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2021 The ROMS/TOMS Group                           :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# WRF Compiling CSH Script                                              :::
#                                                                       :::
# Script to compile WRF where source code files are kept separate       :::
# from the application configuration and build objects.                 :::
#                                                                       :::
# Q: How/why does this script work?                                     :::
#                                                                       :::
# A:                                                                    :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    ./build_wrf.csh [options]                                          :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    -j [N]      Compile in parallel using N CPUs                       :::
#                  omit argument for all available CPUs                 :::
#                                                                       :::
#    -move       Move compiled objects to build directory               :::
#                                                                       :::
#    -noclean    Do not run clean -a script                             :::
#                                                                       :::
#    -noconfig   Do not run configure compilation script                :::
#                                                                       :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

setenv which_MPI openmpi                        # default, overwriten below

# Initialize.

set separator = `perl -e "print ':' x 100;"`

set parallel = 0
set clean = 1
set config = 1
set move = 0

setenv CPLFLAG ''
setenv MY_CPP_FLAGS ''

while ( ($#argv) > 0 )
  switch ($1)
    case "-j"
      shift
      set parallel = 1
      if (`echo $1 | grep '^[0-9]\+$'` != "" ) then
        set NCPUS = $1
        shift
      else
        set NCPUS = 2
      endif
    breaksw

    case "-noclean"
      shift
      set clean = 0
    breaksw

    case "-noconfig"
      shift
      set config = 0
    breaksw

    case "-move"
      shift
      set move = 1
    breaksw

    case "-*":
      echo ""
      echo "$0 : Unknown option [ $1 ]"
      echo ""
      echo "Available Options:"
      echo ""
      echo "-j [N]        Compile in parallel using N CPUs"
      echo ""
      echo "-move         Move compiled objects to build directory"
      echo ""
      echo "-noclean      Do not run clean script"
      echo ""
      echo "-noconfig     Do not run configure compilation script"
      echo ""
     exit 1
    breaksw

  endsw
end

#--------------------------------------------------------------------------
# Set a local environmental variable to define the root path to the
# directories where ROMS and WRF source files are located.  The ROMS
# source directory is needed for replacing several WRF files for
# ESMF/NUOPC Coupling (see below).
#--------------------------------------------------------------------------

 setenv ROMS_SRC_DIR         ${HOME}/ocean/repository/coupling

 setenv WRF_ROOT_DIR         ${HOME}/ocean/repository/WRF
 setenv WRF_SRC_DIR          ${WRF_ROOT_DIR}

#--------------------------------------------------------------------------
# Set a local environmental variable to define the path of the working
# application directory where all this project's files are kept.
#--------------------------------------------------------------------------

 setenv MY_PROJECT_DIR       ${PWD}
 setenv MY_PROJECT_DATA      `dirname ${PWD}`/Data

#--------------------------------------------------------------------------
# COAMPS configuration CPP options.
#--------------------------------------------------------------------------

# Sometimes it is desirable to activate one or more CPP options to
# configure a particular application. If it is the case, specify each
# option here using the -D syntax. Notice also that you need to use
# shell's quoting syntax to enclose the definition. Both single or
# double quotes work. For example,
#
#setenv MY_CPP_FLAGS         "${MY_CPP_FLAGS} -DHDF5"
#
# can be used to read input data from HDF5 file instead of flat files.
# Notice that you can have as many definitions as you want by appending
# values.

 setenv MY_CPP_FLAGS         "${MY_CPP_FLAGS} -DHDF5"

#--------------------------------------------------------------------------
# Set Fortran compiler and MPI library to use.
#--------------------------------------------------------------------------

# Set path of the directory containing makefile configuration (*.mk) files.
# The user has the option to specify a customized version of these files
# in a different directory than the one distributed with the source code,
# ${COAMPS_ROOT_DIR}/compilers. If this is the case, the you need to keep
# these configurations files up-to-date.

 setenv USE_MPI              on          # distributed-memory parallelism
 setenv USE_MPIF90           on          # compile with mpif90 script
#setenv which_MPI            mpich       # compile with MPICH library
#setenv which_MPI            mpich2      # compile with MPICH2 library
 setenv which_MPI            openmpi     # compile with OpenMPI library

 setenv FORT                 ifort
#setenv FORT                 gfortran
#setenv FORT                 pgi

#setenv USE_REAL_DOUBLE      on          # use real double precision (-r8)
#setenv USE_DEBUG            on          # use Fortran debugging flags
 setenv USE_HDF5             on          # compile with HDF5 library
 setenv USE_NETCDF           on          # compile with NetCDF
#setenv USE_NETCDF4          on          # compile with NetCDF-4 library
                                         # (Must also set USE_NETCDF)

#--------------------------------------------------------------------------
# Use my specified library paths. It is not needed but it is added for
# for checking in the future.
#--------------------------------------------------------------------------

 setenv USE_MY_LIBS no           # use system default library paths
#setenv USE_MY_LIBS yes          # use my customized library paths

#set MY_PATHS = ${ROMS_SRC_DIR}/Compilers/my_build_paths.csh
 set MY_PATHS = ${HOME}/Compilers/ROMS/my_build_paths.csh

if ($USE_MY_LIBS == 'yes') then
  source ${MY_PATHS} ${MY_PATHS}
endif

#--------------------------------------------------------------------------
# WRF build and executable directory.
#--------------------------------------------------------------------------

# Put the *.a, .f and .f90 files in a project specific Build directory to
# avoid conflict with other projects.

if ($?USE_DEBUG) then
  setenv  WRF_BUILD_DIR   ${MY_PROJECT_DIR}/Build_wrfG
else
  setenv  WRF_BUILD_DIR   ${MY_PROJECT_DIR}/Build_wrf
endif

# Put WRF executables in the following directory.

setenv  WRF_BIN_DIR       ${WRF_BUILD_DIR}/Bin

# Go to the users source directory to compile. The options set above will
# pick up the application-specific code from the appropriate place.

 cd ${WRF_ROOT_DIR}

#--------------------------------------------------------------------------
# Configure. It creates configure.wrf script used for compilation.
#
#   -d    build with debugging information and no optimization
#   -D    build with -d AND floating traps, traceback, uninitialized variables
#   -r8   build with 8 byte reals
#
# During configuration the WRF/arch/Config.pl perl script is executed and
# we need to interactively select the combination of compiler and parallel
# comunications option. For example, for Darwin operating system we get:
#
# Please select from among the following Darwin ARCH options:
#
#  1. (serial)   2. (smpar)   3. (dmpar)   4. (dm+sm)   PGI (pgf90/pgcc)
#  5. (serial)   6. (smpar)   7. (dmpar)   8. (dm+sm)   INTEL (ifort/icc)
#  9. (serial)  10. (smpar)  11. (dmpar)  12. (dm+sm)   INTEL (ifort/clang)
# 13. (serial)               14. (dmpar)                GNU (g95/gcc)
# 15. (serial)  16. (smpar)  17. (dmpar)  18. (dm+sm)   GNU (gfortran/gcc)
# 19. (serial)  20. (smpar)  21. (dmpar)  22. (dm+sm)   GNU (gfortran/clang)
# 23. (serial)               24. (dmpar)                IBM (xlf90_r/cc)
# 25. (serial)  26. (smpar)  27. (dmpar)  28. (dm+sm)   PGI (pgf90/pgcc): -f90=pgf90
# 29. (serial)  30. (smpar)  31. (dmpar)  32. (dm+sm)   INTEL (ifort/icc): Open MPI
# 33. (serial)  34. (smpar)  35. (dmpar)  36. (dm+sm)   INTEL/GNU (ifort/gcc): Open MPI
# 37. (serial)  38. (smpar)  39. (dmpar)  40. (dm+sm)   GNU (gfortran/gcc): Open MPI
#
# Enter selection [1-40] : ??
#
# For coupling with ESMF/NUOPC, we need select an option from the (dmpar)
# column for distributed-memory configuration.
#--------------------------------------------------------------------------

# Clean source code and remove build directory.

if ( $clean == 1 ) then
  echo ""
  echo "${separator}"
  echo "Cleaning WRF source code:  ${WRF_ROOT_DIR}/clean -a"
  echo "${separator}"
  echo ""
  ${WRF_ROOT_DIR}/clean -a            # clean source code
  /bin/rm -rf ${WRF_BUILD_DIR}        # remove existing build directories
endif

if ($?USE_DEBUG) then
#  set DEBUG_FLAG = -d
   set DEBUG_FLAG = -D
endif

setenv CONFIG_FLAGS ''

if ( $config == 1 ) then
  if ($?USE_DEBUG && $?USE_REAL_DOUBLE) then
    setenv CONFIG_FLAGS   "${DEBUG_FLAG} -r8"
  else if ($?USE_DEBUG) then
    setenv CONFIG_FLAGS   ${DEBUG_FLAG}
  else if ($?USE_REAL_DOUBLE) then
    setenv CONFIG_FLAGS   -r8
  endif

  set CHECK_STRING = 'WRF-ROMS ESMF-NUOPC Coupling'
  echo ""
  echo "${separator}"
  echo "If applicable, replacing several WRF files for WRF ESMF/NUOPC Coupling"
  echo "${separator}"
  echo ""

# Reworking linking NetCDF4 library dependencies

  if (! `grep -c "${CHECK_STRING}" ${WRF_ROOT_DIR}/configure`) then
    mv -v ${WRF_ROOT_DIR}/configure ${WRF_ROOT_DIR}/configure.orig
    cp -fv ${ROMS_SRC_DIR}/ESM/wrf_configure ${WRF_ROOT_DIR}/configure
  else
    echo "   No need to replace: ${WRF_ROOT_DIR}/configure"
  endif

  if (! `grep -c "${CHECK_STRING}" ${WRF_ROOT_DIR}/Makefile`) then
    mv -v  ${WRF_ROOT_DIR}/Makefile ${WRF_ROOT_DIR}/Makefile.orig
    cp -fv ${ROMS_SRC_DIR}/ESM/wrf_Makefile  ${WRF_ROOT_DIR}/Makefile
  else
    echo "   No need to replace: ${WRF_ROOT_DIR}/Makefile"
  endif

  if (! `grep -c "${CHECK_STRING}" ${WRF_ROOT_DIR}/arch/postamble`) then
    mv -v  ${WRF_ROOT_DIR}/arch/postamble ${WRF_ROOT_DIR}/arch/postamble.orig
    cp -fv ${ROMS_SRC_DIR}/ESM/wrf_postamble ${WRF_ROOT_DIR}/arch/postamble
  else
    echo "   No need to replace: ${WRF_ROOT_DIR}/arch/postamble"
  endif

# Changing -openmp to -qopenmp, renaming ESMF/esmf to MYESMF/myesmf, adding
# Intel/GNU with OpenMPI

  if (! `grep -c "${CHECK_STRING}" ${WRF_ROOT_DIR}/arch/configure.defaults`) then
    mv -v  ${WRF_ROOT_DIR}/arch/configure.defaults ${WRF_ROOT_DIR}/arch/configure.defaults.orig
    cp -fv ${ROMS_SRC_DIR}/ESM/wrf_configure.defaults ${WRF_ROOT_DIR}/arch/configure.defaults
  else
    echo "   No need to replace: ${WRF_ROOT_DIR}/arch/configure.defaults"
  endif

# Renaming ESMF/esmf to MYESMF/myesmf

  if (! `grep -c "${CHECK_STRING}" ${WRF_ROOT_DIR}/arch/Config.pl`) then
    mv -v  ${WRF_ROOT_DIR}/arch/Config.pl ${WRF_ROOT_DIR}/arch/Config.pl.orig
    cp -fv ${ROMS_SRC_DIR}/ESM/wrf_Config.pl ${WRF_ROOT_DIR}/arch/Config.pl
  else
    echo "   No need to replace: ${WRF_ROOT_DIR}/arch/Config.pl"
  endif

  if (! `grep -c "${CHECK_STRING}" ${WRF_ROOT_DIR}/external/esmf_time_f90/Makefile`) then
    mv -v  ${WRF_ROOT_DIR}/external/esmf_time_f90/Makefile ${WRF_ROOT_DIR}/external/esmf_time_f90/Makefile.orig
    cp -fv ${ROMS_SRC_DIR}/ESM/wrf_Makefile.esmf ${WRF_ROOT_DIR}/external/esmf_time_f90/Makefile
  else
    echo "   No need to replace: ${WRF_ROOT_DIR}/external/esmf_time_f90/Makefile"
  endif

# Correcting optional argument from defaultCalendar to defaultCalKind to
# ESMF_Initialize call

  if (! `grep -c "${CHECK_STRING}" ${WRF_ROOT_DIR}/external/esmf_time_f90/Test1.F90`) then
    mv -v  ${WRF_ROOT_DIR}/external/esmf_time_f90/Test1.F90 ${WRF_ROOT_DIR}/external/esmf_time_f90/Test1.F90.orig
    cp -fv ${ROMS_SRC_DIR}/ESM/wrf_Test1.F90 ${WRF_ROOT_DIR}/external/esmf_time_f90/Test1.F90
  else
    echo "   No need to replace: ${WRF_ROOT_DIR}/external/esmf_time_f90/Test1.F90"
  endif

  echo ""
  echo "${separator}"
  echo "Configuring WRF code:  ${WRF_ROOT_DIR}/configure ${CONFIG_FLAGS}"
  echo "${separator}"
  echo ""

  ${WRF_ROOT_DIR}/configure ${CONFIG_FLAGS}

#  Custom CPP Macros for renaming ESMF/esmf to MYESMF/myesmf to avoid
#  conflict with newer versions of the ESMF/NUOPC libraries

  cat ${ROMS_SRC_DIR}/ESM/wrf_add_configure.wrf >> ${WRF_ROOT_DIR}/configure.wrf

endif

#--------------------------------------------------------------------------
# WRF Compile script:
#
# Usage:
#
#     compile [-j n] wrf   compile wrf in run dir (NOTE: no real.exe,
#                          ndown.exe, or ideal.exe generated)
#
#   or choose a test case (see README_test_cases for details) :
#
#     compile [-j n] em_b_wave
#     compile [-j n] em_convrad
#     compile [-j n] em_esmf_exp
#     compile [-j n] em_fire
#     compile [-j n] em_grav2d_x
#     compile [-j n] em_heldsuarez
#     compile [-j n] em_hill2d_x
#     compile [-j n] em_les
#     compile [-j n] em_quarter_ss
#     compile [-j n] em_real
#     compile [-j n] em_scm_xy
#     compile [-j n] em_seabreeze2d_x
#     compile [-j n] em_squall2d_x
#     compile [-j n] em_squall2d_y
#     compile [-j n] em_tropical_cyclone
#     compile [-j n] nmm_real
#     compile [-j n] nmm_tropical_cyclone
#
#     compile -j n         parallel make using n tasks if supported
#                          (default 2)
#
#     compile -h           help message
#--------------------------------------------------------------------------

setenv WRF_DA_CORE 0             # no Data Assimilation core
setenv WRF_EM_CORE 1             # Eurelian Mass-coordinate core
setenv WRF_NMM_CORE 0            # Nonhydrostatic Mesoscale Model core

# Compile (the binary will go to BINDIR set above).

#set WRF_CASE = wrf
 set WRF_CASE = em_real

if ( $parallel == 1 ) then
  setenv J "-j ${NCPUS}"
else
  setenv J "-j 2"
endif

echo ""
echo "${separator}"
echo "Compiling WRF using  ${MY_PROJECT_DIR}/${0}:"
echo ""
echo "   ${WRF_ROOT_DIR}/compile ${WRF_CASE}"
echo "        WRF_DA_CORE = ${WRF_DA_CORE},    Data Assimilation core"
echo "        WRF_EM_CORE = ${WRF_EM_CORE},    Eurelian Mass-coordinate core"
echo "        WRF_NMM_CORE = ${WRF_NMM_CORE}, Nonhydrostatic Mesoscale Model core"
echo "        J = ${J},          number of compiling CPUs"
echo "${separator}"
echo ""

${WRF_ROOT_DIR}/compile ${WRF_CASE}

#--------------------------------------------------------------------------
# Move WRF objects and executables.
#--------------------------------------------------------------------------

if ( $move == 1 ) then

  echo ""
  echo "${separator}"
  echo "Moving WRF objects to Build directory  ${WRF_BUILD_DIR}:"
  echo "${separator}"
  echo ""

  if ( ! -d ${WRF_BUILD_DIR} ) then
    /bin/mkdir -pv ${WRF_BUILD_DIR}
    /bin/mkdir -pv ${WRF_BUILD_DIR}/Bin
    echo ""
  endif

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
# compact "find" function for these file extensions.

  foreach DIR ( frame chem share dyn_em dyn_exp dyn_nmm phys main tools wrftladj )
    if ( -d $DIR ) then
      /bin/mv -fv ${DIR}/*.f90 ${WRF_BUILD_DIR}
      /bin/mv -fv ${DIR}/*.F90 ${WRF_BUILD_DIR}
      /bin/mv -fv ${DIR}/*.inc ${WRF_BUILD_DIR}
    endif
  end

  find ${WRF_ROOT_DIR} -type f -name "*.mod" -exec /bin/mv -fv {} ${WRF_BUILD_DIR} \;
  find ${WRF_ROOT_DIR} -type f -name "*.o"   -exec /bin/mv -fv {} ${WRF_BUILD_DIR} \;
  find ${WRF_ROOT_DIR} -type f -name "*.obj" -exec /bin/mv -fv {} ${WRF_BUILD_DIR} \;
  find ${WRF_ROOT_DIR} -type f -name "*.a"   -exec /bin/mv -fv {} ${WRF_BUILD_DIR} \;

  find ${WRF_ROOT_DIR} -type f -name "*.exe" -exec /bin/mv -fv {} ${WRF_BIN_DIR} \;
  find ${WRF_ROOT_DIR} -type l -name "*.exe" -exec /bin/rm -fv {} \;

  /bin/mv -fv external/esmf_time_f90/*.f ${WRF_BUILD_DIR}
  /bin/mv -fv external/esmf_time_f90/MYESMF*.inc ${WRF_BUILD_DIR}

  /bin/mv -fv external/io_int/diffwrf ${WRF_BIN_DIR}/diffwrf_int
  /bin/mv -fv external/io_int/test_io_idx ${WRF_BIN_DIR}
  /bin/mv -fv external/io_netcdf/diffwrf ${WRF_BIN_DIR}/diffwrf_nc

  /bin/mv -fv tools/fseeko_test ${WRF_BIN_DIR}
  /bin/mv -fv tools/registry ${WRF_BIN_DIR}

endif

#--------------------------------------------------------------------------
# Create WRF data links.
#--------------------------------------------------------------------------

if ( $WRF_CASE == "em_real" ) then

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
  ln -sfv ${WRF_ROOT_DIR}/run/p3_lookup_table_1.dat-v2.8.2 .
  ln -sfv ${WRF_ROOT_DIR}/run/p3_lookup_table_2.dat-v2.8.2 .

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

# Remove links in WRF/test/em_real sub-directory

  echo ""
  echo "${separator}"
  echo "Removing WRF data links from  ${WRF_ROOT_DIR}/test/em_real:"
  echo "${separator}"
  echo ""

  find ${WRF_ROOT_DIR}/test/em_real -type l -exec /bin/rm -fv {} \;

endif
