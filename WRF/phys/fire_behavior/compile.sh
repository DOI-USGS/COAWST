#!/bin/bash
set -e

# usage instructions
usage () {
  printf "Usage: $0 [OPTIONS]...\n"
  printf "\n"
  printf "OPTIONS\n"
  printf "  --env-auto\n"
  printf "      load preconfigured environment based on system\n"
  printf "  --env-file=ENV_FILE\n"
  printf "      load environment from file\n"
  printf "  --build-dir=BUILD_DIR\n"
  printf "      build directory\n"
  printf "  --build-type=BUILD_TYPE\n"
  printf "      build type; valid options are 'debug', 'release',\n"
  printf "      'relWithDebInfo'\n"
  printf "  --build-jobs=BUILD_JOBS\n"
  printf "      number of jobs used for building esmx and components\n"
  printf "  --mpi-off\n"
  printf "      build without mpi library\n"
  printf "  --nuopc, -n\n"
  printf "      build NUOPC library and module\n"
  printf "  --esmx, -x\n"
  printf "      build ESMX application (includes NUOPC)\n"
  printf "  --openmp-on\n"
  printf "      enable OpenMP parallelization\n"
  printf "  --prefix=INSTALL_PREFIX\n"
  printf "      installation prefix\n"
  printf "  --verbose, -v\n"
  printf "      build with verbose output\n"
  printf "  --test[=TEST_NAME], -t[=TEST_NAME]\n"
  printf "      run tests\n"
  printf "  --clean\n"
  printf "      delete build and install directories\n"
  printf "\n"
}

# print settings
settings () {
  printf "#######################################################\n"
  printf "Settings:\n\n"
  printf "\tSYSTEM=${SYSTEM}\n"
  printf "\tENV_AUTO=${ENV_AUTO}\n"
  printf "\tENV_FILE=${ENV_FILE}\n"
  printf "\tBUILD_DIR=${BUILD_DIR}\n"
  printf "\tBUILD_TYPE=${BUILD_TYPE}\n"
  printf "\tBUILD_JOBS=${BUILD_JOBS}\n"
  printf "\tMPI=${MPI}\n"
  printf "\tNUOPC=${NUOPC}\n"
  printf "\tESMX=${ESMX}\n"
  printf "\tINSTALL_PREFIX=${INSTALL_PREFIX}\n"
  printf "\tVERBOSE=${VERBOSE}\n"
  printf "\tTEST=${TEST}\n"
  printf "\tTEST_NAME=${TEST_NAME}\n"
  printf "\tCLEAN=${CLEAN}\n"
  printf "#######################################################\n"
}

# find system name
find_system () {
    local sysname=`hostname`
    sysname="${sysname//[[:digit:]]/}"
    echo "$sysname"
}

#------------------------------------------------------------------------------

# default settings
FIRE_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
SYSTEM=""
ENV_AUTO=false
ENV_DIR="${FIRE_DIR}/env"
ENV_FILE=""
BUILD_DIR="${FIRE_DIR}/build"
BUILD_TYPE="release"
BUILD_JOBS=""
INSTALL_PREFIX="${FIRE_DIR}/install"
MPI=true
NUOPC=false
ESMX=false
VERBOSE=false
TEST=false
TEST_NAME=""
CLEAN=false
OPENMP=false

#------------------------------------------------------------------------------

# process arguments
POSITIONAL_ARGS=()
while :; do
  case $1 in
    --help|-h) usage; exit 0 ;;
    --system=?*) SYSTEM=${1#*=} ;;
    --system) printf "ERROR: $1 requires an argument.\n"; usage; exit 1 ;;
    --system=) printf "ERROR: $1 requires an argument.\n"; usage; exit 1 ;;
    --env-auto) ENV_AUTO=true ;;
    --env-auto=?*) printf "ERROR: $1 argument ignored.\n"; usage; exit 1 ;;
    --env-auto=) printf "ERROR: $1 argument ignored.\n"; usage; exit 1 ;;
    --env-dir=?*) ENV_DIR=${1#*=} ;;
    --env-dir) printf "ERROR: $1 requires an argument.\n"; usage; exit 1 ;;
    --env-dir=) printf "ERROR: $1 requires an argument.\n"; usage; exit 1 ;;
    --env-file=?*) ENV_FILE=${1#*=} ;;
    --env-file) printf "ERROR: $1 requires an argument.\n"; usage; exit 1 ;;
    --env-file=) printf "ERROR: $1 requires an argument.\n"; usage; exit 1 ;;
    --build-dir=?*) BUILD_DIR=${1#*=} ;;
    --build-dir) printf "ERROR: $1 requires an argument.\n"; usage; exit 1 ;;
    --build-dir=) printf "ERROR: $1 requires an argument.\n"; usage; exit 1 ;;
    --build-type=?*) BUILD_TYPE=${1#*=} ;;
    --build-type) printf "ERROR: $1 requires an argument.\n"; usage; exit 1 ;;
    --build-type=) printf "ERROR: $1 requires an argument.\n"; usage; exit 1 ;;
    --build-jobs=?*) BUILD_JOBS=${1#*=} ;;
    --build-jobs) printf "ERROR: $1 requires an argument.\n"; usage; exit 1 ;;
    --build-jobs=) printf "ERROR: $1 requires an argument.\n"; usage; exit 1 ;;
    --prefix=?*) INSTALL_PREFIX=${1#*=} ;;
    --prefix) printf "ERROR: $1 requires an argument.\n"; usage; exit 1 ;;
    --prefix=) printf "ERROR: $1 requires an argument.\n"; usage; exit 1 ;;
    --mpi-off) MPI=false ;;
    --mpi-off=?*) printf "ERROR: $1 argument ignored.\n"; usage; exit 1 ;;
    --mpi-off=) printf "ERROR: $1 argument ignored.\n"; usage; exit 1 ;;
    --nuopc|-n) NUOPC=true ;;
    --nuopc=?*) printf "ERROR: $1 argument ignored.\n"; usage; exit 1 ;;
    --nuopc=) printf "ERROR: $1 argument ignored.\n"; usage; exit 1 ;;
    --esmx|-x) NUOPC=true; ESMX=true ;;
    --esmx=?*) printf "ERROR: $1 argument ignored.\n"; usage; exit 1 ;;
    --esmx=) printf "ERROR: $1 argument ignored.\n"; usage; exit 1 ;;
    --verbose|-v) VERBOSE=true ;;
    --verbose=?*) printf "ERROR: $1 argument ignored.\n"; usage; exit 1 ;;
    --verbose=) printf "ERROR: $1 argument ignored.\n"; usage; exit 1 ;;
    --test|-t) TEST=true ;;
    --test=?*|-t=?*) TEST=true; TEST_NAME=${1#*=} ;;
    --test=) TEST=true ;;
    --clean) CLEAN=true ;;
    --openmp-on) OPENMP=true ;;
    --clean=?*) printf "ERROR: $1 argument ignored.\n"; usage; exit 1 ;;
    --clean=) printf "ERROR: $1 argument ignored.\n"; usage; exit 1 ;;
    -?*) printf "ERROR: Unknown option $1\n"; usage; exit 1 ;;
    ?*) POSITIONAL_ARGS+=("${1}") ;;
    *) break ;;
  esac
  shift
done
set -- "${POSITIONAL_ARGS[@]}"

if [[ $# -ge 1 ]]; then
  printf "ERROR: Unknown argument $1\n"
  usage
  exit 1
fi

#------------------------------------------------------------------------------

# automatically determine system
if [ -z "${SYSTEM}" ] ; then
  SYSTEM=$(find_system)
fi

# print settings
if [ "${VERBOSE}" = true ] ; then
  settings
fi

# auto environment
if [ "${ENV_AUTO}" = true ] ; then
  case ${SYSTEM} in
    cheyenne) AUTOFILE="${ENV_DIR}/cheyenne/gnu-10.1.0";;
    derecho) AUTOFILE="${ENV_DIR}/derecho/gnu-12.2.0";;
    *) printf "ERROR: unspecified --env-auto for ${SYSTEM}\n"; exit 1 ;;
  esac
  if [ -f "${AUTOFILE}" ]; then
    source ${AUTOFILE}
  else
    printf "ERROR: ${AUTOFILE} does not exist\n"
    exit 1
  fi
fi

# user environment
if [ ! -z "${ENV_FILE}" ]; then
  if [ -f "${ENV_FILE}" ]; then
    source ${ENV_FILE}
  else
    printf "ERROR: ${ENV_FILE} does not exist\n"
    exit 1
  fi
fi

set -u

# clean
if [ "${CLEAN}" = true ]; then
  rm -rf ${BUILD_DIR}
fi

# generate
CMAKE_SETTINGS=("")
if [ ! -z "${BUILD_TYPE}" ]; then
  CMAKE_SETTINGS+=("-DCMAKE_BUILD_TYPE=${BUILD_TYPE}")
fi
if [ ! -z "${INSTALL_PREFIX}" ]; then
  CMAKE_SETTINGS+=("-DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}")
  CMAKE_SETTINGS+=("-DCMAKE_PREFIX_PATH=${INSTALL_PREFIX}")
fi
if [ "${MPI}" = true ]; then
  CMAKE_SETTINGS+=("-DDM_PARALLEL=ON")
else
  CMAKE_SETTINGS+=("-DDM_PARALLEL=OFF")
fi
if [ "${NUOPC}" = true ]; then
  CMAKE_SETTINGS+=("-DNUOPC=ON")
else
  CMAKE_SETTINGS+=("-DNUOPC=OFF")
fi
if [ "${ESMX}" = true ]; then
  CMAKE_SETTINGS+=("-DESMX=ON")
else
  CMAKE_SETTINGS+=("-DESMX=OFF")
fi
if [ "${OPENMP}" = true ]; then
  CMAKE_SETTINGS+=("-DOPENMP=ON")
else
  CMAKE_SETTINGS+=("-DOPENMP=OFF")
fi
cmake -S${FIRE_DIR} -B${BUILD_DIR} ${CMAKE_SETTINGS[@]}
if [ "$?" !=  "0" ]; then
  echo "$0 Failed: (cmake)"
  exit -1
fi

# build
BUILD_SETTINGS=("")
if [ "${VERBOSE}" = true ]; then
  BUILD_SETTINGS+=("-v")
fi
if [ ! -z "${BUILD_JOBS}" ]; then
  BUILD_SETTINGS+=("-j ${BUILD_JOBS}")
fi
cmake --build ${BUILD_DIR} ${BUILD_SETTINGS[@]}
if [ "$?" !=  "0" ]; then
  echo "$0 Failed: (cmake --build)"
  exit -2
fi

# install
INSTALL_SETTINGS=("")
cmake --install ${BUILD_DIR} ${INSTALL_SETTINGS[@]}
if [ "$?" !=  "0" ]; then
  echo "$0 Failed: (cmake --install)"
  exit -3
fi

# build and install esmx
if [ "${ESMX}" = true ]; then
  cmake --build ${BUILD_DIR} ${BUILD_SETTINGS[@]} --target esmx
  if [ "$?" !=  "0" ]; then
    echo "$0 Failed: (cmake --build)"
    exit -4
  fi
  cmake --install ${BUILD_DIR} ${INSTALL_SETTINGS[@]} --component esmx
  if [ "$?" !=  "0" ]; then
    echo "$0 Failed: (cmake --install)"
    exit -5
  fi
fi

# test
TEST_SETTINGS=("")
if [ "${VERBOSE}" = true ]; then
  TEST_SETTINGS=("--output-on-failure")
  TEST_SETTINGS=("--verbose")
fi
if [ ! -z "${TEST_NAME}" ]; then
  TEST_SETTINGS+=("--tests-regex ${TEST_NAME}")
fi 
if [ "${TEST}" = true ]; then
  ctest --test-dir ${BUILD_DIR}/tests ${TEST_SETTINGS[@]} 
  if [ "$?" !=  "0" ]; then
    echo "$0 Failed: (ctest)"
    exit -6
  fi
fi
