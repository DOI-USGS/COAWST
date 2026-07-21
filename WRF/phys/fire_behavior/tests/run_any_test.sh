#!/bin/bash
#
#################################################
#
# Purpose: Check code reproduces results of any test case
#
#################################################


#################################################
# Functions

# usage instructions
usage () {
  printf "Usage: $0 --test=[NUMBER] [OPTIONS]...\n"
  printf "\n"
  printf "OPTIONS\n"
  printf "  --esmx\n"
  printf "      run tests with ESMX instead of standalone executable\n"
  printf "  --no-purge\n"
  printf "      turn off purging of test file\n"
  printf "  -p1=0\n"
  printf "      turn off check of fire area\n"
  printf "  -p2=0\n"
  printf "      turn off check of heat output\n"
  printf "  -p3=0\n"
  printf "      turn off check of latent heat output\n"
  printf "  -p4=0\n"
  printf "      turn off check of max heat flux\n"
  printf "  -p5=0\n"
  printf "      turn off check of max latent heat flux\n"
  printf "  -p6=0\n"
  printf "      turn off check of (only in test1)\n"
  printf "\n"
}

# -----------------------------------------------

# find system name
find_system () {
    local sysname=`hostname`
    sysname="${sysname//[[:digit:]]/}"
    echo "$sysname"
}

# -----------------------------------------------

testp1to5 () {
    myvar=$1
    mycol=$2
    
    #echo "Testing for var=${myvar} col=${mycol}"

    n_tests=$(expr $n_tests + 1)
    [[ -f file1.dat ]] && rm file1.dat
    [[ -f file2.dat ]] && rm file2.dat

    grep "$myvar" $file_output  | awk '{print $2, $"${mycol}"}' | tr -d " " > ./file1.dat
    grep "$myvar" $file_wrf     | awk '{print $2, $"${mycol}"}' | tr -d " " > ./file2.dat

    test=$(diff ./file1.dat ./file2.dat | wc -l)
    if [ $test -eq 0 ]
    then
	echo "  Test ${myvar} PASSED"
	n_test_passed=$(expr $n_test_passed + 1)
    else
	echo "  Test ${myvar} FAILS"
	#exit 1
	diff ./file1.dat ./file2.dat
    fi

}

# -----------------------------------------------

testp6 () {
    
    echo "Testing p6 "

    n_tests=$(expr $n_tests + 1)

    [[ -f file1.dat ]] && rm file1.dat
    [[ -f file2.dat ]] && rm file2.dat

    head -n 1 ./fort.34 | tr -d " 045" > ./file1.dat # offline
    head -n 1 ./th_qv_tend.dat | tr -d " 045" > ./file2.dat # online

    test=$(diff ./file1.dat ./file2.dat | wc -l)
    if [ $test -eq 0 ]
    then
	echo "  Test p6 PASSED"
	n_test_passed=$(expr $n_test_passed + 1)
    else
	echo "  Test p6 FAILS"
    fi

}

#################################################
# Main Code
#################################################

# -----------------------------------------------
# Read command line arguments

((!$#)) && echo 'No arguments supplied!' && usage && exit 1

while [ $# -gt 0 ] 
do
    key="${1}"

    case "${key}" in
	--help|-h) usage; exit 0 ;;

	# -----------------------------------------------
	# Read test number
	-t=*|--test=*)
	    thistest="${key#*=}"
	    echo "TestAny: running test ${thistest}"
	    shift ## shift past key and value
	    ;;

	# Optional Arguments:

	# -----------------------------------------------
	# Test using esmx (default is 0, i.e. use standalone)
	--esmx)
	    useesmx=1
	    echo "Using esmx"
	    shift ## shift past key and value
	    ;;

	# -----------------------------------------------
	--no-purge)
	    purge_output=0
	    echo "Purge test files is off"
	    shift ## shift past key and value
	    ;;

	# -----------------------------------------------
	# Specify test part (default is 1)
	-p1=*)
	    testp1="${key#*=}"
	    echo "Test part1: ${testp1}"
	    shift ## shift past key and value
	    ;;
	-p2=*)
	    testp2="${key#*=}"
	    echo "Test part2: ${testp2}"
	    shift ## shift past key and value
	    ;;
	-p3=*)
	    testp3="${key#*=}"
	    echo "Test part3: ${testp3}"
	    shift ## shift past key and value
	    ;;
	-p4=*)
	    testp4="${key#*=}"
	    echo "Test part4: ${testp4}"
	    shift ## shift past key and value
	    ;;
	-p5=*)
	    testp5="${key#*=}"
	    echo "Test part5: ${testp5}"
	    shift ## shift past key and value
	    ;;
	# -----------------------------------------------

	-?*) printf "ERROR: Unknown option $1\n"; usage; exit 1 ;;

	*)
    esac
done

# -----------------------------------------------
# Identify test to run

thistest="${thistest:?Invalid test argument. Provide arg for test: -t=test1 or -t=1}"

if [[ ${thistest} != test* ]]
then
    if [[ "${#thistest}" -eq 1 ]]
    then
	echo "Setting test name to test${thistest}"
	thistest=test${thistest}
    else
	echo "Something is not right with the argument test"
	exit 1
    fi
fi

#################################################
# Handle System Modules

SYSTEM=""
ENV_DIR="../env"
ENV_FILE=""

# automatically determine system
if [ -z "${SYSTEM}" ] ; then
  SYSTEM=$(find_system)
fi

# automatically determine module file
if [ -z "${ENV_FILE}" ] ; then
    if [ "${SYSTEM}" == "cheyenne" ] ; then
	ENV_FILE="${SYSTEM}/19.1.1"
    else
	ENV_FILE="unknown"
    fi
fi

# load environment 
if [ ! -f "${ENV_DIR}/${ENV_FILE}" ]; then
  printf "ERROR: ${ENV_FILE} does not exist in ${ENV_DIR}.\n"
  printf "\n"
  exit 1
fi
source ${ENV_DIR}/${ENV_FILE}

#################################################
# defaults

useesmx=${useesmx:=0}
purge_output=${purge_output:=1} # 0) No, 1) yes
testp1=${testp1:=1} # Check fire area
testp2=${testp2:=1} # Check heat output
testp3=${testp3:=1} # Check latent heat output
testp4=${testp4:=1} # Check Max heat flux
testp5=${testp5:=1} # Check Max latent heat flux
testp6=${testp6:=1} # Check tendency t,qv first time step, test1 only

esmx_exe="../../build/esmx"
standalone_exe="../../install/bin/fire_behavior.exe"

#################################################
# Tasks for standalone and esmx

file_wrf=test_solution.txt
file_output=${thistest}_output.txt

TEST_DIR="$PWD"
cd ${thistest}

# clean up and link files for this test
[[ -f "${file_output}" ]] && rm -f ${file_output}

#################################################
# Run standalone executable

if [[ ${useesmx} -eq 0 ]]
then
    # -----------------------------------------------
    # clean old links
    [[ -L "fire_behavior.exe" ]] && rm fire_behavior.exe
    
    # -----------------------------------------------
    # check executable is present, link it, run it
    if [[ -f ${standalone_exe} ]]
    then
	echo "Running standalone code"
	ln -sf ${standalone_exe} ./fire_behavior.exe
	./fire_behavior.exe > ${file_output}
    else
	echo 'Please compile ${standalone_exe}'
	exit 1
    fi

#################################################
# Run ESMX

elif  [[ ${useesmx} -eq 1 ]]
then

    # -----------------------------------------------
    # Check for PBS_ACCOUNT
    if [[ $(env | grep PBS_ACCOUNT=) != PBS_ACCOUNT* ]] 
    then
	echo "Set a PBS account in your shell environment: export PBS_ACCOUNT=XYZ000K"
	exit 1
    # else
    # 	echo "Using $(env | grep PBS_ACCOUNT=)"
    fi

    # -----------------------------------------------
    # clean old esmx log files
    for oldfile in PET*.ESMF_LogFile 
    do
	[[ -f "${oldfile}" ]] && rm ${oldfile}
    done

    # -----------------------------------------------
    # clean old links
    [[ -L "esmx" ]] && rm esmx
    [[ -L "esmxRun.config" ]] && rm esmxRun.config

    # -----------------------------------------------
    # check executable is present & link it
    if [[ ! -f "${esmx_exe}" ]] 
    then
	echo "Executable build/esmx is not present. Is it compiled?"
	exit 1
    fi 
    ln -sf ${esmx_exe} .
    ln -sf ${TEST_DIR}/esmxRun.config .
    ln -sf ${TEST_DIR}/fd_nems.yaml .

    # -----------------------------------------------
    # Submit job with qcmd
    if [[ -L esmx && -L esmxRun.config ]]
    then
	echo "Running ESMX code"

	cmd="./esmx > ${file_output}"
	qcmd -- mpirun -np 1 ${cmd}
	if [[ $? -ne 0 ]]
	then
	    echo "error running test with qcmd"
	    exit 1
	fi
    else
	echo "esmx or esmxRun.config are not linked."
	exit 1
    fi
fi

#################################################
# Check results

n_tests=0
n_test_passed=0
myvar=""
mycol=""

echo "Results for ${thistest}:"

[[ ${testp1} -eq 1 ]] && var="Fire area"            && col=6 && testp1to5 "${var}" "${col}"
[[ ${testp2} -eq 1 ]] && var="Heat output"          && col=6 && testp1to5 "${var}" "${col}"	      
[[ ${testp3} -eq 1 ]] && var="Latent heat output"   && col=7 && testp1to5 "${var}" "${col}" 
[[ ${testp4} -eq 1 ]] && var="Max heat flux"        && col=7 && testp1to5 "${var}" "${col}"	      
[[ ${testp5} -eq 1 ]] && var="Max latent heat flux" && col=8 && testp1to5 "${var}" "${col}"
[[ ${thistest} == test1  && ${testp6} -eq 1 ]] && testp6 

#################################################
# Purge

if [[ $purge_output -eq 1 ]]
then

    # remove files
    for myfile in ${file_output} fort.34 file1.dat file2.dat namelist.fire.output
    do
	[[ -f ${myfile} ]] && rm ${myfile}
    done
    # remove links
    # for mylink in fire_behavior.exe esmxRun.config esmx
    # do
    # 	[[ -L ${mylink} ]] && rm ${mylink}
    # done

fi

cd ${TEST_DIR}

#################################################
# Print summary of Test

if [ $n_test_passed -eq $n_tests ]
then
  echo "SUCCESS: $n_test_passed PASSED of $n_tests"
  echo ''
else
  echo "FAILED: $n_test_passed PASSED of $n_tests"
  echo ''
  exit 1
fi


