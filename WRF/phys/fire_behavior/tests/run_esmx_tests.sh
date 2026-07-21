#!/bin/bash

TEST_DIR="$PWD"
MAIN_DIR="${TEST_DIR}/../"

TEST1="test1" 
TEST2="test2" 
TEST3="test3" 

echo "Test dir: ${TEST_DIR}"
echo "Main dir: ${MAIN_DIR}"

# module use ${MAIN_DIR}/modules
# module load cheyenne
# #module list

# #------------------------------------------------------------------------------
# # set ESMF_ESMXDIR using ESMFMKFILE
# if [ ! -f "${ESMFMKFILE}" ]; then
#   echo "ERROR: ESMFMKFILE does not exists."
#   exit 1
# fi
# ESMF_ESMXDIR=`grep "ESMF_ESMXDIR" ${ESMFMKFILE}`
# export ESMF_ESMXDIR=${ESMF_ESMXDIR#*=}
#------------------------------------------------------------------------------

echo "Starting tests..."

for mytest in ${TEST2} ${TEST1} ${TEST3}
do 
    echo "------------------------------------------------------------------------------" 
    echo "Running ${mytest}"
    
    ./run_any_test.sh -t=${mytest} --esmx
    # cat ${mytest}/PET0.ESMF_LogFile
    # what needs to be checked here?

done    

exit 0
