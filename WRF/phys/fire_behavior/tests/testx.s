#!/bin/bash
#
#################################################
#
# Purpose: Test for a artificial wildland fire ingesting data from ESMX_Data
#
#################################################
#
purge_output=1 # 0) No, 1) yes
plot=0 # 0) No, 1) yes
#
#################################################
#
testxp1=1 # Check fire area
testxp2=1 # Check heat output
testxp3=1 # Check latent heat output
testxp4=1 # Check Max heat flux
testxp5=1 # Check Max latent heat flux
#
#################################################
#
rm -f ATMD_final_export.nc
rm -f ATMD_final_import.nc

file_baseline=./testx/test_solution.txt
file_exe=../install/bin/esmx_fire
file_output=testx_output.txt

cp ./testx/namelist.fire .
cp ./testx/geo_em.d01.nc .
cp ./testx/testx.yaml .

rm -f  ./$file_output
if [ -f $file_exe ]
then
  $file_exe testx.yaml > ./$file_output
else
  echo 'Please compile the code first'
  exit 1
fi

n_tests=0
n_test_passed=0
pass=false

#
# ----------------------------------------
#

echo "TEST X:"

if [ $testxp1 -eq 1 ]
then
  n_tests=$(expr $n_tests + 1)
  rm -f ./file1.dat ./file2.dat
  var="Fire area"
  grep "$var" $file_output   | awk '{print $2, $6}' > ./file1.dat
  grep "$var" $file_baseline | awk '{print $2, $6}' > ./file2.dat

  test=$(diff ./file1.dat ./file2.dat | wc -l)
  if [ $test -eq 0 ]
  then
    echo '  TestX.1 PASSED'
    n_test_passed=$(expr $n_test_passed + 1)
  else
    echo '  TestX.1 FAILS'
  fi

fi

#
# ----------------------------------------
#

if [ $testxp2 -eq 1 ]
then
  n_tests=$(expr $n_tests + 1)
  rm -f ./file1.dat ./file2.dat
  var="Heat output"   # 6
  grep "$var" $file_output   | awk '{print $2, $6}' > ./file1.dat
  grep "$var" $file_baseline | awk '{print $2, $6}' > ./file2.dat

  test=$(diff ./file1.dat ./file2.dat | wc -l)
    # Here we allow one difference since we are not expecting bit4bit results
  if [ $test -eq 0 ]
  then
    echo '  TestX.2 PASSED'
    n_test_passed=$(expr $n_test_passed + 1)
  else
    echo '  TestX.2 FAILS'
  fi
fi

#
# ----------------------------------------
#

if [ $testxp3 -eq 1 ]
then
  n_tests=$(expr $n_tests + 1)
  rm -f ./file1.dat ./file2.dat
  var="Latent heat output" # 7
  grep "$var" $file_output   | awk '{print $2, $7}' > ./file1.dat
  grep "$var" $file_baseline | awk '{print $2, $7}' > ./file2.dat

  test=$(diff ./file1.dat ./file2.dat | wc -l)
  if [ $test -eq 0 ]
  then
    echo '  TestX.3 PASSED'
    n_test_passed=$(expr $n_test_passed + 1)
  else
    echo '  TestX.3 FAILS'
  fi
fi

#
# ----------------------------------------
#

if [ $testxp4 -eq 1 ]
then
  n_tests=$(expr $n_tests + 1)
  rm -f ./file1.dat ./file2.dat
  var="Max heat flux" # 7
  grep "$var" $file_output   | awk '{print $2, $7}' > ./file1.dat
  grep "$var" $file_baseline | awk '{print $2, $7}' > ./file2.dat

  test=$(diff ./file1.dat ./file2.dat | wc -l)
    # Here we allow one difference since we are not expecting bit4bit results
  if [ $test -eq 0 ]
  then
    echo '  TestX.4 PASSED'
    n_test_passed=$(expr $n_test_passed + 1)
  else
    echo '  TestX.4 FAILS'
  fi
fi

#
# ----------------------------------------
#

if [ $testxp5 -eq 1 ]
then
  n_tests=$(expr $n_tests + 1)
  rm -f ./file1.dat ./file2.dat
  var="Max latent heat flux" # 8
  grep "$var" $file_output   | awk '{print $2, $8}' > ./file1.dat
  grep "$var" $file_baseline | awk '{print $2, $8}' > ./file2.dat

  test=$(diff ./file1.dat ./file2.dat | wc -l)
    # Here we allow one difference since we are not expecting bit4bit results
  if [ $test -lt 4 ]
  then
    echo '  TestX.5 PASSED'
    echo '    Ignore this difference:'
    diff ./file1.dat ./file2.dat
    n_test_passed=$(expr $n_test_passed + 1)
  else
    echo '  TestX.5 FAILS'
  fi
fi

#
# ----------------------------------------
#

  # Purge
rm -f ./namelist.fire.output ./file1.dat ./file2.dat ./geo_em.d01.nc ./namelist.fire testx.yaml
if [ $purge_output -eq 1 ]
then
  rm -rf ./$file_output
  rm -f ./fire_output_2012-06-25_18:00:??.nc
fi

  # Print summary of Test 8
if [ $n_test_passed -eq $n_tests ]
then
  echo "SUCCESS: $n_test_passed PASSED of $n_tests"
  echo ''
  pass=true
else
  echo "FAILED: $n_test_passed PASSED of $n_tests"
  echo ''
  pass=false
fi

rm -f ./latlons.dat ./latlons_c.dat ./wrf_latlons_atm.dat ./wrf_latlons_fire.dat
rm -f PET0.ESMF_LogFile ATMD_final_export.nc ATMD_final_import.nc

if [ $pass = true ]
then
  exit 0
else
  exit 1
fi
