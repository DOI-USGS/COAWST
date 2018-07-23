#!/bin/bash
#
# This routine will extract the Wave Watch 3 tarfile and 
# copy necessary routines for coupling into COAWST.
#
# jcwarner 13July2018
#
# 1) You need to register to get the Wave Watch code. Go here:
# http://polar.ncep.noaa.gov/waves/wavewatch/
# and register to get the code. They will send to you a 
# username and pwd. Follow their instructions and get the tar file.
#
# 2) Copy that tar file to this directory:
# COAWST/WW3
#
# 3) Then you can run this bash file.
# ./add_ww3_to_coawst
#
# End of user input.
#
# First we extract the main tar file.
#
tar -xvf wwatch3.v5.16.tar.gz
#
# use an updatetd version of install_ww3_tar
#
cp coawst_ww3_files/install_ww3_tar .
./install_ww3_tar
#
# now copy files needed by COAWST
#
cp -r coawst_ww3_files/coawst_compile_ww3 .
cp -r coawst_ww3_files/bin/* bin/
cp -r coawst_ww3_files/ftn/* ftn/
cp -r coawst_ww3_files/inp/* inp/
cp -r coawst_ww3_files/regtests/bin/* regtests/bin/
#
cd work
rm -rf ww3*
cd ..
cp -r inp/ww3* work
cp -r coawst_ww3_files/work/* work


