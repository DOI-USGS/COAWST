#!/bin/csh
# @ job_type		= parallel
# @ environment		= COPY_ALL;MP_EUILIB=us
# @ job_name		= real.$(jobid)
# @ output		= real.$(jobid).out
# @ error		= real.$(jobid).err
# @ network.MPI		= csss,shared,us
# @ node_usage		= shared
# @ checkpoint		= no
# @ wall_clock_limit	= 0:30:00
# @ node		= 1
# @ total_tasks		= 1
# @ class		= share
# @ ja_report		= yes
# @ queue

# This is a script to run WRF real.exe for a simple MPI-only em_real case 
# on bluesky.  jan00 case is used.  
# All settings are hard-coded.  
# To use it, type "llsubmit real.csh" from the 
# test/em_real subdirectory.  

unalias cd cp rm ls pushd popd mv

#####################################################################

# Set up info for bluesky
set Num_Procs		= 1     # best to match with "#@total_tasks" above...  
setenv MP_PROCS  $Num_Procs
setenv MP_RMPOOL 1
set MPIRUNCOMMAND       =  poe 
setenv MP_SHARED_MEMORY yes
setenv OMP_NUM_THREADS 1
setenv XLSMPOPTS "parthds=1"

# Modify namelist if needed...  
cp -f namelist.input.jan00 namelist.input

# Link in SI data sets for em_real jan00 case
set WRFREGDATAEM = /mmm/users/gill/WRF-data-EM
set thedataem = ${WRFREGDATAEM}/jan00
ln -sf $thedataem/* .

# Build "real" input data files for WRF from SI data sets
real.exe >! real.exe.out

