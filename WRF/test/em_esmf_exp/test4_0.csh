#!/bin/csh
# @ job_type		= parallel
# @ environment		= COPY_ALL;MP_EUILIB=us
# @ job_name		= test4_0.$(jobid)
# @ output		= test4_0.$(jobid).out
# @ error		= test4_0.$(jobid).err
# @ network.MPI		= csss,shared,us
# @ node_usage		= shared
# @ checkpoint		= no
# @ wall_clock_limit	= 0:30:00
# @ node		= 1
# @ total_tasks		= 4
# @ class		= share
# @ ja_report		= yes
# @ queue

# This is a script to run WRF for a simple MPI-only em_real case 
# on bluesky.  jan00 case is used.  
# All settings are hard-coded.  
# To use it, type "llsubmit test4_0.csh" from the 
# test/em_real subdirectory.  

unalias cd cp rm ls pushd popd mv

#####################################################################

# Set up info for bluesky
set Num_Procs		= 4     # best to match with "#@total_tasks" above...  
setenv MP_PROCS  $Num_Procs
setenv MP_RMPOOL 1
set MPIRUNCOMMAND       =  poe 
setenv MP_SHARED_MEMORY yes
setenv OMP_NUM_THREADS 1
setenv XLSMPOPTS "parthds=1"

# Modify namelist if needed...  
cp -f namelist.input.jan00 namelist.input

# Run WRF
$MPIRUNCOMMAND wrf.exe

