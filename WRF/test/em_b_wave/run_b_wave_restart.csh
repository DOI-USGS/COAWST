#!/bin/csh

#	IBM only

if ( `uname` != AIX ) then
	echo Sorry, only on the IBM right now.
	exit ( 1 )
endif

#	Are we in the right place?

set dir = `pwd`
set tail = $dir:t
if ( $tail != em_b_wave ) then
	echo This script needs to be run from the em_b_wave directory
	exit ( 2 )
endif

#	Is everything ready?

if ( ( ! -e ideal.exe ) || ( ! -e wrf.exe ) ) then
	echo Build the WRF code for em_b_wave, optimized, with either RSL or RSL_LITE
	exit ( 3 )
endif

if ( ! -e wrfinput_d01 ) then
	echo We need to have the ideal.exe run already
	exit ( 4 )
endif

#	Make our runtime scripts for the load leveler

m4 -DDIR=`pwd` template_run_9p.csh >! run_9p.csh
m4 -DDIR=`pwd` template_run_1p.csh >! run_1p.csh

#	1 proc restart job

echo running 1p b_wave restart test
llsubmit run_1p.csh >&! llsub.out

set ok = 0
set in_already = 0
while ( $ok == 0 )
	sleep 10 ; llq -u $USER >&! llq.report
	grep `cat llsub.out | grep '"bs1' | cut -d\" -f2` llq.report >& /dev/null
	set ok = $status
	if ( ( $ok == 0 ) && ( $in_already == 0 ) ) then
		set in_already = 1
		set joe_id = `cat llsub.out | grep '"bs1' | cut -d\" -f2 | cut -d. -f2`
	endif
end

rm fort.* >& /dev/null
../../external/io_netcdf/diffwrf 12h_1p/wrfout_d01_0001-01-01_12:00:00 next6h_1p/wrfout_d01_0001-01-01_12:00:00 >& /dev/null
if ( -e fort.88 ) then
	echo PROBS differences in the 1p b_wave restart results
	exit ( 1 ) 
else
	echo 1p b_wave restarts are bit for bit identical
endif

#	9 proc restart job

echo running 9p b_wave restart test
llsubmit run_9p.csh >&! llsub.out

set ok = 0
set in_already = 0
while ( $ok == 0 )
	sleep 10 ; llq -u $USER >&! llq.report
	grep `cat llsub.out | grep '"bs1' | cut -d\" -f2` llq.report >& /dev/null
	set ok = $status
	if ( ( $ok == 0 ) && ( $in_already == 0 ) ) then
		set in_already = 1
		set joe_id = `cat llsub.out | grep '"bs1' | cut -d\" -f2 | cut -d. -f2`
	endif
end

rm fort.* >& /dev/null
../../external/io_netcdf/diffwrf 12h_9p/wrfout_d01_0001-01-01_12:00:00 next6h_9p/wrfout_d01_0001-01-01_12:00:00 >& /dev/null
if ( -e fort.88 ) then
	echo PROBS differences in the 9p b_wave restart results
	exit ( 2 ) 
else
	echo 9p b_wave restarts are bit for bit identical
endif

#	... and 1p vs 9p

rm fort.* >& /dev/null
../../external/io_netcdf/diffwrf 12h_1p/wrfout_d01_0001-01-01_12:00:00 next6h_9p/wrfout_d01_0001-01-01_12:00:00 >& /dev/null
if ( -e fort.88 ) then
	echo PROBS differences in the 1p vs 9p b_wave restart results
	exit ( 3 ) 
else
	echo 1p vs 9p b_wave restarts are bit for bit identical
        echo " "
        echo SUCCESS
endif
