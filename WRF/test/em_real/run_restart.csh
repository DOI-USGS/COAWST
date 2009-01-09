#!/bin/csh

#	Run the WRF restart with 2 domains.

#	This shell assumes that you have already built the WRF code for nesting
#	and that the real.exe and wrf.exe files exist.  It assumes that the input
#	data for the Jan00 case is in a "normal" location.  This code is executed
#	from the ./WRFV2/test/em_real directory.

#	Steps in the shell:
#	1. make coarse+fine grid forecast - 2 CG timesteps
#	2. restart coarse+fine grid forecast - 1 CG timestep
#	3. diff orig vs restart final fcst time

##################################################################################

#	1.
#	Make the coarse + fine grid namelist: 6 min forecast, output every 3 min, 
#	dump restart at 3 min intervals

echo 1.1
sed -e '/^ run_hours/s/12/0/' \
    -e '/^ run_minutes/s/0/6/' \
    -e '/^ frames_per_outfile/s/1000, 1000,/1, 1,/' \
    -e '/^ history_interval/s/180,  60,/3, 3,/' \
    -e '/^ restart_interval/s/5000/3/' \
    -e '/^ sf_surface_physics/s/1/2/g' \
    -e '/^ num_soil_layers/s/5/4/' \
    -e '/^ max_dom/s/1/2/' \
    -e '/^ input_from_file/s/false/true/' \
    namelist.input.jan00 >! namelist.input

#	Get the coarse grid input data from the SI

echo 1.2
if      ( ( `uname` == AIX ) && ( `hostname | cut -c 1-2` == bs ) ) then
	ln -sf /mmm/users/gill/WRF-data-EM/jan00/wrf_real* .
else if ( ( `uname` == OSF1 ) && ( `hostname | cut -c 1-6` == joshua ) ) then
	ln -sf /users/gill/WRF-data-EM/jan00/wrf_real* .
else if ( ( `uname` == Linux ) && ( `hostname` == bay-mmm ) ) then
	ln -sf /users/gill/WRF-data-EM/jan00/wrf_real* .
endif

#	Run coarse + fine grid real

echo 1.3
if      ( ( `uname` == AIX ) && ( `hostname | cut -c 1-2` == bs ) ) then
	setenv MP_PROCS 1
	setenv MP_RMPOOL 1
	real.exe
else if ( ( `uname` == OSF1 ) && ( `hostname | cut -c 1-6` == joshua ) ) then
	cat >! machfile << EOF
`hostname`
`hostname`
`hostname`
`hostname`
EOF
	mpirun -np 1 -machinefile ./machfile real.exe
else if ( ( `uname` == Linux ) && ( `hostname` == bay-mmm ) ) then
	cat >! machfile << EOF
`hostname`
`hostname`
`hostname`
`hostname`
EOF
	mpirun -np 1 -machinefile ./machfile real.exe
endif

#	Well, did the coarse grid ic/bc work?

echo 1.4
set ic_size = `ls -ls wrfinput_d01 | awk ' { print $6 } '`
set bc_times = `ncdump -h wrfbdy_d01 | grep "Time = UNLIMITED" | cut -d"(" -f2 | cut -dc -f1`
if ( ( $ic_size > 10000 ) && ( $bc_times == 4 ) ) then
	rm rsl*
else
	echo coarse grid ic bc wrong size
	exit ( 2 ) 
endif

#	Run coarse+fine grid wrf

echo 1.5
if      ( ( `uname` == AIX ) && ( `hostname | cut -c 1-2` == bs ) ) then
	setenv MP_PROCS 1
	wrf.exe
else if ( ( `uname` == OSF1 ) && ( `hostname | cut -c 1-6` == joshua ) ) then
	mpirun -np 2 -machinefile ./machfile wrf.exe
else if ( ( `uname` == Linux ) && ( `hostname` == bay-mmm ) ) then
	mpirun -np 2 -machinefile ./machfile wrf.exe
endif

#	Well, did the coarse+fine grid wrf work?

echo 1.6
set num_files = `ls -ls wrfout* | wc -l`
if ( $num_files == 6 ) then
	echo first wrf run worked, generated restart file
	rm rsl*
	if ( ! -d 1 ) mkdir 1
	mv wrfo* 1
else
	echo coarse grid wrf output wrong size
	exit ( 3 ) 
endif

##################################################################################

#	2.
#	Make the CG/FG forecast, starting at the first restart time.

echo 2.1
sed -e '/^ run_hours/s/12/0/' \
    -e '/^ run_minutes/s/0/3/' \
    -e '/^ start_minute/s/00,   00,/3, 3, /' \
    -e '/^ frames_per_outfile/s/1000, 1000,/1, 1,/' \
    -e '/^ history_interval/s/180,  60,/3, 3,/' \
    -e '/^ restart/s/.false.,/.true.,/' \
    -e '/^ sf_surface_physics/s/1/2/g' \
    -e '/^ num_soil_layers/s/5/4/' \
    -e '/^ max_dom/s/1/2/' \
    -e '/^ input_from_file/s/false/true/' \
    namelist.input.jan00 >! namelist.input

#	Run coarse+fine grid wrf

echo 2.2
if      ( ( `uname` == AIX ) && ( `hostname | cut -c 1-2` == bs ) ) then
	setenv MP_PROCS 1
	wrf.exe
else if ( ( `uname` == OSF1 ) && ( `hostname | cut -c 1-6` == joshua ) ) then
	mpirun -np 2 -machinefile ./machfile wrf.exe
else if ( ( `uname` == Linux ) && ( `hostname` == bay-mmm ) ) then
	mpirun -np 2 -machinefile ./machfile wrf.exe
endif

#	Well, did the coarse+fine grid wrf work?

echo 2.3
set num_files = `ls -ls wrfout* | wc -l`
if ( $num_files == 4 ) then
	echo second wrf run worked, used restart file
	rm rsl*
	if ( ! -d 2 ) mkdir 2
	mv wrfo* 2
else
	echo coarse grid wrf output wrong size
	exit ( 3 ) 
endif

##################################################################################

#	3.
#	Diff the output from the first fcst and the second.

echo 3.1
echo CG: orig fcst vs restart fcst
if ( -e fort.88 ) rm fort.88
if ( -e fort.98 ) rm fort.98
../../external/io_netcdf/diffwrf 1/wrfout_d01_2000-01-24_12:06:00 2/wrfout_d01_2000-01-24_12:06:00 >& /dev/null
if ( ( -e fort.88 ) || ( -e fort.98 ) ) then
	echo FAIL restart bit for bit test for CG
else
	echo PASS restart bit for bit test for CG
endif

echo 3.2
echo FG: orig fcst vs restart fcst 
if ( -e fort.88 ) rm fort.88
if ( -e fort.98 ) rm fort.98
../../external/io_netcdf/diffwrf 1/wrfout_d02_2000-01-24_12:06:00 2/wrfout_d02_2000-01-24_12:06:00 >& /dev/null
if ( ( -e fort.88 ) || ( -e fort.98 ) ) then
	echo FAIL restart bit for bit test for FG
else
	echo PASS restart bit for bit test for FG
endif
