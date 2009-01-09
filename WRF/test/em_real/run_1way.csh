#!/bin/csh

#	Run the WRF 1-way using the ndown program to generate a nested forecast.
#	This shell assumes that you have already built the WRF code and that the
#	real.exe, ndown.exe, and wrf.exe files exist.  It assumes that the input
#	data for the Jan00 case is in a "normal" location.  This code is executed
#	from the ./WRFV2/test/em_real directory.

#	Steps in the shell:
#	1. make coarse grid forecast - 24 h Jan 00, 30 km 
#	2. make fine grid ic, only a single time period
#	3. use ndown to make fine grid input, use hourly BC
#	4. make fine grid forecast, 24 h, 10 km nested domain
#	5. a last check of whether it all worked

#goto number_one
#goto number_two
#goto number_three
#goto number_four
#goto number_five

number_one:
#	1.
#	Make the coarse grid namelist: 24 h, output every hour, and use the Noah LSM

echo 1.1
sed -e '/^ run_hours/s/12/24/' \
    -e '/^ history_interval/s/180/60/' \
    -e '/^ sf_surface_physics/s/1/2/' \
    -e '/^ num_soil_layers/s/5/4/' \
    namelist.input.jan00 >! namelist.input

#	Get the coarse grid input data from the SI

if      ( ( `uname` == AIX ) && ( `hostname | cut -c 1-2` == bs ) ) then
	ln -sf /mmm/users/gill/WRF-data-EM/jan00/*d01* .
else if ( ( `uname` == OSF1 ) && ( `hostname | cut -c 1-6` == joshua ) ) then
	ln -sf /users/gill/WRF-data-EM/jan00/*d01* .
else if ( ( `uname` == Linux ) && ( `hostname` == bay-mmm ) ) then
	ln -sf /users/gill/WRF-data-EM/jan00/*d01* .
endif

#	Run coarse grid real

echo 1.2
if      ( ( `uname` == AIX ) && ( `hostname | cut -c 1-2` == bs ) ) then
	setenv MP_PROCS 1
	setenv MP_RMPOOL 1
	real.exe
else if ( ( `uname` == OSF1 ) && ( `hostname | cut -c 1-6` == joshua ) ) then
	real.exe
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

echo 1.3
set ic_size = `ls -ls wrfinput_d01 | awk ' { print $6 } '`
set bc_times = `ncdump -h wrfbdy_d01 | grep "Time = UNLIMITED" | cut -d"(" -f2 | cut -dc -f1`
if ( ( $ic_size > 10000 ) && ( $bc_times == 4 ) ) then
	rm rsl*
else
	echo coarse grid ic bc wrong size
	exit ( 2 ) 
endif

#	Run coarse grid wrf

echo 1.4
if      ( ( `uname` == AIX ) && ( `hostname | cut -c 1-2` == bs ) ) then
	setenv MP_PROCS 2
	wrf.exe
else if ( ( `uname` == OSF1 ) && ( `hostname | cut -c 1-6` == joshua ) ) then
	wrf.exe
else if ( ( `uname` == Linux ) && ( `hostname` == bay-mmm ) ) then
	mpirun -np 2 -machinefile ./machfile wrf.exe
endif

#	Well, did the coarse grid wrf work?

echo 1.5
set m_times = `ncdump -h wrfout_d01_2000-01-24_12:00:00 | grep "Time = UNLIMITED" | cut -d"(" -f2 | cut -dc -f1`
if ( $m_times == 25 ) then
	rm rsl*
else
	echo coarse grid wrf output wrong size
	exit ( 3 ) 
endif

##################################################################################

number_two:
#	2.
#	Make the fine grid namelist for a 1-time run for real

echo 2.1
sed -e '/^ end_day/s/25/24/' \
    -e '/^ sf_surface_physics/s/1/2/' \
    -e '/^ num_soil_layers/s/5/4/' \
    -e '/^ e_we/s/74,//' \
    -e '/^ e_sn/s/61,//' \
    -e '/^ dx/s/30000/10000/' \
    -e '/^ dy/s/30000/10000/' \
    namelist.input.jan00 >! namelist.input

#	Get the fine grid input data from the SI

echo 2.2
rm -rf wrf_real*d01* met_*d01*
rm -rf wrfi* wrfb*
if      ( ( `uname` == AIX ) && ( `hostname | cut -c 1-2` == bs ) ) then
	rm rsl*
	ln -sf /mmm/users/gill/WRF-data-EM/jan00/met_*d02*24_12:00:00.nc met_em.d01.2000-01-24_12:00:00.nc
else if ( ( `uname` == OSF1 ) && ( `hostname | cut -c 1-6` == joshua ) ) then
	ln -sf /users/gill/WRF-data-EM/jan00/met*d02*24_12:00:00.nc met_em.d01.2000-01-24_12:00:00.nc
else if ( ( `uname` == Linux ) && ( `hostname` == bay-mmm ) ) then
	ln -sf /users/gill/WRF-data-EM/jan00/met*d02*24_12:00:00.nc met_em.d01.2000-01-24_12:00:00.nc
endif

#	Run fine grid real

echo 2.3
if      ( ( `uname` == AIX ) && ( `hostname | cut -c 1-2` == bs ) ) then
	setenv MP_PROCS 1
	setenv MP_RMPOOL 1
	real.exe
else if ( ( `uname` == OSF1 ) && ( `hostname | cut -c 1-6` == joshua ) ) then
	real.exe
else if ( ( `uname` == Linux ) && ( `hostname` == bay-mmm ) ) then
	mpirun -np 1 -machinefile ./machfile real.exe
endif

#	Well, did the fine grid ic (remember, 1-time, no bc) work?

echo 2.4
set ic_size = `ls -ls wrfinput_d01 | awk ' { print $6 } '`
if ( $ic_size > 10000 ) then
	rm rsl*
	rm met_em*
	mv wrfinput_d01 wrfndi_d02 # tricky, I know
else
	echo fine grid ic 1 time wrong size
	exit ( 12 ) 
endif


##################################################################################

number_three:
#	3.
#	Make the ndown namelist: 24 h, output every hour, and use the Noah LSM

echo 3.1
sed -e '/^ interval_seconds/s/21600/3600/' \
    -e '/^ max_dom/s/1/2/' \
    -e '/^ sf_surface_physics/s/1/2/' \
    -e '/^ sf_surface_physics/s/1/2/' \
    -e '/^ num_soil_layers/s/5/4/' \
    namelist.input.jan00 >! namelist.input

#	Run ndown

echo 3.2
if      ( ( `uname` == AIX ) && ( `hostname | cut -c 1-2` == bs ) ) then
	setenv MP_PROCS 1
	setenv MP_RMPOOL 1
	ndown.exe
else if ( ( `uname` == OSF1 ) && ( `hostname | cut -c 1-6` == joshua ) ) then
	ndown.exe
else if ( ( `uname` == Linux ) && ( `hostname` == bay-mmm ) ) then
	mpirun -np 1 -machinefile ./machfile ndown.exe
endif

#	Well, did ndown work?

echo 3.3
set ic_size = `ls -ls wrfinput_d02 | awk ' { print $6 } '`
set bc_times = `ncdump -h wrfbdy_d02 | grep "Time = UNLIMITED" | cut -d"(" -f2 | cut -dc -f1`
if ( ( $ic_size > 10000 ) && ( $bc_times == 24 ) ) then
	rm rsl*
	mv wrfinput_d02 wrfinput_d01
	mv wrfbdy_d02 wrfbdy_d01
	if ( ! -d hold ) mkdir hold
	mv wrfo* wrfndi_d02 hold
else
	echo ndown output wrong size
	exit ( 22 ) 
endif

##################################################################################

number_four:
#	4.
#	Make the fine grid namelist for the full forecast run

echo 4.1
sed -e '/^ run_hours/s/12/24/' \
    -e '/^ interval_seconds/s/21600/3600/' \
    -e '/^ history_interval/s/180/60/' \
    -e '/^ time_step/s/180/60/' \
    -e '/^ e_we/s/74,//' \
    -e '/^ e_sn/s/61,//' \
    -e '/^ dx/s/30000,//' \
    -e '/^ dy/s/30000,//' \
    -e '/^ sf_surface_physics/s/1/2/' \
    -e '/^ num_soil_layers/s/5/4/' \
    namelist.input.jan00 >! namelist.input

#	Run fine grid wrf

echo 4.2
if      ( ( `uname` == AIX ) && ( `hostname | cut -c 1-2` == bs ) ) then
	setenv MP_PROCS 2
	#wrf.exe
	echo submitting wrf job in batch queue
	m4 -DPWD=`pwd` run.csh.template >! run.csh
	llsubmit run.csh
	echo need to manually see if it worked
	echo 'set m_times = `ncdump -h wrfout_d01_2000-01-24_12:00:00 | grep "Time = UNLIMITED" | cut -d"(" -f2 | cut -dc -f1`'
	echo 'if $m_times == 25, you are OK'
	exit ( 0 ) 
else if ( ( `uname` == OSF1 ) && ( `hostname | cut -c 1-6` == joshua ) ) then
	wrf.exe
else if ( ( `uname` == Linux ) && ( `hostname` == bay-mmm ) ) then
	mpirun -np 2 -machinefile ./machfile wrf.exe
endif

##################################################################################

number_five:
#	5.
#	Well, did the fine grid wrf work?

echo 5.1
set m_times = `ncdump -h wrfout_d01_2000-01-24_12:00:00 | grep "Time = UNLIMITED" | cut -d"(" -f2 | cut -dc -f1`
if ( $m_times == 25 ) then
	rm rsl*
else
	echo fine grid wrf output wrong size
	exit ( 43 ) 
endif
