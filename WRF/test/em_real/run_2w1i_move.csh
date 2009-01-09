#!/bin/csh

#	Run the WRF 2-way using the WRF program to generate the fine grid info
#	from the coarse grid.

#	This shell assumes that you have already built the WRF code and that the
#	real.exe, and wrf.exe files exist.  It assumes that the input
#	data for the Jan00 case is in a "normal" location.  This code is executed
#	from the ./WRFV2/test/em_real directory.

#	Are we allowed to have a moving nest?  Check the build from configure.wrf.  There
#	should be a flag -DMOVE_NESTS on the ARCHFLAGS.

if ( -e ../../configure.wrf ) then
	grep DMOVE_NESTS ../../configure.wrf | grep -v "#" >& /dev/null
	set ok = $status
	if ( $ok == 0 ) then
		#echo things are ok to run moving nest
	else
		echo $ok
		echo "probably did not build for a moving nest ..."
		exit ( 1 )
	endif
endif

#	Make the coarse+fine grid namelist: 24 h, output every hour, and use the Noah LSM

echo 1.1
cat >! move_stuff << EOF
 num_moves                           = 4
 move_id                             = 2 , 2 ,  2 , 2
 move_interval                       = 7, 22 , 37 , 52
 move_cd_x                           = 1 , 1 ,  1 , 1
 move_cd_y                           = 1 , 1 ,  1 , 1
EOF

sed -e '/^ run_hours/s/12/1/' \
    -e '/^ history_interval/s/180/60/' \
    -e '/^ sf_surface_physics/s/1/2/' \
    -e '/^ sf_surface_physics/s/1/2/' \
    -e '/^ num_soil_layers/s/5/4/' \
    -e '/^ smooth/r move_stuff' \
    namelist.input.jan00 >! namelist.input

#	Get the coarse grid input data from the SI

echo 1.2
if      ( ( `uname` == AIX ) && ( `hostname | cut -c 1-2` == bs ) ) then
	ln -sf /mmm/users/gill/WRF-data-EM/jan00/wrf_real*d01* .
else if ( ( `uname` == OSF1 ) && ( `hostname | cut -c 1-6` == joshua ) ) then
	ln -sf /users/gill/WRF-data-EM/jan00/wrf_real*d01* .
else if ( ( `uname` == Linux ) && ( `hostname` == bay-mmm ) ) then
	ln -sf /users/gill/WRF-data-EM/jan00/wrf_real*d01* .
endif

#	Run coarse grid real

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

#	Coarse + fine grid namelist (add max_dom=2).

echo 1.5
sed -e '/^ run_hours/s/12/1/' \
    -e '/^ history_interval/s/180/60/' \
    -e '/^ sf_surface_physics/s/1/2/' \
    -e '/^ sf_surface_physics/s/1/2/' \
    -e '/^ num_soil_layers/s/5/4/' \
    -e '/^ max_dom/s/1/2/' \
    -e '/^ smooth/r move_stuff' \
    namelist.input.jan00 >! namelist.input
rm move_stuff

#	Run coarse+fine grid wrf

echo 1.6
if      ( ( `uname` == AIX ) && ( `hostname | cut -c 1-2` == bs ) ) then
	setenv MP_PROCS 1
	wrf.exe
else if ( ( `uname` == OSF1 ) && ( `hostname | cut -c 1-6` == joshua ) ) then
	mpirun -np 2 -machinefile ./machfile wrf.exe
else if ( ( `uname` == Linux ) && ( `hostname` == bay-mmm ) ) then
	mpirun -np 2 -machinefile ./machfile wrf.exe
endif

#	Well, did the coarse grid wrf work?

echo 1.7
set c_times = `ncdump -h wrfout_d01_2000-01-24_12:00:00 | grep "Time = UNLIMITED" | cut -d"(" -f2 | cut -dc -f1`
set f_times = `ncdump -h wrfout_d02_2000-01-24_12:00:00 | grep "Time = UNLIMITED" | cut -d"(" -f2 | cut -dc -f1`
if ( ( $c_times == 2 ) && ( $f_times == 2 ) ) then
	echo 2w1i_move worked 
	rm rsl*
else
	echo coarse grid wrf output wrong size
	exit ( 3 ) 
endif
