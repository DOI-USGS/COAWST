#!/bin/csh

if ( ( ! -e ../../external/io_netcdf/diffwrf ) || ( ! -e ./wrf.exe ) ) then
	echo missing executables
	echo "missing executables" >>! results
	exit ( 2 )
endif

set dir = `pwd`
set tail_dir = $dir:t
if ( ( $tail_dir == em_real ) || ( $tail_dir == em_quarter_ss ) || ( $tail_dir == em_b_wave ) ) then
# no op
else
	echo this needs to be run from one of the em directories
	echo "this needs to be run from one of the em directories" >> results
	exit ( 4 ) 
endif

set which_em_case = $1

if      ( $which_em_case == em_real ) then
	set cu = $2
	set mp = $3
	set ra = $4
	set sf = $5
	set bl = $6
	set ls = $7
	set ns = $8
	set datecase =  2000-01-24_12:00:00
else if ( $which_em_case == em_quarter_ss ) then
	./run_me_first.csh  >& /dev/null
	set mp = $2
	set di = $3
	set km = $4
	set da = $5
	set px = $6
	set ox = $7
	set py = $8
	set oy = $9
	set datecase =  0001-01-01_00:00:00
else if ( $which_em_case == em_b_wave ) then
	./run_me_first.csh  >& /dev/null
	set mp = $2
	set nh = $3
	set datecase =  0001-01-01_00:00:00
endif

setenv MP_RMPOOL 1

cat >! machfile << EOF
`hostname`
`hostname`
`hostname`
`hostname`
EOF

if ( -e namelist.input ) rm namelist.input

if      ( $which_em_case == em_real ) then
	m4 -DCU_OPTION=$cu -DMP_OPTION=$mp -DRA_OPTION=$ra -DSF_OPTION=$sf -DBL_OPTION=$bl \
	   -DLS_OPTION=$ls -DNS_OPTION=$ns \
	   namelist.input.template.$which_em_case >! namelist.input
else if ( $which_em_case == em_quarter_ss ) then
	m4 -DMP_OPTION=$mp -DDIFF_OPT=$di -DKM_OPT=$km -DDAMP_OPT=$da \
	   -DPERIODIC_X=$px -DOPEN_XS=$ox -DOPEN_XE=$ox \
	   -DPERIODIC_Y=$py -DOPEN_YS=$oy -DOPEN_YE=$oy \
	   namelist.input.template.$which_em_case >! namelist.input
else if ( $which_em_case == em_b_wave ) then
	m4 -DMP_OPTION=$mp -DNH_OPTION=$nh \
	   namelist.input.template.$which_em_case >! namelist.input
endif

rm rsl* wrfo* >& /dev/null
			
echo " "
echo " " >> results
set foo = ( `date` )
echo $foo
echo "$foo" >> results

set procs = ( 1 2 )
set procs = ( 1 4 )

foreach p ( $procs )

	if      ( $which_em_case == em_real ) then 
		set opts =  (       cu=${cu} mp=${mp} ra=${ra} sf=${sf} bl=${bl} ls=${ls} )
		set dir  =    p${p}_cu=${cu}_mp=${mp}_ra=${ra}_sf=${sf}_bl=${bl}_ls=${ls}
		set dir2 =          cu=${cu}_mp=${mp}_ra=${ra}_sf=${sf}_bl=${bl}_ls=${ls}
	else if ( $which_em_case == em_quarter_ss ) then
		set opts = (        mp=${mp} di=${di} km=${km} da=${da} px=${px} ox=${ox} py=${py} oy=${oy} )
		set dir  =    p${p}_mp=${mp}_di=${di}_km=${km}_da=${da}_px=${px}_ox=${ox}_py=${py}_oy=${oy}
		set dir2 =          mp=${mp}_di=${di}_km=${km}_da=${da}_px=${px}_ox=${ox}_py=${py}_oy=${oy}
	else if ( $which_em_case == em_b_wave ) then
		set opts = (        mp=${mp} nh=${nh} )
		set dir  =    p${p}_mp=${mp}_nh=${nh}
		set dir2 =          mp=${mp}_nh=${nh}
	endif

	if ( `uname` == AIX ) then 
		set RUNDIFFWRFCOMMAND = ../../external/io_netcdf/diffwrf
		set RUNWRFCOMMAND     = ./wrf.exe
		if      ( $which_em_case == em_real ) then
			set RUNINITCOMMAND    = ./real.exe
		else if ( $which_em_case == em_quarter_ss ) then
			set RUNINITCOMMAND    = ./ideal.exe
		else if ( $which_em_case == em_b_wave ) then
			set RUNINITCOMMAND    = ./ideal.exe
		endif
        else if ( `uname` == Linux ) then
                set RUNDIFFWRFCOMMAND = ( /usr/local/mpich/bin/mpirun -np 1 ../../external/io_netcdf/diffwrf )
                set RUNWRFCOMMAND     = ( /usr/local/mpich/bin/mpirun -np $p ./wrf.exe )
                if      ( $which_em_case == em_real ) then
                        set RUNINITCOMMAND    = ( /usr/local/mpich/bin/mpirun -np $p ./real.exe )
                else if ( $which_em_case == em_quarter_ss ) then
                        set RUNINITCOMMAND    = ( /usr/local/mpich/bin/mpirun -np 1  ./ideal.exe )
                else if ( $which_em_case == em_b_wave ) then
                        set RUNINITCOMMAND    = ( /usr/local/mpich/bin/mpirun -np 1  ./ideal.exe )
                endif
	else if ( `uname` != AIX ) then 
		set RUNDIFFWRFCOMMAND = ( /usr/local/mpich/bin/mpirun -np 1 -machinefile ./machfile ../../external/io_netcdf/diffwrf )
		set RUNWRFCOMMAND     = ( /usr/local/mpich/bin/mpirun -np $p -machinefile ./machfile ./wrf.exe )
		if      ( $which_em_case == em_real ) then
			set RUNINITCOMMAND    = ( /usr/local/mpich/bin/mpirun -np $p -machinefile ./machfile ./real.exe )
		else if ( $which_em_case == em_quarter_ss ) then
			set RUNINITCOMMAND    = ( /usr/local/mpich/bin/mpirun -np 1  -machinefile ./machfile ./ideal.exe )
		else if ( $which_em_case == em_b_wave ) then
			set RUNINITCOMMAND    = ( /usr/local/mpich/bin/mpirun -np 1  -machinefile ./machfile ./ideal.exe )
		endif
	endif

	setenv MP_PROCS $p
	echo running p=$p $opts
	echo "running p=$p $opts" >> results
	if ( -d $dir ) rm -rf $dir
	mkdir $dir

	if      ( $which_em_case == em_real   ) then
		setenv MP_PROCS $p
	else if ( $which_em_case != em_real   ) then
		setenv MP_PROCS 1
	endif

	$RUNINITCOMMAND >& /dev/null
	rm rsl* >& /dev/null

	setenv MP_PROCS $p
	$RUNWRFCOMMAND >& /dev/null

	if ( `uname` == AIX ) then
		mv rsl* wrfo* $dir
	else
		mv      wrfo* $dir
	endif
end
			
if ( -e fort.88 ) rm fort.88
if ( -e fort.98 ) rm fort.98

set file1 =  p$procs[1]_$dir2/wrfout_d01_$datecase
set file2 =  p$procs[2]_$dir2/wrfout_d01_$datecase
			
setenv MP_PROCS 1
$RUNDIFFWRFCOMMAND $file1 $file2 >& /dev/null

set size1 = `ls -ls $file1 | awk ' { print $6 } '`
set size2 = `ls -ls $file2 | awk ' { print $6 } '`
				
if ( ( -e fort.88 ) || ( $size1 != $size2 ) || ( $size1 <= 10000 ) ) then	
	echo probs on 1 vs 2 on coarse $opts
	echo "probs on 1 vs 2 on coarse $opts" >> results
else
	echo OK on coarse $opts
	echo "OK on coarse $opts" >> results
endif
			
if ( -e fort.88 ) rm fort.88
if ( -e fort.98 ) rm fort.98

set file1 =  p${procs[1]}_$dir2/wrfout_d02_$datecase
set file2 =  p${procs[2]}_$dir2/wrfout_d02_$datecase

$RUNDIFFWRFCOMMAND $file1 $file2 >& /dev/null

set size1 = `ls -ls $file1 | awk ' { print $6 } '`
set size2 = `ls -ls $file2 | awk ' { print $6 } '`
				
if ( ( -e fort.88 ) || ( $size1 != $size2 ) || ( $size1 <= 10000 ) ) then	
	echo probs on 1 vs 2 on nest $opts
	echo "probs on 1 vs 2 on $opts" >> results
	exit ( 1 )
else
	echo OK on fine $opts
	echo "OK on fine $opts" >> results
	exit ( 0 ) 
endif
