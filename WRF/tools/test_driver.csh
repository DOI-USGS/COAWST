#!/bin/csh

#PBS -V -A wrfijet
#PBS -lnodes=4:comp -l walltime=3600

if ( ( `uname` == Linux ) && ( `hostname | cut -d. -f 2-` == fsl.noaa.gov ) ) then
	source /usr/local/bin/setup-mpi.csh
	cd /misc/p30/wrfijet/users/michalak/WRFV1/test/em_real
endif

if ( -e results ) rm results
touch  results

#	real-data cases have additional options

set method = LONG
set method = QUICK 
set method = MEDIUM

#	Make sure we are running the right case for the 
#	directory we are in.

set dir = `pwd`
set tail_dir = $dir:t
if      ( $tail_dir == em_real ) then
	set which_em_case = em_real
else if ( $tail_dir == em_quarter_ss ) then
	set which_em_case = em_quarter_ss
else if ( $tail_dir == em_b_wave ) then
	set which_em_case = em_b_wave
else
	echo this needs to be run from one of the em directories, with WRF already built
	echo "this needs to be run from one of the em directories, with WRF already built" >> results
	exit ( 1 ) 
endif

#	More nice print outs.

echo "$tail_dir for `uname` on `date`"
echo "$tail_dir for `uname` on `date`" >> results

#	Three options: real, quarter_ss, b_wave.

if      ( $which_em_case == em_real ) then
	if ( $method == LONG ) then
		foreach mp ( 0 1 2 3 4 5 6 )
			foreach cu ( 1 2 3 )
				foreach ra ( 1 2 ) 
					foreach sf ( 1 2 ) 
						foreach bl ( 1 2 ) 
							foreach ls ( 1 2 3 ) 
								if ( $ls == 1 ) set ns = 5
								if ( $ls == 2 ) set ns = 4
								if ( $ls == 3 ) then
									set ns = 6 ; set bl = 2 ; set sf = 2
								endif
								if ( ( ( $sf == 1 ) && ( $bl == 2 ) ) || \
								     ( ( $sf == 2 ) && ( $bl == 1 ) ) ) then
		#							skip invalid phys option combination
								else
									./test_inner.csh $which_em_case $cu $mp $ra $sf $bl $ls $ns
									set ok = $status
									if ( $ok != 0 ) then
										exit ( 2 )
									endif
								endif
							end
						end
					end
				end
			end
		end
	
	else if ( $method == MEDIUM ) then
		set mp = 1
		set cu = 1
		set ra = 1
		set sf = 1
		set bl = 1
		set ls = 1
		set ns = 5
		foreach mp ( 0 1 2 3 4 5 6 )
			./test_inner.csh $which_em_case $cu $mp $ra $sf $bl $ls $ns
			set ok = $status
			if ( $ok != 0 ) then
				exit ( 2 )
			endif
		end
	
		set mp = 1
		set cu = 1
		set ra = 1
		set sf = 1
		set bl = 1
		set ls = 1
		set ns = 5
		foreach cu ( 1 2 3 )
			./test_inner.csh $which_em_case $cu $mp $ra $sf $bl $ls $ns
			set ok = $status
			if ( $ok != 0 ) then
				exit ( 2 )
			endif
		end
	
		set mp = 1
		set cu = 1
		set ra = 1
		set sf = 1
		set bl = 1
		set ls = 1
		set ns = 5
		foreach ra ( 1 2 )
			./test_inner.csh $which_em_case $cu $mp $ra $sf $bl $ls $ns
			set ok = $status
			if ( $ok != 0 ) then
				exit ( 2 )
			endif
		end
	
		set mp = 1
		set cu = 1
		set ra = 1
		set sf = 1
		set bl = 1
		set ls = 1
		set ns = 5
		foreach sf ( 1 2 )
			if ( $sf == 2 ) set bl = 2
			./test_inner.csh $which_em_case $cu $mp $ra $sf $bl $ls $ns
			set ok = $status
			if ( $ok != 0 ) then
				exit ( 2 )
			endif
		end
	
		set mp = 1
		set cu = 1
		set ra = 1
		set sf = 1
		set bl = 1
		set ls = 1
		set ns = 5
		foreach bl ( 1 2 )
			if ( $bl == 2 ) set sf = 2
			./test_inner.csh $which_em_case $cu $mp $ra $sf $bl $ls $ns
			set ok = $status
			if ( $ok != 0 ) then
				exit ( 2 )
			endif
		end
	
		set mp = 1
		set cu = 1
		set ra = 1
		set sf = 1
		set bl = 1
		set ls = 1
		set ns = 5
		foreach ls ( 1 2 3 )
			if ( $ls == 1 ) set ns = 5
			if ( $ls == 2 ) set ns = 4
			if ( $ls == 3 ) then
				set ns = 6 ; set sf = 2 ; set bl = 2
			endif
			./test_inner.csh $which_em_case $cu $mp $ra $sf $bl $ls $ns
			set ok = $status
			if ( $ok != 0 ) then
				exit ( 2 )
			endif
		end
	
	else if ( $method == QUICK ) then
		set mp = 1
		set cu = 1
		set ra = 1
		set sf = 1
		set bl = 1
		set ls = 1
		set ns = 5
		./test_inner.csh $which_em_case $cu $mp $ra $sf $bl $ls $ns
		set ok = $status
		if ( $ok != 0 ) then
			exit ( 2 )
		endif
	endif

else if ( $which_em_case == em_quarter_ss ) then
	set mp = 0
	set di = 1
	set km = 1
	set da = 0
	set px = F
	set ox = T
	set py = F
	set oy = T
	./test_inner.csh $which_em_case $mp $di $km $da $px $ox $py $oy
	set ok = $status
	if ( $ok != 0 ) then
		exit ( 2 )
	endif

	set mp = 1
	set di = 2
	set km = 2
	set da = 1
	set px = F
	set ox = T
	set py = F
	set oy = T
	./test_inner.csh $which_em_case $mp $di $km $da $px $ox $py $oy
	set ok = $status
	if ( $ok != 0 ) then
		exit ( 2 )
	endif

	set mp = 2
	set di = 2
	set km = 3
	set da = 1
	set px = T
	set ox = F
	set py = T
	set oy = F
	./test_inner.csh $which_em_case $mp $di $km $da $px $ox $py $oy
	set ok = $status
	if ( $ok != 0 ) then
		exit ( 2 )
	endif

else if ( $which_em_case == em_b_wave ) then
	set nh = T
	foreach mp ( 0 1 2 ) 
		./test_inner.csh $which_em_case $mp $nh
		set ok = $status
		if ( $ok != 0 ) then
			exit ( 2 )
		endif
	end

	set mp = 0
	set nh = F
	./test_inner.csh $which_em_case $mp $nh
	set ok = $status
	if ( $ok != 0 ) then
		exit ( 2 )
	endif

endif

set foo = ( `date` )
echo " "
echo " " >> results
echo END $foo
echo "END $foo" >> results
