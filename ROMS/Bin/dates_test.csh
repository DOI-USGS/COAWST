#!/bin/csh -f
#
# svn $Id: dates_test.csh 1054 2021-03-06 19:47:12Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2021 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# Script to test and show usage of the "dates" perl script.             :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set today = `date +"%d-%m-%Y %r"`

set dn1 = `dates datenum`
set ds1 = `dates numdate $dn1`
set yd1 = `dates yday $ds1`

set dn2 = `dates datenum 1900-01-01`
set ds2 = `dates numdate $dn2`
set yd2 = `dates yday $ds2`
set d21 = `dates daysdiff $ds2 $ds1`

set dn2 = `dates datenum 19000101`
set ds2 = `dates numdate $dn2`
set yd2 = `dates yday $ds2`
set d21 = `dates daysdiff $ds2 $ds1`

set dn3 = `dates datenum 1968-05-23`
set ds3 = `dates numdate $dn3`
set yd3 = `dates yday $ds3`
set d31 = `dates daysdiff $ds3 $ds1`

echo
echo "Testing 'dates' Perl Script on $today"
echo
echo "Today's Date:               $ds1"
echo "Today's Date Number:        $dn1"
echo "Today's Day-of-the-year:    $yd1"
echo
echo "Reference Date:             $ds2"
echo "Reference Date Number:      $dn2"
echo "Reference Day-of-the-year:  $yd2"
echo "Days since Reference date:  $d21"
echo
echo "Truncated Date:             $ds3"
echo "Truncated Date Number:      $dn3"
echo "Truncated Day-of-the-year:  $yd3"
echo "Days since TRuncated date:  $d31"
echo
