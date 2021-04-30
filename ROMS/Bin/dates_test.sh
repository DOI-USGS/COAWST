#!/bin/bash
#
# svn $Id: dates_test.sh 1054 2021-03-06 19:47:12Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2021 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# Script to test and show usage of the "dates" perl script.             :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

today=`date +"%d-%m-%Y %r"`

dn1=`dates datenum`
ds1=`dates numdate $dn1`
yd1=`dates yday $ds1`

dn2=`dates datenum 1900-01-01`
ds2=`dates numdate $dn2`
yd2=`dates yday $ds2`
d21=`dates daysdiff $ds2 $ds1`

dn2=`dates datenum 19000101`
ds2=`dates numdate $dn2`
yd2=`dates yday $ds2`
d21=`dates daysdiff $ds2 $ds1`

dn3=`dates datenum 1968-05-23`
ds3=`dates numdate $dn3`
yd3=`dates yday $ds3`
d31=`dates daysdiff $ds3 $ds1`

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
