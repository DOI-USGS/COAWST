#!/bin/bash
#
# git $Id$
# svn $Id: check_openfiles.sh 1151 2023-02-09 03:08:53Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# ROMS BASH Script to check open files                                  :::
#                                                                       :::
# In the UNIX environment, there is a limit to the number of open files :::
# during program execution. Use the commands to check such limit:       :::
#                                                                       :::
#   limit                                                               :::
#   ulimit -a                                                           :::
#   ulimit -S -n                                                        :::
#                                                                       :::
# Usually, 256 files can be openned by default. If the number of open   :::
# files is exceeded, you will get the 'Too many open files' error.      :::
#                                                                       :::
# For example, in Linux we can change the default number:               :::
#                                                                       :::
#   limit descriptors 2048   or any other value                         :::
#                                                                       :::
# The C-preprocessing option CHECK_OPEN_FILES in ROMS can be used to    :::
# report the number of files created, opened, and closed for an         :::
# application. The report is written to Fortran file "fort.1000".       :::
#                                                                       :::
# Script to update the copyright information on ROMS source files.      :::
# This script replaces the copyright string in the source files and     :::
# updates the copyright svn property. This script must be executed      :::
# from top level of the ROMS source code.                               :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    ./ROMS/Bin/check_openfiles.sh                                      :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Set report file.

report="fort.1000"

# Get the number of files created, opened, and closed.

CREATED=`grep CREATE ${report} | wc -l`
OPENED=`grep OPEN  ${report} | wc -l`
CLOSED=`grep CLOSE ${report} | wc -l`

# Report ROMS I/O NetCDF files.

echo " "
grep CREATE ${report}
echo " "
grep _obs   ${report}
echo " "
grep _ini   ${report}
echo " "
grep _swrad ${report}
echo " "
grep _bry   ${report}
echo " "
grep _clm   ${report}

echo " "
grep _adj ${report}
grep _avg ${report}
grep _dai ${report}
grep _dia ${report}
grep _flt ${report}
grep _fwd ${report}
grep _grd ${report}
grep _gri ${report}
grep _his ${report}
grep _irp ${report}
grep _itl ${report}
grep _mod ${report}
grep _qck ${report}
grep _rst ${report}
grep _sta ${report}
grep _tlf ${report}
grep _tlm ${report}

echo " "
echo "*** Number of opened  files = ${OPENED}"
echo "*** Number of closed  files = ${CLOSED}"
echo "*** Number of created files = ${CREATED}"
