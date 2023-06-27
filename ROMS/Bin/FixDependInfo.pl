# git $Id$
# svn $Id: FixDependInfo.pl 1151 2023-02-09 03:08:53Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# PERL function to remove .F files the compiler dependency information
# when using CMake. The dependency generator is easily confused by CPP
# if-directives. Only the preprocesed .f90 files are used to determine
# the dependencies.

use strict;
use warnings;
use File::Copy;

open my $in, '<', $ARGV[0] or die "Cannot open $ARGV[0] : $!";

# Create a temporary file to write the lines which don't match pattern.

open my $tmp, ">", "tmp" or die "Unable to open file tmp: $!\n";

while(my $line = <$in>)
{
    next if $line =~ m/\.F"/;   # ignore lines which matches pattern
    print $tmp $line;           # write the rest lines in temp file
}

close $in;
close $tmp;
move("tmp", "$ARGV[0]");        # move temp file into original file.
