#!/usr/bin/perl
#---------------------------------------------------
#wrf_rename.plx
#rename wrf out files from :00:00 to _00_00
@files=<*>;
foreach $file (@files){

if ($file =~ /:/) {
#	print "$file \n";
	($newfile = $file) =~s/:/_/g;
#	print "$newfile \n";
}
if (-e $newfile){ #if file already exists
	warn "cannot rename $file to $newfile: $newfile exists\n";
}elsif (rename $file, $newfile){ #rename file
	##success, do nothing
}else{ #could not rename the file
	warn "rename $file to $newfile failed: $!\n";
}
}
