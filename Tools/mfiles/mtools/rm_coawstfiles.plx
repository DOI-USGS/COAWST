#!/urs/bin/perl
#------------------------------------------------------
#remove files when restarting coawst/roms systems
#removes PRINT, Errfile, .mat and .dat files
#------------------------------------------------------

@matfile=<*.mat*>;
#print@matfile;
foreach $file(@matfile){
unlink($file);
}

@datfile=<*.dat*>;
#print @datfile;
foreach $file(@datfile){
unlink($file);
}

@Errfile=<Err*>;
#print @Errfile;
foreach $file(@Errfile){
unlink($file);
}

@Prifile=<PRI*>;
#print @Prifile;
foreach $file(@Prifile){
unlink($file);
}