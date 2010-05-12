# --- this script slightly modify LaTeX documentation for SWAN
#     to make appropriate for HTML documentation
@filecontent=<>;
$k=0;
$m=0;
foreach $line (@filecontent) {
  $k++;
  chomp $line;
  if ("$line" eq "\\cleardoublepage") {
     $m=$k+4;
  } elsif ($k>$m || $m==0) {
     if ( "$line" ne "\\tableofcontents" ) {
        print "$line\n";
     }
  }
}
