if (open(P,"uname -s 2>&1 |") && (@tmp = <P>) && close(P))
{
   chop($os = $tmp[$#tmp]);
}
else
{
   $os = $^O;
}
@tmp = ();
close P;

if (open(P,"uname -m 2>&1 |") && (@tmp = <P>) && close(P))
{
   chop($cpu = $tmp[$#tmp]);
}
else
{
   $cpu = $^O;
}
@tmp = ();
close P;

open(OUTFILE,">macros.inc");

if ($os =~ /IRIX64/i) {
  print OUTFILE "##############################################################################\n";
  print OUTFILE "# MIPS_SGI:		SGI Origin 3000 with IRIX using SGI compiler.\n";
  print OUTFILE "##############################################################################\n";
  print OUTFILE "F90_SER = f90\n";
  print OUTFILE "F90_OMP = f90\n";
  print OUTFILE "F90_MPI = f90\n";
  print OUTFILE "FLAGS_OPT = -Ofast=IP35 -mips4 -r12000\n";
  print OUTFILE "FLAGS_MSC = -64\n";
  print OUTFILE "FLAGS_SER =\n";
  print OUTFILE "FLAGS_OMP = -mp\n";
  print OUTFILE "FLAGS_MPI =\n";
  print OUTFILE "INCS_SER = \n";
  print OUTFILE "INCS_OMP =\n";
  print OUTFILE "INCS_MPI =\n";
  print OUTFILE "LIBS_SER =\n";
  print OUTFILE "LIBS_OMP =\n";
  print OUTFILE "LIBS_MPI = -lmpi\n";
  print OUTFILE "OUT = -o \n";
  print OUTFILE "EXTO = o\n";
  print OUTFILE "MAKE = make\n";
  print OUTFILE "RM = rm -f\n";
  print OUTFILE "swch = -unix -f95 -timg -sgi\n";
}
elsif ($os =~ /AIX/i) {
  print OUTFILE "##############################################################################\n";
  print OUTFILE "# RS6000_IBM:		IBM SP with AIX using IBM compiler.\n";
  print OUTFILE "##############################################################################\n";
  print OUTFILE "F90_SER = xlf90_r\n";
  print OUTFILE "F90_OMP = xlf90_r\n";
  print OUTFILE "F90_MPI = mpxlf90_r\n";
  print OUTFILE "FLAGS_OPT = -O3 -qstrict -qarch=auto -qtune=auto -qnohot -qcache=auto \\\n";
  print OUTFILE "            -qunroll -qalign=4k -qfloat=hsflt\n";
  print OUTFILE "FLAGS_MSC = -w -qfixed -qnosave -q64\n";
  print OUTFILE "FLAGS_SER =\n";
  print OUTFILE "FLAGS_OMP = -qsmp=omp\n";
  print OUTFILE "FLAGS_MPI =\n";
  print OUTFILE "INCS_SER = \n";
  print OUTFILE "INCS_OMP =\n";
  print OUTFILE "INCS_MPI =\n";
  print OUTFILE "LIBS_SER = -lessl -lmass\n";
  print OUTFILE "LIBS_OMP = -lessl -lmass\n";
  print OUTFILE "LIBS_MPI = -lessl -lmass\n";
  print OUTFILE "OUT = -o \n";
  print OUTFILE "EXTO = o\n";
  print OUTFILE "MAKE = make\n";
  print OUTFILE "RM = rm -f\n";
  print OUTFILE "swch = -unix -timg\n";
}
elsif ($os =~ /OSF1/i) {
  print OUTFILE "##############################################################################\n";
  print OUTFILE "# Alpha_Compaq:		Compaq True 64 Alpha with OSF1 using Compaq compiler.\n";
  print OUTFILE "##############################################################################\n";
  print OUTFILE "F90_SER = f90\n";
  print OUTFILE "F90_OMP = f90\n";
  print OUTFILE "F90_MPI = f90\n";
  print OUTFILE "FLAGS_OPT = -fast\n";
  print OUTFILE "FLAGS_MSC = -w -fixed\n";
  print OUTFILE "FLAGS_SER =\n";
  print OUTFILE "FLAGS_OMP = -Wp,-C -omp\n";
  print OUTFILE "FLAGS_MPI =\n";
  print OUTFILE "INCS_SER =\n";
  print OUTFILE "INCS_OMP =\n";
  print OUTFILE "INCS_MPI =\n";
  print OUTFILE "LIBS_SER =\n";
  print OUTFILE "LIBS_OMP =\n";
  print OUTFILE "LIBS_MPI = -lfmpi -lmpi -lelan\n";
  print OUTFILE "OUT = -o \n";
  print OUTFILE "EXTO = o\n";
  print OUTFILE "MAKE = make\n";
  print OUTFILE "RM = rm -f\n";
  print OUTFILE "swch = -unix -timg\n";
}
elsif ($os =~ /SunOS/i) {
  print OUTFILE "##############################################################################\n";
  print OUTFILE "# SPARC_Sun:		Sun SPARC with Solaris using Sun compiler.\n";
  print OUTFILE "##############################################################################\n";
  print OUTFILE "F90_SER = f90\n";
  print OUTFILE "F90_OMP = f90\n";
  print OUTFILE "F90_MPI = mpf90\n";
  print OUTFILE "FLAGS_OPT = -xO3 -xtarget=native -fsimple=1 -depend -libmil -xlibmopt -xlic_lib=sunperf\n";
  print OUTFILE "FLAGS_MSC = -w -silent\n";
  print OUTFILE "FLAGS_OMP = -openmp\n";
  print OUTFILE "INCS_SER = \n";
  print OUTFILE "INCS_OMP = \n";
  print OUTFILE "INCS_MPI =\n";
  print OUTFILE "LIBS_SER =\n";
  print OUTFILE "LIBS_OMP =\n";
  print OUTFILE "LIBS_MPI = -lmpi\n";
  print OUTFILE "OUT = -o \n";
  print OUTFILE "EXTO = o\n";
  print OUTFILE "MAKE = make\n";
  print OUTFILE "RM = rm -f\n";
  print OUTFILE "swch = -unix -timg\n";
}
elsif ($os =~ /HP-UX/i) {
  print OUTFILE "##############################################################################\n";
  print OUTFILE "# RISC_HPUX:		PA-RISC with HP-UX 11 using HP Fortran compiler.\n";
  print OUTFILE "##############################################################################\n";
  print OUTFILE "F90_SER = f90\n";
  print OUTFILE "F90_OMP = f90\n";
  print OUTFILE "F90_MPI =\n";
  print OUTFILE "FLAGS_OPT = +O2 +Onolimit\n";
  print OUTFILE "FLAGS_MSC =\n";
  print OUTFILE "FLAGS_SER = +Onoopenmp\n";
  print OUTFILE "FLAGS_OMP = +Oopenmp\n";
  print OUTFILE "FLAGS_MPI =\n";
  print OUTFILE "INCS_SER =\n";
  print OUTFILE "INCS_OMP =\n";
  print OUTFILE "INCS_MPI =\n";
  print OUTFILE "LIBS_SER =\n";
  print OUTFILE "LIBS_OMP =\n";
  print OUTFILE "LIBS_MPI =\n";
  print OUTFILE "OUT = -o \n";
  print OUTFILE "EXTO = o\n";
  print OUTFILE "MAKE = make\n";
  print OUTFILE "RM = rm -f\n";
  print OUTFILE "swch = -unix -f95 -timg\n";
}
elsif ($os =~ /Linux/i) {
  system 'sh ./getcmpl';
  if ( -f "ifort" )
  {
    system 'rm ifort';
    print OUTFILE "##############################################################################\n";
    print OUTFILE "# IA32_Intel/x86-64_Intel:	Intel Pentium with Linux using Intel compiler 9.1.\n";
    print OUTFILE "##############################################################################\n";
    print OUTFILE "F90_SER = ifort\n";
    print OUTFILE "F90_OMP = ifort\n";
    print OUTFILE "F90_MPI = mpif90\n";
    if ($cpu =~ /i686/i) {
      print OUTFILE "FLAGS_OPT = -O2 -xN -mp1\n";
    }
    elsif ($cpu =~ /x86_64/i) {
      print OUTFILE "FLAGS_OPT = -O2 -ipo -xW -mp1\n";
    }
    else {
      print OUTFILE "FLAGS_OPT = -O2\n";
    }
    print OUTFILE "FLAGS_MSC = -W0 -assume byterecl -traceback\n";
    print OUTFILE "FLAGS_SER =\n";
    print OUTFILE "FLAGS_OMP = -openmp -assume cc_omp -fpp2\n";
    print OUTFILE "FLAGS_MPI =\n";
    print OUTFILE "INCS_SER =\n";
    print OUTFILE "INCS_OMP =\n";
    print OUTFILE "INCS_MPI =\n";
    print OUTFILE "LIBS_SER =\n";
    print OUTFILE "LIBS_OMP =\n";
    print OUTFILE "LIBS_MPI =\n";
    print OUTFILE "OUT = -o \n";
    print OUTFILE "EXTO = o\n";
    print OUTFILE "MAKE = make\n";
    print OUTFILE "RM = rm -f\n";
    print OUTFILE "swch = -unix -f95 -timg -impi\n";
  }
  elsif ( -f "ifc" )
  {
    system 'rm ifc';
    print OUTFILE "##############################################################################\n";
    print OUTFILE "# IA32_Intel:		Intel Pentium with Linux using Intel compiler 7.0.\n";
    print OUTFILE "# Note: -tpp6 is for Pentium III & -tpp7 is for Pentium IV\n";
    print OUTFILE "##############################################################################\n";
    print OUTFILE "F90_SER = ifc\n";
    print OUTFILE "F90_OMP = ifc\n";
    print OUTFILE "F90_MPI = mpif90\n";
    print OUTFILE "FLAGS_OPT = -O2 -tpp7\n";
    print OUTFILE "FLAGS_MSC = -W0 -auto\n";
    print OUTFILE "FLAGS_SER =\n";
    print OUTFILE "FLAGS_OMP = -openmp -fpp2\n";
    print OUTFILE "FLAGS_MPI =\n";
    print OUTFILE "INCS_SER =\n";
    print OUTFILE "INCS_OMP =\n";
    print OUTFILE "INCS_MPI =\n";
    print OUTFILE "LIBS_SER =\n";
    print OUTFILE "LIBS_OMP =\n";
    print OUTFILE "LIBS_MPI =\n";
    print OUTFILE "OUT = -o \n";
    print OUTFILE "EXTO = o\n";
    print OUTFILE "MAKE = make\n";
    print OUTFILE "RM = rm -f\n";
    print OUTFILE "swch = -unix -f95 -timg\n";
  }
  elsif ( -f "efc" )
  {
    system 'rm efc';
    print OUTFILE "##############################################################################\n";
    print OUTFILE "# IA64_Intel:		Intel Itanium with Linux using Intel compiler.\n";
    print OUTFILE "# Note: -tpp1 is for Itanium I & -tpp2 is for Itanium II\n";
    print OUTFILE "##############################################################################\n";
    print OUTFILE "F90_SER = efc\n";
    print OUTFILE "F90_OMP = efc\n";
    print OUTFILE "F90_MPI = mpif90\n";
    print OUTFILE "FLAGS_OPT = -O2 -tpp1\n";
    print OUTFILE "FLAGS_MSC = -W0 -auto\n";
    print OUTFILE "FLAGS_SER =\n";
    print OUTFILE "FLAGS_OMP = -openmp -fpp2\n";
    print OUTFILE "FLAGS_MPI =\n";
    print OUTFILE "INCS_SER =\n";
    print OUTFILE "INCS_OMP =\n";
    print OUTFILE "INCS_MPI =\n";
    print OUTFILE "LIBS_SER =\n";
    print OUTFILE "LIBS_OMP =\n";
    print OUTFILE "LIBS_MPI =\n";
    print OUTFILE "OUT = -o \n";
    print OUTFILE "EXTO = o\n";
    print OUTFILE "MAKE = make\n";
    print OUTFILE "RM = rm -f\n";
    print OUTFILE "swch = -unix -f95 -timg\n";
  }
  elsif ( -f "pgf90" )
  {
    system 'rm pgf90';
    print OUTFILE "##############################################################################\n";
    print OUTFILE "# IA32_PGF:		Intel Pentium with Linux using Portland Group compiler\n";
    print OUTFILE "# Note: -lpgc is required for OpenMP.\n";
    print OUTFILE "##############################################################################\n";
    print OUTFILE "F90_SER = pgf90\n";
    print OUTFILE "F90_OMP = pgf90\n";
    print OUTFILE "F90_MPI = mpif90\n";
    print OUTFILE "FLAGS_OPT =\n";
    print OUTFILE "FLAGS_MSC = -mcmodel=medium\n";
    print OUTFILE "FLAGS_SER = -fast\n";
    print OUTFILE "FLAGS_OMP = -mp\n";
    print OUTFILE "FLAGS_MPI = -fast\n";
    print OUTFILE "INCS_SER =\n";
    print OUTFILE "INCS_OMP =\n";
    print OUTFILE "INCS_MPI =\n";
    print OUTFILE "LIBS_SER =\n";
    print OUTFILE "LIBS_OMP = -lpgc\n";
    print OUTFILE "LIBS_MPI = -lmpi\n";
    print OUTFILE "OUT = -o \n";
    print OUTFILE "EXTO = o\n";
    print OUTFILE "MAKE = make\n";
    print OUTFILE "RM = rm -f\n";
    print OUTFILE "swch = -unix -impi -timg\n";
  }
  elsif ( -f "lf95" )
  {
    system 'rm lf95';
    print OUTFILE "##############################################################################\n";
    print OUTFILE "# IA32_Lahey:		Intel Pentium with Linux using Lahey compiler.\n";
    print OUTFILE "##############################################################################\n";
    print OUTFILE "F90_SER = lf95\n";
    print OUTFILE "F90_OMP =\n";
    print OUTFILE "F90_MPI = mpif90\n";
    print OUTFILE "FLAGS_OPT = -O --tpp\n";
    print OUTFILE "FLAGS_MSC = --staticlink --nwo\n";
    print OUTFILE "FLAGS_SER =\n";
    print OUTFILE "FLAGS_OMP =\n";
    print OUTFILE "FLAGS_MPI =\n";
    print OUTFILE "INCS_SER =\n";
    print OUTFILE "INCS_OMP =\n";
    print OUTFILE "INCS_MPI =\n";
    print OUTFILE "LIBS_SER =\n";
    print OUTFILE "LIBS_OMP =\n";
    print OUTFILE "LIBS_MPI =\n";
    print OUTFILE "OUT = -o \n";
    print OUTFILE "EXTO = o\n";
    print OUTFILE "MAKE = make\n";
    print OUTFILE "RM = rm -f\n";
    print OUTFILE "swch = -unix -f95 -timg\n";
  }
  elsif ( -f "gfortran" )
  {
    system 'rm gfortran';
    print OUTFILE "##############################################################################\n";
    print OUTFILE "# IA32_GNU:		Intel Pentium with Linux using GNU compiler gfortran.\n";
    print OUTFILE "##############################################################################\n";
    print OUTFILE "F90_SER = gfortran\n";
    print OUTFILE "F90_OMP = \n";
    print OUTFILE "F90_MPI = mpif90\n";
    print OUTFILE "FLAGS_OPT = -O\n";
    print OUTFILE "FLAGS_MSC = -w\n";
    print OUTFILE "FLAGS_SER =\n";
    print OUTFILE "FLAGS_OMP =\n";
    print OUTFILE "FLAGS_MPI =\n";
    print OUTFILE "INCS_SER =\n";
    print OUTFILE "INCS_OMP =\n";
    print OUTFILE "INCS_MPI =\n";
    print OUTFILE "LIBS_SER =\n";
    print OUTFILE "LIBS_OMP =\n";
    print OUTFILE "LIBS_MPI =\n";
    print OUTFILE "OUT = -o \n";
    print OUTFILE "EXTO = o\n";
    print OUTFILE "MAKE = make\n";
    print OUTFILE "RM = rm -f\n";
    print OUTFILE "swch = -unix -f95 -timg\n";
  }
  elsif ( -f "g95" )
  {
    system 'rm g95';
    print OUTFILE "##############################################################################\n";
    print OUTFILE "# IA32_GNU:		Intel Pentium with Linux using GNU compiler g95.\n";
    print OUTFILE "##############################################################################\n";
    print OUTFILE "F90_SER = g95\n";
    print OUTFILE "F90_OMP = \n";
    print OUTFILE "F90_MPI = mpif90\n";
    print OUTFILE "FLAGS_OPT = -O\n";
    print OUTFILE "FLAGS_MSC = -w\n";
    print OUTFILE "FLAGS_SER =\n";
    print OUTFILE "FLAGS_OMP =\n";
    print OUTFILE "FLAGS_MPI =\n";
    print OUTFILE "INCS_SER =\n";
    print OUTFILE "INCS_OMP =\n";
    print OUTFILE "INCS_MPI =\n";
    print OUTFILE "LIBS_SER =\n";
    print OUTFILE "LIBS_OMP =\n";
    print OUTFILE "LIBS_MPI =\n";
    print OUTFILE "OUT = -o \n";
    print OUTFILE "EXTO = o\n";
    print OUTFILE "MAKE = make\n";
    print OUTFILE "RM = rm -f\n";
    print OUTFILE "swch = -unix -f95 -timg\n";
  }
  else
  {
    die "Current Fortran compiler not supported.... \n";
  }
}
elsif ($os =~ /WindowsNT/i || $os =~ /MSWin32/i) {
  print OUTFILE "##############################################################################\n";
  print OUTFILE "# IA32_Intel/EM64T_Intel:	Intel Pentium with MS Windows using Intel compiler 9.1.\n";
  print OUTFILE "##############################################################################\n";
  print OUTFILE "F90_SER = ifort\n";
  print OUTFILE "F90_OMP = ifort\n";
  print OUTFILE "F90_MPI = ifort\n";
  print OUTFILE "FLAGS_OPT = /optimize:2\n";
  print OUTFILE "FLAGS_MSC = /assume:byterecl /traceback /nowarn /nologo\n";
  print OUTFILE "FLAGS_SER =\n";
  print OUTFILE "FLAGS_OMP = /Qopenmp /assume:cc_omp /fpp\n";
  print OUTFILE "FLAGS_MPI =\n";
  print OUTFILE "INCS_SER =\n";
  print OUTFILE "INCS_OMP =\n";
  print OUTFILE "INCS_MPI = /include:\"c:\\\progra~1\\\MPICH\\\SDK\\\include\"\n";
  print OUTFILE "LIBS_SER =\n";
  print OUTFILE "LIBS_OMP =\n";
  print OUTFILE "LIBS_MPI = c:\\\progra~1\\\MPICH\\\SDK\\\llib\\\mpe.lib c:\\\progra~1\\\MPICH\\\SDK\\\llib\\\mpich.lib\n";
  print OUTFILE "OUT = /exe:\n";
  print OUTFILE "EXTO = obj\n";
  print OUTFILE "MAKE = nmake\n";
  print OUTFILE "RM = del\n";
  print OUTFILE "swch = -dos -f95 -impi -cvis -timg\n";
}
elsif ($os =~ /Darwin/i) {
  print OUTFILE "##############################################################################\n";
  print OUTFILE "# MAC_IBM:		MAC OS X Apple with IBM Fortran Compiler.\n";
  print OUTFILE "##############################################################################\n";
  print OUTFILE "F90_SER = xlf90\n";
  print OUTFILE "F90_OMP =\n";
  print OUTFILE "F90_MPI =\n";
  print OUTFILE "FLAGS_OPT = -O3 -qstrict -qtune=auto -qcache=auto -qalign=4k\n";
  print OUTFILE "FLAGS_MSC = -w -qfixed\n";
  print OUTFILE "FLAGS_SER =\n";
  print OUTFILE "FLAGS_OMP =\n";
  print OUTFILE "FLAGS_MPI =\n";
  print OUTFILE "INCS_SER = \n";
  print OUTFILE "INCS_OMP =\n";
  print OUTFILE "INCS_MPI =\n";
  print OUTFILE "LIBS_SER =\n";
  print OUTFILE "LIBS_OMP =\n";
  print OUTFILE "LIBS_MPI =\n";
  print OUTFILE "OUT = -o \n";
  print OUTFILE "EXTO = o\n";
  print OUTFILE "MAKE = make\n";
  print OUTFILE "RM = rm -f\n";
  print OUTFILE "swch = -unix -timg\n";
}
else
{
   die "Current operating system not supported.... \n";
}

close(OUTFILE);
