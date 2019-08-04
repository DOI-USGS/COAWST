my $os  = `uname -s`;
   $os  = $^O unless chomp($os);

my $cpu = `uname -m`;
   $cpu = $^O unless chomp($cpu);

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
  print OUTFILE "FLAGS90_MSC = \$(FLAGS_MSC)\n";
  print OUTFILE "FLAGS_DYN =\n";
  print OUTFILE "FLAGS_SER =\n";
  print OUTFILE "FLAGS_OMP = -mp\n";
  print OUTFILE "FLAGS_MPI =\n";
  print OUTFILE "NETCDFROOT =\n";
  print OUTFILE "ifneq (\$(NETCDFROOT),)\n";
  print OUTFILE "  INCS_SER = -I\$(NETCDFROOT)/include\n";
  print OUTFILE "  INCS_OMP = -I\$(NETCDFROOT)/include\n";
  print OUTFILE "  INCS_MPI = -I\$(NETCDFROOT)/include\n";
  print OUTFILE "  LIBS_SER = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
  print OUTFILE "  LIBS_OMP = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
  print OUTFILE "  LIBS_MPI = -lmpi -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
  print OUTFILE "  NCF_OBJS = nctablemd.o agioncmd.o swn_outnc.o\n";
  print OUTFILE "else\n";
  print OUTFILE "  INCS_SER =\n";
  print OUTFILE "  INCS_OMP =\n";
  print OUTFILE "  INCS_MPI =\n";
  print OUTFILE "  LIBS_SER =\n";
  print OUTFILE "  LIBS_OMP =\n";
  print OUTFILE "  LIBS_MPI = -lmpi\n";
  print OUTFILE "  NCF_OBJS =\n";
  print OUTFILE "endif\n";
  print OUTFILE "O_DIR = ../work/odir4/\n";
  print OUTFILE "OUT = -o \n";
  print OUTFILE "EXTO = o\n";
  print OUTFILE "MAKE = make\n";
  print OUTFILE "RM = rm -f\n";
  print OUTFILE "ifneq (\$(NETCDFROOT),)\n";
  print OUTFILE "  swch = -unix -f95 -timg -sgi -netcdf\n";
  print OUTFILE "else\n";
  print OUTFILE "  swch = -unix -f95 -timg -sgi\n";
  print OUTFILE "endif\n";
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
  print OUTFILE "FLAGS90_MSC = -w -qfree=f90 -qnosave -q64\n";
  print OUTFILE "FLAGS_DYN =\n";
  print OUTFILE "FLAGS_SER =\n";
  print OUTFILE "FLAGS_OMP = -qsmp=omp\n";
  print OUTFILE "FLAGS_MPI =\n";
  print OUTFILE "NETCDFROOT =\n";
  print OUTFILE "ifneq (\$(NETCDFROOT),)\n";
  print OUTFILE "  INCS_SER = -I\$(NETCDFROOT)/include\n";
  print OUTFILE "  INCS_OMP = -I\$(NETCDFROOT)/include\n";
  print OUTFILE "  INCS_MPI = -I\$(NETCDFROOT)/include\n";
  print OUTFILE "  LIBS_SER = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff -lessl -lmass\n";
  print OUTFILE "  LIBS_OMP = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff -lessl -lmass\n";
  print OUTFILE "  LIBS_MPI = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff -lessl -lmass\n";
  print OUTFILE "  NCF_OBJS = nctablemd.o agioncmd.o swn_outnc.o\n";
  print OUTFILE "else\n";
  print OUTFILE "  INCS_SER =\n";
  print OUTFILE "  INCS_OMP =\n";
  print OUTFILE "  INCS_MPI =\n";
  print OUTFILE "  LIBS_SER = -lessl -lmass\n";
  print OUTFILE "  LIBS_OMP = -lessl -lmass\n";
  print OUTFILE "  LIBS_MPI = -lessl -lmass\n";
  print OUTFILE "  NCF_OBJS =\n";
  print OUTFILE "endif\n";
  print OUTFILE "O_DIR = ../work/odir4/\n";
  print OUTFILE "OUT = -o \n";
  print OUTFILE "EXTO = o\n";
  print OUTFILE "MAKE = make\n";
  print OUTFILE "RM = rm -f\n";
  print OUTFILE "ifneq (\$(NETCDFROOT),)\n";
  print OUTFILE "  swch = -unix -f95 -timg -netcdf\n";
  print OUTFILE "else\n";
  print OUTFILE "  swch = -unix -f95 -timg\n";
  print OUTFILE "endif\n";
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
  print OUTFILE "FLAGS90_MSC = -w -free\n";
  print OUTFILE "FLAGS_DYN =\n";
  print OUTFILE "FLAGS_SER =\n";
  print OUTFILE "FLAGS_OMP = -Wp,-C -omp\n";
  print OUTFILE "FLAGS_MPI =\n";
  print OUTFILE "NETCDFROOT =\n";
  print OUTFILE "INCS_SER =\n";
  print OUTFILE "INCS_OMP =\n";
  print OUTFILE "INCS_MPI =\n";
  print OUTFILE "LIBS_SER =\n";
  print OUTFILE "LIBS_OMP =\n";
  print OUTFILE "LIBS_MPI = -lfmpi -lmpi -lelan\n";
  print OUTFILE "O_DIR = ../work/odir4/\n";
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
  print OUTFILE "FLAGS90_MSC = \$(FLAGS_MSC)\n";
  print OUTFILE "FLAGS_DYN =\n";
  print OUTFILE "FLAGS_SER =\n";
  print OUTFILE "FLAGS_OMP = -openmp\n";
  print OUTFILE "FLAGS_MPI =\n";
  print OUTFILE "NETCDFROOT =\n";
  print OUTFILE "ifneq (\$(NETCDFROOT),)\n";
  print OUTFILE "  INCS_SER = -I\$(NETCDFROOT)/include\n";
  print OUTFILE "  INCS_OMP = -I\$(NETCDFROOT)/include\n";
  print OUTFILE "  INCS_MPI = -I\$(NETCDFROOT)/include\n";
  print OUTFILE "  LIBS_SER = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
  print OUTFILE "  LIBS_OMP = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
  print OUTFILE "  LIBS_MPI = -lmpi -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
  print OUTFILE "  NCF_OBJS = nctablemd.o agioncmd.o swn_outnc.o\n";
  print OUTFILE "else\n";
  print OUTFILE "  INCS_SER =\n";
  print OUTFILE "  INCS_OMP =\n";
  print OUTFILE "  INCS_MPI =\n";
  print OUTFILE "  LIBS_SER =\n";
  print OUTFILE "  LIBS_OMP =\n";
  print OUTFILE "  LIBS_MPI = -lmpi\n";
  print OUTFILE "  NCF_OBJS =\n";
  print OUTFILE "endif\n";
  print OUTFILE "O_DIR = ../work/odir4/\n";
  print OUTFILE "OUT = -o \n";
  print OUTFILE "EXTO = o\n";
  print OUTFILE "MAKE = make\n";
  print OUTFILE "RM = rm -f\n";
  print OUTFILE "ifneq (\$(NETCDFROOT),)\n";
  print OUTFILE "  swch = -unix -timg -netcdf\n";
  print OUTFILE "else\n";
  print OUTFILE "  swch = -unix -timg\n";
  print OUTFILE "endif\n";
}
elsif ($os =~ /HP-UX/i) {
  print OUTFILE "##############################################################################\n";
  print OUTFILE "# RISC_HPUX:		PA-RISC with HP-UX 11 using HP compiler.\n";
  print OUTFILE "##############################################################################\n";
  print OUTFILE "F90_SER = f90\n";
  print OUTFILE "F90_OMP = f90\n";
  print OUTFILE "F90_MPI =\n";
  print OUTFILE "FLAGS_OPT = +O2 +Onolimit\n";
  print OUTFILE "FLAGS_MSC =\n";
  print OUTFILE "FLAGS90_MSC =\n";
  print OUTFILE "FLAGS_DYN =\n";
  print OUTFILE "FLAGS_SER = +Onoopenmp\n";
  print OUTFILE "FLAGS_OMP = +Oopenmp\n";
  print OUTFILE "FLAGS_MPI =\n";
  print OUTFILE "NETCDFROOT =\n";
  print OUTFILE "ifneq (\$(NETCDFROOT),)\n";
  print OUTFILE "  INCS_SER = -I\$(NETCDFROOT)/include\n";
  print OUTFILE "  INCS_OMP = -I\$(NETCDFROOT)/include\n";
  print OUTFILE "  INCS_MPI = -I\$(NETCDFROOT)/include\n";
  print OUTFILE "  LIBS_SER = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
  print OUTFILE "  LIBS_OMP = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
  print OUTFILE "  LIBS_MPI = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
  print OUTFILE "  NCF_OBJS = nctablemd.o agioncmd.o swn_outnc.o\n";
  print OUTFILE "else\n";
  print OUTFILE "  INCS_SER =\n";
  print OUTFILE "  INCS_OMP =\n";
  print OUTFILE "  INCS_MPI =\n";
  print OUTFILE "  LIBS_SER =\n";
  print OUTFILE "  LIBS_OMP =\n";
  print OUTFILE "  LIBS_MPI =\n";
  print OUTFILE "  NCF_OBJS =\n";
  print OUTFILE "endif\n";
  print OUTFILE "O_DIR =\n";
  print OUTFILE "OUT = -o \n";
  print OUTFILE "EXTO = o\n";
  print OUTFILE "MAKE = make\n";
  print OUTFILE "RM = rm -f\n";
  print OUTFILE "ifneq (\$(NETCDFROOT),)\n";
  print OUTFILE "  swch = -unix -f95 -timg -netcdf\n";
  print OUTFILE "else\n";
  print OUTFILE "  swch = -unix -f95 -timg\n";
  print OUTFILE "endif\n";
}
elsif ($os =~ /Linux/i) {
  my $compiler = getcmpl();
  if ( $compiler eq "ifort" )
  {
    print OUTFILE "##############################################################################\n";
    print OUTFILE "# IA32_Intel/x86-64_Intel:	Intel Pentium with Linux using Intel compiler 17.\n";
    print OUTFILE "##############################################################################\n";
    print OUTFILE "F90_SER = ifort\n";
    print OUTFILE "F90_OMP = ifort\n";
    print OUTFILE "# if appropriate, use mpiifort of Intel instead\n";
    print OUTFILE "F90_MPI = mpif90\n";
#    if ($cpu =~ /i686/i) {
#      print OUTFILE "FLAGS_OPT = -O2 -xN -mp1\n";
#    }
#    elsif ($cpu =~ /x86_64/i) {
#      print OUTFILE "FLAGS_OPT = -O2 -ipo -xW -mp1\n";
#    }
#    else {
      print OUTFILE "FLAGS_OPT = -O2\n";
#    }
#    print OUTFILE "FLAGS_MSC = -W0 -assume byterecl -traceback -diag-disable remark\n";
    print OUTFILE "FLAGS_MSC = -W0 -assume byterecl -traceback -diag-disable 8290 -diag-disable 8291 -diag-disable 8293\n";
    print OUTFILE "FLAGS90_MSC = \$(FLAGS_MSC)\n";
    print OUTFILE "FLAGS_DYN = -fPIC\n";
    print OUTFILE "FLAGS_SER =\n";
    print OUTFILE "FLAGS_OMP = -qopenmp\n";
    print OUTFILE "FLAGS_MPI =\n";
    print OUTFILE "NETCDFROOT =\n";
    print OUTFILE "ifneq (\$(NETCDFROOT),)\n";
    print OUTFILE "  INCS_SER = -I\$(NETCDFROOT)/include\n";
    print OUTFILE "  INCS_OMP = -I\$(NETCDFROOT)/include\n";
    print OUTFILE "  INCS_MPI = -I\$(NETCDFROOT)/include\n";
    print OUTFILE "  LIBS_SER = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
    print OUTFILE "  LIBS_OMP = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
    print OUTFILE "  LIBS_MPI = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
    print OUTFILE "  NCF_OBJS = nctablemd.o agioncmd.o swn_outnc.o\n";
    print OUTFILE "else\n";
    print OUTFILE "  INCS_SER =\n";
    print OUTFILE "  INCS_OMP =\n";
    print OUTFILE "  INCS_MPI =\n";
    print OUTFILE "  LIBS_SER =\n";
    print OUTFILE "  LIBS_OMP =\n";
    print OUTFILE "  LIBS_MPI =\n";
    print OUTFILE "  NCF_OBJS =\n";
    print OUTFILE "endif\n";
    print OUTFILE "O_DIR = ../work/odir4/\n";
    print OUTFILE "OUT = -o \n";
    print OUTFILE "EXTO = o\n";
    print OUTFILE "MAKE = make\n";
    print OUTFILE "RM = rm -f\n";
    print OUTFILE "ifneq (\$(NETCDFROOT),)\n";
    print OUTFILE "  swch = -unix -impi -netcdf\n";
    print OUTFILE "else\n";
    print OUTFILE "  swch = -unix -impi\n";
    print OUTFILE "endif\n";
  }
  elsif ( $compiler eq "ifc" )
  {
    print OUTFILE "##############################################################################\n";
    print OUTFILE "# IA32_Intel:		Intel Pentium with Linux using Intel compiler 7.0.\n";
    print OUTFILE "# Note: -tpp6 is for Pentium III & -tpp7 is for Pentium IV\n";
    print OUTFILE "##############################################################################\n";
    print OUTFILE "F90_SER = ifc\n";
    print OUTFILE "F90_OMP = ifc\n";
    print OUTFILE "F90_MPI = mpif90\n";
    print OUTFILE "FLAGS_OPT = -O2 -tpp7\n";
    print OUTFILE "FLAGS_MSC = -W0 -auto\n";
    print OUTFILE "FLAGS90_MSC = \$(FLAGS_MSC)\n";
    print OUTFILE "FLAGS_DYN =\n";
    print OUTFILE "FLAGS_SER =\n";
    print OUTFILE "FLAGS_OMP = -openmp\n";
    print OUTFILE "FLAGS_MPI =\n";
    print OUTFILE "NETCDFROOT =\n";
    print OUTFILE "INCS_SER =\n";
    print OUTFILE "INCS_OMP =\n";
    print OUTFILE "INCS_MPI =\n";
    print OUTFILE "LIBS_SER =\n";
    print OUTFILE "LIBS_OMP =\n";
    print OUTFILE "LIBS_MPI =\n";
    print OUTFILE "O_DIR = ../work/odir4/\n";
    print OUTFILE "OUT = -o \n";
    print OUTFILE "EXTO = o\n";
    print OUTFILE "MAKE = make\n";
    print OUTFILE "RM = rm -f\n";
    print OUTFILE "swch = -unix -f95 -timg\n";
  }
  elsif ( $compiler eq "efc" )
  {
    print OUTFILE "##############################################################################\n";
    print OUTFILE "# IA64_Intel:		Intel Itanium with Linux using Intel compiler.\n";
    print OUTFILE "# Note: -tpp1 is for Itanium I & -tpp2 is for Itanium II\n";
    print OUTFILE "##############################################################################\n";
    print OUTFILE "F90_SER = efc\n";
    print OUTFILE "F90_OMP = efc\n";
    print OUTFILE "F90_MPI = mpif90\n";
    print OUTFILE "FLAGS_OPT = -O2 -tpp1\n";
    print OUTFILE "FLAGS_MSC = -W0 -auto\n";
    print OUTFILE "FLAGS90_MSC = \$(FLAGS_MSC)\n";
    print OUTFILE "FLAGS_DYN =\n";
    print OUTFILE "FLAGS_SER =\n";
    print OUTFILE "FLAGS_OMP = -openmp\n";
    print OUTFILE "FLAGS_MPI =\n";
    print OUTFILE "NETCDFROOT =\n";
    print OUTFILE "INCS_SER =\n";
    print OUTFILE "INCS_OMP =\n";
    print OUTFILE "INCS_MPI =\n";
    print OUTFILE "LIBS_SER =\n";
    print OUTFILE "LIBS_OMP =\n";
    print OUTFILE "LIBS_MPI =\n";
    print OUTFILE "O_DIR = ../work/odir4/\n";
    print OUTFILE "OUT = -o \n";
    print OUTFILE "EXTO = o\n";
    print OUTFILE "MAKE = make\n";
    print OUTFILE "RM = rm -f\n";
    print OUTFILE "swch = -unix -f95 -timg\n";
  }
  elsif ( $compiler eq "pgf90" )
  {
    print OUTFILE "##############################################################################\n";
    print OUTFILE "# IA32_PGF:		Intel Pentium with Linux using Portland Group compiler\n";
    print OUTFILE "##############################################################################\n";
    print OUTFILE "F90_SER = pgf90\n";
    print OUTFILE "F90_OMP = pgf90\n";
    print OUTFILE "F90_MPI = mpif90\n";
    print OUTFILE "FLAGS_OPT = -fast\n";
    print OUTFILE "FLAGS_MSC = -Mfixed\n";
    print OUTFILE "FLAGS90_MSC = -Mfree\n";
    print OUTFILE "FLAGS_DYN =\n";
    print OUTFILE "FLAGS_SER =\n";
    print OUTFILE "FLAGS_OMP = -mp\n";
    print OUTFILE "FLAGS_MPI = -tp barcelona-64\n";
    print OUTFILE "NETCDFROOT =\n";
    print OUTFILE "ifneq (\$(NETCDFROOT),)\n";
    print OUTFILE "  INCS_SER = -I\$(NETCDFROOT)/include\n";
    print OUTFILE "  INCS_OMP = -I\$(NETCDFROOT)/include\n";
    print OUTFILE "  INCS_MPI = -I\$(NETCDFROOT)/include\n";
    print OUTFILE "  LIBS_SER = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
    print OUTFILE "  LIBS_OMP = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
    print OUTFILE "  LIBS_MPI = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
    print OUTFILE "  NCF_OBJS = nctablemd.o agioncmd.o swn_outnc.o\n";
    print OUTFILE "else\n";
    print OUTFILE "  INCS_SER =\n";
    print OUTFILE "  INCS_OMP =\n";
    print OUTFILE "  INCS_MPI =\n";
    print OUTFILE "  LIBS_SER =\n";
    print OUTFILE "  LIBS_OMP =\n";
    print OUTFILE "  LIBS_MPI =\n";
    print OUTFILE "  NCF_OBJS =\n";
    print OUTFILE "endif\n";
    print OUTFILE "O_DIR = ../work/odir4/\n";
    print OUTFILE "OUT = -o \n";
    print OUTFILE "EXTO = o\n";
    print OUTFILE "MAKE = make\n";
    print OUTFILE "RM = rm -f\n";
    print OUTFILE "ifneq (\$(NETCDFROOT),)\n";
    print OUTFILE "  swch = -unix -impi -netcdf\n";
    print OUTFILE "else\n";
    print OUTFILE "  swch = -unix -impi\n";
    print OUTFILE "endif\n";
  }
  elsif ( $compiler eq "lf95" )
  {
    print OUTFILE "##############################################################################\n";
    print OUTFILE "# IA32_Lahey:		Intel Pentium with Linux using Lahey compiler.\n";
    print OUTFILE "##############################################################################\n";
    print OUTFILE "F90_SER = lf95\n";
    print OUTFILE "F90_OMP =\n";
    print OUTFILE "F90_MPI = mpif90\n";
    print OUTFILE "FLAGS_OPT = -O --tpp\n";
    print OUTFILE "FLAGS_MSC = --staticlink --nwo\n";
    print OUTFILE "FLAGS90_MSC = \$(FLAGS_MSC)\n";
    print OUTFILE "FLAGS_DYN =\n";
    print OUTFILE "FLAGS_SER =\n";
    print OUTFILE "FLAGS_OMP =\n";
    print OUTFILE "FLAGS_MPI =\n";
    print OUTFILE "NETCDFROOT =\n";
    print OUTFILE "INCS_SER =\n";
    print OUTFILE "INCS_OMP =\n";
    print OUTFILE "INCS_MPI =\n";
    print OUTFILE "LIBS_SER =\n";
    print OUTFILE "LIBS_OMP =\n";
    print OUTFILE "LIBS_MPI =\n";
    print OUTFILE "O_DIR =\n";
    print OUTFILE "OUT = -o \n";
    print OUTFILE "EXTO = o\n";
    print OUTFILE "MAKE = make\n";
    print OUTFILE "RM = rm -f\n";
    print OUTFILE "swch = -unix\n";
  }
  elsif ( $compiler eq "gfortran" )
  {
    print OUTFILE "##############################################################################\n";
    print OUTFILE "# IA32_GNU:		Intel Pentium with Linux using GNU compiler gfortran.\n";
    print OUTFILE "##############################################################################\n";
    print OUTFILE "F90_SER = gfortran\n";
    print OUTFILE "F90_OMP = gfortran\n";
    print OUTFILE "F90_MPI = mpif90\n";
    print OUTFILE "FLAGS_OPT = -O\n";
    print OUTFILE "FLAGS_MSC = -w -fno-second-underscore\n";
    print OUTFILE "FLAGS90_MSC = \$(FLAGS_MSC) -ffree-line-length-none\n";
    print OUTFILE "FLAGS_DYN =\n";
    print OUTFILE "FLAGS_SER =\n";
    print OUTFILE "FLAGS_OMP = -fopenmp\n";
    print OUTFILE "FLAGS_MPI =\n";
    print OUTFILE "NETCDFROOT =\n";
    print OUTFILE "ifneq (\$(NETCDFROOT),)\n";
    print OUTFILE "  INCS_SER = -I\$(NETCDFROOT)/include\n";
    print OUTFILE "  INCS_OMP = -I\$(NETCDFROOT)/include\n";
    print OUTFILE "  INCS_MPI = -I\$(NETCDFROOT)/include\n";
    print OUTFILE "  LIBS_SER = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
    print OUTFILE "  LIBS_OMP = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff -static-libgcc\n";
    print OUTFILE "  LIBS_MPI = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
    print OUTFILE "  NCF_OBJS = nctablemd.o agioncmd.o swn_outnc.o\n";
    print OUTFILE "else\n";
    print OUTFILE "  INCS_SER =\n";
    print OUTFILE "  INCS_OMP =\n";
    print OUTFILE "  INCS_MPI =\n";
    print OUTFILE "  LIBS_SER =\n";
    print OUTFILE "  LIBS_OMP = -static-libgcc\n";
    print OUTFILE "  LIBS_MPI =\n";
    print OUTFILE "  NCF_OBJS =\n";
    print OUTFILE "endif\n";
    print OUTFILE "O_DIR = ../work/odir4/\n";
    print OUTFILE "OUT = -o \n";
    print OUTFILE "EXTO = o\n";
    print OUTFILE "MAKE = make\n";
    print OUTFILE "RM = rm -f\n";
    print OUTFILE "ifneq (\$(NETCDFROOT),)\n";
    print OUTFILE "  swch = -unix -netcdf\n";
    print OUTFILE "else\n";
    print OUTFILE "  swch = -unix\n";
    print OUTFILE "endif\n";
  }
  elsif ( $compiler eq "g95" )
  {
    print OUTFILE "##############################################################################\n";
    print OUTFILE "# IA32_GNU:		Intel Pentium with Linux using GNU compiler g95.\n";
    print OUTFILE "##############################################################################\n";
    print OUTFILE "F90_SER = g95\n";
    print OUTFILE "F90_OMP = NO_OPENMP_WITH_G95\n";
    print OUTFILE "F90_MPI = mpif90\n";
    print OUTFILE "FLAGS_OPT = -O\n";
    print OUTFILE "FLAGS_MSC = -fno-second-underscore\n";
    print OUTFILE "FLAGS90_MSC = \$(FLAGS_MSC) -ffree-line-length-huge\n";
    print OUTFILE "FLAGS_DYN =\n";
    print OUTFILE "FLAGS_SER =\n";
    print OUTFILE "FLAGS_OMP =\n";
    print OUTFILE "FLAGS_MPI =\n";
    print OUTFILE "NETCDFROOT =\n";
    print OUTFILE "ifneq (\$(NETCDFROOT),)\n";
    print OUTFILE "  INCS_SER = -I\$(NETCDFROOT)/include\n";
    print OUTFILE "  INCS_OMP = -I\$(NETCDFROOT)/include\n";
    print OUTFILE "  INCS_MPI = -I\$(NETCDFROOT)/include\n";
    print OUTFILE "  LIBS_SER = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
    print OUTFILE "  LIBS_OMP = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
    print OUTFILE "  LIBS_MPI = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
    print OUTFILE "  NCF_OBJS = nctablemd.o agioncmd.o swn_outnc.o\n";
    print OUTFILE "else\n";
    print OUTFILE "  INCS_SER =\n";
    print OUTFILE "  INCS_OMP =\n";
    print OUTFILE "  INCS_MPI =\n";
    print OUTFILE "  LIBS_SER =\n";
    print OUTFILE "  LIBS_OMP =\n";
    print OUTFILE "  LIBS_MPI =\n";
    print OUTFILE "  NCF_OBJS =\n";
    print OUTFILE "endif\n";
    print OUTFILE "O_DIR = ../work/odir4/\n";
    print OUTFILE "OUT = -o \n";
    print OUTFILE "EXTO = o\n";
    print OUTFILE "MAKE = make\n";
    print OUTFILE "RM = rm -f\n";
    print OUTFILE "ifneq (\$(NETCDFROOT),)\n";
    print OUTFILE "  swch = -unix -netcdf\n";
    print OUTFILE "else\n";
    print OUTFILE "  swch = -unix\n";
    print OUTFILE "endif\n";
  }
  elsif ( $compiler eq "xlf90" )
  {
    print OUTFILE "##############################################################################\n";
    print OUTFILE "# PPC64_IBM:		IBM Power6 with Linux using IBM Compiler.\n";
    print OUTFILE "##############################################################################\n";
    print OUTFILE "F90_SER = xlf90\n";
    print OUTFILE "F90_OMP =\n";
    print OUTFILE "F90_MPI = mpfort\n";
    print OUTFILE "FLAGS_OPT = -O3 -qstrict -qarch=auto -qtune=auto\n";
    print OUTFILE "FLAGS_MSC = -qfixed -qzerosize -qwarn64\n";
    print OUTFILE "FLAGS90_MSC = -qfree=f90 -qzerosize -qwarn64\n";
    print OUTFILE "FLAGS_DYN =\n";
    print OUTFILE "FLAGS_SER =\n";
    print OUTFILE "FLAGS_OMP =\n";
    print OUTFILE "FLAGS_MPI =\n";
    print OUTFILE "NETCDFROOT =\n";
    print OUTFILE "ifneq (\$(NETCDFROOT),)\n";
    print OUTFILE "  INCS_SER = -I\$(NETCDFROOT)/include\n";
    print OUTFILE "  INCS_OMP = -I\$(NETCDFROOT)/include\n";
    print OUTFILE "  INCS_MPI = -I\$(NETCDFROOT)/include\n";
    print OUTFILE "  LIBS_SER = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
    print OUTFILE "  LIBS_OMP = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
    print OUTFILE "  LIBS_MPI = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
    print OUTFILE "  NCF_OBJS = nctablemd.o agioncmd.o swn_outnc.o\n";
    print OUTFILE "else\n";
    print OUTFILE "  INCS_SER =\n";
    print OUTFILE "  INCS_OMP =\n";
    print OUTFILE "  INCS_MPI =\n";
    print OUTFILE "  LIBS_SER =\n";
    print OUTFILE "  LIBS_OMP =\n";
    print OUTFILE "  LIBS_MPI =\n";
    print OUTFILE "  NCF_OBJS =\n";
    print OUTFILE "endif\n";
    print OUTFILE "O_DIR = ../work/odir4/\n";
    print OUTFILE "OUT = -o \n";
    print OUTFILE "EXTO = o\n";
    print OUTFILE "MAKE = make\n";
    print OUTFILE "RM = rm -f\n";
    print OUTFILE "ifneq (\$(NETCDFROOT),)\n";
    print OUTFILE "  swch = -unix -matl4 -netcdf\n";
    print OUTFILE "else\n";
    print OUTFILE "  swch = -unix -matl4\n";
    print OUTFILE "endif\n";
  }
  else
  {
    die "Current Fortran compiler '$compiler' not supported.... \n";
  }
}
elsif ($os =~ /WindowsNT/i || $os =~ /MSWin32/i || $os =~ /CYGWIN/i) {
  my $compiler = getcmpl();
  if ( $compiler eq "ifort" )
  {
    print OUTFILE "##############################################################################\n";
    print OUTFILE "# IA32_Intel/EM64T_Intel:	Intel Pentium with MS Windows using Intel compiler 17.\n";
    print OUTFILE "##############################################################################\n";
    print OUTFILE "F90_SER = ifort\n";
    print OUTFILE "F90_OMP = ifort\n";
    print OUTFILE "F90_MPI = ifort\n";
    print OUTFILE "FLAGS_OPT = /O2\n";
#    print OUTFILE "FLAGS_MSC = /assume:byterecl /traceback /nowarn /nologo /Qdiag-disable:remark\n";
    print OUTFILE "FLAGS_MSC = /assume:byterecl /traceback /nowarn /nologo /Qdiag-disable:8290 /Qdiag-disable:8291 /Qdiag-disable:8293\n";
    print OUTFILE "FLAGS90_MSC = \$(FLAGS_MSC)\n";
    print OUTFILE "FLAGS_DYN =\n";
    print OUTFILE "FLAGS_SER =\n";
    print OUTFILE "FLAGS_OMP = /Qopenmp\n";
    print OUTFILE "FLAGS_MPI =\n";
#    print OUTFILE "#NETCDFROOT = C:\\\PROGRA~2\\\lnetcdf\n";
#    print OUTFILE "!IF DEFINED(NETCDFROOT)\n";
#    print OUTFILE "INCS_SER = /include:\$(NETCDFROOT)\n";
#    print OUTFILE "INCS_OMP = /include:\$(NETCDFROOT)\n";
#    print OUTFILE "INCS_MPI = /include:\"C:\\\PROGRA~1\\\MPICH2\\\include\" /include:\$(NETCDFROOT)\n";
#    print OUTFILE "LIBSC = /link /NODEFAULTLIB:MSVCRT.lib /NODEFAULTLIB:LIBC.lib\n";
#    print OUTFILE "LIBS_SER = \$(NETCDFROOT)\\\lf90_netcdf.lib \$(NETCDFROOT)\\\lf77_netcdf.lib \$(NETCDFROOT)\\\lnetcdf.lib \$(LIBSC)\n";
#    print OUTFILE "LIBS_OMP = \$(NETCDFROOT)\\\lf90_netcdf.lib \$(NETCDFROOT)\\\lf77_netcdf.lib \$(NETCDFROOT)\\\lnetcdf.lib \$(LIBSC)\n";
#    print OUTFILE "LIBS_MPI = C:\\\PROGRA~1\\\MPICH2\\\llib\\\lfmpich2.lib \$(NETCDFROOT)\\\lf90_netcdf.lib \$(NETCDFROOT)\\\lf77_netcdf.lib \$(NETCDFROOT)\\\lnetcdf.lib \$(LIBSC)\n";
#    print OUTFILE "NCF_OBJS = nctablemd.obj agioncmd.obj swn_outnc.obj\n";
#    print OUTFILE "!ELSE\n";
    print OUTFILE "INCS_SER =\n";
    print OUTFILE "INCS_OMP =\n";
    print OUTFILE "INCS_MPI = /include:\"C:\\\PROGRA~1\\\MPICH2\\\include\"\n";
    print OUTFILE "LIBS_SER =\n";
    print OUTFILE "LIBS_OMP =\n";
    print OUTFILE "LIBS_MPI = C:\\\PROGRA~1\\\MPICH2\\\llib\\\lfmpich2.lib\n";
#    print OUTFILE "NCF_OBJS =\n";
#    print OUTFILE "!ENDIF\n";
    print OUTFILE "O_DIR =\n";
    print OUTFILE "OUT = /exe:\n";
    print OUTFILE "EXTO = obj\n";
    print OUTFILE "MAKE = nmake\n";
    print OUTFILE "RM = del\n";
#    print OUTFILE "!IF DEFINED(NETCDFROOT)\n";
#    print OUTFILE "swch = -dos -impi -cvis -netcdf\n";
#    print OUTFILE "!ELSE\n";
    print OUTFILE "swch = -dos -impi -cvis\n";
#    print OUTFILE "!ENDIF\n";
  }
  elsif ( $compiler eq "f90" )
  {
    print OUTFILE "##############################################################################\n";
    print OUTFILE "# IA32_Intel/EM64T_Intel:	Intel Pentium with MS Windows using Compaq Visual Fortran 6.6x.\n";
    print OUTFILE "##############################################################################\n";
    print OUTFILE "F90_SER = f90\n";
    print OUTFILE "F90_OMP = f90\n";
    print OUTFILE "F90_MPI = f90\n";
    print OUTFILE "FLAGS_OPT = /optimize:2\n";
    print OUTFILE "FLAGS_MSC = /assume:byterecl /traceback /names:lowercase /nowarn /nologo\n";
    print OUTFILE "FLAGS90_MSC = \$(FLAGS_MSC)\n";
    print OUTFILE "FLAGS_DYN =\n";
    print OUTFILE "FLAGS_SER =\n";
    print OUTFILE "FLAGS_OMP = /Qopenmp /Qopenmp-link:static\n";
    print OUTFILE "FLAGS_MPI =\n";
    print OUTFILE "NETCDFROOT:\n";
    print OUTFILE "!IF DEFINED(NETCDFROOT)\n";
    print OUTFILE "INCS_SER =\n";
    print OUTFILE "INCS_OMP =\n";
    print OUTFILE "INCS_MPI = /include:\"C:\\\PROGRA~1\\\MPICH2\\\include\"\n";
    print OUTFILE "LIBS_SER =\n";
    print OUTFILE "LIBS_OMP =\n";
    print OUTFILE "LIBS_MPI = C:\\\PROGRA~1\\\MPICH2\\\llib\\\lfmpich2.lib\n";
    print OUTFILE "NCF_OBJS = nctablemd.obj agioncmd.obj swn_outnc.obj\n";
    print OUTFILE "!ELSE\n";
    print OUTFILE "INCS_SER =\n";
    print OUTFILE "INCS_OMP =\n";
    print OUTFILE "INCS_MPI = /include:\"C:\\\PROGRA~1\\\MPICH2\\\include\"\n";
    print OUTFILE "LIBS_SER =\n";
    print OUTFILE "LIBS_OMP =\n";
    print OUTFILE "LIBS_MPI = C:\\\PROGRA~1\\\MPICH2\\\llib\\\lfmpich2.lib\n";
    print OUTFILE "NCF_OBJS =\n";
    print OUTFILE "!ENDIF\n";
    print OUTFILE "O_DIR =\n";
    print OUTFILE "OUT = /exe:\n";
    print OUTFILE "EXTO = obj\n";
    print OUTFILE "MAKE = nmake\n";
    print OUTFILE "RM = del\n";
    print OUTFILE "!IF DEFINED(NETCDFROOT)\n";
    print OUTFILE "swch = -dos -impi -cvis -netcdf\n";
    print OUTFILE "!ELSE\n";
    print OUTFILE "swch = -dos -impi -cvis\n";
    print OUTFILE "!ENDIF\n";
  }
  else
  {
    die "Current Fortran compiler '$compiler' not supported.... \n";
  }
}
elsif ($os =~ /Darwin/i) {
  my $compiler = getcmpl();
  if ( $compiler eq "gfortran" )
  {
    print OUTFILE "##############################################################################\n";
    print OUTFILE "# MAC_GNU:		macOS with GNU compiler gfortran.\n";
    print OUTFILE "##############################################################################\n";
    print OUTFILE "F90_SER = gfortran\n";
    print OUTFILE "F90_OMP = gfortran\n";
    print OUTFILE "F90_MPI = mpif90\n";
    print OUTFILE "FLAGS_OPT = -O\n";
    print OUTFILE "FLAGS_MSC = -w -fno-second-underscore\n";
    print OUTFILE "FLAGS90_MSC = \$(FLAGS_MSC) -ffree-line-length-none\n";
    print OUTFILE "FLAGS_DYN =\n";
    print OUTFILE "FLAGS_SER =\n";
    print OUTFILE "FLAGS_OMP = -fopenmp\n";
    print OUTFILE "FLAGS_MPI =\n";
    print OUTFILE "NETCDFROOT =\n";
    print OUTFILE "ifneq (\$(NETCDFROOT),)\n";
    print OUTFILE "  INCS_SER = -I\$(NETCDFROOT)/include\n";
    print OUTFILE "  INCS_OMP = -I\$(NETCDFROOT)/include\n";
    print OUTFILE "  INCS_MPI = -I\$(NETCDFROOT)/include\n";
    print OUTFILE "  LIBS_SER = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
    print OUTFILE "  LIBS_OMP = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff -static-libgcc\n";
    print OUTFILE "  LIBS_MPI = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
    print OUTFILE "  NCF_OBJS = nctablemd.o agioncmd.o swn_outnc.o\n";
    print OUTFILE "else\n";
    print OUTFILE "  INCS_SER =\n";
    print OUTFILE "  INCS_OMP =\n";
    print OUTFILE "  INCS_MPI =\n";
    print OUTFILE "  LIBS_SER =\n";
    print OUTFILE "  LIBS_OMP = -static-libgcc\n";
    print OUTFILE "  LIBS_MPI =\n";
    print OUTFILE "  NCF_OBJS =\n";
    print OUTFILE "endif\n";
    print OUTFILE "O_DIR = ../work/odir4/\n";
    print OUTFILE "OUT = -o \n";
    print OUTFILE "EXTO = o\n";
    print OUTFILE "MAKE = make\n";
    print OUTFILE "RM = rm -f\n";
    print OUTFILE "ifneq (\$(NETCDFROOT),)\n";
    print OUTFILE "  swch = -unix -netcdf\n";
    print OUTFILE "else\n";
    print OUTFILE "  swch = -unix\n";
    print OUTFILE "endif\n";
  }
  elsif ( $compiler eq "ifort" )
  {
    print OUTFILE "##############################################################################\n";
    print OUTFILE "# MAC_Intel:		macOS with Intel compiler 17.\n";
    print OUTFILE "##############################################################################\n";
    print OUTFILE "F90_SER = ifort\n";
    print OUTFILE "F90_OMP = ifort\n";
    print OUTFILE "# if appropriate, use mpiifort of Intel instead\n";
    print OUTFILE "F90_MPI = mpif90\n";
    print OUTFILE "FLAGS_OPT = -O2\n";
    print OUTFILE "FLAGS_MSC = -W0 -assume byterecl -traceback -diag-disable 8290 -diag-disable 8291 -diag-disable 8293\n";
    print OUTFILE "FLAGS90_MSC = \$(FLAGS_MSC)\n";
    print OUTFILE "FLAGS_DYN = -fPIC\n";
    print OUTFILE "FLAGS_SER =\n";
    print OUTFILE "FLAGS_OMP = -qopenmp\n";
    print OUTFILE "FLAGS_MPI =\n";
    print OUTFILE "NETCDFROOT =\n";
    print OUTFILE "ifneq (\$(NETCDFROOT),)\n";
    print OUTFILE "  INCS_SER = -I\$(NETCDFROOT)/include\n";
    print OUTFILE "  INCS_OMP = -I\$(NETCDFROOT)/include\n";
    print OUTFILE "  INCS_MPI = -I\$(NETCDFROOT)/include\n";
    print OUTFILE "  LIBS_SER = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
    print OUTFILE "  LIBS_OMP = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
    print OUTFILE "  LIBS_MPI = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
    print OUTFILE "  NCF_OBJS = nctablemd.o agioncmd.o swn_outnc.o\n";
    print OUTFILE "else\n";
    print OUTFILE "  INCS_SER =\n";
    print OUTFILE "  INCS_OMP =\n";
    print OUTFILE "  INCS_MPI =\n";
    print OUTFILE "  LIBS_SER =\n";
    print OUTFILE "  LIBS_OMP =\n";
    print OUTFILE "  LIBS_MPI =\n";
    print OUTFILE "  NCF_OBJS =\n";
    print OUTFILE "endif\n";
    print OUTFILE "O_DIR = ../work/odir4/\n";
    print OUTFILE "OUT = -o \n";
    print OUTFILE "EXTO = o\n";
    print OUTFILE "MAKE = make\n";
    print OUTFILE "RM = rm -f\n";
    print OUTFILE "ifneq (\$(NETCDFROOT),)\n";
    print OUTFILE "  swch = -unix -impi -netcdf\n";
    print OUTFILE "else\n";
    print OUTFILE "  swch = -unix -impi\n";
    print OUTFILE "endif\n";
  }
  elsif ( $compiler eq "xlf90" )
  {
    print OUTFILE "##############################################################################\n";
    print OUTFILE "# MAC_IBM:		Mac OS X with IBM Compiler.\n";
    print OUTFILE "##############################################################################\n";
    print OUTFILE "F90_SER = xlf90\n";
    print OUTFILE "F90_OMP =\n";
    print OUTFILE "F90_MPI =\n";
    print OUTFILE "FLAGS_OPT = -O3 -qstrict -qtune=auto -qcache=auto -qalign=4k\n";
    print OUTFILE "FLAGS_MSC = -w -qfixed\n";
    print OUTFILE "FLAGS90_MSC = -w -qfree=f90\n";
    print OUTFILE "FLAGS_DYN =\n";
    print OUTFILE "FLAGS_SER =\n";
    print OUTFILE "FLAGS_OMP =\n";
    print OUTFILE "FLAGS_MPI =\n";
    print OUTFILE "NETCDFROOT =\n";
    print OUTFILE "ifneq (\$(NETCDFROOT),)\n";
    print OUTFILE "  INCS_SER = -I\$(NETCDFROOT)/include\n";
    print OUTFILE "  INCS_OMP = -I\$(NETCDFROOT)/include\n";
    print OUTFILE "  INCS_MPI = -I\$(NETCDFROOT)/include\n";
    print OUTFILE "  LIBS_SER = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
    print OUTFILE "  LIBS_OMP = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
    print OUTFILE "  LIBS_MPI = -L\$(NETCDFROOT)/lib -lnetcdf -lnetcdff\n";
    print OUTFILE "  NCF_OBJS = nctablemd.o agioncmd.o swn_outnc.o\n";
    print OUTFILE "else\n";
    print OUTFILE "  INCS_SER =\n";
    print OUTFILE "  INCS_OMP =\n";
    print OUTFILE "  INCS_MPI =\n";
    print OUTFILE "  LIBS_SER =\n";
    print OUTFILE "  LIBS_OMP =\n";
    print OUTFILE "  LIBS_MPI =\n";
    print OUTFILE "  NCF_OBJS =\n";
    print OUTFILE "endif\n";
    print OUTFILE "O_DIR =\n";
    print OUTFILE "OUT = -o \n";
    print OUTFILE "EXTO = o\n";
    print OUTFILE "MAKE = make\n";
    print OUTFILE "RM = rm -f\n";
    print OUTFILE "ifneq (\$(NETCDFROOT),)\n";
    print OUTFILE "  swch = -unix -netcdf\n";
    print OUTFILE "else\n";
    print OUTFILE "  swch = -unix\n";
    print OUTFILE "endif\n";
  }
  else
  {
    die "Current Fortran compiler '$compiler' not supported.... \n";
  }
}
else
{
   die "Current operating system '$os' not supported.... \n";
}

close(OUTFILE);

# =============================================================================

sub getcmpl {

   my $compiler = $ENV{'FC'};

   unless ( $compiler ) {
      foreach ('ifort','gfortran','f90','ifc','efc','pgf90','xlf90', 'lf95','g95') {
         $compiler = $_;
         my $path  = `which $compiler`;
         last if $path;
      }
   }

   return $compiler;
}

# =============================================================================
