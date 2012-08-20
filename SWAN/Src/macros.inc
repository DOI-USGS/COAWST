##############################################################################
# IA32_PGF:		Intel Pentium with Linux using Portland Group compiler
##############################################################################
F90_SER = pgf90
F90_OMP = pgf90
F90_MPI = mpif90
FLAGS_OPT = -fast
FLAGS_MSC = -Mfixed
FLAGS90_MSC = -Mfree
FLAGS_SER =
FLAGS_OMP = -mp
FLAGS_MPI = -tp barcelona-64
INCS_SER =
INCS_OMP =
INCS_MPI =
LIBS_SER =
LIBS_OMP =
LIBS_MPI =
OUT = -o 
EXTO = o
MAKE = make
RM = rm -f
swch = -unix -impi -mpi -f95 -cdate14
