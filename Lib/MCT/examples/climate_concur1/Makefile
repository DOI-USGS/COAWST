
SHELL           = /bin/sh

# SOURCE FILES

SRCS_F90	= master.F90 coupler.F90 model.F90

OBJS_ALL	= $(SRCS_F90:.F90=.o)

# MACHINE AND COMPILER FLAGS

include ../../Makefile.conf

# ADDITIONAL FLAGS SPECIFIC FOR UTMCT COMPILATION

MCTLIBS    = -L$(MPEUPATH) -L$(MCTPATH) -lmct -lmpeu 
UTLDFLAGS  = $(REAL8)
UTCMPFLAGS = $(REAL8) $(INCFLAG)$(MPEUPATH) $(INCFLAG)$(MCTPATH)

# TARGETS

all: climate

climate: $(OBJS_ALL) 
	$(FC) -o $@ $(OBJS_ALL) $(FCFLAGS) $(UTLDFLAGS) $(MCTLIBS) $(MPILIBS)

# RULES

.SUFFIXES:
.SUFFIXES: .F90 .o

$(F90RULE):
	$(FC) -c $(INCPATH) $(DEFS) $(FCFLAGS) $(F90FLAGS) $(UTCMPFLAGS) $*.F90

$(F90RULECPP):
	$(FPP) $(DEFS) $(FPPFLAGS) $*.F90 $*.f90
	$(FC) -c $(INCPATH) $(FCFLAGS) $(F90FLAGS) $(UTCMPFLAGS) $*.f90
	$(RM) $*.f90

clean:
	${RM} *.o *.mod climate

# DEPENDENCIES:

$(OBJS_ALL): $(MCTPATH)/libmct.a











