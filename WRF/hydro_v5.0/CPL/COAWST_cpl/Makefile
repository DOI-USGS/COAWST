# Makefile 
#
.SUFFIXES:
.SUFFIXES: .o .F


include ../../macros

MODFLAG = -I./ -I$(MCT_PARAMS_DIR)

OBJS = hydro_coupler.o

CPPFLAGS_COAWST  = -I../../../CPL/COAWST_cpl
CPPFLAGS_COAWST += -D$(ROMS_APPLICATION)
CPPFLAGS_COAWST += -D'PROJECT_HEADER="$(MY_HEADER_DIR)/$(HEADER)"'
CPPFLAGS_COAWST += -D'GLOBAL_DEFS="$(GLOBAL_DEFS_FILE)"'

all:	$(OBJS) 

.F.o:
	@echo ""
	$(CPP) $(CPPFLAGS) $(CPPFLAGS_COAWST) -I$(NETCDFINC) $(*).F > $(*).f
	$(COMPILER90) -o $(@) $(F90FLAGS) $(MODFLAG) -I../../Land_models/NoahMP/phys $(*).f
#	$(RMD) $(*).f
	@echo ""
	ar -r ../../lib/libHYDRO_COAWST.a $(@)

#
# Dependencies:
#

hydro_coupler.o:

clean:
	rm -f *.o *.mod *.stb *~ *.f Noah_hrldas_beta
