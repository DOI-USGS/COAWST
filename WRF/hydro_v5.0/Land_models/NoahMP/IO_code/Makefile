# Makefile 
#
.SUFFIXES:
.SUFFIXES: .o .F

include ../user_build_options

CMODFLAG = -I./ -I$(MCT_PARAMS_DIR) -I../phys -I../Utility_routines

OBJS_NoahMP = module_NoahMP_hrldas_driver.o

OBJS =	../../../CPL/COAWST_cpl/hydro_coupler.o \
	main_hrldas_driver.o \
	module_hrldas_netcdf_io.o

CPPHRLDAS = -D_HRLDAS_OFFLINE_ $(MOD_OPT)

CPPFLAGS_COAWST  = -I../../../CPL/COAWST_cpl
CPPFLAGS_COAWST += -D$(ROMS_APPLICATION)
CPPFLAGS_COAWST += -D'PROJECT_HEADER="$(MY_HEADER_DIR)/$(HEADER)"'
CPPFLAGS_COAWST += -D'GLOBAL_DEFS="$(GLOBAL_DEFS_FILE)"'

all:	$(OBJS_NoahMP) $(OBJS)

NoahMP : $(OBJS_NoahMP) $(OBJS)

module_NoahMP_hrldas_driver.o: module_NoahMP_hrldas_driver.F ../../../HYDRO_drv/module_HYDRO_drv.o  ../../../Data_Rec/module_namelist.o ../../../Data_Rec/module_RT_data.o
	@echo ""
	$(RM) $(*).f90
	$(CPP) $(CPPFLAGS) $(CPPHRLDAS) $(*).F > $(*).f90
	$(COMPILERF90) -o $(@) -c $(F90FLAGS) $(FREESOURCE) $(MODFLAG) -I. \
	-I../phys -I../Utility_routines -I../../../mod -I ../../../CPL/COAWST_cpl/ $(NETCDFMOD) $(*).f90
	@echo ""

../../../CPL/COAWST_cpl/hydro_coupler.o: ../../../CPL/COAWST_cpl/hydro_coupler.F
	@echo ""
	$(RM) $(*).f90
	$(CPP) $(CPPFLAGS) $(CPPFLAGS_COAWST) $(*).F > $(*).f90
	$(COMPILER90) -o $(@) -c $(F90FLAGS) $(CMODFLAG) $(MODFLAG) $(NETCDFMOD) $(*).f90
	@echo ""
#	ar -r ../../lib/libHYDRO_COAWST.a $(@)

main_hrldas_driver.o: main_hrldas_driver.F
	@echo ""
	$(RM) $(*).f90
	$(CPP) $(CPPFLAGS) $(CPPHRLDAS) $(*).F > $(*).f90
	$(COMPILERF90) -o $(@) -c $(F90FLAGS) $(LDFLAGS) $(FREESOURCE) -I ../MPP -I. \
    -I../../../CPL/COAWST_cpl -I$(MCT_PARAMS_DIR) \
	-I../phys -I../Utility_routines -I../../../mod -I../../../MPP $(NETCDFMOD) $(*).f90
#	$(COMPILERF90) -o $(@) -c $(F90FLAGS) $(FREESOURCE) -I ../MPP -I. \
#	-I../phys -I../Utility_routines $(NETCDFMOD) $(*).f90
	@echo ""

module_hrldas_netcdf_io.o: module_hrldas_netcdf_io.F
	@echo ""
	$(RM) $(*).f90
	$(CPP) $(CPPFLAGS) $(CPPHRLDAS) $(NETCDFMOD) $(*).F > $(*).f90
	$(COMPILERF90) -o $(@) -c $(F90FLAGS) $(FREESOURCE) -I ../MPP -I../Utility_routines $(NETCDFMOD)  $(*).f90
	@echo ""

.F.o:
	@echo ""
	$(RM) $(*).f90
	$(CPP) $(CPPFLAGS) $(CPPHRLDAS) $(*).F > $(*).f90
	$(COMPILERF90) -o $(@) -c $(F90FLAGS) $(FREESOURCE)  -I$(MCT_INCDIR) -I../Utility_routines  $(NETCDFMOD)  -I../Noah $(*).f90
	@echo ""


#
# Dependencies:
#

main_hrldas_driver.o:	$(OBJS_NoahMP) \
#		        ../../../nudging/module_nudging_utils.o \
#		        ../../../nudging/module_stream_nudging.o \
			../../../MPP/mpp_land.o \
                        ../../../Routing/module_NWM_io.o
$(OBJS_NoahMP):	module_hrldas_netcdf_io.o

# This command cleans up object files, etc.
clean:
	$(RM) *.o *.mod *.stb *~ *.f90
