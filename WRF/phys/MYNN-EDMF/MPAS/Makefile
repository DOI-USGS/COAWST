.SUFFIXES: .F90 .o

all: dummy mynnedmf

dummy:
	echo "****** compiling MYNN-EDMF PBL Scheme ******"

OBJS = \
	module_bl_mynnedmf_common.o	\
	module_bl_mynnedmf.o		\
	module_bl_mynnedmf_driver.o

mynnedmf: $(OBJS)
	ar -ru ./../../libphys.a $(OBJS)

# DEPENDENCIES:
module_bl_mynnedmf_common.o: \
	../../mpas_atmphys_constants.o

module_bl_mynnedmf.o: \
	module_bl_mynnedmf_common.o

module_bl_mynnedmf_driver.o: \
	module_bl_mynnedmf.o

clean:
	$(RM) *.f90 *.o *.mod
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F90.o:
ifeq "$(GEN_F90)" "true"
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../.. -I../../physics_wrf -I../../physics_mmm -I../../../../framework -I../../../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F90 $(CPPINCLUDES) $(FCINCLUDES) -I../.. -I../../physics_wrf -I../../physics_mmm -I../../../../framework -I../../../../external/esmf_time_f90
endif
