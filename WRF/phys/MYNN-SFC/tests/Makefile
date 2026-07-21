FC     = gfortran
FCFLAGS = --coverage -Dstandalone -fprofile-arcs -ftest-coverage -fdefault-real-8 -fdefault-double-8 -g -fbacktrace

SRC=machine.F90 ccpp_kind_types.F90 physcons.F90 ../module_sf_mynnsfc.F90  module_sf_mynnsfc_ccpp_tests.F90 module_state_description.F90 module_model_constants.F90 ../WRF/module_sf_mynnsfc_common.F90 ../module_sf_mynnsfc_land.F90 ../module_sf_mynnsfc_water.F90 ../module_sf_mynnsfc_ice.F90 ../WRF/module_sf_mynnsfc_driver.F90  module_sf_mynnsfc_wrf_tests.F90 driver.F90
OBJ=machine.o ccpp_kind_types.o physcons.o module_sf_mynnsfc.o module_sf_mynnsfc_ccpp_tests.o module_state_description.o module_model_constants.o module_sf_mynnsfc_common.o module_sf_mynnsfc_land.o module_sf_mynnsfc_water.o module_sf_mynnsfc_ice.o module_sf_mynnsfc_driver.o  module_sf_mynnsfc_wrf_tests.o driver.o

#OBJ=$(SRC:.F90=.o)
%.o: %.F90
	$(FC) $(FCFLAGS) -c $(SRC) -lgcov

driver: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) -lgcov

clean:
	@rm -f *.mod *.o driver *.gcda *.gcno *.gcov
