FC     = gfortran
FCFLAGS = --coverage -fprofile-arcs -ftest-coverage \
	  -Dstandalone -fdefault-real-8 -fdefault-double-8 -g \
          -fbacktrace \
	  -fsanitize=address\
          -fcheck=all \
          -Wall -Wextra -pedantic \
          -Wuninitialized -Wmaybe-uninitialized \
          -ffpe-trap=invalid,zero,overflow \
  	  -finit-real=snan \
	  -finit-integer=-999 \
          -finit-logical=false \
          -finit-character=65 \
          -ffpe-summary=all \
          -I$(NFHOME)/include
LDFLAGS = -L$(NFHOME)/lib -lnetcdff -lnetcdf -fsanitize=address

SRC = module_model_constants.F90 ccpp_kind_types.F90 ../WRF/module_bl_mynnedmf_common.F90 ../module_bl_mynnedmf.F90 ../WRF/module_bl_mynnedmf_driver.F90 module_bl_mynnedmf_wrf_tests.F90 driver.F90
OBJ = module_model_constants.o ccpp_kind_types.o module_bl_mynnedmf_common.o module_bl_mynnedmf.o module_bl_mynnedmf_driver.o module_bl_mynnedmf_wrf_tests.o driver.o

%.o: %.F90
	$(FC) $(FCFLAGS) -c $(SRC) -lgcov

driver: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LDFLAGS) -lgcov

clean:
	@rm -f *.mod *.o driver *.gcda *.gcno *.gcov
