# Makefile for configuring, building and installing SWAN

.PHONY: allclean build clobber config help install uninstall
.SILENT:

CMAKE_FLAGS = -G Ninja
ifneq ($(fc),)
    CMAKE_FLAGS += -DCMAKE_Fortran_COMPILER=$(fc)
endif
ifneq ($(prefix),)
    CMAKE_FLAGS += -DCMAKE_INSTALL_PREFIX=$(prefix)
endif
ifneq ($(openmp),)
    CMAKE_FLAGS += -DOPENMP=$(openmp)
endif
ifneq ($(mpi),)
    CMAKE_FLAGS += -DMPI=$(mpi)
endif
ifneq ($(metis),)
    CMAKE_FLAGS += -DMETIS=$(metis)
endif
ifneq ($(netcdf),)
    CMAKE_FLAGS += -DNETCDF=$(netcdf)
endif

build:
	cmake --build build

help:
	@echo "make config    - configures the SWAN build"
	@echo "make           - builds the SWAN software"
	@echo "make install   - installs the SWAN package"
	@echo "make clobber   - removes the build files"
	@echo "make uninstall - uninstalls the SWAN package"

allclean:
	cmake -P clobber.cmake

clobber:
	cmake -P clobber.cmake

config: clobber
	mkdir build
	cd build && cmake $(CURDIR) $(CMAKE_FLAGS)

install:
	cmake --install build

uninstall:
	cat build/install_manifest.txt | xargs rm -vf
