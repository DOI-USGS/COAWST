#!/bin/csh
#
#  Makefile for scrip coawst
#  This make will look in the ../../Compilers folder and
#  use those files as make control options.
#
#  1. You need to set: FORT
#

FORT = ifort
#FORT = gfortran
#FORT  = pgi

# 2. If you want to run in MPI, then set USE_MPI=on,
# otherwise just keep it blank USE_MPI=  [just keep it blank].
# The USE_MPF90 will activate the call to mpif90.

USE_MPI=on               # distributed-memory parallelism
USE_MPIF90=on            # compile with mpif90 script
USE_NETCDF4=on

#
#   3. At the command prompt, type-->      make
#
#   4. YOU  SHOULD NOT NEED TO CHANGE ANYTHING BELOW HERE
#

USE_SCRIP := on

COMPILERS = ../../Compilers

#export HDF5_DISABLE_VERSION_CHECK=1

# This determine the operating system.
OS := $(shell uname -s | sed 's/[\/ ]/-/g')
OS := $(patsubst CYGWIN_%,CYGWIN,$(OS))
OS := $(patsubst MINGW%,MINGW,$(OS))
OS := $(patsubst sn%,UNICOS-sn,$(OS))

include $(COMPILERS)/$(OS)-$(strip $(FORT)).mk

SRCDIR  = .
EXEDIR  = .
     NETCDF_INCFILE := netcdf.inc
#    NETCDF_MODFILE := netcdf.mod
# TYPESIZES_MODFILE := typesizes.mod
FLAGS = $(FFLAGS)
#FLAGS += -g -fbounds-check -fbacktrace -fcheck=all

INCLUDE =

OBJSET  = \
	kinds_mod.o \
	constants.o \
	iounits.o \
	netcdf.o \
	scripwrap_mod.o \
	create_fullgrid.o \
	grids.o \
	read_swan.o \
	read_ww3.o \
	read_roms.o \
	read_wrf.o \
	read_hydro.o \
	remap_vars.o \
	remap_distwgt.o \
	remap_conserv.o \
	remap_bilinear.o \
	remap_bicubic.o \
	timers_scrip.o \
	remap_write.o \
	scrip.o \
	create_masks.o \
	scrip_coawst.o

#all: $(SRCDIR)/$(NETCDF_MODFILE) $(SRCDIR)/$(TYPESIZES_MODFILE)
all: $(SRCDIR)/$(NETCDF_INCFILE)
all: $(EXEDIR)/scrip_coawst

$(EXEDIR)/scrip_coawst: $(OBJSET)
	$(FC) $(FLAGS) $(OBJSET)  -o $(EXEDIR)/scrip_coawst $(LIBS)
##	$(FC) $(FLAGS) $(OBJSET) $(LIB) -o $(EXEDIR)/scrip

kinds_mod.o: $(SRCDIR)/kinds_mod.F $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/kinds_mod.F

constants.o: $(SRCDIR)/constants.F kinds_mod.o $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/constants.F

iounits.o: $(SRCDIR)/iounits.F kinds_mod.o $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/iounits.F

netcdf.o: $(SRCDIR)/netcdf.F kinds_mod.o constants.o $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/netcdf.F

grids.o: $(SRCDIR)/grids.F kinds_mod.o constants.o iounits.o netcdf.o \
	$(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/grids.F

scripwrap_mod.o: $(SRCDIR)/scripwrap_mod.F kinds_mod.o \
	$(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/scripwrap_mod.F

create_fullgrid.o: $(SRCDIR)/create_fullgrid.F kinds_mod.o $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/create_fullgrid.F

read_swan.o: $(SRCDIR)/read_swan.F kinds_mod.o scripwrap_mod.o \
	               create_fullgrid.o     $(INCLUDE)
	$(CPP) $(CPPFLAGS)  $(SRCDIR)/read_swan.F  $(SRCDIR)/read_swan.fpp
	$(FC) $(FLAGS) -c $(SRCDIR)/read_swan.fpp

read_ww3.o: $(SRCDIR)/read_ww3.F kinds_mod.o scripwrap_mod.o \
	               create_fullgrid.o     $(INCLUDE)
	$(CPP) $(CPPFLAGS)  $(SRCDIR)/read_ww3.F  $(SRCDIR)/read_ww3.fpp
	$(FC) $(FLAGS) -c $(SRCDIR)/read_ww3.fpp

read_roms.o: $(SRCDIR)/read_roms.F kinds_mod.o scripwrap_mod.o  \
                        netcdf.o create_fullgrid.o $(INCLUDE)
	$(CPP) $(CPPFLAGS)  $(SRCDIR)/read_roms.F  $(SRCDIR)/read_roms.fpp
	$(FC) $(FLAGS) -c $(SRCDIR)/read_roms.fpp

read_wrf.o: $(SRCDIR)/read_wrf.f kinds_mod.o scripwrap_mod.o  \
                          netcdf.o  create_fullgrid.o $(INCLUDE)
	$(CPP) $(CPPFLAGS)  $(SRCDIR)/read_wrf.F  $(SRCDIR)/read_wrf.fpp
	$(FC) $(FLAGS) -c $(SRCDIR)/read_wrf.fpp

read_hydro.o: $(SRCDIR)/read_hydro.F kinds_mod.o scripwrap_mod.o  \
                        netcdf.o create_fullgrid.o $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/read_hydro.F

remap_vars.o: $(SRCDIR)/remap_vars.F kinds_mod.o constants.o grids.o \
                            $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/remap_vars.F

remap_conserv.o: $(SRCDIR)/remap_conserv.F kinds_mod.o constants.o \
		timers_scrip.o remap_vars.o grids.o $(INCLUDE)
	$(CPP) $(CPPFLAGS)  $(SRCDIR)/remap_conserv.F  $(SRCDIR)/remap_conserv.fpp
	$(FC) $(FLAGS) -c $(SRCDIR)/remap_conserv.fpp

remap_distwgt.o: $(SRCDIR)/remap_distwgt.F kinds_mod.o constants.o \
		remap_vars.o grids.o $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/remap_distwgt.F

remap_bilinear.o: $(SRCDIR)/remap_bilinear.F kinds_mod.o constants.o \
		remap_vars.o grids.o timers_scrip.o $(INCLUDE)
	$(CPP) $(CPPFLAGS)  $(SRCDIR)/remap_bilinear.F  $(SRCDIR)/remap_bilinear.fpp
	$(FC) $(FLAGS) -c $(SRCDIR)/remap_bilinear.fpp

remap_bicubic.o: $(SRCDIR)/remap_bicubic.F kinds_mod.o constants.o \
		remap_vars.o grids.o $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/remap_bicubic.F

timers_scrip.o: $(SRCDIR)/timers_scrip.F kinds_mod.o constants.o $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/timers_scrip.F

remap_write.o: $(SRCDIR)/remap_write.F kinds_mod.o constants.o \
		netcdf.o remap_vars.o grids.o $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/remap_write.F

remap_read.o: $(SRCDIR)/remap_read.F kinds_mod.o constants.o netcdf.o \
		remap_vars.o grids.o $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/remap_read.F

remap.o: $(SRCDIR)/remap.F kinds_mod.o constants.o 
	$(FC) $(FLAGS) -c $(SRCDIR)/remap.F

scrip.o: $(SRCDIR)/scrip.F kinds_mod.o constants.o iounits.o timers_scrip.o \
		remap_vars.o grids.o remap_conserv.o remap_distwgt.o \
		remap_bilinear.o remap_bicubic.o remap_write.o \
		$(INCLUDE)
	$(CPP) $(CPPFLAGS)  $(SRCDIR)/scrip.F  $(SRCDIR)/scrip.fpp
	$(FC) $(FLAGS) -c $(SRCDIR)/scrip.fpp

create_masks.o: $(SRCDIR)/create_masks.F kinds_mod.o scripwrap_mod.o  \
                          scrip.o   $(INCLUDE)
	$(CPP) $(CPPFLAGS)  $(SRCDIR)/create_masks.F  $(SRCDIR)/create_masks.fpp
	$(FC) $(FLAGS) -c $(SRCDIR)/create_masks.fpp

scrip_coawst.o: $(SRCDIR)/scrip_coawst.F kinds_mod.o constants.o iounits.o \
		scripwrap_mod.o create_fullgrid.o read_swan.o read_roms.o  \
		read_wrf.o read_hydro.o read_ww3.o scrip.o create_masks.o \
		$(INCLUDE)
	$(CPP) $(CPPFLAGS)  $(SRCDIR)/scrip_coawst.F  $(SRCDIR)/scrip_coawst.fpp
	$(FC) $(FLAGS) -c $(SRCDIR)/scrip_coawst.fpp


clean: 
	/bin/rm *.o *.mod


$(SRCDIR)/$(NETCDF_INCFILE): | $(SRCDIR)
	cp -f $(NETCDF_INCDIR)/$(NETCDF_INCFILE) $(SRCDIR)
#$(SRCDIR)/$(NETCDF_MODFILE): | $(SRCDIR)
#	cp -f $(NETCDF_INCDIR)/$(NETCDF_MODFILE) $(SRCDIR)
#$(SRCDIR)/$(TYPESIZES_MODFILE): | $(SRCDIR)
#	cp -f $(NETCDF_INCDIR)/$(TYPESIZES_MODFILE) $(SRCDIR)


