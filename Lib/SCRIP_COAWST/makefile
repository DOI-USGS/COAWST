#!/bin/csh
#
#  Makefile for scrip coawst
#  You need to set: FORT
#

#FORT = ifort
#FORT = gfortran
FORT  = pgi

#
#   YOU  SHOULD NOT NEED TO CHANGE ANYTHING BELOW HERE
#

COMPILERS = ../../Compilers

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

INCLUDE =

OBJSET  = \
	kinds_mod.o \
	constants.o \
	iounits.o \
	netcdf.o \
	scripwrap_mod.o \
	create_fullgrid.o \
	read_swan.o \
	read_roms.o \
	read_wrf.o \
	grids.o \
	remap_vars.o \
	remap_distwgt.o \
	remap_conserv.o \
	remap_bilinear.o \
	remap_bicubic.o \
	timers.o \
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

kinds_mod.o: $(SRCDIR)/kinds_mod.f $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/kinds_mod.f

constants.o: $(SRCDIR)/constants.f kinds_mod.o $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/constants.f

iounits.o: $(SRCDIR)/iounits.f kinds_mod.o $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/iounits.f

netcdf.o: $(SRCDIR)/netcdf.f kinds_mod.o constants.o $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/netcdf.f

grids.o: $(SRCDIR)/grids.f kinds_mod.o constants.o iounits.o netcdf.o \
	$(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/grids.f

scripwrap_mod.o: $(SRCDIR)/scripwrap_mod.f kinds_mod.o \
	$(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/scripwrap_mod.f

create_fullgrid.o: $(SRCDIR)/create_fullgrid.f kinds_mod.o $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/create_fullgrid.f

read_swan.o: $(SRCDIR)/read_swan.f kinds_mod.o scripwrap_mod.o \
	               create_fullgrid.o     $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/read_swan.f

read_roms.o: $(SRCDIR)/read_roms.f kinds_mod.o scripwrap_mod.o  \
                        netcdf.o create_fullgrid.o $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/read_roms.f

read_wrf.o: $(SRCDIR)/read_wrf.f kinds_mod.o scripwrap_mod.o  \
                          netcdf.o  create_fullgrid.o $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/read_wrf.f

remap_vars.o: $(SRCDIR)/remap_vars.f kinds_mod.o constants.o grids.o \
                            $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/remap_vars.f

remap_conserv.o: $(SRCDIR)/remap_conserv.f kinds_mod.o constants.o \
		timers.o remap_vars.o grids.o $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/remap_conserv.f

remap_distwgt.o: $(SRCDIR)/remap_distwgt.f kinds_mod.o constants.o \
		remap_vars.o grids.o $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/remap_distwgt.f

remap_bilinear.o: $(SRCDIR)/remap_bilinear.f kinds_mod.o constants.o \
		remap_vars.o grids.o timers.o $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/remap_bilinear.f

remap_bicubic.o: $(SRCDIR)/remap_bicubic.f kinds_mod.o constants.o \
		remap_vars.o grids.o $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/remap_bicubic.f

timers.o: $(SRCDIR)/timers.f kinds_mod.o constants.o $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/timers.f

remap_write.o: $(SRCDIR)/remap_write.f kinds_mod.o constants.o \
		netcdf.o remap_vars.o grids.o $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/remap_write.f

remap_read.o: $(SRCDIR)/remap_read.f kinds_mod.o constants.o netcdf.o \
		remap_vars.o grids.o $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/remap_read.f

remap.o: $(SRCDIR)/remap.f kinds_mod.o constants.o 
	$(FC) $(FLAGS) -c $(SRCDIR)/remap.f

scrip.o: $(SRCDIR)/scrip.f kinds_mod.o constants.o iounits.o timers.o \
		remap_vars.o grids.o remap_conserv.o remap_distwgt.o \
		remap_bilinear.o remap_bicubic.o remap_write.o \
		$(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/scrip.f

create_masks.o: $(SRCDIR)/create_masks.f kinds_mod.o scripwrap_mod.o  \
                          scrip.o   $(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/create_masks.f

scrip_coawst.o: $(SRCDIR)/scrip_coawst.f kinds_mod.o constants.o iounits.o \
		scripwrap_mod.o create_fullgrid.o read_swan.o read_roms.o  \
		read_wrf.o scrip.o create_masks.o \
		$(INCLUDE)
	$(FC) $(FLAGS) -c $(SRCDIR)/scrip_coawst.f


clean: 
	/bin/rm *.o *.mod


$(SRCDIR)/$(NETCDF_INCFILE): | $(SRCDIR)
	cp -f $(NETCDF_INCDIR)/$(NETCDF_INCFILE) $(SRCDIR)
#$(SRCDIR)/$(NETCDF_MODFILE): | $(SRCDIR)
#	cp -f $(NETCDF_INCDIR)/$(NETCDF_MODFILE) $(SRCDIR)
#$(SRCDIR)/$(TYPESIZES_MODFILE): | $(SRCDIR)
#	cp -f $(NETCDF_INCDIR)/$(TYPESIZES_MODFILE) $(SRCDIR)


