# $Id: makefile 787 2008-10-14 00:37:50Z jcwarner $
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
# Copyright (c) 2002-2008 The ROMS/TOMS Group             Kate Hedstrom :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                                                                       :::
#  ROMS/TOMS Framework Master Makefile                                  :::
#                                                                       :::
#  This makefile is designed to work only with GNU Make version 3.80 or :::
#  higher. It can be used in any architecture provided that there is a  :::
#  machine/compiler rules file in the  "Compilers"  subdirectory.  You  :::
#  may need to modify the rules file to specify the  correct path  for  :::
#  the NetCDF and ARPACK libraries. The ARPACK library is only used in  :::
#  the Generalized Stability Theory analysis and Laczos algorithm.      :::
#                                                                       :::
#  If appropriate,  the USER needs to modify the  macro definitions in  :::
#  in user-defined section below.  To activate an option set the macro  :::
#  to "on". For example, if you want to compile with debugging options  :::
#  set:                                                                 :::
#                                                                       :::
#      USE_DEBUG := on                                                  :::
#                                                                       :::
#  Otherwise, leave macro definition blank.                             :::
#                                                                       :::
#  The USER needs to provide a value for the  macro FORT.  Choose  the  :::
#  appropriate value from the list below.                               :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

NEED_VERSION := 3.80 3.81
$(if $(filter $(MAKE_VERSION),$(NEED_VERSION)),,        \
 $(error This makefile requires one of GNU make version $(NEED_VERSION).))

#--------------------------------------------------------------------------
#  Initialize some things.
#--------------------------------------------------------------------------

  sources    := 
  libraries  :=

#==========================================================================
#  Start of user-defined options. In some macro definitions below: "on" or
#  any other string means TRUE while blank (or spaces) is FALSE.
#==========================================================================
#
#  The CPP option defining a particular application is specified below.
#  See header file "ROMS/Include/cppdefs.h" for all available idealized
#  and realistic applications CPP flags. For example, to activate the
#  upwelling test case (UPWELLING) set:
#
#    ROMS_APPLICATION ?= UPWELLING
#
#  Notice that this makefile will include the associated application header
#  file, which is located either in the "ROMS/Include" or MY_HEADER_DIR
#  directory.  This makefile is designed to search in both directories.
#  The only constrain is that the application CPP option must be unique
#  and header file name is the lowercase value of ROMS_APPLICATION with
#  the .h extension. For example, the upwelling application includes the
#  "upwelling.h" header file.  

ROMS_APPLICATION ?= INLET_TEST

#  If application header files is not located in "ROMS/Include",
#  provide an alternate directory FULL PATH.

#MY_HEADER_DIR ?= /raid1/jcwarner/Models/COAWST/Projects/Inlet_test/Refined
#MY_HEADER_DIR ?= /cygdrive/d/data/models/COAWST/Projects/Inlet_test/Refined
MY_HEADER_DIR ?= /cygdrive/d/data/models/COAWST/Projects/Inlet_test/Coupled

#  If your application requires analytical expressions and they are not
#  located in "ROMS/Functionals", provide an alternate directory.
#  Notice that a set analytical expressions templates can be found in
#  "User/Functionals".

#MY_ANALYTICAL_DIR ?= /raid1/jcwarner/Models/COAWST/Projects/Inlet_test/Refined
#MY_ANALYTICAL_DIR ?= /cygdrive/d/data/models/COAWST/Projects/Inlet_test/Refined
MY_ANALYTICAL_DIR ?= /cygdrive/d/data/models/COAWST/Projects/Inlet_test/Coupled

#  Sometimes it is desirable to activate one or more CPP options to
#  run different variants of the same application without modifying
#  its header file. If this is the case, specify such options here
#  using the -D syntax.  For example, to write time-averaged fields
#  set:
#
#    MY_CPP_FLAGS ?= -DAVERAGES
#

MY_CPP_FLAGS ?=

#  Set number of ROMS nested and/or composed grid.  Currently, only
#  one grid is supported.  This option will be available in the near
#  future.

 NestedGrids ?= 1

#  Activate debugging compiler options:

   USE_DEBUG ?= on

#  If parallel applications, use at most one of these definitions
#  (leave both definitions blank in serial applications):

     USE_MPI ?= on
  USE_OpenMP ?=

#  If distributed-memory, turn on compilation via the script "mpif90".
#  This is needed in some Linux operating systems. In some systems with
#  native MPI libraries the compilation does not require MPICH type
#  scripts. This macro is also convient when there are several fortran
#  compiliers (ifort, pgf90, pathf90) in the system that use mpif90.
#  In this, case the user need to select the desired compiler below and
#  turn on both USE_MPI and USE_MPIF90 macros.

  USE_MPIF90 ?= on

#  If applicable, activate 64-bit compilation:

   USE_LARGE ?=

#  If applicable, link with NetCDF-4 library. Notice that the NetCDF-4
#  library needs both the HDF5 and MPI libraries.

 USE_NETCDF4 ?=

#--------------------------------------------------------------------------
#  We are going to include a file with all the settings that depend on
#  the system and the compiler. We are going to build up the name of the
#  include file using information on both. Set your compiler here from
#  the following list:
#
#  Operating System        Compiler(s)
#
#     AIX:                    xlf
#     ALPHA:                  f90
#     CYGWIN:                 g95, df, ifort
#     Darwin:                 f90, xlf
#     IRIX:                   f90
#     Linux:                  ftn, ifc, ifort, pgi, path, g95, gfortran
#     SunOS:                  f95
#     UNICOS-mp:              ftn
#     SunOS/Linux:            ftn (Cray cross-compiler)
#
#  Feel free to send us additional rule files to include! Also, be sure
#  to check the appropriate file to make sure it has the right paths to
#  NetCDF and so on.
#--------------------------------------------------------------------------

        FORT ?= ifort
#       FORT ?= pgi

#--------------------------------------------------------------------------
#  Set directory for executable.
#--------------------------------------------------------------------------

      BINDIR ?= .

#==========================================================================
#  End of user-defined options. See also the machine-dependent include
#  file being used above.
#==========================================================================

#--------------------------------------------------------------------------
#  Set directory for temporary objects.
#--------------------------------------------------------------------------

SCRATCH_DIR ?= Build
 clean_list := core *.ipo $(SCRATCH_DIR)

ifeq "$(strip $(SCRATCH_DIR))" "."
  clean_list := core *.o *.oo *.mod *.f90 lib*.a *.bak
  clean_list += $(CURDIR)/*.ipo
endif
ifeq "$(strip $(SCRATCH_DIR))" "./"
  clean_list := core *.o *.oo *.ipo *.mod *.f90 lib*.a *.bak
  clean_list += $(CURDIR)/*.ipo
endif

#--------------------------------------------------------------------------
#  Set Pattern rules.
#--------------------------------------------------------------------------

%.o: %.F

%.o: %.f90
	cd $(SCRATCH_DIR); $(FC) -c $(FFLAGS) $(notdir $<)

%.f90: %.F
	$(CPP) $(CPPFLAGS) $(MY_CPP_FLAGS) $< > $*.f90
	$(CLEAN) $*.f90

CLEAN := ROMS/Bin/cpp_clean

#--------------------------------------------------------------------------
#  Set C-preprocessing flags associated with ROMS application. They are
#  used in "ROMS/Include/cppdefs.h" to include the appropriate application
#  header file.
#--------------------------------------------------------------------------

ifdef ROMS_APPLICATION
        HEADER := $(addsuffix .h, \
			$(shell echo ${ROMS_APPLICATION} | tr [A-Z] [a-z]))
 ROMS_CPPFLAGS := -D$(ROMS_APPLICATION)
 ROMS_CPPFLAGS += -D'HEADER="$(HEADER)"'
 ifdef MY_HEADER_DIR
  ROMS_CPPFLAGS += -D'ROMS_HEADER="$(MY_HEADER_DIR)/$(HEADER)"'
 else
  ROMS_CPPFLAGS += -D'ROMS_HEADER="$(HEADER)"'
 endif
 ifdef MY_CPP_FLAGS
  ROMS_CPPFLAGS += $(MY_CPP_FLAGS)
 endif
endif

#--------------------------------------------------------------------------
#  Internal macro definitions used to select the code to compile and
#  additional libraries to link. It uses the CPP activated in the
#  header file ROMS/Include/cppdefs.h to determine macro definitions.
#--------------------------------------------------------------------------

MAKE_MACROS := Compilers/make_macros.mk

ifneq "$(MAKECMDGOALS)" "clean"
 MACROS := $(shell cpp -P $(ROMS_CPPFLAGS) Compilers/make_macros.h > \
		$(MAKE_MACROS); $(CLEAN) $(MAKE_MACROS))

 GET_MACROS := $(wildcard $(SCRATCH_DIR)/make_macros.*)

 ifdef GET_MACROS
  include $(SCRATCH_DIR)/make_macros.mk
  $(if ,, $(warning INCLUDING FILE $(SCRATCH_DIR)/make_macros.mk \
                    WHICH CONTAINS APPLICATION-DEPENDENT MAKE DEFINITIONS))
 else
  include $(MAKE_MACROS)
  $(if ,, $(warning INCLUDING FILE $(MAKE_MACROS) \
                   WHICH CONTAINS APPLICATION-DEPENDENT MAKE DEFINITIONS))
 endif
endif

clean_list += $(MAKE_MACROS)

#--------------------------------------------------------------------------
#  Make functions for putting the temporary files in $(SCRATCH_DIR)
#  DO NOT modify this section; spaces and blank lineas are needed.
#--------------------------------------------------------------------------

# $(call source-dir-to-binary-dir, directory-list)
source-dir-to-binary-dir = $(addprefix $(SCRATCH_DIR)/, $(notdir $1))

# $(call source-to-object, source-file-list)
source-to-object = $(call source-dir-to-bindary-dir,   \
                   $(subst .F,.o,$1))

# $(call make-library, library-name, source-file-list)
define make-library
   libraries += $(SCRATCH_DIR)/$1
   sources   += $2

   $(SCRATCH_DIR)/$1: $(call source-dir-to-binary-dir,    \
                      $(subst .F,.o,$2))
	$(AR) $(ARFLAGS) $$@ $$^
	$(RANLIB) $$@
endef

# $(call f90-source, source-file-list)
f90-source = $(call source-dir-to-binary-dir,     \
                   $(subst .F,.f90,$1))

# $(compile-rules)
define compile-rules
  $(foreach f, $(local_src),       \
    $(call one-compile-rule,$(call source-to-object,$f), \
    $(call f90-source,$f),$f))
endef

# $(call one-compile-rule, binary-file, f90-file, source-files)
define one-compile-rule
  $1: $2 $3
	cd $$(SCRATCH_DIR); $$(FC) -c $$(FFLAGS) $(notdir $2)

  $2: $3
	$$(CPP) $$(CPPFLAGS) $$(MY_CPP_FLAGS) $$< > $$@
	$$(CLEAN) $$@

endef

#--------------------------------------------------------------------------
#  Set ROMS/TOMS executable file name.
#--------------------------------------------------------------------------

BIN := $(BINDIR)/coawstS
ifdef USE_DEBUG
  BIN := $(BINDIR)/coawstG
else
 ifdef USE_MPI
   BIN := $(BINDIR)/coawstM
 endif
 ifdef USE_OpenMP
   BIN := $(BINDIR)/coawstO
 endif
endif

#--------------------------------------------------------------------------
#  Set name of module files for netCDF F90 interface. On some platforms
#  these will need to be overridden in the machine-dependent include file.
#--------------------------------------------------------------------------

   NETCDF_MODFILE := netcdf.mod
TYPESIZES_MODFILE := typesizes.mod

#--------------------------------------------------------------------------
#  "uname -s" should return the OS or kernel name and "uname -m" should
#  return the CPU or hardware name. In practice the results can be pretty
#  flaky. Run the results through sed to convert "/" and " " to "-",
#  then apply platform-specific conversions.
#--------------------------------------------------------------------------

OS := $(shell uname -s | sed 's/[\/ ]/-/g')
OS := $(patsubst CYGWIN_%,CYGWIN,$(OS))
OS := $(patsubst MINGW%,MINGW,$(OS))
OS := $(patsubst sn%,UNICOS-sn,$(OS))

CPU := $(shell uname -m | sed 's/[\/ ]/-/g')

SVNREV ?= $(shell svnversion -n .)

ROOTDIR := $(shell pwd)

COMPILERS := ./Compilers

ifndef FORT
  $(error Variable FORT not set)
endif

ifneq "$(MAKECMDGOALS)" "clean"
  include $(COMPILERS)/$(OS)-$(strip $(FORT)).mk
endif

ifdef USE_MPI
 ifdef USE_OpenMP
  $(error You cannot activate USE_MPI and USE_OpenMP at the same time!)
 endif
endif

#--------------------------------------------------------------------------
#  Pass the platform variables to the preprocessor as macros. Convert to
#  valid, upper-case identifiers. Attach ROMS application  CPP options.
#--------------------------------------------------------------------------

CPPFLAGS += -D$(shell echo ${OS} | tr "-" "_" | tr [a-z] [A-Z])
CPPFLAGS += -D$(shell echo ${CPU} | tr "-" "_" | tr [a-z] [A-Z])
CPPFLAGS += -D$(shell echo ${FORT} | tr "-" "_" | tr [a-z] [A-Z])

CPPFLAGS += -D'ROOT_DIR="$(ROOTDIR)"'
ifdef ROMS_APPLICATION
  CPPFLAGS  += $(ROMS_CPPFLAGS)
  CPPFLAGS  += -DNestedGrids=$(NestedGrids)
  MDEPFLAGS += -DROMS_HEADER="$(HEADER)"
endif

ifndef MY_ANALYTICAL_DIR
  MY_ANALYTICAL_DIR := $(ROOTDIR)/ROMS/Functionals
endif
ifeq (,$(findstring ROMS/Functionals,$(MY_ANALYTICAL_DIR)))
  MY_ANALYTICAL := on
endif
CPPFLAGS += -D'ANALYTICAL_DIR="$(MY_ANALYTICAL_DIR)"'

ifdef MY_ANALYTICAL
  CPPFLAGS += -D'MY_ANALYTICAL="$(MY_ANALYTICAL)"'
endif

ifdef SVNREV
  CPPFLAGS += -D'SVN_REV="$(SVNREV)"'
else
  SVNREV := $(shell grep Revision ./ROMS/Version | sed 's/.* \([0-9]*\) .*/\1/')
  CPPFLAGS += -D'SVN_REV="$(SVNREV)"'
endif  

#--------------------------------------------------------------------------
#  Build target directories.
#--------------------------------------------------------------------------

.PHONY: all

all: $(SCRATCH_DIR) $(SCRATCH_DIR)/MakeDepend $(BIN) rm_macros

 modules  :=
ifdef USE_ADJOINT
 modules  +=	ROMS/Adjoint
endif
ifdef USE_REPRESENTER
 modules  +=	ROMS/Representer
endif
ifdef USE_SEAICE
 modules  +=	ROMS/SeaIce
endif
ifdef USE_TANGENT
 modules  +=	ROMS/Tangent
endif
 modules  +=	ROMS/Nonlinear \
		ROMS/Functionals \
		ROMS/Utility \
		ROMS/Modules

 includes :=	ROMS/Include
ifdef USE_ADJOINT
 includes +=	ROMS/Adjoint
endif
ifdef USE_REPRESENTER
 includes +=	ROMS/Representer
endif
ifdef USE_SEAICE
 includes +=	ROMS/SeaIce
endif
ifdef USE_TANGENT
 includes +=	ROMS/Tangent
endif
 includes +=	ROMS/Nonlinear \
		ROMS/Utility \
		ROMS/Drivers

ifdef MY_ANALYTICAL
 includes +=	$(MY_ANALYTICAL_DIR)
endif
 includes +=    ROMS/Functionals

ifdef MY_HEADER_DIR
 includes +=	$(MY_HEADER_DIR)
endif

ifdef USE_SWAN
 modules  +=	SWAN/Src
 includes +=	SWAN/Src
endif

ifdef USE_REFDIF
 modules  +=	REFDIF
 includes +=	REFDIF
endif

 modules  +=	Master
 includes +=	Master Compilers

vpath %.F $(modules)
vpath %.h $(includes)
vpath %.f90 $(SCRATCH_DIR)
vpath %.o $(SCRATCH_DIR)

include $(addsuffix /Module.mk,$(modules))

MDEPFLAGS += $(patsubst %,-I %,$(includes)) --silent --moddir $(SCRATCH_DIR)

CPPFLAGS += $(patsubst %,-I%,$(includes))

ifdef MY_HEADER_DIR
  CPPFLAGS += -D'HEADER_DIR="$(MY_HEADER_DIR)"'
else
  CPPFLAGS += -D'HEADER_DIR="./ROMS/Include"'
endif

$(SCRATCH_DIR):
	$(shell $(TEST) -d $(SCRATCH_DIR) || $(MKDIR) $(SCRATCH_DIR) )

#--------------------------------------------------------------------------
#  Add profiling.
#--------------------------------------------------------------------------

# FFLAGS += -check bounds                 # ifort
# FFLAGS += -C                            # pgi
# FFLAGS += -xpg                          # Sun
# FFLAGS += -pg                           # g95
# FFLAGS += -qp                           # ifort
# FFLAGS += -Mprof=func,lines             # pgi
# FFLAGS += -Mprof=mpi,lines              # pgi
# FFLAGS += -Mprof=mpi,hwcts              # pgi
# FFLAGS += -Mprof=func                   # pgi

#--------------------------------------------------------------------------
#  Special CPP macros for mod_strings.F
#--------------------------------------------------------------------------

$(SCRATCH_DIR)/mod_strings.f90: CPPFLAGS += -DMY_OS='"$(OS)"' \
              -DMY_CPU='"$(CPU)"' -DMY_FORT='"$(FORT)"' \
              -DMY_FC='"$(FC)"' -DMY_FFLAGS='"$(FFLAGS)"'

#--------------------------------------------------------------------------
#  ROMS/TOMS libraries.
#--------------------------------------------------------------------------

MYLIB := libocean.a

.PHONY: libraries

libraries: $(libraries)

#--------------------------------------------------------------------------
#  Target to create ROMS/TOMS dependecies.
#--------------------------------------------------------------------------

$(SCRATCH_DIR)/$(NETCDF_MODFILE): | $(SCRATCH_DIR)
	cp -f $(NETCDF_INCDIR)/$(NETCDF_MODFILE) $(SCRATCH_DIR)

$(SCRATCH_DIR)/$(TYPESIZES_MODFILE): | $(SCRATCH_DIR)
	cp -f $(NETCDF_INCDIR)/$(TYPESIZES_MODFILE) $(SCRATCH_DIR)

$(SCRATCH_DIR)/MakeDepend: makefile \
                           $(SCRATCH_DIR)/$(NETCDF_MODFILE) \
                           $(SCRATCH_DIR)/$(TYPESIZES_MODFILE) \
                           | $(SCRATCH_DIR)
	$(SFMAKEDEPEND) $(MDEPFLAGS) $(sources) > $(SCRATCH_DIR)/MakeDepend
	cp -p $(CURDIR)/$(MAKE_MACROS) $(SCRATCH_DIR)

.PHONY: depend

SFMAKEDEPEND := ./ROMS/Bin/sfmakedepend

depend: $(SCRATCH_DIR)
	$(SFMAKEDEPEND) $(MDEPFLAGS) $(sources) > $(SCRATCH_DIR)/MakeDepend 

ifneq "$(MAKECMDGOALS)" "clean"
  -include $(SCRATCH_DIR)/MakeDepend
endif

#--------------------------------------------------------------------------
#  Target to create ROMS/TOMS tar file.
#--------------------------------------------------------------------------

.PHONY: tarfile

tarfile:
#		tar --exclude=".svn" --exclude Output -cvf coawst_v1.1.tar *
		tar --exclude=".svn" -cvf coawst_v1.0.tar run_* *.TBL *.tbl RRTM* makefile ROMS/ SWAN/ WRF/ Master/ Tools/ Compilers/ Projects/JOE_TC Projects/Rip_current Projects/wetdry Projects/Visser Projects/Hole7 Projects/Griz_Bay Projects/Dogbone Projects/Inlet_test

.PHONY: zipfile

zipfile:
		zip -r roms-3_0.zip *

.PHONY: gzipfile

gzipfile:
		gzip -v roms-3_0.gzip *

#--------------------------------------------------------------------------
#  Cleaning targets.
#--------------------------------------------------------------------------

.PHONY: clean

clean:
	$(RM) -r $(clean_list)

.PHONY: rm_macros

rm_macros:
	$(RM) -r $(CURDIR)/$(MAKE_MACROS)

#--------------------------------------------------------------------------
#  A handy debugging target. This will allow to print the value of any
#  makefile defined macro (see http://tinyurl.com/8ax3j). For example,
#  to find the value of CPPFLAGS execute:
#
#        gmake print-CPPFLAGS
#  or
#        make print-CPPFLAGS
#--------------------------------------------------------------------------

.PHONY: print-%

print-%:
	@echo $* = $($*)
# DO NOT DELETE THIS LINE - used by make depend
esmf_roms.o: cppdefs.h globaldefs.h
esmf_roms.o: Build/distribute.o Build/mod_coupler.o Build/mod_forces.o
esmf_roms.o: Build/mod_grid.o Build/mod_iounits.o Build/mod_ncparam.o
esmf_roms.o: Build/mod_ocean.o Build/mod_parallel.o Build/mod_param.o
esmf_roms.o: Build/mod_scalars.o Build/mod_stepping.o Build/ocean_control.o
esmf_roms.o: Build/roms_import.o

master.o: esmf_coupler.h ocean.h mct_coupler.h cppdefs.h globaldefs.h
master.o: Build/esmf_roms.o Build/mod_coupler.o Build/mod_iounits.o
master.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o
master.o: Build/ocean_control.o Build/ocean_coupler.o Build/waves_control.o

ocean_control.o: ad_ocean.h w4dpsas_ocean.h tl_ocean.h correlation.h
ocean_control.o: s4dvar_ocean.h w4dvar_ocean.h afte_ocean.h fsv_ocean.h
ocean_control.o: grad_ocean.h so_semi_ocean.h is4dvar_lanczos_ocean.h cppdefs.h
ocean_control.o: globaldefs.h convolution.h is4dvar_ocean.h pert_ocean.h
ocean_control.o: nl_ocean.h optobs_ocean.h fte_ocean.h op_ocean.h
ocean_control.o: tlcheck_ocean.h obs_sen_ocean.h picard_ocean.h adsen_ocean.h
ocean_control.o: is4dvar_ocean_old.h rp_ocean.h symmetry.h
ocean_control.o: Build/analytical.o Build/back_cost.o Build/back_cov.o
ocean_control.o: Build/back_step.o Build/cgradient.o Build/cost_grad.o
ocean_control.o: Build/cost_norm.o Build/descent.o Build/distribute.o
ocean_control.o: Build/dotproduct.o Build/downhill.o Build/impulse.o
ocean_control.o: Build/ini_adjust.o Build/ini_fields.o Build/mod_forces.o
ocean_control.o: Build/mod_fourdvar.o Build/mod_grid.o Build/mod_iounits.o
ocean_control.o: Build/mod_ncparam.o Build/mod_netcdf.o Build/mod_ocean.o
ocean_control.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o
ocean_control.o: Build/mod_stepping.o Build/mod_storage.o Build/normalization.o
ocean_control.o: Build/ocean_coupler.o Build/packing.o Build/propagator.o

ocean_coupler.o: mct_roms_wrf.h tile.h mct_roms_refdif.h mct_roms_swan.h
ocean_coupler.o: cppdefs.h globaldefs.h
ocean_coupler.o: Build/bc_2d.o Build/exchange_2d.o Build/mod_coupler.o
ocean_coupler.o: Build/mod_forces.o Build/mod_grid.o Build/mod_iounits.o
ocean_coupler.o: Build/mod_kinds.o Build/mod_ocean.o Build/mod_parallel.o
ocean_coupler.o: Build/mod_param.o Build/mod_scalars.o Build/mod_sediment.o
ocean_coupler.o: Build/mod_stepping.o Build/mp_exchange.o

propagator.o: propagator_fsv.h propagator_fte.h propagator_so_semi.h
propagator.o: propagator_op.h propagator_afte.h cppdefs.h globaldefs.h
propagator.o: Build/dotproduct.o Build/ini_adjust.o Build/mod_coupling.o
propagator.o: Build/mod_iounits.o Build/mod_ocean.o Build/mod_parallel.o
propagator.o: Build/mod_param.o Build/mod_scalars.o Build/mod_stepping.o
propagator.o: Build/packing.o Build/set_depth.o

roms_export.o: set_bounds.h cppdefs.h globaldefs.h
roms_export.o: Build/distribute.o Build/mod_kinds.o Build/mod_ncparam.o
roms_export.o: Build/mod_param.o

roms_import.o: cppdefs.h globaldefs.h
roms_import.o: Build/distribute.o Build/exchange_2d.o Build/mod_kinds.o
roms_import.o: Build/mod_ncparam.o Build/mod_param.o Build/mp_exchange.o

analytical.o: ana_ssh.h set_bounds.h tile.h ana_tobc.h ana_fsobc.h ana_bmflux.h
analytical.o: ana_passive.h ana_stflux.h ana_spinning.h ana_m2clima.h
analytical.o: ana_psource.h ana_rain.h ana_cloud.h ana_specir.h ana_scope.h
analytical.o: ana_humid.h ana_initial.h ana_sst.h ana_grid.h ana_nudgcoef.h
analytical.o: ana_wwave.h ana_m3obc.h ana_btflux.h ana_hmixcoef.h ana_sss.h
analytical.o: ana_tclima.h ana_diag.h ana_mask.h ana_pair.h ana_m2obc.h
analytical.o: ana_winds.h ana_biology.h ana_m3clima.h ana_smflux.h
analytical.o: ana_perturb.h ana_vmix.h ana_tair.h cppdefs.h globaldefs.h
analytical.o: ana_srflux.h ana_sediment.h
analytical.o: Build/distribute.o Build/exchange_2d.o Build/exchange_3d.o
analytical.o: Build/mod_biology.o Build/mod_boundary.o Build/mod_clima.o
analytical.o: Build/mod_eclight.o Build/mod_forces.o Build/mod_grid.o
analytical.o: Build/mod_iounits.o Build/mod_mixing.o Build/mod_ncparam.o
analytical.o: Build/mod_ocean.o Build/mod_parallel.o Build/mod_param.o
analytical.o: Build/mod_scalars.o Build/mod_sediment.o Build/mod_sources.o
analytical.o: Build/mod_stepping.o Build/mp_exchange.o

mod_arrays.o: cppdefs.h globaldefs.h
mod_arrays.o: Build/mod_average.o Build/mod_bbl.o Build/mod_boundary.o
mod_arrays.o: Build/mod_clima.o Build/mod_coupling.o Build/mod_diags.o
mod_arrays.o: Build/mod_forces.o Build/mod_grid.o Build/mod_mixing.o
mod_arrays.o: Build/mod_obs.o Build/mod_ocean.o Build/mod_parallel.o
mod_arrays.o: Build/mod_param.o Build/mod_sources.o Build/mod_tides.o

mod_average.o: tile.h cppdefs.h globaldefs.h
mod_average.o: Build/mod_kinds.o Build/mod_param.o Build/mod_sediment.o

mod_bbl.o: tile.h cppdefs.h globaldefs.h
mod_bbl.o: Build/mod_kinds.o Build/mod_param.o

mod_biology.o: cppdefs.h globaldefs.h
mod_biology.o: Build/mod_eclight.o Build/mod_param.o Build/mod_scalars.o

mod_boundary.o: tile.h cppdefs.h globaldefs.h
mod_boundary.o: Build/mod_kinds.o Build/mod_param.o

mod_clima.o: tile.h cppdefs.h globaldefs.h
mod_clima.o: Build/mod_kinds.o Build/mod_param.o

mod_coupler.o: cppdefs.h globaldefs.h
mod_coupler.o: Build/mod_iounits.o Build/mod_ncparam.o Build/mod_parallel.o
mod_coupler.o: Build/mod_param.o Build/mod_scalars.o

mod_coupling.o: tile.h cppdefs.h globaldefs.h
mod_coupling.o: Build/mod_kinds.o Build/mod_param.o

mod_diags.o: tile.h cppdefs.h globaldefs.h
mod_diags.o: Build/mod_kinds.o Build/mod_param.o

mod_eclight.o: cppdefs.h globaldefs.h
mod_eclight.o: Build/mod_param.o

mod_eoscoef.o: cppdefs.h globaldefs.h
mod_eoscoef.o: Build/mod_kinds.o

mod_floats.o: cppdefs.h globaldefs.h
mod_floats.o: Build/mod_param.o Build/mod_scalars.o

mod_forces.o: tile.h cppdefs.h globaldefs.h
mod_forces.o: Build/mod_kinds.o Build/mod_param.o Build/mod_scalars.o

mod_fourdvar.o: cppdefs.h globaldefs.h
mod_fourdvar.o: Build/mod_iounits.o Build/mod_ncparam.o Build/mod_netcdf.o
mod_fourdvar.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o

mod_grid.o: tile.h cppdefs.h globaldefs.h
mod_grid.o: Build/mod_kinds.o Build/mod_param.o

mod_iounits.o: cppdefs.h globaldefs.h
mod_iounits.o: Build/mod_param.o


mod_mixing.o: tile.h cppdefs.h globaldefs.h
mod_mixing.o: Build/mod_kinds.o Build/mod_param.o Build/mod_scalars.o

mod_ncparam.o: cppdefs.h globaldefs.h
mod_ncparam.o: Build/mod_iounits.o Build/mod_parallel.o Build/mod_param.o
mod_ncparam.o: Build/mod_scalars.o Build/mod_sediment.o

mod_nesting.o: cppdefs.h globaldefs.h
mod_nesting.o: Build/mod_kinds.o

mod_netcdf.o: cppdefs.h globaldefs.h
mod_netcdf.o: Build/distribute.o Build/mod_iounits.o Build/mod_kinds.o
mod_netcdf.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o

mod_obs.o: tile.h cppdefs.h globaldefs.h
mod_obs.o: Build/mod_kinds.o Build/mod_param.o

mod_ocean.o: tile.h cppdefs.h globaldefs.h
mod_ocean.o: Build/mod_kinds.o Build/mod_param.o Build/mod_sediment.o

mod_parallel.o: cppdefs.h globaldefs.h
mod_parallel.o: Build/mod_iounits.o Build/mod_param.o Build/mod_scalars.o
mod_parallel.o: Build/mod_strings.o

mod_param.o: cppdefs.h globaldefs.h
mod_param.o: Build/mod_kinds.o

mod_scalars.o: cppdefs.h globaldefs.h
mod_scalars.o: Build/mod_param.o

mod_sediment.o: cppdefs.h globaldefs.h
mod_sediment.o: Build/mod_param.o

mod_sources.o: cppdefs.h globaldefs.h
mod_sources.o: Build/distribute.o Build/mod_iounits.o Build/mod_kinds.o
mod_sources.o: Build/mod_ncparam.o Build/mod_parallel.o Build/mod_param.o

mod_stepping.o: cppdefs.h globaldefs.h
mod_stepping.o: Build/mod_param.o

mod_storage.o: cppdefs.h globaldefs.h
mod_storage.o: Build/mod_param.o Build/mod_scalars.o

mod_strings.o: cppdefs.h globaldefs.h

mod_tides.o: tile.h cppdefs.h globaldefs.h
mod_tides.o: Build/distribute.o Build/mod_iounits.o Build/mod_kinds.o
mod_tides.o: Build/mod_ncparam.o Build/mod_parallel.o Build/mod_param.o
mod_tides.o: Build/mod_scalars.o

bbl.o: mb_bbl.h set_bounds.h tile.h sg_bbl.h ssw_bbl.h cppdefs.h globaldefs.h
bbl.o: Build/bc_2d.o Build/mod_bbl.o Build/mod_forces.o Build/mod_grid.o
bbl.o: Build/mod_ocean.o Build/mod_parallel.o Build/mod_param.o
bbl.o: Build/mod_scalars.o Build/mod_sediment.o Build/mod_stepping.o
bbl.o: Build/mp_exchange.o

bc_2d.o: set_bounds.h cppdefs.h globaldefs.h
bc_2d.o: Build/exchange_2d.o Build/mod_grid.o Build/mod_param.o
bc_2d.o: Build/mod_scalars.o

bc_3d.o: set_bounds.h cppdefs.h globaldefs.h
bc_3d.o: Build/exchange_3d.o Build/mod_grid.o Build/mod_param.o
bc_3d.o: Build/mod_scalars.o

bed.o: set_bounds.h tile.h cppdefs.h globaldefs.h
bed.o: Build/bc_3d.o Build/exchange_2d.o Build/mod_bbl.o Build/mod_forces.o
bed.o: Build/mod_grid.o Build/mod_ocean.o Build/mod_param.o Build/mod_scalars.o
bed.o: Build/mod_sediment.o Build/mod_stepping.o Build/mp_exchange.o

bed_cohesive.o: set_bounds.h tile.h cppdefs.h globaldefs.h
bed_cohesive.o: Build/bc_3d.o Build/exchange_2d.o Build/mod_bbl.o
bed_cohesive.o: Build/mod_forces.o Build/mod_grid.o Build/mod_ocean.o
bed_cohesive.o: Build/mod_param.o Build/mod_scalars.o Build/mod_sediment.o
bed_cohesive.o: Build/mod_stepping.o Build/mp_exchange.o

bedbiodiff.o: set_bounds.h tile.h cppdefs.h globaldefs.h
bedbiodiff.o: Build/bc_3d.o Build/exchange_2d.o Build/mod_bbl.o
bedbiodiff.o: Build/mod_forces.o Build/mod_grid.o Build/mod_ocean.o
bedbiodiff.o: Build/mod_param.o Build/mod_scalars.o Build/mod_sediment.o
bedbiodiff.o: Build/mod_stepping.o Build/mp_exchange.o

bedload.o: set_bounds.h tile.h cppdefs.h globaldefs.h
bedload.o: Build/bc_3d.o Build/exchange_2d.o Build/mod_bbl.o Build/mod_forces.o
bedload.o: Build/mod_grid.o Build/mod_ocean.o Build/mod_param.o
bedload.o: Build/mod_scalars.o Build/mod_sediment.o Build/mod_stepping.o
bedload.o: Build/mp_exchange.o

biology.o: ecosim.h set_bounds.h tile.h nemuro.h fasham.h npzd_Franks.h
biology.o: npzd_Powell.h cppdefs.h globaldefs.h
biology.o: Build/mod_biology.o Build/mod_diags.o Build/mod_forces.o
biology.o: Build/mod_grid.o Build/mod_iounits.o Build/mod_kinds.o
biology.o: Build/mod_ncparam.o Build/mod_ocean.o Build/mod_param.o
biology.o: Build/mod_scalars.o Build/mod_stepping.o

bottom.o: set_bounds.h tile.h cppdefs.h globaldefs.h
bottom.o: Build/bc_3d.o Build/exchange_2d.o Build/mod_ocean.o Build/mod_param.o
bottom.o: Build/mod_scalars.o Build/mod_sediment.o Build/mod_stepping.o
bottom.o: Build/mp_exchange.o

bulk_flux.o: set_bounds.h tile.h cppdefs.h globaldefs.h
bulk_flux.o: Build/exchange_2d.o Build/mod_forces.o Build/mod_grid.o
bulk_flux.o: Build/mod_kinds.o Build/mod_mixing.o Build/mod_ocean.o
bulk_flux.o: Build/mod_param.o Build/mod_scalars.o Build/mod_stepping.o
bulk_flux.o: Build/mp_exchange.o

bvf_mix.o: set_bounds.h tile.h cppdefs.h globaldefs.h
bvf_mix.o: Build/exchange_3d.o Build/mod_mixing.o Build/mod_param.o
bvf_mix.o: Build/mod_scalars.o Build/mp_exchange.o

conv_2d.o: set_bounds.h cppdefs.h globaldefs.h
conv_2d.o: Build/bc_2d.o Build/mod_param.o Build/mp_exchange.o

conv_3d.o: set_bounds.h cppdefs.h globaldefs.h
conv_3d.o: Build/bc_3d.o Build/mod_param.o Build/mp_exchange.o

dep_ero.o: set_bounds.h tile.h cppdefs.h globaldefs.h
dep_ero.o: Build/mod_bbl.o Build/mod_forces.o Build/mod_grid.o
dep_ero.o: Build/mod_ocean.o Build/mod_param.o Build/mod_scalars.o
dep_ero.o: Build/mod_sediment.o Build/mod_stepping.o

diag.o: set_bounds.h tile.h cppdefs.h globaldefs.h
diag.o: Build/analytical.o Build/distribute.o Build/mod_grid.o
diag.o: Build/mod_iounits.o Build/mod_ocean.o Build/mod_parallel.o
diag.o: Build/mod_param.o Build/mod_scalars.o Build/mod_stepping.o

exchange_2d.o: set_bounds.h cppdefs.h globaldefs.h
exchange_2d.o: Build/mod_param.o

exchange_3d.o: set_bounds.h cppdefs.h globaldefs.h
exchange_3d.o: Build/mod_param.o

forcing.o: set_bounds.h tile.h cppdefs.h globaldefs.h
forcing.o: Build/mod_coupling.o Build/mod_iounits.o Build/mod_ocean.o
forcing.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o

get_data.o: cppdefs.h globaldefs.h
get_data.o: Build/mod_boundary.o Build/mod_clima.o Build/mod_forces.o
get_data.o: Build/mod_grid.o Build/mod_iounits.o Build/mod_ncparam.o
get_data.o: Build/mod_obs.o Build/mod_ocean.o Build/mod_param.o
get_data.o: Build/mod_scalars.o Build/mod_sources.o Build/mod_stepping.o

get_idata.o: cppdefs.h globaldefs.h
get_idata.o: Build/distribute.o Build/mod_grid.o Build/mod_iounits.o
get_idata.o: Build/mod_ncparam.o Build/mod_netcdf.o Build/mod_parallel.o
get_idata.o: Build/mod_param.o Build/mod_scalars.o Build/mod_sources.o
get_idata.o: Build/mod_stepping.o Build/mod_tides.o

gls_corstep.o: set_bounds.h tile.h cppdefs.h globaldefs.h
gls_corstep.o: Build/exchange_3d.o Build/mod_forces.o Build/mod_grid.o
gls_corstep.o: Build/mod_mixing.o Build/mod_ocean.o Build/mod_param.o
gls_corstep.o: Build/mod_scalars.o Build/mod_stepping.o Build/mp_exchange.o
gls_corstep.o: Build/tkebc_im.o

gls_prestep.o: set_bounds.h tile.h cppdefs.h globaldefs.h
gls_prestep.o: Build/exchange_3d.o Build/mod_grid.o Build/mod_mixing.o
gls_prestep.o: Build/mod_ocean.o Build/mod_param.o Build/mod_scalars.o
gls_prestep.o: Build/mod_stepping.o Build/mp_exchange.o Build/tkebc_im.o

hmixing.o: set_bounds.h tile.h cppdefs.h globaldefs.h
hmixing.o: Build/exchange_3d.o Build/mod_grid.o Build/mod_mixing.o
hmixing.o: Build/mod_ocean.o Build/mod_param.o Build/mod_scalars.o
hmixing.o: Build/mod_stepping.o Build/mp_exchange.o

ini_fields.o: set_bounds.h tile.h cppdefs.h globaldefs.h
ini_fields.o: Build/exchange_2d.o Build/exchange_3d.o Build/mod_coupling.o
ini_fields.o: Build/mod_grid.o Build/mod_ncparam.o Build/mod_nesting.o
ini_fields.o: Build/mod_ocean.o Build/mod_param.o Build/mod_scalars.o
ini_fields.o: Build/mod_sediment.o Build/mod_stepping.o Build/mp_exchange.o
ini_fields.o: Build/nesting.o Build/set_depth.o Build/t3dbc_im.o
ini_fields.o: Build/u2dbc_im.o Build/u3dbc_im.o Build/v2dbc_im.o
ini_fields.o: Build/v3dbc_im.o Build/zetabc.o

initial.o: cppdefs.h globaldefs.h
initial.o: Build/analytical.o Build/distribute.o Build/get_2dparent_data.o
initial.o: Build/get_3dparent_data.o Build/ini_adjust.o
initial.o: Build/init_child_hindices.o Build/init_parent_hindices.o
initial.o: Build/metrics.o Build/mod_bbl.o Build/mod_fourdvar.o
initial.o: Build/mod_grid.o Build/mod_iounits.o Build/mod_ncparam.o
initial.o: Build/mod_ocean.o Build/mod_parallel.o Build/mod_param.o
initial.o: Build/mod_scalars.o Build/mod_stepping.o Build/ocean_coupler.o
initial.o: Build/omega.o Build/rho_eos.o Build/set_depth.o Build/set_massflux.o
initial.o: Build/stiffness.o Build/wpoints.o

interp_floats.o: cppdefs.h globaldefs.h
interp_floats.o: Build/mod_ncparam.o Build/mod_param.o Build/mod_scalars.o

lmd_bkpp.o: set_bounds.h tile.h cppdefs.h globaldefs.h
lmd_bkpp.o: Build/bc_2d.o Build/mod_forces.o Build/mod_grid.o
lmd_bkpp.o: Build/mod_mixing.o Build/mod_ocean.o Build/mod_param.o
lmd_bkpp.o: Build/mod_scalars.o Build/mod_stepping.o Build/mp_exchange.o
lmd_bkpp.o: Build/shapiro.o

lmd_skpp.o: set_bounds.h tile.h cppdefs.h globaldefs.h
lmd_skpp.o: Build/bc_2d.o Build/mod_forces.o Build/mod_grid.o
lmd_skpp.o: Build/mod_mixing.o Build/mod_ocean.o Build/mod_param.o
lmd_skpp.o: Build/mod_scalars.o Build/mod_stepping.o Build/mp_exchange.o
lmd_skpp.o: Build/shapiro.o

lmd_swfrac.o: set_bounds.h cppdefs.h globaldefs.h
lmd_swfrac.o: Build/mod_mixing.o Build/mod_param.o Build/mod_scalars.o

lmd_vmix.o: set_bounds.h tile.h cppdefs.h globaldefs.h
lmd_vmix.o: Build/bc_3d.o Build/lmd_bkpp.o Build/lmd_skpp.o Build/mod_grid.o
lmd_vmix.o: Build/mod_mixing.o Build/mod_ocean.o Build/mod_param.o
lmd_vmix.o: Build/mod_scalars.o Build/mod_stepping.o Build/mp_exchange.o

main2d.o: cppdefs.h globaldefs.h
main2d.o: Build/diag.o Build/dotproduct.o Build/forcing.o Build/frc_adjust.o
main2d.o: Build/get_2dparent_data.o Build/ini_fields.o Build/mod_coupler.o
main2d.o: Build/mod_iounits.o Build/mod_parallel.o Build/mod_param.o
main2d.o: Build/mod_scalars.o Build/mod_stepping.o Build/nesting.o
main2d.o: Build/ocean_coupler.o Build/oi_update.o Build/radiation_stress.o
main2d.o: Build/set_2dchild_data.o Build/set_2dparent_data.o Build/set_avg.o
main2d.o: Build/set_tides.o Build/set_vbc.o Build/step2d.o Build/step_floats.o

main3d.o: cppdefs.h globaldefs.h
main3d.o: Build/analytical.o Build/bbl.o Build/biology.o Build/bulk_flux.o
main3d.o: Build/bvf_mix.o Build/diag.o Build/dotproduct.o Build/forcing.o
main3d.o: Build/frc_adjust.o Build/get_2dparent_data.o
main3d.o: Build/get_3dparent_data.o Build/gls_corstep.o Build/gls_prestep.o
main3d.o: Build/hmixing.o Build/ini_fields.o Build/lmd_vmix.o
main3d.o: Build/mod_coupler.o Build/mod_iounits.o Build/mod_parallel.o
main3d.o: Build/mod_param.o Build/mod_scalars.o Build/mod_stepping.o
main3d.o: Build/my25_corstep.o Build/my25_prestep.o Build/nesting.o
main3d.o: Build/ocean_coupler.o Build/oi_update.o Build/omega.o
main3d.o: Build/radiation_stress.o Build/rho_eos.o Build/rhs3d.o
main3d.o: Build/sediment.o Build/set_2dchild_data.o Build/set_2dparent_data.o
main3d.o: Build/set_3dchild_data.o Build/set_3dparent_data.o Build/set_avg.o
main3d.o: Build/set_depth.o Build/set_massflux.o Build/set_tides.o
main3d.o: Build/set_vbc.o Build/set_zeta.o Build/step2d.o Build/step3d_t.o
main3d.o: Build/step3d_uv.o Build/step_floats.o Build/wvelocity.o

mpdata_adiff.o: set_bounds.h cppdefs.h globaldefs.h
mpdata_adiff.o: Build/bc_3d.o Build/mod_param.o Build/mod_scalars.o

my25_corstep.o: set_bounds.h tile.h cppdefs.h globaldefs.h
my25_corstep.o: Build/exchange_3d.o Build/mod_forces.o Build/mod_grid.o
my25_corstep.o: Build/mod_mixing.o Build/mod_ocean.o Build/mod_param.o
my25_corstep.o: Build/mod_scalars.o Build/mod_stepping.o Build/mp_exchange.o
my25_corstep.o: Build/tkebc_im.o

my25_prestep.o: set_bounds.h tile.h cppdefs.h globaldefs.h
my25_prestep.o: Build/exchange_3d.o Build/mod_grid.o Build/mod_mixing.o
my25_prestep.o: Build/mod_ocean.o Build/mod_param.o Build/mod_scalars.o
my25_prestep.o: Build/mod_stepping.o Build/mp_exchange.o Build/tkebc_im.o

obc_volcons.o: set_bounds.h tile.h cppdefs.h globaldefs.h
obc_volcons.o: Build/distribute.o Build/mod_grid.o Build/mod_ocean.o
obc_volcons.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o
obc_volcons.o: Build/mp_exchange.o

omega.o: set_bounds.h tile.h cppdefs.h globaldefs.h
omega.o: Build/bc_3d.o Build/mod_grid.o Build/mod_ocean.o Build/mod_param.o
omega.o: Build/mod_scalars.o Build/mod_sources.o Build/mod_stepping.o
omega.o: Build/mp_exchange.o

output.o: cppdefs.h globaldefs.h
output.o: Build/mod_floats.o Build/mod_fourdvar.o Build/mod_iounits.o
output.o: Build/mod_ncparam.o Build/mod_netcdf.o Build/mod_parallel.o
output.o: Build/mod_param.o Build/mod_scalars.o

pre_step3d.o: set_bounds.h tile.h cppdefs.h globaldefs.h
pre_step3d.o: Build/exchange_3d.o Build/mod_diags.o Build/mod_forces.o
pre_step3d.o: Build/mod_grid.o Build/mod_mixing.o Build/mod_ocean.o
pre_step3d.o: Build/mod_param.o Build/mod_scalars.o Build/mod_sources.o
pre_step3d.o: Build/mod_stepping.o Build/mp_exchange.o Build/t3dbc_im.o

prsgrd.o: prsgrd32.h set_bounds.h tile.h prsgrd31.h prsgrd42.h prsgrd40.h
prsgrd.o: prsgrd44.h cppdefs.h globaldefs.h
prsgrd.o: Build/mod_diags.o Build/mod_forces.o Build/mod_grid.o
prsgrd.o: Build/mod_ocean.o Build/mod_param.o Build/mod_scalars.o
prsgrd.o: Build/mod_stepping.o

radiation_stress.o: set_bounds.h tile.h cppdefs.h globaldefs.h
radiation_stress.o: Build/bc_2d.o Build/bc_3d.o Build/exchange_2d.o
radiation_stress.o: Build/exchange_3d.o Build/mod_coupling.o Build/mod_diags.o
radiation_stress.o: Build/mod_forces.o Build/mod_grid.o Build/mod_mixing.o
radiation_stress.o: Build/mod_ocean.o Build/mod_param.o Build/mod_scalars.o
radiation_stress.o: Build/mod_stepping.o Build/mp_exchange.o

rho_eos.o: set_bounds.h tile.h cppdefs.h globaldefs.h
rho_eos.o: Build/exchange_2d.o Build/exchange_3d.o Build/mod_coupling.o
rho_eos.o: Build/mod_eoscoef.o Build/mod_grid.o Build/mod_mixing.o
rho_eos.o: Build/mod_ocean.o Build/mod_param.o Build/mod_scalars.o
rho_eos.o: Build/mod_sediment.o Build/mod_stepping.o Build/mp_exchange.o

rhs3d.o: set_bounds.h tile.h cppdefs.h globaldefs.h
rhs3d.o: Build/exchange_2d.o Build/mod_clima.o Build/mod_coupling.o
rhs3d.o: Build/mod_diags.o Build/mod_forces.o Build/mod_grid.o
rhs3d.o: Build/mod_mixing.o Build/mod_obs.o Build/mod_ocean.o Build/mod_param.o
rhs3d.o: Build/mod_scalars.o Build/mod_stepping.o Build/mp_exchange.o
rhs3d.o: Build/pre_step3d.o Build/prsgrd.o Build/t3dmix.o Build/uv3dmix.o

sediment.o: cppdefs.h globaldefs.h
sediment.o: Build/bed.o Build/bed_cohesive.o Build/bedbiodiff.o Build/bedload.o
sediment.o: Build/bottom.o Build/dep_ero.o Build/settling.o

set_avg.o: set_bounds.h tile.h cppdefs.h globaldefs.h
set_avg.o: Build/mod_average.o Build/mod_coupling.o Build/mod_forces.o
set_avg.o: Build/mod_grid.o Build/mod_mixing.o Build/mod_ocean.o
set_avg.o: Build/mod_param.o Build/mod_scalars.o Build/mod_stepping.o
set_avg.o: Build/mod_tides.o

set_data.o: set_bounds.h tile.h cppdefs.h globaldefs.h
set_data.o: Build/analytical.o Build/exchange_2d.o Build/exchange_3d.o
set_data.o: Build/frc_adjust.o Build/mod_boundary.o Build/mod_clima.o
set_data.o: Build/mod_forces.o Build/mod_grid.o Build/mod_mixing.o
set_data.o: Build/mod_ncparam.o Build/mod_obs.o Build/mod_ocean.o
set_data.o: Build/mod_param.o Build/mod_scalars.o Build/mod_sources.o
set_data.o: Build/mod_stepping.o Build/mp_exchange.o Build/set_2dfld.o
set_data.o: Build/set_3dfld.o

set_depth.o: set_bounds.h tile.h cppdefs.h globaldefs.h
set_depth.o: Build/exchange_2d.o Build/exchange_3d.o Build/mod_coupling.o
set_depth.o: Build/mod_grid.o Build/mod_ocean.o Build/mod_param.o
set_depth.o: Build/mod_scalars.o Build/mod_stepping.o Build/mp_exchange.o

set_massflux.o: set_bounds.h tile.h cppdefs.h globaldefs.h
set_massflux.o: Build/exchange_3d.o Build/mod_coupling.o Build/mod_grid.o
set_massflux.o: Build/mod_ocean.o Build/mod_param.o Build/mod_scalars.o
set_massflux.o: Build/mod_stepping.o Build/mp_exchange.o

set_tides.o: set_bounds.h tile.h cppdefs.h globaldefs.h
set_tides.o: Build/distribute.o Build/exchange_2d.o Build/mod_boundary.o
set_tides.o: Build/mod_clima.o Build/mod_grid.o Build/mod_param.o
set_tides.o: Build/mod_scalars.o Build/mod_stepping.o Build/mod_tides.o
set_tides.o: Build/mp_exchange.o

set_vbc.o: set_bounds.h tile.h cppdefs.h globaldefs.h
set_vbc.o: Build/bc_2d.o Build/mod_forces.o Build/mod_grid.o Build/mod_ocean.o
set_vbc.o: Build/mod_param.o Build/mod_scalars.o Build/mod_stepping.o
set_vbc.o: Build/mp_exchange.o

set_zeta.o: set_bounds.h tile.h cppdefs.h globaldefs.h
set_zeta.o: Build/exchange_2d.o Build/mod_coupling.o Build/mod_ocean.o
set_zeta.o: Build/mod_param.o Build/mod_scalars.o Build/mp_exchange.o

settling.o: set_bounds.h tile.h cppdefs.h globaldefs.h
settling.o: Build/mod_forces.o Build/mod_grid.o Build/mod_ocean.o
settling.o: Build/mod_param.o Build/mod_scalars.o Build/mod_sediment.o
settling.o: Build/mod_stepping.o

step2d.o: step2d_LF_AM3.h set_bounds.h tile.h cppdefs.h globaldefs.h
step2d.o: Build/exchange_2d.o Build/mod_clima.o Build/mod_coupling.o
step2d.o: Build/mod_diags.o Build/mod_forces.o Build/mod_grid.o
step2d.o: Build/mod_mixing.o Build/mod_ocean.o Build/mod_param.o
step2d.o: Build/mod_scalars.o Build/mod_sediment.o Build/mod_sources.o
step2d.o: Build/mod_stepping.o Build/mp_exchange.o Build/obc_volcons.o
step2d.o: Build/set_depth.o Build/u2dbc_im.o Build/v2dbc_im.o Build/wetdry.o
step2d.o: Build/zetabc.o

step3d_t.o: set_bounds.h tile.h cppdefs.h globaldefs.h
step3d_t.o: Build/exchange_3d.o Build/mod_clima.o Build/mod_diags.o
step3d_t.o: Build/mod_grid.o Build/mod_mixing.o Build/mod_ncparam.o
step3d_t.o: Build/mod_obs.o Build/mod_ocean.o Build/mod_param.o
step3d_t.o: Build/mod_scalars.o Build/mod_sources.o Build/mod_stepping.o
step3d_t.o: Build/mp_exchange.o Build/mpdata_adiff.o Build/t3dbc_im.o

step3d_uv.o: set_bounds.h tile.h cppdefs.h globaldefs.h
step3d_uv.o: Build/exchange_2d.o Build/exchange_3d.o Build/mod_coupling.o
step3d_uv.o: Build/mod_diags.o Build/mod_forces.o Build/mod_grid.o
step3d_uv.o: Build/mod_mixing.o Build/mod_ocean.o Build/mod_param.o
step3d_uv.o: Build/mod_scalars.o Build/mod_sources.o Build/mod_stepping.o
step3d_uv.o: Build/mp_exchange.o Build/u3dbc_im.o Build/v3dbc_im.o

step_floats.o: cppdefs.h globaldefs.h
step_floats.o: Build/distribute.o Build/mod_floats.o Build/mod_grid.o
step_floats.o: Build/mod_iounits.o Build/mod_mixing.o Build/mod_ncparam.o
step_floats.o: Build/mod_ocean.o Build/mod_parallel.o Build/mod_param.o
step_floats.o: Build/mod_scalars.o Build/mod_stepping.o Build/utility.o

t3dbc_im.o: set_bounds.h tile.h cppdefs.h globaldefs.h
t3dbc_im.o: Build/mod_boundary.o Build/mod_grid.o Build/mod_ocean.o
t3dbc_im.o: Build/mod_param.o Build/mod_scalars.o Build/mod_stepping.o

t3dmix.o: t3dmix4_geo.h set_bounds.h tile.h t3dmix2_geo.h t3dmix2_s.h
t3dmix.o: t3dmix4_iso.h t3dmix4_s.h t3dmix2_iso.h cppdefs.h globaldefs.h
t3dmix.o: Build/mod_diags.o Build/mod_grid.o Build/mod_mixing.o
t3dmix.o: Build/mod_ocean.o Build/mod_param.o Build/mod_scalars.o
t3dmix.o: Build/mod_stepping.o

tkebc_im.o: set_bounds.h tile.h cppdefs.h globaldefs.h
tkebc_im.o: Build/mod_grid.o Build/mod_mixing.o Build/mod_param.o
tkebc_im.o: Build/mod_scalars.o Build/mod_stepping.o

u2dbc_im.o: set_bounds.h tile.h cppdefs.h globaldefs.h
u2dbc_im.o: Build/mod_boundary.o Build/mod_forces.o Build/mod_grid.o
u2dbc_im.o: Build/mod_ocean.o Build/mod_param.o Build/mod_scalars.o
u2dbc_im.o: Build/mod_stepping.o

u3dbc_im.o: set_bounds.h tile.h cppdefs.h globaldefs.h
u3dbc_im.o: Build/mod_boundary.o Build/mod_grid.o Build/mod_ocean.o
u3dbc_im.o: Build/mod_param.o Build/mod_scalars.o Build/mod_stepping.o

uv3dmix.o: uv3dmix2_s.h set_bounds.h tile.h uv3dmix2_geo.h uv3dmix4_s.h
uv3dmix.o: uv3dmix4_geo.h cppdefs.h globaldefs.h
uv3dmix.o: Build/mod_coupling.o Build/mod_diags.o Build/mod_grid.o
uv3dmix.o: Build/mod_mixing.o Build/mod_ocean.o Build/mod_param.o
uv3dmix.o: Build/mod_scalars.o Build/mod_stepping.o

v2dbc_im.o: set_bounds.h tile.h cppdefs.h globaldefs.h
v2dbc_im.o: Build/mod_boundary.o Build/mod_forces.o Build/mod_grid.o
v2dbc_im.o: Build/mod_ocean.o Build/mod_param.o Build/mod_scalars.o
v2dbc_im.o: Build/mod_stepping.o

v3dbc_im.o: set_bounds.h tile.h cppdefs.h globaldefs.h
v3dbc_im.o: Build/mod_boundary.o Build/mod_grid.o Build/mod_ocean.o
v3dbc_im.o: Build/mod_param.o Build/mod_scalars.o Build/mod_stepping.o

wetdry.o: set_bounds.h cppdefs.h globaldefs.h
wetdry.o: Build/exchange_2d.o Build/mod_param.o Build/mod_scalars.o
wetdry.o: Build/mp_exchange.o

wvelocity.o: set_bounds.h tile.h cppdefs.h globaldefs.h
wvelocity.o: Build/bc_3d.o Build/exchange_2d.o Build/mod_coupling.o
wvelocity.o: Build/mod_grid.o Build/mod_ocean.o Build/mod_param.o
wvelocity.o: Build/mod_scalars.o Build/mod_stepping.o Build/mp_exchange.o

zetabc.o: set_bounds.h tile.h cppdefs.h globaldefs.h
zetabc.o: Build/mod_boundary.o Build/mod_grid.o Build/mod_ocean.o
zetabc.o: Build/mod_param.o Build/mod_scalars.o Build/mod_stepping.o

abort.o: cppdefs.h globaldefs.h
abort.o: Build/ocean_control.o

back_cost.o: set_bounds.h tile.h cppdefs.h globaldefs.h
back_cost.o: Build/distribute.o Build/mod_forces.o Build/mod_fourdvar.o
back_cost.o: Build/mod_ncparam.o Build/mod_ocean.o Build/mod_parallel.o
back_cost.o: Build/mod_param.o Build/mod_scalars.o

back_cov.o: set_bounds.h tile.h cppdefs.h globaldefs.h
back_cov.o: Build/mod_coupling.o Build/mod_fourdvar.o Build/mod_grid.o
back_cov.o: Build/mod_iounits.o Build/mod_ncparam.o Build/mod_ocean.o
back_cov.o: Build/mod_param.o Build/mod_scalars.o Build/mod_stepping.o
back_cov.o: Build/set_depth.o

back_step.o: set_bounds.h tile.h cppdefs.h globaldefs.h
back_step.o: Build/distribute.o Build/mod_fourdvar.o Build/mod_ocean.o
back_step.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o

cgradient.o: cgradient_lanczos.h set_bounds.h tile.h cgradient.h cppdefs.h
cgradient.o: globaldefs.h
cgradient.o: Build/distribute.o Build/mod_coupling.o Build/mod_forces.o
cgradient.o: Build/mod_fourdvar.o Build/mod_grid.o Build/mod_iounits.o
cgradient.o: Build/mod_ncparam.o Build/mod_netcdf.o Build/mod_ocean.o
cgradient.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o
cgradient.o: Build/mod_stepping.o Build/state_addition.o Build/state_copy.o
cgradient.o: Build/state_dotprod.o Build/state_initialize.o Build/state_scale.o

checkdefs.o: cppdefs.h globaldefs.h
checkdefs.o: Build/mod_iounits.o Build/mod_parallel.o Build/mod_param.o
checkdefs.o: Build/mod_scalars.o Build/mod_strings.o

checkerror.o: cppdefs.h globaldefs.h
checkerror.o: Build/mod_iounits.o Build/mod_parallel.o Build/mod_param.o
checkerror.o: Build/mod_scalars.o

checkvars.o: cppdefs.h globaldefs.h
checkvars.o: Build/mod_iounits.o Build/mod_ncparam.o Build/mod_netcdf.o
checkvars.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o
checkvars.o: Build/mod_sediment.o

close_io.o: cppdefs.h globaldefs.h
close_io.o: Build/mod_iounits.o Build/mod_ncparam.o Build/mod_netcdf.o
close_io.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o

congrad.o: cppdefs.h globaldefs.h
congrad.o: Build/distribute.o Build/mod_fourdvar.o Build/mod_iounits.o
congrad.o: Build/mod_ncparam.o Build/mod_netcdf.o Build/mod_parallel.o
congrad.o: Build/mod_param.o Build/mod_scalars.o

cost_grad.o: set_bounds.h tile.h cppdefs.h globaldefs.h
cost_grad.o: Build/mod_forces.o Build/mod_ocean.o Build/mod_param.o
cost_grad.o: Build/mod_scalars.o

cost_norm.o: set_bounds.h tile.h cppdefs.h globaldefs.h
cost_norm.o: Build/distribute.o Build/mod_forces.o Build/mod_fourdvar.o
cost_norm.o: Build/mod_ncparam.o Build/mod_ocean.o Build/mod_parallel.o
cost_norm.o: Build/mod_param.o Build/mod_scalars.o

def_avg.o: cppdefs.h globaldefs.h
def_avg.o: Build/distribute.o Build/mod_fourdvar.o Build/mod_iounits.o
def_avg.o: Build/mod_ncparam.o Build/mod_netcdf.o Build/mod_parallel.o
def_avg.o: Build/mod_param.o Build/mod_scalars.o Build/mod_sediment.o

def_diags.o: cppdefs.h globaldefs.h
def_diags.o: Build/distribute.o Build/mod_fourdvar.o Build/mod_iounits.o
def_diags.o: Build/mod_ncparam.o Build/mod_netcdf.o Build/mod_parallel.o
def_diags.o: Build/mod_param.o Build/mod_scalars.o Build/mod_sediment.o

def_dim.o: cppdefs.h globaldefs.h
def_dim.o: Build/distribute.o Build/mod_iounits.o Build/mod_netcdf.o
def_dim.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o

def_floats.o: cppdefs.h globaldefs.h
def_floats.o: Build/distribute.o Build/mod_floats.o Build/mod_fourdvar.o
def_floats.o: Build/mod_grid.o Build/mod_iounits.o Build/mod_ncparam.o
def_floats.o: Build/mod_netcdf.o Build/mod_parallel.o Build/mod_param.o
def_floats.o: Build/mod_scalars.o Build/mod_sediment.o

def_gst.o: cppdefs.h globaldefs.h
def_gst.o: Build/mod_iounits.o Build/mod_ncparam.o Build/mod_netcdf.o
def_gst.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o
def_gst.o: Build/mod_storage.o

def_hessian.o: cppdefs.h globaldefs.h
def_hessian.o: Build/distribute.o Build/mod_fourdvar.o Build/mod_iounits.o
def_hessian.o: Build/mod_ncparam.o Build/mod_netcdf.o Build/mod_parallel.o
def_hessian.o: Build/mod_param.o Build/mod_scalars.o Build/mod_sediment.o

def_his.o: cppdefs.h globaldefs.h
def_his.o: Build/distribute.o Build/mod_fourdvar.o Build/mod_iounits.o
def_his.o: Build/mod_ncparam.o Build/mod_netcdf.o Build/mod_parallel.o
def_his.o: Build/mod_param.o Build/mod_scalars.o Build/mod_sediment.o

def_impulse.o: cppdefs.h globaldefs.h
def_impulse.o: Build/distribute.o Build/mod_fourdvar.o Build/mod_iounits.o
def_impulse.o: Build/mod_ncparam.o Build/mod_netcdf.o Build/mod_parallel.o
def_impulse.o: Build/mod_param.o Build/mod_scalars.o Build/mod_sediment.o

def_info.o: cppdefs.h globaldefs.h
def_info.o: Build/mod_grid.o Build/mod_iounits.o Build/mod_ncparam.o
def_info.o: Build/mod_netcdf.o Build/mod_param.o Build/mod_scalars.o
def_info.o: Build/mod_strings.o

def_ini.o: cppdefs.h globaldefs.h
def_ini.o: Build/distribute.o Build/mod_iounits.o Build/mod_ncparam.o
def_ini.o: Build/mod_netcdf.o Build/mod_parallel.o Build/mod_param.o
def_ini.o: Build/mod_scalars.o

def_lanczos.o: cppdefs.h globaldefs.h
def_lanczos.o: Build/distribute.o Build/mod_fourdvar.o Build/mod_iounits.o
def_lanczos.o: Build/mod_ncparam.o Build/mod_netcdf.o Build/mod_parallel.o
def_lanczos.o: Build/mod_param.o Build/mod_scalars.o Build/mod_sediment.o

def_mod.o: cppdefs.h globaldefs.h
def_mod.o: Build/mod_fourdvar.o Build/mod_iounits.o Build/mod_ncparam.o
def_mod.o: Build/mod_netcdf.o Build/mod_parallel.o Build/mod_param.o
def_mod.o: Build/mod_scalars.o Build/mod_strings.o

def_norm.o: cppdefs.h globaldefs.h
def_norm.o: Build/distribute.o Build/mod_fourdvar.o Build/mod_iounits.o
def_norm.o: Build/mod_ncparam.o Build/mod_netcdf.o Build/mod_parallel.o
def_norm.o: Build/mod_param.o Build/mod_scalars.o Build/mod_sediment.o

def_rst.o: cppdefs.h globaldefs.h
def_rst.o: Build/distribute.o Build/mod_fourdvar.o Build/mod_iounits.o
def_rst.o: Build/mod_ncparam.o Build/mod_netcdf.o Build/mod_parallel.o
def_rst.o: Build/mod_param.o Build/mod_scalars.o Build/mod_sediment.o

def_station.o: cppdefs.h globaldefs.h
def_station.o: Build/distribute.o Build/mod_fourdvar.o Build/mod_iounits.o
def_station.o: Build/mod_ncparam.o Build/mod_netcdf.o Build/mod_parallel.o
def_station.o: Build/mod_param.o Build/mod_scalars.o Build/mod_sediment.o

def_tides.o: cppdefs.h globaldefs.h
def_tides.o: Build/distribute.o Build/mod_iounits.o Build/mod_ncparam.o
def_tides.o: Build/mod_netcdf.o Build/mod_parallel.o Build/mod_param.o
def_tides.o: Build/mod_scalars.o Build/mod_stepping.o Build/mod_tides.o

def_var.o: cppdefs.h globaldefs.h
def_var.o: Build/mod_iounits.o Build/mod_ncparam.o Build/mod_netcdf.o
def_var.o: Build/mod_param.o Build/mod_scalars.o

descent.o: set_bounds.h tile.h cppdefs.h globaldefs.h
descent.o: Build/distribute.o Build/mod_coupling.o Build/mod_fourdvar.o
descent.o: Build/mod_grid.o Build/mod_ocean.o Build/mod_parallel.o
descent.o: Build/mod_param.o Build/mod_scalars.o Build/mod_stepping.o

distribute.o: set_bounds.h cppdefs.h globaldefs.h
distribute.o: Build/mod_iounits.o Build/mod_ncparam.o Build/mod_parallel.o
distribute.o: Build/mod_param.o Build/mod_scalars.o

dotproduct.o: set_bounds.h tile.h cppdefs.h globaldefs.h
dotproduct.o: Build/distribute.o Build/mod_fourdvar.o Build/mod_grid.o
dotproduct.o: Build/mod_iounits.o Build/mod_ncparam.o Build/mod_ocean.o
dotproduct.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o
dotproduct.o: Build/mod_stepping.o

downhill.o: set_bounds.h tile.h cppdefs.h globaldefs.h
downhill.o: Build/distribute.o Build/mod_coupling.o Build/mod_fourdvar.o
downhill.o: Build/mod_grid.o Build/mod_ocean.o Build/mod_parallel.o
downhill.o: Build/mod_param.o Build/mod_scalars.o Build/mod_stepping.o

extract_obs.o: cppdefs.h globaldefs.h
extract_obs.o: Build/mod_kinds.o Build/mod_param.o

extract_sta.o: cppdefs.h globaldefs.h
extract_sta.o: Build/distribute.o Build/mod_grid.o Build/mod_ncparam.o
extract_sta.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o

find_string.o: cppdefs.h globaldefs.h

frc_adjust.o: set_bounds.h tile.h cppdefs.h globaldefs.h
frc_adjust.o: Build/mod_forces.o Build/mod_grid.o Build/mod_param.o
frc_adjust.o: Build/mod_scalars.o

gasdev.o: cppdefs.h globaldefs.h
gasdev.o: Build/mod_kinds.o Build/nrutil.o

get_2dfld.o: cppdefs.h globaldefs.h
get_2dfld.o: Build/distribute.o Build/mod_iounits.o Build/mod_ncparam.o
get_2dfld.o: Build/mod_netcdf.o Build/mod_parallel.o Build/mod_param.o
get_2dfld.o: Build/mod_scalars.o

get_2dfldr.o: cppdefs.h globaldefs.h
get_2dfldr.o: Build/distribute.o Build/mod_iounits.o Build/mod_ncparam.o
get_2dfldr.o: Build/mod_netcdf.o Build/mod_parallel.o Build/mod_param.o
get_2dfldr.o: Build/mod_scalars.o

get_2dparent_data.o: set_bounds.h tile.h cppdefs.h globaldefs.h
get_2dparent_data.o: Build/distribute.o Build/interpolate.o
get_2dparent_data.o: Build/mod_boundary.o Build/mod_coupling.o Build/mod_grid.o
get_2dparent_data.o: Build/mod_iounits.o Build/mod_ncparam.o Build/mod_ocean.o
get_2dparent_data.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o
get_2dparent_data.o: Build/mod_stepping.o

get_3dfld.o: cppdefs.h globaldefs.h
get_3dfld.o: Build/distribute.o Build/mod_iounits.o Build/mod_ncparam.o
get_3dfld.o: Build/mod_netcdf.o Build/mod_parallel.o Build/mod_param.o
get_3dfld.o: Build/mod_scalars.o

get_3dfldr.o: cppdefs.h globaldefs.h
get_3dfldr.o: Build/distribute.o Build/mod_iounits.o Build/mod_ncparam.o
get_3dfldr.o: Build/mod_netcdf.o Build/mod_parallel.o Build/mod_param.o
get_3dfldr.o: Build/mod_scalars.o

get_3dparent_data.o: set_bounds.h tile.h cppdefs.h globaldefs.h
get_3dparent_data.o: Build/distribute.o Build/interpolate.o
get_3dparent_data.o: Build/mod_boundary.o Build/mod_grid.o Build/mod_iounits.o
get_3dparent_data.o: Build/mod_ncparam.o Build/mod_ocean.o Build/mod_parallel.o
get_3dparent_data.o: Build/mod_param.o Build/mod_scalars.o Build/mod_stepping.o

get_bounds.o: cppdefs.h globaldefs.h
get_bounds.o: Build/mod_ncparam.o Build/mod_parallel.o Build/mod_param.o

get_cycle.o: cppdefs.h globaldefs.h
get_cycle.o: Build/mod_iounits.o Build/mod_netcdf.o Build/mod_param.o
get_cycle.o: Build/mod_scalars.o

get_date.o: cppdefs.h globaldefs.h
get_date.o: Build/mod_kinds.o

get_grid.o: cppdefs.h globaldefs.h
get_grid.o: Build/distribute.o Build/exchange_2d.o Build/mod_grid.o
get_grid.o: Build/mod_iounits.o Build/mod_ncparam.o Build/mod_netcdf.o
get_grid.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o
get_grid.o: Build/mp_exchange.o

get_gridcoords.o: cppdefs.h globaldefs.h
get_gridcoords.o: Build/mod_grid.o Build/mod_iounits.o Build/mod_ncparam.o
get_gridcoords.o: Build/mod_netcdf.o Build/mod_parallel.o Build/mod_param.o
get_gridcoords.o: Build/mod_scalars.o

get_gst.o: cppdefs.h globaldefs.h
get_gst.o: Build/distribute.o Build/mod_iounits.o Build/mod_ncparam.o
get_gst.o: Build/mod_netcdf.o Build/mod_parallel.o Build/mod_param.o
get_gst.o: Build/mod_scalars.o Build/mod_storage.o

get_ngfld.o: cppdefs.h globaldefs.h
get_ngfld.o: Build/distribute.o Build/mod_iounits.o Build/mod_ncparam.o
get_ngfld.o: Build/mod_netcdf.o Build/mod_parallel.o Build/mod_param.o
get_ngfld.o: Build/mod_scalars.o

get_ngfldr.o: cppdefs.h globaldefs.h
get_ngfldr.o: Build/distribute.o Build/mod_iounits.o Build/mod_ncparam.o
get_ngfldr.o: Build/mod_netcdf.o Build/mod_parallel.o Build/mod_param.o
get_ngfldr.o: Build/mod_scalars.o

get_sparse_matrix.o: cppdefs.h globaldefs.h
get_sparse_matrix.o: Build/mod_coupler.o Build/mod_iounits.o
get_sparse_matrix.o: Build/mod_ncparam.o Build/mod_netcdf.o Build/mod_scalars.o

get_state.o: cppdefs.h globaldefs.h
get_state.o: Build/distribute.o Build/mod_forces.o Build/mod_fourdvar.o
get_state.o: Build/mod_grid.o Build/mod_iounits.o Build/mod_mixing.o
get_state.o: Build/mod_ncparam.o Build/mod_netcdf.o Build/mod_ocean.o
get_state.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o
get_state.o: Build/mod_sediment.o Build/mod_stepping.o Build/mod_strings.o

get_varcoords.o: cppdefs.h globaldefs.h
get_varcoords.o: Build/mod_grid.o Build/mod_iounits.o Build/mod_netcdf.o
get_varcoords.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o

grid_coords.o: cppdefs.h globaldefs.h
grid_coords.o: Build/distribute.o Build/interpolate.o Build/mod_floats.o
grid_coords.o: Build/mod_grid.o Build/mod_parallel.o Build/mod_param.o
grid_coords.o: Build/mod_scalars.o

impulse.o: set_bounds.h cppdefs.h globaldefs.h
impulse.o: Build/distribute.o Build/mod_grid.o Build/mod_iounits.o
impulse.o: Build/mod_ncparam.o Build/mod_netcdf.o Build/mod_ocean.o
impulse.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o

ini_adjust.o: set_bounds.h tile.h cppdefs.h globaldefs.h
ini_adjust.o: Build/distribute.o Build/exchange_2d.o Build/exchange_3d.o
ini_adjust.o: Build/mod_coupling.o Build/mod_forces.o Build/mod_fourdvar.o
ini_adjust.o: Build/mod_grid.o Build/mod_iounits.o Build/mod_ncparam.o
ini_adjust.o: Build/mod_ocean.o Build/mod_parallel.o Build/mod_param.o
ini_adjust.o: Build/mod_scalars.o Build/mod_stepping.o Build/mp_exchange.o
ini_adjust.o: Build/set_depth.o Build/t3dbc_im.o Build/u2dbc_im.o
ini_adjust.o: Build/u3dbc_im.o Build/v2dbc_im.o Build/v3dbc_im.o Build/zetabc.o

ini_lanczos.o: set_bounds.h tile.h cppdefs.h globaldefs.h
ini_lanczos.o: Build/distribute.o Build/mod_fourdvar.o Build/mod_grid.o
ini_lanczos.o: Build/mod_iounits.o Build/mod_ncparam.o Build/mod_netcdf.o
ini_lanczos.o: Build/mod_ocean.o Build/mod_parallel.o Build/mod_param.o
ini_lanczos.o: Build/mod_scalars.o Build/state_addition.o Build/state_dotprod.o
ini_lanczos.o: Build/state_initialize.o Build/state_scale.o

init_child_hindices.o: set_bounds.h tile.h cppdefs.h globaldefs.h
init_child_hindices.o: Build/interpolate.o Build/mod_boundary.o
init_child_hindices.o: Build/mod_grid.o Build/mod_iounits.o Build/mod_ncparam.o
init_child_hindices.o: Build/mod_ocean.o Build/mod_parallel.o Build/mod_param.o
init_child_hindices.o: Build/mod_scalars.o Build/mod_stepping.o

init_parent_hindices.o: set_bounds.h tile.h cppdefs.h globaldefs.h
init_parent_hindices.o: Build/interpolate.o Build/mod_boundary.o
init_parent_hindices.o: Build/mod_grid.o Build/mod_iounits.o
init_parent_hindices.o: Build/mod_ncparam.o Build/mod_ocean.o
init_parent_hindices.o: Build/mod_parallel.o Build/mod_param.o
init_parent_hindices.o: Build/mod_scalars.o Build/mod_stepping.o

inp_par.o: cppdefs.h globaldefs.h
inp_par.o: Build/distribute.o Build/mod_biology.o Build/mod_coupler.o
inp_par.o: Build/mod_floats.o Build/mod_fourdvar.o Build/mod_iounits.o
inp_par.o: Build/mod_kinds.o Build/mod_ncparam.o Build/mod_parallel.o
inp_par.o: Build/mod_param.o Build/mod_scalars.o Build/mod_sediment.o
inp_par.o: Build/mod_storage.o Build/mod_strings.o Build/ran_state.o

interpolate.o: cppdefs.h globaldefs.h
interpolate.o: Build/mod_grid.o Build/mod_kinds.o Build/mod_param.o
interpolate.o: Build/mod_scalars.o

lubksb.o: cppdefs.h globaldefs.h
lubksb.o: Build/mod_kinds.o

ludcmp.o: cppdefs.h globaldefs.h
ludcmp.o: Build/mod_kinds.o

metrics.o: set_bounds.h tile.h cppdefs.h globaldefs.h
metrics.o: Build/distribute.o Build/exchange_2d.o Build/mod_fourdvar.o
metrics.o: Build/mod_grid.o Build/mod_iounits.o Build/mod_mixing.o
metrics.o: Build/mod_ncparam.o Build/mod_ocean.o Build/mod_parallel.o
metrics.o: Build/mod_param.o Build/mod_scalars.o Build/mod_stepping.o
metrics.o: Build/mp_exchange.o Build/set_depth.o

mp_exchange.o: set_bounds.h cppdefs.h globaldefs.h
mp_exchange.o: Build/mod_iounits.o Build/mod_parallel.o Build/mod_param.o
mp_exchange.o: Build/mod_scalars.o

mp_routines.o: cppdefs.h globaldefs.h
mp_routines.o: Build/mod_kinds.o

nesting.o: cppdefs.h globaldefs.h
nesting.o: Build/interpolate.o Build/mod_coupling.o Build/mod_forces.o
nesting.o: Build/mod_grid.o Build/mod_mixing.o Build/mod_ncparam.o
nesting.o: Build/mod_nesting.o Build/mod_ocean.o Build/mod_param.o
nesting.o: Build/mod_scalars.o Build/mod_stepping.o

nf_fread2d.o: cppdefs.h globaldefs.h
nf_fread2d.o: Build/distribute.o Build/mod_grid.o Build/mod_ncparam.o
nf_fread2d.o: Build/mod_netcdf.o Build/mod_parallel.o Build/mod_param.o
nf_fread2d.o: Build/mod_scalars.o

nf_fread3d.o: cppdefs.h globaldefs.h
nf_fread3d.o: Build/distribute.o Build/mod_ncparam.o Build/mod_netcdf.o
nf_fread3d.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o

nf_fread4d.o: cppdefs.h globaldefs.h
nf_fread4d.o: Build/distribute.o Build/mod_ncparam.o Build/mod_netcdf.o
nf_fread4d.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o

nf_fwrite2d.o: cppdefs.h globaldefs.h
nf_fwrite2d.o: Build/distribute.o Build/mod_ncparam.o Build/mod_netcdf.o
nf_fwrite2d.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o

nf_fwrite3d.o: cppdefs.h globaldefs.h
nf_fwrite3d.o: Build/distribute.o Build/mod_ncparam.o Build/mod_netcdf.o
nf_fwrite3d.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o

nf_fwrite4d.o: cppdefs.h globaldefs.h
nf_fwrite4d.o: Build/distribute.o Build/mod_ncparam.o Build/mod_netcdf.o
nf_fwrite4d.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o

normalization.o: set_bounds.h tile.h cppdefs.h globaldefs.h
normalization.o: Build/bc_2d.o Build/bc_3d.o Build/distribute.o
normalization.o: Build/mod_forces.o Build/mod_fourdvar.o Build/mod_grid.o
normalization.o: Build/mod_iounits.o Build/mod_kinds.o Build/mod_mixing.o
normalization.o: Build/mod_ncparam.o Build/mod_netcdf.o Build/mod_ocean.o
normalization.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o
normalization.o: Build/mod_stepping.o Build/mp_exchange.o Build/set_depth.o
normalization.o: Build/white_noise.o

nrutil.o: Build/mod_kinds.o

obs_cost.o: cppdefs.h globaldefs.h
obs_cost.o: Build/mod_fourdvar.o Build/mod_parallel.o Build/mod_param.o
obs_cost.o: Build/mod_scalars.o

obs_depth.o: set_bounds.h cppdefs.h globaldefs.h
obs_depth.o: Build/mod_fourdvar.o Build/mod_ncparam.o Build/mod_param.o
obs_depth.o: Build/mod_scalars.o

obs_initial.o: cppdefs.h globaldefs.h
obs_initial.o: Build/distribute.o Build/mod_fourdvar.o Build/mod_iounits.o
obs_initial.o: Build/mod_ncparam.o Build/mod_netcdf.o Build/mod_parallel.o
obs_initial.o: Build/mod_param.o Build/mod_scalars.o

obs_read.o: cppdefs.h globaldefs.h
obs_read.o: Build/distribute.o Build/mod_fourdvar.o Build/mod_iounits.o
obs_read.o: Build/mod_ncparam.o Build/mod_netcdf.o Build/mod_parallel.o
obs_read.o: Build/mod_param.o Build/mod_scalars.o

obs_scale.o: cppdefs.h globaldefs.h
obs_scale.o: Build/distribute.o Build/mod_fourdvar.o Build/mod_grid.o
obs_scale.o: Build/mod_ncparam.o Build/mod_param.o Build/mod_scalars.o

obs_write.o: set_bounds.h cppdefs.h globaldefs.h
obs_write.o: Build/distribute.o Build/extract_obs.o Build/mod_fourdvar.o
obs_write.o: Build/mod_grid.o Build/mod_iounits.o Build/mod_ncparam.o
obs_write.o: Build/mod_netcdf.o Build/mod_ocean.o Build/mod_parallel.o
obs_write.o: Build/mod_param.o Build/mod_scalars.o Build/mod_stepping.o

oi_update.o: set_bounds.h tile.h cppdefs.h globaldefs.h
oi_update.o: Build/mod_grid.o Build/mod_iounits.o Build/mod_obs.o
oi_update.o: Build/mod_ocean.o Build/mod_parallel.o Build/mod_param.o
oi_update.o: Build/mod_scalars.o Build/mod_stepping.o

opencdf.o: cppdefs.h globaldefs.h
opencdf.o: Build/mod_iounits.o Build/mod_ncparam.o Build/mod_netcdf.o
opencdf.o: Build/mod_param.o Build/mod_scalars.o

packing.o: set_bounds.h tile.h cppdefs.h globaldefs.h
packing.o: Build/distribute.o Build/mod_forces.o Build/mod_grid.o
packing.o: Build/mod_iounits.o Build/mod_ncparam.o Build/mod_ocean.o
packing.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o
packing.o: Build/mod_stepping.o Build/mod_storage.o

ran1.o: cppdefs.h globaldefs.h
ran1.o: Build/mod_kinds.o Build/ran_state.o

ran_state.o: Build/mod_kinds.o Build/nrutil.o

regrid.o: cppdefs.h globaldefs.h
regrid.o: Build/distribute.o Build/interpolate.o Build/mod_iounits.o
regrid.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o

set_2dchild_data.o: set_bounds.h tile.h cppdefs.h globaldefs.h
set_2dchild_data.o: Build/distribute.o Build/exchange_2d.o Build/interpolate.o
set_2dchild_data.o: Build/mod_boundary.o Build/mod_coupling.o Build/mod_grid.o
set_2dchild_data.o: Build/mod_iounits.o Build/mod_ncparam.o Build/mod_ocean.o
set_2dchild_data.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o
set_2dchild_data.o: Build/mod_stepping.o Build/mp_exchange.o

set_2dfld.o: set_bounds.h cppdefs.h globaldefs.h
set_2dfld.o: Build/exchange_2d.o Build/mod_iounits.o Build/mod_ncparam.o
set_2dfld.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o
set_2dfld.o: Build/mp_exchange.o

set_2dfldr.o: set_bounds.h cppdefs.h globaldefs.h
set_2dfldr.o: Build/exchange_2d.o Build/mod_iounits.o Build/mod_ncparam.o
set_2dfldr.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o
set_2dfldr.o: Build/mp_exchange.o

set_2dparent_data.o: set_bounds.h tile.h cppdefs.h globaldefs.h
set_2dparent_data.o: Build/distribute.o Build/exchange_2d.o Build/interpolate.o
set_2dparent_data.o: Build/mod_boundary.o Build/mod_coupling.o Build/mod_grid.o
set_2dparent_data.o: Build/mod_iounits.o Build/mod_ncparam.o Build/mod_ocean.o
set_2dparent_data.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o
set_2dparent_data.o: Build/mod_stepping.o Build/mp_exchange.o

set_3dchild_data.o: set_bounds.h tile.h cppdefs.h globaldefs.h
set_3dchild_data.o: Build/distribute.o Build/exchange_2d.o Build/interpolate.o
set_3dchild_data.o: Build/mod_boundary.o Build/mod_coupling.o Build/mod_grid.o
set_3dchild_data.o: Build/mod_iounits.o Build/mod_ncparam.o Build/mod_ocean.o
set_3dchild_data.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o
set_3dchild_data.o: Build/mod_stepping.o Build/mp_exchange.o

set_3dfld.o: set_bounds.h cppdefs.h globaldefs.h
set_3dfld.o: Build/exchange_3d.o Build/mod_iounits.o Build/mod_ncparam.o
set_3dfld.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o
set_3dfld.o: Build/mp_exchange.o

set_3dfldr.o: set_bounds.h cppdefs.h globaldefs.h
set_3dfldr.o: Build/exchange_3d.o Build/mod_iounits.o Build/mod_ncparam.o
set_3dfldr.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o
set_3dfldr.o: Build/mp_exchange.o

set_3dparent_data.o: set_bounds.h tile.h cppdefs.h globaldefs.h
set_3dparent_data.o: Build/distribute.o Build/interpolate.o
set_3dparent_data.o: Build/mod_boundary.o Build/mod_grid.o Build/mod_iounits.o
set_3dparent_data.o: Build/mod_ncparam.o Build/mod_ocean.o Build/mod_parallel.o
set_3dparent_data.o: Build/mod_param.o Build/mod_scalars.o Build/mod_stepping.o
set_3dparent_data.o: Build/mp_exchange.o

set_diags.o: set_bounds.h tile.h cppdefs.h globaldefs.h
set_diags.o: Build/bc_2d.o Build/bc_3d.o Build/mod_diags.o Build/mod_grid.o
set_diags.o: Build/mod_param.o Build/mod_scalars.o Build/mod_stepping.o
set_diags.o: Build/mp_exchange.o

set_ngfld.o: cppdefs.h globaldefs.h
set_ngfld.o: Build/mod_iounits.o Build/mod_ncparam.o Build/mod_parallel.o
set_ngfld.o: Build/mod_param.o Build/mod_scalars.o

set_ngfldr.o: cppdefs.h globaldefs.h
set_ngfldr.o: Build/mod_iounits.o Build/mod_ncparam.o Build/mod_parallel.o
set_ngfldr.o: Build/mod_param.o Build/mod_scalars.o

set_scoord.o: cppdefs.h globaldefs.h
set_scoord.o: Build/mod_fourdvar.o Build/mod_grid.o Build/mod_iounits.o
set_scoord.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o

set_weights.o: cppdefs.h globaldefs.h
set_weights.o: Build/mod_iounits.o Build/mod_parallel.o Build/mod_param.o
set_weights.o: Build/mod_scalars.o

shapiro.o: set_bounds.h cppdefs.h globaldefs.h
shapiro.o: Build/mod_param.o

state_addition.o: set_bounds.h cppdefs.h globaldefs.h
state_addition.o: Build/mod_param.o Build/mod_scalars.o

state_copy.o: set_bounds.h cppdefs.h globaldefs.h
state_copy.o: Build/mod_param.o Build/mod_scalars.o

state_dotprod.o: set_bounds.h cppdefs.h globaldefs.h
state_dotprod.o: Build/distribute.o Build/mod_ncparam.o Build/mod_parallel.o
state_dotprod.o: Build/mod_param.o Build/mod_scalars.o

state_initialize.o: set_bounds.h cppdefs.h globaldefs.h
state_initialize.o: Build/mod_param.o Build/mod_scalars.o

state_scale.o: set_bounds.h cppdefs.h globaldefs.h
state_scale.o: Build/mod_param.o Build/mod_scalars.o

stats_modobs.o: cppdefs.h globaldefs.h
stats_modobs.o: Build/mod_fourdvar.o Build/mod_iounits.o Build/mod_ncparam.o
stats_modobs.o: Build/mod_netcdf.o Build/mod_parallel.o Build/mod_param.o
stats_modobs.o: Build/mod_scalars.o

stiffness.o: set_bounds.h tile.h cppdefs.h globaldefs.h
stiffness.o: Build/distribute.o Build/mod_grid.o Build/mod_iounits.o
stiffness.o: Build/mod_ocean.o Build/mod_parallel.o Build/mod_param.o
stiffness.o: Build/mod_scalars.o

timers.o: cppdefs.h globaldefs.h
timers.o: Build/distribute.o Build/mod_iounits.o Build/mod_parallel.o
timers.o: Build/mod_param.o Build/mod_strings.o

utility.o: cppdefs.h globaldefs.h
utility.o: Build/mod_kinds.o

white_noise.o: cppdefs.h globaldefs.h
white_noise.o: Build/distribute.o Build/mod_kinds.o Build/mod_parallel.o
white_noise.o: Build/mod_param.o Build/mod_scalars.o Build/nrutil.o

wpoints.o: set_bounds.h tile.h cppdefs.h globaldefs.h
wpoints.o: Build/distribute.o Build/mod_grid.o Build/mod_iounits.o
wpoints.o: Build/mod_ncparam.o Build/mod_parallel.o Build/mod_param.o
wpoints.o: Build/mod_scalars.o Build/mod_storage.o

wrt_avg.o: cppdefs.h globaldefs.h
wrt_avg.o: Build/mod_average.o Build/mod_forces.o Build/mod_grid.o
wrt_avg.o: Build/mod_iounits.o Build/mod_ncparam.o Build/mod_netcdf.o
wrt_avg.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o
wrt_avg.o: Build/mod_sediment.o Build/mod_tides.o

wrt_diags.o: cppdefs.h globaldefs.h
wrt_diags.o: Build/mod_diags.o Build/mod_grid.o Build/mod_iounits.o
wrt_diags.o: Build/mod_ncparam.o Build/mod_netcdf.o Build/mod_parallel.o
wrt_diags.o: Build/mod_param.o Build/mod_scalars.o

wrt_floats.o: cppdefs.h globaldefs.h
wrt_floats.o: Build/mod_floats.o Build/mod_iounits.o Build/mod_ncparam.o
wrt_floats.o: Build/mod_netcdf.o Build/mod_parallel.o Build/mod_param.o
wrt_floats.o: Build/mod_scalars.o Build/mod_stepping.o

wrt_gst.o: cppdefs.h globaldefs.h
wrt_gst.o: Build/distribute.o Build/mod_iounits.o Build/mod_ncparam.o
wrt_gst.o: Build/mod_netcdf.o Build/mod_parallel.o Build/mod_param.o
wrt_gst.o: Build/mod_scalars.o Build/mod_storage.o

wrt_hessian.o: cppdefs.h globaldefs.h
wrt_hessian.o: Build/mod_forces.o Build/mod_grid.o Build/mod_iounits.o
wrt_hessian.o: Build/mod_mixing.o Build/mod_ncparam.o Build/mod_netcdf.o
wrt_hessian.o: Build/mod_ocean.o Build/mod_parallel.o Build/mod_param.o
wrt_hessian.o: Build/mod_scalars.o Build/mod_stepping.o

wrt_his.o: cppdefs.h globaldefs.h
wrt_his.o: Build/mod_bbl.o Build/mod_coupling.o Build/mod_forces.o
wrt_his.o: Build/mod_grid.o Build/mod_iounits.o Build/mod_mixing.o
wrt_his.o: Build/mod_ncparam.o Build/mod_netcdf.o Build/mod_ocean.o
wrt_his.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o
wrt_his.o: Build/mod_sediment.o Build/mod_stepping.o Build/omega.o

wrt_info.o: cppdefs.h globaldefs.h
wrt_info.o: Build/extract_sta.o Build/mod_biology.o Build/mod_fourdvar.o
wrt_info.o: Build/mod_grid.o Build/mod_iounits.o Build/mod_ncparam.o
wrt_info.o: Build/mod_netcdf.o Build/mod_parallel.o Build/mod_param.o
wrt_info.o: Build/mod_scalars.o

wrt_ini.o: cppdefs.h globaldefs.h
wrt_ini.o: Build/distribute.o Build/mod_forces.o Build/mod_fourdvar.o
wrt_ini.o: Build/mod_grid.o Build/mod_iounits.o Build/mod_mixing.o
wrt_ini.o: Build/mod_ncparam.o Build/mod_netcdf.o Build/mod_ocean.o
wrt_ini.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o
wrt_ini.o: Build/mod_sediment.o Build/mod_stepping.o

wrt_rst.o: cppdefs.h globaldefs.h
wrt_rst.o: Build/mod_grid.o Build/mod_iounits.o Build/mod_mixing.o
wrt_rst.o: Build/mod_ncparam.o Build/mod_netcdf.o Build/mod_ocean.o
wrt_rst.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o
wrt_rst.o: Build/mod_sediment.o Build/mod_stepping.o

wrt_station.o: cppdefs.h globaldefs.h
wrt_station.o: Build/extract_sta.o Build/mod_bbl.o Build/mod_forces.o
wrt_station.o: Build/mod_grid.o Build/mod_iounits.o Build/mod_mixing.o
wrt_station.o: Build/mod_ncparam.o Build/mod_netcdf.o Build/mod_ocean.o
wrt_station.o: Build/mod_parallel.o Build/mod_param.o Build/mod_scalars.o
wrt_station.o: Build/mod_sediment.o Build/mod_stepping.o

wrt_tides.o: cppdefs.h globaldefs.h
wrt_tides.o: Build/mod_grid.o Build/mod_iounits.o Build/mod_ncparam.o
wrt_tides.o: Build/mod_netcdf.o Build/mod_parallel.o Build/mod_param.o
wrt_tides.o: Build/mod_scalars.o Build/mod_stepping.o Build/mod_tides.o

get_sparse_waves_matrix.o: cppdefs.h globaldefs.h
get_sparse_waves_matrix.o: Build/mod_coupler.o Build/mod_iounits.o
get_sparse_waves_matrix.o: Build/mod_ncparam.o Build/mod_netcdf.o
get_sparse_waves_matrix.o: Build/mod_scalars.o

m_constants.o: Build/swmod1.o

m_fileio.o: Build/swmod2.o

mod_xnl4v5.o: Build/m_constants.o Build/m_fileio.o Build/serv_xnl4v5.o

ocpcre.o: swancpp.h cppdefs.h globaldefs.h
ocpcre.o: Build/swmod1.o Build/swmod1.o Build/swmod1.o Build/swmod1.o

ocpids.o: swancpp.h cppdefs.h globaldefs.h
ocpids.o: Build/swmod1.o Build/swmod1.o Build/swmod1.o Build/swmod1.o
ocpids.o: Build/swmod2.o

ocpmix.o: swancpp.h cppdefs.h globaldefs.h
ocpmix.o: Build/swmod1.o Build/swmod1.o Build/swmod1.o Build/swmod1.o
ocpmix.o: Build/swmod2.o Build/swmod2.o


swancom1.o: swancpp.h cppdefs.h globaldefs.h
swancom1.o: Build/m_constants.o Build/m_fileio.o Build/mod_xnl4v5.o
swancom1.o: Build/swmod1.o Build/swmod1.o Build/swmod1.o Build/swmod1.o
swancom1.o: Build/swmod1.o Build/swmod1.o Build/swmod1.o Build/swmod1.o
swancom1.o: Build/swmod1.o Build/swmod2.o Build/swmod2.o Build/swmod2.o
swancom1.o: Build/swmod2.o Build/swmod2.o

swancom2.o: swancpp.h cppdefs.h globaldefs.h
swancom2.o: Build/swmod1.o Build/swmod1.o Build/swmod1.o Build/swmod1.o
swancom2.o: Build/swmod1.o Build/swmod2.o

swancom3.o: swancpp.h cppdefs.h globaldefs.h
swancom3.o: Build/swmod1.o Build/swmod1.o Build/swmod1.o Build/swmod1.o
swancom3.o: Build/swmod1.o Build/swmod1.o Build/swmod1.o Build/swmod1.o

swancom4.o: swancpp.h cppdefs.h globaldefs.h
swancom4.o: Build/mod_xnl4v5.o Build/serv_xnl4v5.o Build/swmod1.o
swancom4.o: Build/swmod1.o Build/swmod1.o Build/swmod2.o Build/swmod2.o

swancom5.o: swancpp.h cppdefs.h globaldefs.h
swancom5.o: Build/swmod1.o Build/swmod1.o Build/swmod1.o Build/swmod1.o
swancom5.o: Build/swmod1.o Build/swmod1.o Build/swmod2.o Build/swmod2.o

swanmain.o: swancpp.h cppdefs.h globaldefs.h
swanmain.o: Build/mod_param.o Build/swmod1.o Build/swmod1.o Build/swmod1.o
swanmain.o: Build/swmod1.o Build/swmod1.o Build/swmod1.o Build/swmod1.o
swanmain.o: Build/swmod1.o Build/swmod1.o Build/swmod2.o Build/swmod2.o
swanmain.o: Build/swmod2.o Build/swmod2.o Build/swmod2.o Build/swmod2.o
swanmain.o: Build/swmod2.o Build/swmod2.o Build/swmod2.o Build/swmod2.o
swanmain.o: Build/swmod2.o Build/swpoint.o Build/waves_coupler.o

swanout1.o: swancpp.h cppdefs.h globaldefs.h
swanout1.o: Build/mod_param.o Build/swmod1.o Build/swmod1.o Build/swmod1.o
swanout1.o: Build/swmod1.o Build/swmod1.o Build/swmod1.o Build/swmod1.o
swanout1.o: Build/swmod2.o Build/swmod2.o Build/swmod2.o Build/swmod2.o
swanout1.o: Build/waves_coupler.o

swanout2.o: swancpp.h cppdefs.h globaldefs.h
swanout2.o: Build/swmod1.o Build/swmod1.o Build/swmod1.o Build/swmod1.o
swanout2.o: Build/swmod1.o Build/swmod1.o Build/swmod1.o Build/swmod2.o
swanout2.o: Build/swmod2.o

swanparll.o: swancpp.h cppdefs.h globaldefs.h
swanparll.o: Build/swmod1.o Build/swmod1.o Build/swmod1.o Build/swmod1.o
swanparll.o: Build/swmod1.o Build/swmod1.o Build/swmod2.o Build/swmod2.o
swanparll.o: Build/swmod2.o Build/swmod2.o

swanpre1.o: swancpp.h cppdefs.h globaldefs.h
swanpre1.o: Build/swmod1.o Build/swmod1.o Build/swmod1.o Build/swmod1.o
swanpre1.o: Build/swmod1.o Build/swmod1.o Build/swmod1.o Build/swmod1.o
swanpre1.o: Build/swmod1.o Build/swmod2.o Build/swmod2.o Build/swmod2.o
swanpre1.o: Build/swmod2.o Build/swmod2.o Build/swmod2.o

swanpre2.o: swancpp.h cppdefs.h globaldefs.h
swanpre2.o: Build/swmod1.o Build/swmod1.o Build/swmod1.o Build/swmod1.o
swanpre2.o: Build/swmod1.o Build/swmod1.o Build/swmod1.o Build/swmod1.o
swanpre2.o: Build/swmod2.o Build/swmod2.o Build/swmod2.o

swanser.o: swancpp.h cppdefs.h globaldefs.h
swanser.o: Build/swmod1.o Build/swmod1.o Build/swmod1.o Build/swmod1.o
swanser.o: Build/swmod1.o Build/swmod1.o Build/swmod1.o Build/swmod2.o
swanser.o: Build/swmod2.o Build/swmod2.o Build/swmod2.o Build/swmod2.o

swmod1.o: swancpp.h cppdefs.h globaldefs.h

swmod2.o: swancpp.h cppdefs.h globaldefs.h
swmod2.o: Build/swmod1.o

swpoint.o: cppdefs.h globaldefs.h
swpoint.o: Build/swmod1.o Build/swmod1.o Build/swmod1.o Build/swmod1.o
swpoint.o: Build/swmod1.o Build/swmod1.o Build/swmod1.o Build/swmod2.o
swpoint.o: Build/swmod2.o Build/swmod2.o Build/swmod2.o Build/swmod2.o
swpoint.o: Build/swmod2.o Build/swmod2.o

waves_control.o: cppdefs.h globaldefs.h
waves_control.o: Build/mod_coupler.o Build/mod_fourdvar.o Build/mod_iounits.o
waves_control.o: Build/mod_ncparam.o Build/mod_parallel.o Build/mod_param.o
waves_control.o: Build/mod_scalars.o Build/swmod1.o Build/swmod1.o
waves_control.o: Build/swmod2.o Build/swmod2.o Build/swmod2.o Build/swmod2.o

waves_coupler.o: cppdefs.h globaldefs.h
waves_coupler.o: Build/mod_coupler.o Build/mod_iounits.o Build/mod_scalars.o
waves_coupler.o: Build/swmod1.o Build/swmod1.o Build/swmod1.o Build/swmod1.o
waves_coupler.o: Build/swmod2.o Build/swmod2.o Build/swmod2.o Build/swmod2.o

analytical_mod.mod: analytical.o
back_cost_mod.mod: back_cost.o
back_cov_mod.mod: back_cov.o
back_step_mod.mod: back_step.o
bbl_mod.mod: bbl.o
bc_2d_mod.mod: bc_2d.o
bc_3d_mod.mod: bc_3d.o
bed_cohesive_mod.mod: bed_cohesive.o
bed_mod.mod: bed.o
bedbiodiff_mod.mod: bedbiodiff.o
bedload_mod.mod: bedload.o
biology_mod.mod: biology.o
bottom_mod.mod: bottom.o
bulk_flux_mod.mod: bulk_flux.o
bvf_mix_mod.mod: bvf_mix.o
cgradient_mod.mod: cgradient.o
conv_2d_mod.mod: conv_2d.o
conv_3d_mod.mod: conv_3d.o
cost_grad_mod.mod: cost_grad.o
cost_norm_mod.mod: cost_norm.o
dep_ero_mod.mod: dep_ero.o
descent_mod.mod: descent.o
diag_mod.mod: diag.o
distribute_mod.mod: distribute.o
dotproduct_mod.mod: dotproduct.o
downhill_mod.mod: downhill.o
exchange_2d_mod.mod: exchange_2d.o
exchange_3d_mod.mod: exchange_3d.o
extract_obs_mod.mod: extract_obs.o
extract_sta_mod.mod: extract_sta.o
forcing_mod.mod: forcing.o
frc_adjust_mod.mod: frc_adjust.o
get_2dparent_data_mod.mod: get_2dparent_data.o
get_3dparent_data_mod.mod: get_3dparent_data.o
gls_corstep_mod.mod: gls_corstep.o
gls_prestep_mod.mod: gls_prestep.o
hmixing_mod.mod: hmixing.o
impulse_mod.mod: impulse.o
ini_adjust_mod.mod: ini_adjust.o
ini_fields_mod.mod: ini_fields.o
ini_lanczos_mod.mod: ini_lanczos.o
init_child_hindices_mod.mod: init_child_hindices.o
init_parent_hindices_mod.mod: init_parent_hindices.o
interpolate_mod.mod: interpolate.o
lmd_bkpp_mod.mod: lmd_bkpp.o
lmd_skpp_mod.mod: lmd_skpp.o
lmd_vmix_mod.mod: lmd_vmix.o
m_bndspec.mod: swmod2.o
m_coupling.mod: swmod2.o
m_cvmesh.mod: swmod2.o
m_diffr.mod: swmod2.o
m_genarr.mod: swmod2.o
m_mpi.mod: swmod2.o
m_obsta.mod: swmod2.o
m_parall.mod: swmod2.o
m_parall2.mod: swmod2.o
m_snl4.mod: swmod2.o
m_wcap.mod: swmod2.o
m_xnldata.mod: mod_xnl4v5.o
metrics_mod.mod: metrics.o
mp_exchange_mod.mod: mp_exchange.o
mpdata_adiff_mod.mod: mpdata_adiff.o
my25_corstep_mod.mod: my25_corstep.o
my25_prestep_mod.mod: my25_prestep.o
nesting_mod.mod: nesting.o
normalization_mod.mod: normalization.o
obc_volcons_mod.mod: obc_volcons.o
ocean_control_mod.mod: ocean_control.o
ocean_coupler_mod.mod: ocean_coupler.o
ocpcomm1.mod: swmod1.o
ocpcomm2.mod: swmod1.o
ocpcomm3.mod: swmod1.o
ocpcomm4.mod: swmod1.o
oi_update_mod.mod: oi_update.o
omega_mod.mod: omega.o
outp_data.mod: swmod2.o
packing_mod.mod: packing.o
pre_step3d_mod.mod: pre_step3d.o
propagator_mod.mod: propagator.o
prsgrd_mod.mod: prsgrd.o
radiation_stress_mod.mod: radiation_stress.o
rho_eos_mod.mod: rho_eos.o
rhs3d_mod.mod: rhs3d.o
roms_export_mod.mod: roms_export.o
roms_import_mod.mod: roms_import.o
sediment_mod.mod: sediment.o
set_2dchild_data_mod.mod: set_2dchild_data.o
set_2dfld_mod.mod: set_2dfld.o
set_2dfldr_mod.mod: set_2dfldr.o
set_2dparent_data_mod.mod: set_2dparent_data.o
set_3dchild_data_mod.mod: set_3dchild_data.o
set_3dfld_mod.mod: set_3dfld.o
set_3dfldr_mod.mod: set_3dfldr.o
set_3dparent_data_mod.mod: set_3dparent_data.o
set_avg_mod.mod: set_avg.o
set_depth_mod.mod: set_depth.o
set_massflux_mod.mod: set_massflux.o
set_tides_mod.mod: set_tides.o
set_vbc_mod.mod: set_vbc.o
set_zeta_mod.mod: set_zeta.o
settling_mod.mod: settling.o
shapiro_mod.mod: shapiro.o
state_addition_mod.mod: state_addition.o
state_copy_mod.mod: state_copy.o
state_dotprod_mod.mod: state_dotprod.o
state_initialize_mod.mod: state_initialize.o
state_scale_mod.mod: state_scale.o
step2d_mod.mod: step2d.o
step3d_t_mod.mod: step3d_t.o
step3d_uv_mod.mod: step3d_uv.o
step_floats_mod.mod: step_floats.o
stiffness_mod.mod: stiffness.o
swcomm1.mod: swmod1.o
swcomm2.mod: swmod1.o
swcomm3.mod: swmod1.o
swcomm4.mod: swmod1.o
swpoint_mod.mod: swpoint.o
t3dbc_mod.mod: t3dbc_im.o
t3dmix_mod.mod: t3dmix.o
timecomm.mod: swmod1.o
tkebc_mod.mod: tkebc_im.o
u2dbc_mod.mod: u2dbc_im.o
u3dbc_mod.mod: u3dbc_im.o
utility_mod.mod: utility.o
uv3dmix_mod.mod: uv3dmix.o
v2dbc_mod.mod: v2dbc_im.o
v3dbc_mod.mod: v3dbc_im.o
waves_control_mod.mod: waves_control.o
waves_coupler_mod.mod: waves_coupler.o
wetdry_mod.mod: wetdry.o
white_noise_mod.mod: white_noise.o
wpoints_mod.mod: wpoints.o
wvelocity_mod.mod: wvelocity.o
zetabc_mod.mod: zetabc.o
