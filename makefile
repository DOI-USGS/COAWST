# svn $Id: makefile 807 2016-07-09 02:03:55Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
# Copyright (c) 2002-2016 The ROMS/TOMS Group             Kate Hedstrom :::
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

ifneq (3.80,$(firstword $(sort $(MAKE_VERSION) 3.80)))
 $(error This makefile requires GNU make version 3.80 or higher. \
		Your current version is: $(MAKE_VERSION))
endif

#--------------------------------------------------------------------------
#  Initialize some things.
#--------------------------------------------------------------------------

  sources    :=
  libraries  :=
  c_sources  := 

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

ROMS_APPLICATION ?= SHOREFACE

#  If application header files is not located in "ROMS/Include",
#  provide an alternate directory FULL PATH.

MY_HEADER_DIR ?=

#  If your application requires analytical expressions and they are
#  not located in "ROMS/Functionals", provide an alternate directory.
#  Notice that a set analytical expressions templates can be found in
#  "User/Functionals".
#
#  If applicable, also used this directory to place your customized
#  biology model header file (like fennel.h, nemuro.h, ecosim.h, etc).

MY_ANALYTICAL_DIR ?=

# If applicable, where does CICE put its binary files?

MY_CICE_DIR ?= /center/w/kate/CICE/NEP/compile

#  Sometimes it is desirable to activate one or more CPP options to
#  run different variants of the same application without modifying
#  its header file. If this is the case, specify such options here
#  using the -D syntax.  For example, to write time-averaged fields
#  set:
#
#    MY_CPP_FLAGS ?= -DAVERAGES
#

MY_CPP_FLAGS ?=

#  Activate debugging compiler options:

   USE_DEBUG ?=

#  If parallel applications, use at most one of these definitions
#  (leave both definitions blank in serial applications):

     USE_MPI ?=
  USE_OpenMP ?=

#  If distributed-memory, turn on compilation via the script "mpif90".
#  This is needed in some Linux operating systems. In some systems with
#  native MPI libraries the compilation does not require MPICH type
#  scripts. This macro is also convient when there are several fortran
#  compiliers (ifort, pgf90, pathf90) in the system that use mpif90.
#  In this, case the user need to select the desired compiler below and
#  turn on both USE_MPI and USE_MPIF90 macros.

  USE_MPIF90 ?=

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

        FORT ?= pgi

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

ifdef USE_ROMS
#--------------------------------------------------------------------------
#  Notice that the token "libraries" is initialize with the ROMS/Utility
#  library to account for calls to objects in other ROMS libraries or
#  cycling dependencies. These type of dependencies are problematic in
#  some compilers during linking. This library appears twice at linking
#  step (begining and almost the end of ROMS library list).
#--------------------------------------------------------------------------

  libraries  := $(SCRATCH_DIR)/libUTIL.a
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

  COMPILERS ?= $(CURDIR)/Compilers

MAKE_MACROS := $(shell echo ${HOME} | sed 's| |\\ |g')/make_macros.mk

ifneq "$(MAKECMDGOALS)" "clean"
 ifneq "$(MAKECMDGOALS)" "tarfile"
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
endif

clean_list += $(MAKE_MACROS)

#--------------------------------------------------------------------------
#  Make functions for putting the temporary files in $(SCRATCH_DIR)
#  DO NOT modify this section; spaces and blank lines are needed.
#--------------------------------------------------------------------------

# $(call source-dir-to-binary-dir, directory-list)
source-dir-to-binary-dir = $(addprefix $(SCRATCH_DIR)/, $(notdir $1))

# $(call source-to-object, source-file-list)
source-to-object = $(call source-dir-to-binary-dir,   \
                   $(subst .F,.o,$1))

# $(call source-to-object, source-file-list)
c-source-to-object = $(call source-dir-to-binary-dir,       \
                     $(subst .c,.o,$(filter %.c,$1))        \
                     $(subst .cc,.o,$(filter %.cc,$1)))

# $(call make-library, library-name, source-file-list)
define make-library
   libraries += $(SCRATCH_DIR)/$1
   sources   += $2

   $(SCRATCH_DIR)/$1: $(call source-dir-to-binary-dir,    \
                      $(subst .F,.o,$2))
	$(AR) $(ARFLAGS) $$@ $$^
	$(RANLIB) $$@
endef

# $(call make-c-library, library-name, source-file-list)
define make-c-library
   libraries += $(SCRATCH_DIR)/$1
   c_sources += $2

   $(SCRATCH_DIR)/$1: $(call source-dir-to-binary-dir,    \
                      $(subst .c,.o,$(filter %.c,$2))     \
                      $(subst .cc,.o,$(filter %.cc,$2)))
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

# $(c-compile-rules)
define c-compile-rules
  $(foreach f, $(local_c_src),       \
    $(call one-c-compile-rule,$(call c-source-to-object,$f), $f))
endef

# $(call one-compile-rule, binary-file, f90-file, source-file)
define one-compile-rule
  $1: $2 $3
	cd $$(SCRATCH_DIR); $$(FC) -c $$(FFLAGS) $(notdir $2)

  $2: $3
	$$(CPP) $$(CPPFLAGS) $$(MY_CPP_FLAGS) $$< > $$@
	$$(CLEAN) $$@

endef

# $(call one-c-compile-rule, binary-file, source-file)
define one-c-compile-rule
  $1: $2
	cd $$(SCRATCH_DIR); $$(CXX) -c $$(CXXFLAGS) $$<

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

ifndef FORT
  $(error Variable FORT not set)
endif

ifneq "$(MAKECMDGOALS)" "clean"
 ifneq "$(MAKECMDGOALS)" "tarfile"
  include $(COMPILERS)/$(OS)-$(strip $(FORT)).mk
 endif
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

ifdef USE_ROMS
.PHONY: all

ifdef USE_CICE
all: $(SCRATCH_DIR) $(SCRATCH_DIR)/libCICE.a
endif

all: $(SCRATCH_DIR) $(SCRATCH_DIR)/MakeDepend $(BIN) rm_macros
endif

ifdef USE_SWAN
.PHONY: all

all: $(SCRATCH_DIR) $(SCRATCH_DIR)/MakeDepend $(BIN) rm_macros
endif

 modules  :=
ifdef USE_ADJOINT
 modules  +=	ROMS/Adjoint \
		ROMS/Adjoint/Biology
endif
ifdef USE_REPRESENTER
 modules  +=	ROMS/Representer \
		ROMS/Representer/Biology
endif
ifdef USE_TANGENT
 modules  +=	ROMS/Tangent \
		ROMS/Tangent/Biology
endif
ifdef USE_ROMS
 modules  +=	ROMS/Nonlinear \
		ROMS/Nonlinear/Biology \
		ROMS/Nonlinear/Sediment \
		ROMS/Nonlinear/Wec \
		ROMS/Nonlinear/Vegetation \
		ROMS/Functionals \
		ROMS/Utility \
		ROMS/Modules
endif
ifdef USE_SEAICE
 modules  +=	ROMS/SeaIce
endif
ifdef USE_CICE
 modules  +=	SeaIce/Extra
    LIBS  +=    $(SCRATCH_DIR)/libCICE.a
endif

 includes :=	ROMS/Include
ifdef MY_ANALYTICAL
 includes +=	$(MY_ANALYTICAL_DIR)
endif
ifdef USE_ADJOINT
 includes +=	ROMS/Adjoint \
		ROMS/Adjoint/Biology
endif
ifdef USE_REPRESENTER
 includes +=	ROMS/Representer \
		ROMS/Representer/Biology
endif
ifdef USE_SEAICE
 includes +=	ROMS/SeaIce
endif
ifdef USE_TANGENT
 includes +=	ROMS/Tangent \
		ROMS/Tangent/Biology
endif
ifdef USE_ROMS
 includes +=	ROMS/Nonlinear \
		ROMS/Nonlinear/Biology \
		ROMS/Nonlinear/Sediment \
		ROMS/Nonlinear/Wec \
		ROMS/Nonlinear/Vegetation \
		ROMS/Utility \
		ROMS/Drivers \
		ROMS/Functionals
endif
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

ifdef USE_INWAVE
 modules  +=	InWave/Drivers \
             InWave/Action_balance \
             InWave/Boundaries \
             InWave/Modules \
             InWave/Utility
 includes +=	InWave/Drivers \
             InWave/Action_balance \
             InWave/Boundaries \
             InWave/Modules \
             InWave/Utility
endif

ifdef USE_ROMS
 modules  +=	Master
 includes +=	Master Compilers
else ifdef USE_SWAN
 modules  +=	Master
 includes +=	Master Compilers
endif

vpath %.F $(modules)
vpath %.cc $(modules)
vpath %.h $(includes)
vpath %.f90 $(SCRATCH_DIR)
vpath %.o $(SCRATCH_DIR)

include $(addsuffix /Module.mk,$(modules))

MDEPFLAGS += $(patsubst %,-I %,$(includes)) --silent --moddir $(SCRATCH_DIR)

CPPFLAGS  += $(patsubst %,-I%,$(includes))

ifdef MY_HEADER_DIR
  CPPFLAGS += -D'HEADER_DIR="$(MY_HEADER_DIR)"'
else
  CPPFLAGS += -D'HEADER_DIR="$(ROOTDIR)/ROMS/Include"'
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
#  Build WRF.
#--------------------------------------------------------------------------

.PHONY: wrfclean

wrfclean:
ifdef USE_WRF
	cd $(WRF_DIR); ls; ./clean -a;                            \
	./configure;                                              \
	echo " "; echo " ";                                       \
	echo "cleaned wrf";
endif

.PHONY: wrf

wrf:
ifdef USE_WRF
	WRF_MAKE_EXE=NO
 ifndef USE_ROMS
  ifndef USE_SWAN
	WRF_MAKE_EXE=YES
  endif
 endif
	cp -p $(BIN) $(BIN).backup;                               \
	$(RM) -r $(BIN);                                          \
	cd $(WRF_DIR); ls;                                        \
	echo " "; echo " ";                                       \
	echo "Compiling wrf";                                     \
	./compile em_real;                                        \
	echo "";                                                  \
	echo "-------- Finished compiling WRF ------------"
 ifndef USE_ROMS
  ifndef USE_SWAN
	ln -sf WRF/main/wrf.exe coawstM;
  endif
 endif
	echo "";
endif

#--------------------------------------------------------------------------
#  Target to create ROMS/TOMS dependecies.
#--------------------------------------------------------------------------
ifneq "$(MAKECMDGOALS)" "tarfile"
$(SCRATCH_DIR)/$(NETCDF_MODFILE): | $(SCRATCH_DIR)
	cp -f $(NETCDF_INCDIR)/$(NETCDF_MODFILE) $(SCRATCH_DIR)

$(SCRATCH_DIR)/$(TYPESIZES_MODFILE): | $(SCRATCH_DIR)
	cp -f $(NETCDF_INCDIR)/$(TYPESIZES_MODFILE) $(SCRATCH_DIR)

$(SCRATCH_DIR)/libCICE.a: $(MY_CICE_DIR)/libCICE.a
	cp -f $(MY_CICE_DIR)/libCICE.a $(MY_CICE_DIR)/*.mod $(SCRATCH_DIR)

$(MY_CICE_DIR)/libCICE.a:
	SeaIce/comp_ice
ifdef USE_CICE
$(SCRATCH_DIR)/initial.o: $(MY_CICE_DIR)/CICE_InitMod.o
$(SCRATCH_DIR)/ice_fakecpl.o: $(MY_CICE_DIR)/CICE_RunMod.o
$(SCRATCH_DIR)/ice_fakecpl.o: $(MY_CICE_DIR)/ice_blocks.o
$(SCRATCH_DIR)/ice_fakecpl.o: $(MY_CICE_DIR)/ice_broadcast.o
$(SCRATCH_DIR)/ice_fakecpl.o: $(MY_CICE_DIR)/ice_calendar.o
$(SCRATCH_DIR)/ice_fakecpl.o: $(MY_CICE_DIR)/ice_communicate.o
$(SCRATCH_DIR)/ice_fakecpl.o: $(MY_CICE_DIR)/ice_constants.o
$(SCRATCH_DIR)/ice_fakecpl.o: $(MY_CICE_DIR)/ice_domain.o
$(SCRATCH_DIR)/ice_fakecpl.o: $(MY_CICE_DIR)/ice_domain_size.o
$(SCRATCH_DIR)/ice_fakecpl.o: $(MY_CICE_DIR)/ice_fileunits.o
$(SCRATCH_DIR)/ice_fakecpl.o: $(MY_CICE_DIR)/ice_flux.o
$(SCRATCH_DIR)/ice_fakecpl.o: $(MY_CICE_DIR)/ice_gather_scatter.o
$(SCRATCH_DIR)/ice_fakecpl.o: $(MY_CICE_DIR)/ice_grid.o
$(SCRATCH_DIR)/ice_fakecpl.o: $(MY_CICE_DIR)/ice_history.o
$(SCRATCH_DIR)/ice_fakecpl.o: $(MY_CICE_DIR)/ice_init.o
$(SCRATCH_DIR)/ice_fakecpl.o: $(MY_CICE_DIR)/ice_kinds_mod.o
$(SCRATCH_DIR)/ice_fakecpl.o: $(MY_CICE_DIR)/ice_restart.o
$(SCRATCH_DIR)/ice_fakecpl.o: $(MY_CICE_DIR)/ice_restart_shared.o
$(SCRATCH_DIR)/ice_fakecpl.o: $(MY_CICE_DIR)/ice_state.o
$(SCRATCH_DIR)/ice_fakecpl.o: $(MY_CICE_DIR)/ice_timers.o
endif

$(SCRATCH_DIR)/MakeDepend: makefile \
                           $(SCRATCH_DIR)/$(NETCDF_MODFILE) \
                           $(SCRATCH_DIR)/$(TYPESIZES_MODFILE) \
                           | $(SCRATCH_DIR)
	$(SFMAKEDEPEND) $(MDEPFLAGS) $(sources) > $(SCRATCH_DIR)/MakeDepend
	cp -p $(MAKE_MACROS) $(SCRATCH_DIR)

.PHONY: depend

SFMAKEDEPEND := ./ROMS/Bin/sfmakedepend

depend: $(SCRATCH_DIR)
	$(SFMAKEDEPEND) $(MDEPFLAGS) $(sources) > $(SCRATCH_DIR)/MakeDepend
endif

ifneq "$(MAKECMDGOALS)" "clean"
  -include $(SCRATCH_DIR)/MakeDepend
endif

#--------------------------------------------------------------------------
#  Target to create ROMS/TOMS tar file.
#--------------------------------------------------------------------------

.PHONY: tarfile

tarfile:
		tar --exclude=".svn" --exclude Output -cvf coawst_v3.2.tar *

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
	$(RM) -r $(MAKE_MACROS)

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
