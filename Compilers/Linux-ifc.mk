# svn $Id: Linux-ifc.mk 655 2008-07-25 18:57:05Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2016 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Include file for Intel IFC (version 7.x) compiler on Linux
# -------------------------------------------------------------------------
#
# ARPACK_LIBDIR  ARPACK libary directory
# FC             Name of the fortran compiler to use
# FFLAGS         Flags to the fortran compiler
# CPP            Name of the C-preprocessor
# CPPFLAGS       Flags to the C-preprocessor
# CC             Name of the C compiler
# CFLAGS         Flags to the C compiler
# CXX            Name of the C++ compiler
# CXXFLAGS       Flags to the C++ compiler
# CLEAN          Name of cleaning executable after C-preprocessing
# NETCDF_INCDIR  NetCDF include directory
# NETCDF_LIBDIR  NetCDF libary directory
# LD             Program to load the objects into an executable
# LDFLAGS        Flags to the loader
# RANLIB         Name of ranlib command
# MDEPFLAGS      Flags for sfmakedepend  (-s if you keep .f files)
#
# First the defaults
#
               FC := ifc
           FFLAGS :=
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -traditional
               CC := gcc
              CXX := g++
           CFLAGS :=
         CXXFLAGS :=
          LDFLAGS :=
               AR := ar
          ARFLAGS := r
            MKDIR := mkdir -p
               RM := rm -f
           RANLIB := ranlib
             PERL := perl
             TEST := test

        MDEPFLAGS := --cpp --fext=f90 --file=- --objdir=$(SCRATCH_DIR)

#
# Library locations, can be overridden by environment variables.
#

ifdef USE_NETCDF4
        NC_CONFIG ?= nc-config
    NETCDF_INCDIR ?= $(shell $(NC_CONFIG) --prefix)/include
             LIBS := $(shell $(NC_CONFIG) --flibs)
else
    NETCDF_INCDIR ?= /usr/local/include
    NETCDF_LIBDIR ?= /usr/local/lib
             LIBS := -L$(NETCDF_LIBDIR) -lnetcdf -lnetcdff
endif

ifdef USE_ARPACK
 ifdef USE_MPI
   PARPACK_LIBDIR ?= /opt/intelsoft/PARPACK
             LIBS += -L$(PARPACK_LIBDIR) -lparpack
 endif
    ARPACK_LIBDIR ?= /opt/intelsoft/PARPACK
             LIBS += -L$(ARPACK_LIBDIR) -larpack
endif

ifdef USE_MPI
         CPPFLAGS += -DMPI
 ifdef USE_MPIF90
               FC := /opt/intelsoft/mpich2/bin/mpif90
 else
             LIBS += -lfmpi-pgi -lmpi-pgi
 endif
endif

ifdef USE_OpenMP
         CPPFLAGS += -D_OPENMP
           FFLAGS += -openmp
endif

ifdef USE_DEBUG
           FFLAGS += -g -CA -CB -CS -CU
           CFLAGS += -g
         CXXFLAGS += -g
else
           FFLAGS += -ip -O3
           CFLAGS += -g
         CXXFLAGS += -g
endif

ifdef USE_MCT
       MCT_INCDIR ?= /opt/intelsoft/mct/include
       MCT_LIBDIR ?= /opt/intelsoft/mct/lib
           FFLAGS += -I$(MCT_INCDIR)
             LIBS += -L$(MCT_LIBDIR) -lmct -lmpeu
endif

ifdef USE_ESMF
      ESMF_SUBDIR := $(ESMF_OS).$(ESMF_COMPILER).$(ESMF_ABI).$(ESMF_COMM).$(ESMF_SITE)
      ESMF_MK_DIR ?= $(ESMF_DIR)/lib/lib$(ESMF_BOPT)/$(ESMF_SUBDIR)
                     include $(ESMF_MK_DIR)/esmf.mk
           FFLAGS += $(ESMF_F90COMPILEPATHS)
             LIBS += $(ESMF_F90LINKPATHS) -lesmf -lC
endif

ifdef USE_CXX
             LIBS += -lstdc++
endif

ifdef USE_WRF
             FFLAGS += -I$(WRF_DIR)/main -I$(WRF_DIR)/external/esmf_time_f90 -I$(WRF_DIR)/frame -I$(WRF_DIR)/share
             LIBS += WRF/main/module_wrf_top.o
             LIBS += WRF/main/libwrflib.a
             LIBS += WRF/external/fftpack/fftpack5/libfftpack.a
             LIBS += WRF/external/io_grib1/libio_grib1.a
             LIBS += WRF/external/io_grib_share/libio_grib_share.a
             LIBS += WRF/external/io_int/libwrfio_int.a
             LIBS += WRF/external/esmf_time_f90/libesmf_time.a
             LIBS += WRF/external/RSL_LITE/librsl_lite.a
             LIBS += WRF/frame/module_internal_header_util.o
             LIBS += WRF/frame/pack_utils.o
             LIBS += WRF/external/io_netcdf/libwrfio_nf.a
#            LIBS += WRF/external/io_netcdf/wrf_io.o
endif

       clean_list += ifc* work.pc*

#
# Use full path of compiler.
#
               FC := $(shell which ${FC})
               LD := $(FC)

#
# Set free form format in source files to allow long string for
# local directory and compilation flags inside the code.
#

$(SCRATCH_DIR)/mod_ncparam.o: FFLAGS += -free
$(SCRATCH_DIR)/mod_strings.o: FFLAGS += -free
$(SCRATCH_DIR)/analytical.o: FFLAGS += -free
$(SCRATCH_DIR)/biology.o: FFLAGS += -free
ifdef USE_ADJOINT
$(SCRATCH_DIR)/ad_biology.o: FFLAGS += -free
endif
ifdef USE_REPRESENTER
$(SCRATCH_DIR)/rp_biology.o: FFLAGS += -free
endif
ifdef USE_TANGENT
$(SCRATCH_DIR)/tl_biology.o: FFLAGS += -free
endif

#
# Supress free format in SWAN source files since there are comments
# beyond column 72.
#

ifdef USE_SWAN

$(SCRATCH_DIR)/ocpcre.o: FFLAGS += -nofree
$(SCRATCH_DIR)/ocpids.o: FFLAGS += -nofree
$(SCRATCH_DIR)/ocpmix.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swancom1.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swancom2.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swancom3.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swancom4.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swancom5.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swanmain.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swanout1.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swanout2.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swanparll.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swanpre1.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swanpre2.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swanser.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swmod1.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swmod2.o: FFLAGS += -nofree
$(SCRATCH_DIR)/SwanCompdata.o: FFLAGS += -free
$(SCRATCH_DIR)/SwanGriddata.o: FFLAGS += -free
$(SCRATCH_DIR)/m_constants.o:  FFLAGS += -free
$(SCRATCH_DIR)/m_fileio.o:     FFLAGS += -free
$(SCRATCH_DIR)/mod_xnl4v5.o:   FFLAGS += -free
$(SCRATCH_DIR)/serv_xnl4v5.o:  FFLAGS += -free
$(SCRATCH_DIR)/nctablemd.o:    FFLAGS += -free
$(SCRATCH_DIR)/agioncmd.o:     FFLAGS += -free
$(SCRATCH_DIR)/swn_outnc.o:    FFLAGS += -free

endif
