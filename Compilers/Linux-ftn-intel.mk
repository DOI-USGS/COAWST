# svn $Id: Linux-ftn-intel.mk 885 2017-12-27 23:18:30Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2017 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Include file for CRAY FTN cross-compiler with Linux
# -------------------------------------------------------------------------
#
# ARPACK_LIBDIR  ARPACK libary directory
# FC             Name of the fortran compiler to use
# FFLAGS         Flags to the fortran compiler
# CPP            Name of the C-preprocessor
# CPPFLAGS       Flags to the C-preprocessor
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
               FC := ftn
           FFLAGS := -heap-arrays -fp-model precise
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -traditional
          LDFLAGS :=
               AR := ar
          ARFLAGS := -r
            MKDIR := mkdir -p
               RM := rm -f
           RANLIB := touch
             PERL := perl
             TEST := test

        MDEPFLAGS := --cpp --fext=f90 --file=- --objdir=$(SCRATCH_DIR)

#
# Library locations, can be overridden by environment variables.
#

ifdef USE_NETCDF4
        NF_CONFIG ?= nf-config
    NETCDF_INCDIR ?= $(shell $(NF_CONFIG) --prefix)/include
             LIBS := $(shell $(NF_CONFIG) --flibs)
else
    NETCDF_INCDIR ?= /usr/local/include
    NETCDF_LIBDIR ?= /usr/local/lib
             LIBS := -L$(NETCDF_LIBDIR) -lnetcdf
endif

ifdef USE_ARPACK
 ifdef USE_MPI
   PARPACK_LIBDIR ?= /usr/local/lib
             LIBS += -L$(PARPACK_LIBDIR) -lparpack
 endif
    ARPACK_LIBDIR ?= /usr/local/lib
             LIBS += -L$(ARPACK_LIBDIR) -larpack
endif

ifdef USE_MPI
         CPPFLAGS += -DMPI
endif

ifdef USE_OpenMP
         CPPFLAGS += -D_OPENMP
           FFLAGS += -openmp -fpp
endif

ifdef USE_DEBUG
           FFLAGS += -g -check bounds -traceback -check uninit -warn interfaces,nouncalled -gen-interfaces
else
           FFLAGS += -ip -O3
endif

ifdef USE_MCT
       MCT_INCDIR ?= /usr/local/mct/include
       MCT_LIBDIR ?= /usr/local/mct/lib
           FFLAGS += -I$(MCT_INCDIR)
             LIBS += -L$(MCT_LIBDIR) -lmct -lmpeu
endif

ifdef USE_ESMF
          ESMF_OS ?= $(OS)
      ESMF_SUBDIR := $(ESMF_OS).$(ESMF_COMPILER).$(ESMF_ABI).$(ESMF_COMM).$(ESMF_SITE)
      ESMF_MK_DIR ?= $(ESMF_DIR)/lib/lib$(ESMF_BOPT)/$(ESMF_SUBDIR)
                     include $(ESMF_MK_DIR)/esmf.mk
           FFLAGS += $(ESMF_F90COMPILEPATHS)
             LIBS += $(ESMF_F90LINKPATHS) $(ESMF_F90ESMFLINKLIBS)
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
endif

ifdef USE_WW3
             FFLAGS += -I${COAWST_WW3_DIR}/mod_DIST/
             LIBS += WW3/obj/libWW3.a
endif

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
$(SCRATCH_DIR)/m_constants.o: FFLAGS += -free
$(SCRATCH_DIR)/m_fileio.o: FFLAGS += -free
$(SCRATCH_DIR)/mod_xnl4v5.o: FFLAGS += -free
$(SCRATCH_DIR)/serv_xnl4v5.o: FFLAGS += -free

endif
