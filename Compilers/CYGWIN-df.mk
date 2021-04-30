# svn $Id: CYGWIN-df.mk 1054 2021-03-06 19:47:12Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2021 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Include file for Compaq Visual Fortran compiler on Cygwin
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
# LIBS           Required libraries during linking
# NETCDF_INCDIR  NetCDF include directory
# NETCDF_LIBDIR  NetCDF libary directory
# LD             Program to load the objects into an executable
# LDFLAGS        Flags to the loader
# RANLIB         Name of ranlib command
# MDEPFLAGS      Flags for sfmakedepend  (-s if you keep .f files)
#
# First the defaults
#

              BIN := $(BIN).exe

               FC := df
           FFLAGS := /stand:f95
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -DCYGWIN
               CC := gcc
              CXX := g++
           CFLAGS :=
         CXXFLAGS :=
             LIBS :=
ifdef USE_ROMS
             LIBS += $(SCRATCH_DIR)/libNLM.a         # cyclic dependencies
endif
          LDFLAGS := /link /nodefaultlib:libcmt /nodefaultlib:libifcore /stack:67108864
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
# These are be specified in Unix form and will be converted as
# necessary to Windows form for Windows-native commands. The default
# values below assume that Cygwin mounts have been defined pointing to
# the NETCDF and ARPACK library locations.
#

ifdef USE_NETCDF4
        NF_CONFIG ?= nf-config
    NETCDF_INCDIR ?= $(shell $(NF_CONFIG) --prefix)/include
             LIBS += $(shell $(NF_CONFIG) --flibs)
else
    NETCDF_INCDIR ?= /usr/local/include
    NETCDF_LIBDIR ?= /usr/local/lib
             LIBS += -L$(NETCDF_LIBDIR) -lnetcdf
endif

ifdef USE_ARPACK
 ifdef USE_MPI
   PARPACK_LIBDIR ?= /usr/local/arpack-win32/lib
       ARPACK_LIB := $(PARPACK_LIBDIR)/parpack.lib
 endif
    ARPACK_LIBDIR ?= /usr/local/arpack-win32/lib
       ARPACK_LIB := $(ARPACK_LIBDIR)/arpack.lib
endif

ifdef USE_MPI
       MPI_INCDIR ?= c:\\work\\models\\MPICH2\\include
       MPI_LIBDIR ?= c:\\work\\models\\MPICH2\\lib
       LIBS_WIN32 += $(MPI_LIBDIR)\\fmpich2s.lib
           FFLAGS += -I$(MPI_INCDIR)
         CPPFLAGS += -DMPI -I$(MPI_INCDIR)
endif

ifdef USE_ESMF
          ESMF_OS ?= $(OS)
      ESMF_SUBDIR := $(ESMF_OS).$(ESMF_COMPILER).$(ESMF_ABI).$(ESMF_COMM).$(ESMF_SITE)
      ESMF_MK_DIR ?= $(ESMF_DIR)/lib/lib$(ESMF_BOPT)/$(ESMF_SUBDIR)
                     include $(ESMF_MK_DIR)/esmf.mk
           FFLAGS += $(ESMF_F90COMPILEPATHS)
       LIBS_WIN32 += $(ESMF_F90LINKPATHS) $(ESMF_F90ESMFLINKLIBS)
endif

ifndef USE_SCRIP
             LIBS += $(MCT_PARAMS_DIR)/mct_coupler_params.o
endif

ifdef USE_WW3
             FFLAGS += -I${COAWST_WW3_DIR}/mod_MPI
             LIBS += WW3/model/obj_MPI/libWW3.a
endif

ifdef USE_MCT
       MCT_LIBDIR ?= c:\\work\\models\\MCT_v2.2\\mct
      MPEU_LIBDIR ?= c:\\work\\models\\MCT_v2.2\\mpeu
       LIBS_WIN32 += $(MCT_LIBDIR)\\libmct.a $(MPEU_LIBDIR)\\libmpeu.a
         CPPFLAGS += -traditional-cpp
           FFLAGS += -I$(MCT_LIBDIR) -I$(MPEU_LIBDIR)
           FFLAGS += /noextend_source -assume:byterecl
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

#
# Compiler flags
#

ifdef USE_DEBUG
           FFLAGS += /debug:full /traceback /nopdbfile
else
           FFLAGS += /fast
endif

#
# For a Windows compiler, create variables pointing to the Windows
# file names needed when linking. Use of the "=" sign means that
# variables will be evaluated only when needed.
#

         BIN_WIN32 = "$$(cygpath --windows $(BIN))"
        LIBS_WIN32 = "$$(cygpath --windows $(NETCDF_LIB))"
ifdef USE_ARPACK
        LIBS_WIN32 += "$$(cygpath --windows $(ARPACK_LIB))"
endif

        LD_WINDOWS := on

#
# Use full path of compiler.
#
                FC := "$(shell which ${FC})"
                LD := $(FC)

#
# Set free form format in source files to allow long string for
# local directory and compilation flags inside the code.
#

$(SCRATCH_DIR)/mod_ncparam.o: FFLAGS += /free
$(SCRATCH_DIR)/mod_strings.o: FFLAGS += /free
$(SCRATCH_DIR)/analytical.o: FFLAGS += /free
$(SCRATCH_DIR)/biology.o: FFLAGS += /free
ifdef USE_ADJOINT
$(SCRATCH_DIR)/ad_biology.o: FFLAGS += /free
endif
ifdef USE_REPRESENTER
$(SCRATCH_DIR)/rp_biology.o: FFLAGS += /free
endif
ifdef USE_TANGENT
$(SCRATCH_DIR)/tl_biology.o: FFLAGS += /free
endif

#
# Supress free format in SWAN source files since there are comments
# beyond column 72.
#

ifdef USE_SWAN

$(SCRATCH_DIR)/ocpcre.o:   FFLAGS += /fixed
$(SCRATCH_DIR)/ocpids.o:   FFLAGS += /fixed
$(SCRATCH_DIR)/ocpmix.o:   FFLAGS += /fixed
$(SCRATCH_DIR)/swancom1.o: FFLAGS += /fixed
$(SCRATCH_DIR)/swancom2.o: FFLAGS += /fixed
$(SCRATCH_DIR)/swancom3.o: FFLAGS += /fixed
$(SCRATCH_DIR)/swancom4.o: FFLAGS += /fixed
$(SCRATCH_DIR)/swancom5.o: FFLAGS += /fixed
$(SCRATCH_DIR)/swanmain.o: FFLAGS += /fixed
$(SCRATCH_DIR)/swanout1.o: FFLAGS += /fixed
$(SCRATCH_DIR)/swanout2.o: FFLAGS += /fixed
$(SCRATCH_DIR)/swanparll.o:FFLAGS += /fixed
$(SCRATCH_DIR)/swanpre1.o: FFLAGS += /fixed
$(SCRATCH_DIR)/swanpre2.o: FFLAGS += /fixed
$(SCRATCH_DIR)/swanser.o:  FFLAGS += /fixed
$(SCRATCH_DIR)/swmod1.o:   FFLAGS += /fixed
$(SCRATCH_DIR)/swmod2.o:       FFLAGS += /fixed
$(SCRATCH_DIR)/SwanSpectPart.o:FFLAGS += /fixed
$(SCRATCH_DIR)/m_constants.o:  FFLAGS += /free
$(SCRATCH_DIR)/m_fileio.o:     FFLAGS += /free
$(SCRATCH_DIR)/mod_xnl4v5.o:   FFLAGS += /free
$(SCRATCH_DIR)/serv_xnl4v5.o:  FFLAGS += /free
$(SCRATCH_DIR)/nctablemd.o:    FFLAGS += /free
$(SCRATCH_DIR)/agioncmd.o:     FFLAGS += /free
$(SCRATCH_DIR)/swn_outnc.o:    FFLAGS += /free
$(SCRATCH_DIR)/SdsBabanin.o:   FFLAGS += /free
$(SCRATCH_DIR)/SwanBpntlist.o: FFLAGS += /free
$(SCRATCH_DIR)/SwanCheckGrid.o:FFLAGS += /free
$(SCRATCH_DIR)/SwanCompdata.o: FFLAGS += /free
$(SCRATCH_DIR)/SwanCompUnstruc.o:       FFLAGS += /free
$(SCRATCH_DIR)/SwanComputeForce.o:      FFLAGS += /free
$(SCRATCH_DIR)/SwanConvAccur.o:         FFLAGS += /free
$(SCRATCH_DIR)/SwanConvStopc.o:         FFLAGS += /free
$(SCRATCH_DIR)/SwanCreateEdges.o:       FFLAGS += /free
$(SCRATCH_DIR)/SwanCrossObstacle.o:     FFLAGS += /free
$(SCRATCH_DIR)/SwanDiffPar.o:           FFLAGS += /free
$(SCRATCH_DIR)/SwanDispParm.o:          FFLAGS += /free
$(SCRATCH_DIR)/SwanFindObstacles.o:     FFLAGS += /free
$(SCRATCH_DIR)/SwanFindPoint.o:         FFLAGS += /free
$(SCRATCH_DIR)/SwanGridCell.o:          FFLAGS += /free
$(SCRATCH_DIR)/SwanGriddata.o:          FFLAGS += /free
$(SCRATCH_DIR)/SwanGridFace.o:          FFLAGS += /free
$(SCRATCH_DIR)/SwanGridobjects.o:       FFLAGS += /free
$(SCRATCH_DIR)/SwanGridTopology.o:      FFLAGS += /free
$(SCRATCH_DIR)/SwanGridVert.o:          FFLAGS += /free
$(SCRATCH_DIR)/SwanGSECorr.o:           FFLAGS += /free
$(SCRATCH_DIR)/SwanInitCompGrid.o:      FFLAGS += /free
$(SCRATCH_DIR)/SwanInterpolateAc.o:     FFLAGS += /free
$(SCRATCH_DIR)/SwanInterpolateOutput.o: FFLAGS += /free
$(SCRATCH_DIR)/SwanInterpolatePoint.o:  FFLAGS += /free
$(SCRATCH_DIR)/SwanIntgratSpc.o:        FFLAGS += /free
$(SCRATCH_DIR)/SwanPointinMesh.o:       FFLAGS += /free
$(SCRATCH_DIR)/SwanPrepComp.o:          FFLAGS += /free
$(SCRATCH_DIR)/SwanPrintGridInfo.o:     FFLAGS += /free
$(SCRATCH_DIR)/SwanPropvelS.o:          FFLAGS += /free
$(SCRATCH_DIR)/SwanPropvelX.o:          FFLAGS += /free
$(SCRATCH_DIR)/SwanReadADCGrid.o:       FFLAGS += /free
$(SCRATCH_DIR)/SwanReadEasymeshGrid.o:  FFLAGS += /free
$(SCRATCH_DIR)/SwanReadGrid.o:          FFLAGS += /free
$(SCRATCH_DIR)/SwanReadTriangleGrid.o:  FFLAGS += /free
$(SCRATCH_DIR)/SwanSweepSel.o:          FFLAGS += /free
$(SCRATCH_DIR)/SwanThreadBounds.o:      FFLAGS += /free
$(SCRATCH_DIR)/SwanTranspAc.o:          FFLAGS += /free
$(SCRATCH_DIR)/SwanTranspX.o:           FFLAGS += /free
$(SCRATCH_DIR)/SwanVertlist.o:  FFLAGS += /free
$(SCRATCH_DIR)/waves_control.o: FFLAGS += /free
$(SCRATCH_DIR)/waves_coupler.o: FFLAGS += /free

endif

#
# For a Windows compiler, override the compilation rule
#

%.o: %.f90
	cd $(SCRATCH_DIR); $(FC) $(FFLAGS) /compile (notdir $<) /object:$(notdir $@)
