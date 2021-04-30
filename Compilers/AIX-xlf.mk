# svn $Id: AIX-xlf.mk 1054 2021-03-06 19:47:12Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2021 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Include file for IBM xlf95_r Fortran Compiler
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
               FC := xlf95_r
           FFLAGS := -qsuffix=f=f90 -qmaxmem=-1 -qarch=pwr4 -qtune=pwr4
              CPP := /usr/lib/cpp
         CPPFLAGS := -P
               CC := xlc_r
              CXX := xlC_r
           CFLAGS :=
         CXXFLAGS :=
             LIBS :=
ifdef USE_ROMS
             LIBS += $(SCRATCH_DIR)/libNLM.a         # cyclic dependencies
endif
          LDFLAGS :=
               AR := ar
          ARFLAGS := -r
            MKDIR := mkdir -p
               RM := rm -f
           RANLIB := ranlib
             PERL := perl
             TEST := test

        MDEPFLAGS := --cpp --fext=f90 --file=- --objdir=$(SCRATCH_DIR)

#
# Library locations, can be overridden by environment variables.
#

ifdef USE_LARGE
           FFLAGS += -q64 -I.
           CFLAGS += -q64 -I.
         CXXFLAGS += -q64 -I.
          ARFLAGS += -X 64
          LDFLAGS += -bmaxdata:0x200000000
else
          LDFLAGS += -bmaxdata:0x70000000
endif

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
   PARPACK_LIBDIR ?= /usr/local/lib
             LIBS += -L$(PARPACK_LIBDIR) -lparpack
 endif
    ARPACK_LIBDIR ?= /usr/local/lib
             LIBS += -L$(ARPACK_LIBDIR) -larpack
endif

ifdef USE_MPI
         CPPFLAGS += -DMPI
               FC := mpxlf95_r

endif

ifdef USE_OpenMP
         CPPFLAGS += -D_OPENMP
           FFLAGS += -qsmp=omp
endif

ifdef USE_DEBUG
           FFLAGS += -g -qfullpath -qflttrap=enable:zerodivide:invalid
           CFLAGS += -g -qfullpath
         CXXFLAGS += -g -qfullpath
else
           FFLAGS += -O3 -qstrict
           CFLAGS += -O2
         CXXFLAGS += -O2
endif

ifdef USE_ESMF
          ESMF_OS ?= $(OS)
      ESMF_SUBDIR := $(ESMF_OS).$(ESMF_COMPILER).$(ESMF_ABI).$(ESMF_COMM).$(ESMF_SITE)
      ESMF_MK_DIR ?= $(ESMF_DIR)/lib/lib$(ESMF_BOPT)/$(ESMF_SUBDIR)
                     include $(ESMF_MK_DIR)/esmf.mk
           FFLAGS += $(ESMF_F90COMPILEPATHS)
             LIBS += $(ESMF_F90LINKPATHS) $(ESMF_F90ESMFLINKLIBS)
endif

ifdef USE_CXX
             LIBS += -lstdc++
endif

ifndef USE_SCRIP
             LIBS += $(MCT_PARAMS_DIR)/mct_coupler_params.o
endif

ifdef USE_WW3
             FFLAGS += -I${COAWST_WW3_DIR}/mod_MPI
             LIBS += WW3/model/obj_MPI/libWW3.a
endif

ifdef USE_MCT
       MCT_INCDIR ?= /usr/local/pkg/mct/include
       MCT_LIBDIR ?= /usr/local/pkg/mct/lib
           FFLAGS += -I$(MCT_INCDIR)
             LIBS += -L$(MCT_LIBDIR) -lmct -lmpeu
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
# Use full path of compiler.
#
               FC := $(shell which ${FC})
               LD := $(FC)

#
# Set free form format in source files to allow long string for
# local directory and compilation flags inside the code.
#

$(SCRATCH_DIR)/mod_ncparam.o: FFLAGS += -qfree
$(SCRATCH_DIR)/mod_strings.o: FFLAGS += -qfree
$(SCRATCH_DIR)/analytical.o: FFLAGS += -qfree
$(SCRATCH_DIR)/biology.o: FFLAGS += -qfree
ifdef USE_ADJOINT
$(SCRATCH_DIR)/ad_biology.o: FFLAGS += -qfree
endif
ifdef USE_REPRESENTER
$(SCRATCH_DIR)/rp_biology.o: FFLAGS += -qfree
endif
ifdef USE_TANGENT
$(SCRATCH_DIR)/tl_biology.o: FFLAGS += -qfree
endif

#
# Supress free format in SWAN source files since there are comments
# beyond column 72.
#

ifdef USE_SWAN

$(SCRATCH_DIR)/ocpcre.o:   FFLAGS += -qfixed
$(SCRATCH_DIR)/ocpids.o:   FFLAGS += -qfixed
$(SCRATCH_DIR)/ocpmix.o:   FFLAGS += -qfixed
$(SCRATCH_DIR)/swancom1.o: FFLAGS += -qfixed
$(SCRATCH_DIR)/swancom2.o: FFLAGS += -qfixed
$(SCRATCH_DIR)/swancom3.o: FFLAGS += -qfixed
$(SCRATCH_DIR)/swancom4.o: FFLAGS += -qfixed
$(SCRATCH_DIR)/swancom5.o: FFLAGS += -qfixed
$(SCRATCH_DIR)/swanmain.o: FFLAGS += -qfixed
$(SCRATCH_DIR)/swanout1.o: FFLAGS += -qfixed
$(SCRATCH_DIR)/swanout2.o: FFLAGS += -qfixed
$(SCRATCH_DIR)/swanparll.o:FFLAGS += -qfixed
$(SCRATCH_DIR)/swanpre1.o: FFLAGS += -qfixed
$(SCRATCH_DIR)/swanpre2.o: FFLAGS += -qfixed
$(SCRATCH_DIR)/swanser.o:  FFLAGS += -qfixed
$(SCRATCH_DIR)/swmod1.o:   FFLAGS += -qfixed
$(SCRATCH_DIR)/swmod2.o:       FFLAGS += -qfixed
$(SCRATCH_DIR)/SwanSpectPart.o:FFLAGS += -qfixed
$(SCRATCH_DIR)/m_constants.o:  FFLAGS += -qfree
$(SCRATCH_DIR)/m_fileio.o:     FFLAGS += -qfree
$(SCRATCH_DIR)/mod_xnl4v5.o:   FFLAGS += -qfree
$(SCRATCH_DIR)/serv_xnl4v5.o:  FFLAGS += -qfree
$(SCRATCH_DIR)/nctablemd.o:    FFLAGS += -qfree
$(SCRATCH_DIR)/agioncmd.o:     FFLAGS += -qfree
$(SCRATCH_DIR)/swn_outnc.o:    FFLAGS += -qfree
$(SCRATCH_DIR)/SdsBabanin.o:   FFLAGS += -qfree
$(SCRATCH_DIR)/SwanBpntlist.o: FFLAGS += -qfree
$(SCRATCH_DIR)/SwanCheckGrid.o:FFLAGS += -qfree
$(SCRATCH_DIR)/SwanCompdata.o: FFLAGS += -qfree
$(SCRATCH_DIR)/SwanCompUnstruc.o:       FFLAGS += -qfree
$(SCRATCH_DIR)/SwanComputeForce.o:      FFLAGS += -qfree
$(SCRATCH_DIR)/SwanConvAccur.o:         FFLAGS += -qfree
$(SCRATCH_DIR)/SwanConvStopc.o:         FFLAGS += -qfree
$(SCRATCH_DIR)/SwanCreateEdges.o:       FFLAGS += -qfree
$(SCRATCH_DIR)/SwanCrossObstacle.o:     FFLAGS += -qfree
$(SCRATCH_DIR)/SwanDiffPar.o:           FFLAGS += -qfree
$(SCRATCH_DIR)/SwanDispParm.o:          FFLAGS += -qfree
$(SCRATCH_DIR)/SwanFindObstacles.o:     FFLAGS += -qfree
$(SCRATCH_DIR)/SwanFindPoint.o:         FFLAGS += -qfree
$(SCRATCH_DIR)/SwanGridCell.o:          FFLAGS += -qfree
$(SCRATCH_DIR)/SwanGriddata.o:          FFLAGS += -qfree
$(SCRATCH_DIR)/SwanGridFace.o:          FFLAGS += -qfree
$(SCRATCH_DIR)/SwanGridobjects.o:       FFLAGS += -qfree
$(SCRATCH_DIR)/SwanGridTopology.o:      FFLAGS += -qfree
$(SCRATCH_DIR)/SwanGridVert.o:          FFLAGS += -qfree
$(SCRATCH_DIR)/SwanGSECorr.o:           FFLAGS += -qfree
$(SCRATCH_DIR)/SwanInitCompGrid.o:      FFLAGS += -qfree
$(SCRATCH_DIR)/SwanInterpolateAc.o:     FFLAGS += -qfree
$(SCRATCH_DIR)/SwanInterpolateOutput.o: FFLAGS += -qfree
$(SCRATCH_DIR)/SwanInterpolatePoint.o:  FFLAGS += -qfree
$(SCRATCH_DIR)/SwanIntgratSpc.o:        FFLAGS += -qfree
$(SCRATCH_DIR)/SwanPointinMesh.o:       FFLAGS += -qfree
$(SCRATCH_DIR)/SwanPrepComp.o:          FFLAGS += -qfree
$(SCRATCH_DIR)/SwanPrintGridInfo.o:     FFLAGS += -qfree
$(SCRATCH_DIR)/SwanPropvelS.o:          FFLAGS += -qfree
$(SCRATCH_DIR)/SwanPropvelX.o:          FFLAGS += -qfree
$(SCRATCH_DIR)/SwanReadADCGrid.o:       FFLAGS += -qfree
$(SCRATCH_DIR)/SwanReadEasymeshGrid.o:  FFLAGS += -qfree
$(SCRATCH_DIR)/SwanReadGrid.o:          FFLAGS += -qfree
$(SCRATCH_DIR)/SwanReadTriangleGrid.o:  FFLAGS += -qfree
$(SCRATCH_DIR)/SwanSweepSel.o:          FFLAGS += -qfree
$(SCRATCH_DIR)/SwanThreadBounds.o:      FFLAGS += -qfree
$(SCRATCH_DIR)/SwanTranspAc.o:          FFLAGS += -qfree
$(SCRATCH_DIR)/SwanTranspX.o:           FFLAGS += -qfree
$(SCRATCH_DIR)/SwanVertlist.o:  FFLAGS += -qfree
$(SCRATCH_DIR)/waves_control.o: FFLAGS += -qfree
$(SCRATCH_DIR)/waves_coupler.o: FFLAGS += -qfree

endif
