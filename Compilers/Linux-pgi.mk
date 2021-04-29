# svn $Id: Linux-pgi.mk 1054 2021-03-06 19:47:12Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2021 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Include file for PGI Fortran compiler on Linux
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
# HDF5_INCDIR    HDF5 include directory
# HDF5_LIBDIR    HDF5 library directory
# HDF5_LIBS      HDF5 library switches
# LIBS           Required libraries during linking
# NF_CONFIG      NetCDF Fortran configuration script
# NETCDF_INCDIR  NetCDF include directory
# NETCDF_LIBDIR  NetCDF library directory
# NETCDF_LIBS    NetCDF library switches
# LD             Program to load the objects into an executable
# LDFLAGS        Flags to the loader
# RANLIB         Name of ranlib command
# MDEPFLAGS      Flags for sfmakedepend  (-s if you keep .f files)
#
# First the defaults
#
               FC := pgf90
           FFLAGS :=
       FIXEDFLAGS := -Mnofree
        FREEFLAGS := -Mfree
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -traditional -w              # -w turns of warnings
ifdef USE_DEBUG
         CPPFLAGS += -DUSE_DEBUG
endif
               CC := gcc
              CXX := g++
           CFLAGS :=
         CXXFLAGS :=
           INCDIR := /usr/include /usr/local/bin
            SLIBS := -L/usr/local/lib -L/usr/lib
            ULIBS :=
             LIBS :=
ifdef USE_ROMS
             LIBS += $(SCRATCH_DIR)/libNLM.a         # cyclic dependencies
endif
       MOD_SUFFIX := mod
               LD := $(FC)
          LDFLAGS :=
               AR := ar
          ARFLAGS := r
            MKDIR := mkdir -p
               CP := cp -p -v
               RM := rm -f
           RANLIB := ranlib
             PERL := perl
             TEST := test

#--------------------------------------------------------------------------
# Compiling flags for ROMS Applications.
#--------------------------------------------------------------------------
#
# Perform floating-point operations in strict conformance with the
# IEEE standard by using -Kieee. This may slow down computations
# because some optimizations are disabled.  However, we noticed a
# speed-up. The user may want to have this option to allow similar,
# if not identical solutions between different of the PGI compiler.

ifdef USE_ROMS
 ifdef USE_DEBUG
           FFLAGS += -g
           FFLAGS += -C
           FFLAGS += -Mchkstk -Mchkfpstk
 else
#          FFLAGS += -Bstatic
           FFLAGS += -fastsse -Mipa=fast
 endif
           FFLAGS += -Kieee
endif
        MDEPFLAGS := --cpp --fext=f90 --file=- --objdir=$(SCRATCH_DIR)

#--------------------------------------------------------------------------
# Compiling flags for CICE Applications.
#--------------------------------------------------------------------------

ifdef CICE_APPLICATION
          CPPDEFS := -DLINUS $(MY_CPP_FLAGS)
 ifdef USE_DEBUG
           FFLAGS += -g
           FFLAGS += -C
           FFLAGS += -Mchkstk -Mchkfpstk
 else
           FFLAGS += -fastsse -Mipa=fast
 endif
           FFLAGS += -Kieee
endif

#--------------------------------------------------------------------------
# Coupled models.  Notice Linux needs the libraries repeated for
# dependencies for some of the coupled components.
#--------------------------------------------------------------------------

ifdef USE_COAMPS
             LIBS += $(COAMPS_LIB_DIR)/coamps_driver.a
             LIBS += $(COAMPS_LIB_DIR)/libaa.a
             LIBS += $(COAMPS_LIB_DIR)/libam.a
             LIBS += $(COAMPS_LIB_DIR)/libashare.a
             LIBS += $(COAMPS_LIB_DIR)/libcoamps.a
             LIBS += $(COAMPS_LIB_DIR)/libfnoc.a
             LIBS += $(COAMPS_LIB_DIR)/libaa.a
             LIBS += $(COAMPS_LIB_DIR)/libam.a
             LIBS += $(COAMPS_LIB_DIR)/libashare.a
             LIBS += $(COAMPS_LIB_DIR)/libcoamps.a
             LIBS += $(COAMPS_LIB_DIR)/libfnoc.a
             LIBS += $(COAMPS_LIB_DIR)/libfishpak.a
             LIBS += $(COAMPS_LIB_DIR)/libtracer.a
endif

ifdef CICE_APPLICATION
            SLIBS += $(SLIBS) $(LIBS)
endif
# Library locations, can be overridden by environment variables.
#--------------------------------------------------------------------------
#
# According to the PGI manual, the -Bstatic flags initializes
# the symbol table with -Bstatic, which is undefined for the linker.
# An undefined symbol triggers loading of the first member of an
# archive library. The -u flag fails with version 7.x of the compiler
# because it expects an argument.

          LDFLAGS := $(FFLAGS)

ifdef USE_NETCDF4
        NF_CONFIG ?= nf-config
    NETCDF_INCDIR ?= $(shell $(NF_CONFIG) --prefix)/include
             LIBS += $(shell $(NF_CONFIG) --flibs)
           INCDIR += $(NETCDF_INCDIR) $(INCDIR)
else
    NETCDF_INCDIR ?= /opt/pgisoft/serial/netcdf3/include
    NETCDF_LIBDIR ?= /opt/pgisoft/serial/netcdf3/lib
      NETCDF_LIBS ?= -lnetcdf
             LIBS += -L$(NETCDF_LIBDIR) $(NETCDF_LIBS)
           INCDIR += $(NETCDF_INCDIR) $(INCDIR)
endif

ifdef USE_HDF5
      HDF5_INCDIR ?= /opt/pgisoft/serial/hdf5/include
      HDF5_LIBDIR ?= /opt/pgisoft/serial/hdf5/lib
        HDF5_LIBS ?= -lhdf5_fortran -lhdf5hl_fortran -lhdf5 -lz
             LIBS += -L$(HDF5_LIBDIR) $(HDF5_LIBS)
           INCDIR += $(HDF5_INCDIR)
endif

ifdef USE_ARPACK
 ifdef USE_MPI
   PARPACK_LIBDIR ?= /opt/pgisoft/PARPACK
             LIBS += -L$(PARPACK_LIBDIR) -lparpack
 endif
    ARPACK_LIBDIR ?= /opt/pgisoft/ARPACK
             LIBS += -L$(ARPACK_LIBDIR) -larpack
endif

ifdef USE_MPI
         CPPFLAGS += -DMPI
 ifdef USE_MPIF90
               FC := mpif90
 else
             LIBS += -Bdynamic -lfmpi-pgi -lmpi-pgi -Bstatic
 endif
endif

ifdef USE_OpenMP
         CPPFLAGS += -D_OPENMP
           FFLAGS += -mp
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
       MCT_INCDIR ?= /opt/pgisoft/mct/include
       MCT_LIBDIR ?= /opt/pgisoft/mct/lib
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

       clean_list += ifc* work.pc*

#
# Use full path of compiler.
#
               FC := $(shell which ${FC})
               LD := $(FC)

#
# Set free form format in source files to allow long string for
# local directory and compilation flags inside the code.

ifdef USE_ROMS
 $(SCRATCH_DIR)/mod_ncparam.o: FFLAGS += $(FREEFLAGS)
 $(SCRATCH_DIR)/mod_strings.o: FFLAGS += $(FREEFLAGS)
 $(SCRATCH_DIR)/analytical.o: FFLAGS += $(FREEFLAGS)
 $(SCRATCH_DIR)/biology.o: FFLAGS += $(FREEFLAGS)

 ifdef USE_ADJOINT
  $(SCRATCH_DIR)/ad_biology.o: FFLAGS += $(FREEFLAGS)
 endif
 ifdef USE_REPRESENTER
  $(SCRATCH_DIR)/rp_biology.o: FFLAGS += $(FREEFLAGS)
 endif
 ifdef USE_TANGENT
  $(SCRATCH_DIR)/tl_biology.o: FFLAGS += $(FREEFLAGS)
 endif
endif

#--------------------------------------------------------------------------
# Model coupling specific rules.
#--------------------------------------------------------------------------

# Add COAMPS library directory to include path of ESMF coupling files.

ifdef USE_COAMPS
 $(SCRATCH_DIR)/esmf_atm.o: FFLAGS += -I$(COAMPS_LIB_DIR)
 $(SCRATCH_DIR)/esmf_esm.o: FFLAGS += -I$(COAMPS_LIB_DIR)
endif


# Supress free format in SWAN source files since there are comments
# beyond column 72.

ifdef USE_SWAN

$(SCRATCH_DIR)/ocpcre.o:   FFLAGS += -Mnofree
$(SCRATCH_DIR)/ocpids.o:   FFLAGS += -Mnofree
$(SCRATCH_DIR)/ocpmix.o:   FFLAGS += -Mnofree
$(SCRATCH_DIR)/swancom1.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swancom2.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swancom3.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swancom4.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swancom5.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swanmain.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swanout1.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swanout2.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swanparll.o:FFLAGS += -Mnofree
$(SCRATCH_DIR)/swanpre1.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swanpre2.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swanser.o:  FFLAGS += -Mnofree
$(SCRATCH_DIR)/swmod1.o:   FFLAGS += -Mnofree
$(SCRATCH_DIR)/swmod2.o:       FFLAGS += -Mnofree
$(SCRATCH_DIR)/SwanSpectPart.o:FFLAGS += -Mnofree
$(SCRATCH_DIR)/m_constants.o:  FFLAGS += -Mfree
$(SCRATCH_DIR)/m_fileio.o:     FFLAGS += -Mfree
$(SCRATCH_DIR)/mod_xnl4v5.o:   FFLAGS += -Mfree
$(SCRATCH_DIR)/serv_xnl4v5.o:  FFLAGS += -Mfree
$(SCRATCH_DIR)/nctablemd.o:    FFLAGS += -Mfree
$(SCRATCH_DIR)/agioncmd.o:     FFLAGS += -Mfree
$(SCRATCH_DIR)/swn_outnc.o:    FFLAGS += -Mfree
$(SCRATCH_DIR)/SdsBabanin.o:   FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanBpntlist.o: FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanCheckGrid.o:FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanCompdata.o: FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanCompUnstruc.o:       FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanComputeForce.o:      FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanConvAccur.o:         FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanConvStopc.o:         FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanCreateEdges.o:       FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanCrossObstacle.o:     FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanDiffPar.o:           FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanDispParm.o:          FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanFindObstacles.o:     FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanFindPoint.o:         FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanGridCell.o:          FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanGriddata.o:          FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanGridFace.o:          FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanGridobjects.o:       FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanGridTopology.o:      FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanGridVert.o:          FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanGSECorr.o:           FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanInitCompGrid.o:      FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanInterpolateAc.o:     FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanInterpolateOutput.o: FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanInterpolatePoint.o:  FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanIntgratSpc.o:        FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanPointinMesh.o:       FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanPrepComp.o:          FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanPrintGridInfo.o:     FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanPropvelS.o:          FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanPropvelX.o:          FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanReadADCGrid.o:       FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanReadEasymeshGrid.o:  FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanReadGrid.o:          FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanReadTriangleGrid.o:  FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanSweepSel.o:          FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanThreadBounds.o:      FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanTranspAc.o:          FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanTranspX.o:           FFLAGS += -Mfree
$(SCRATCH_DIR)/SwanVertlist.o:  FFLAGS += -Mfree
$(SCRATCH_DIR)/waves_control.o: FFLAGS += -Mfree
$(SCRATCH_DIR)/waves_coupler.o: FFLAGS += -Mfree

endif
