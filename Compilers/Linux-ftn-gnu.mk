# svn $Id: Linux-ftn-gnu.mk 1054 2021-03-06 19:47:12Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2021 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Include file for CRAY ftn compiler with PrgEnv-gnu
# -------------------------------------------------------------------------
#
# ARPACK_LIBDIR  ARPACK libary directory
# FC             Name of the fortran compiler to use
# FFLAGS         Flags to the fortran compiler
# CPP            Name of the C-preprocessor
# CPPFLAGS       Flags to the C-preprocessor
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
               FC := ftn
           FFLAGS := -frepack-arrays
       FIXEDFLAGS := -ffixed-form
        FREEFLAGS := -ffree-form -ffree-line-length-none
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -traditional -w
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
          ARFLAGS := -r
            MKDIR := mkdir -p
               CP := cp -p -v
               RM := rm -f
           RANLIB := ranlib
             PERL := perl
             TEST := test

#--------------------------------------------------------------------------
# Compiling flags for ROMS Applications.
#--------------------------------------------------------------------------

ifdef USE_ROMS
 ifdef USE_DEBUG
           FFLAGS += -g
           FFLAGS += -fbounds-check
           FFLAGS += -fbacktrace
           FFLAGS += -fcheck=all
#          FFLAGS += -fsanitize=address -fsanitize=undefined
           FFLAGS += -finit-real=nan -ffpe-trap=invalid,zero,overflow
 else
           FFLAGS += -O3
           FFLAGS += # -ffast-math ==== bad !!!!
           FFLAGS += -ftree-vectorize -ftree-loop-linear -funroll-loops -w -ffree-form -ffree-line-length-none -frecord-marker=4 -fconvert=big-endian
##                   -fconvert=big-endian
 endif
endif
        MDEPFLAGS := --cpp --fext=f90 --file=- --objdir=$(SCRATCH_DIR)

#--------------------------------------------------------------------------
# Compiling flags for CICE Applications.
#--------------------------------------------------------------------------

ifdef CICE_APPLICATION
          CPPDEFS := -DLINUS $(MY_CPP_FLAGS)
 ifdef USE_DEBUG
           FFLAGS += -g
           FFLAGS += -fbounds-check
 else
           FFLAGS := -O3 -w
 endif
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


#--------------------------------------------------------------------------
# Library locations, can be overridden by environment variables.
#--------------------------------------------------------------------------

          LDFLAGS := $(FFLAGS)

ifdef USE_NETCDF4
        NF_CONFIG ?= nf-config
    NETCDF_INCDIR ?= $(shell $(NF_CONFIG) --prefix)/include
             LIBS += $(shell $(NF_CONFIG) --flibs)
           INCDIR += $(NETCDF_INCDIR) $(INCDIR)
else
    NETCDF_INCDIR ?= /opt/cray/pe/netcdf/default/gnu/7.1/include
    NETCDF_LIBDIR ?= /opt/cray/pe/netcdf/default/gnu/7.1/lib
      NETCDF_LIBS ?= -lnetcdf
             LIBS += -L$(NETCDF_LIBDIR) $(NETCDF_LIBS)
           INCDIR += $(NETCDF_INCDIR) $(INCDIR)
endif

ifdef USE_HDF5
      HDF5_INCDIR ?= /opt/cray/hdf5/default/gnu/7.1/include
      HDF5_LIBDIR ?= /opt/cray/hdf5/default/gnu/7.1/lib
        HDF5_LIBS ?= -lhdf5_fortran -lhdf5hl_fortran -lhdf5 -lz
             LIBS += -L$(HDF5_LIBDIR) $(HDF5_LIBS)
           INCDIR += $(HDF5_INCDIR)
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
           FFLAGS += -fopenmp
endif


ifdef USE_ESMF
          ESMF_OS ?= $(OS)
      ESMF_SUBDIR := $(ESMF_OS).$(ESMF_COMPILER).$(ESMF_ABI).$(ESMF_COMM).$(ESMF_SITE)
      ESMF_MK_DIR ?= $(ESMF_DIR)/lib/lib$(ESMF_BOPT)/$(ESMF_SUBDIR)
                     include $(ESMF_MK_DIR)/esmf.mk
           FFLAGS += $(ESMF_F90COMPILEPATHS)
             LIBS += $(ESMF_F90LINKPATHS) $(ESMF_F90ESMFLINKLIBS)
endif

# Use full path of compiler.

               FC := $(shell which ${FC})
               LD := $(FC)

ifndef USE_SCRIP
             LIBS += $(MCT_PARAMS_DIR)/mct_coupler_params.o
endif

ifdef USE_WW3
             FFLAGS += -I${COAWST_WW3_DIR}/mod_MPI
             LIBS += WW3/model/obj_MPI/libWW3.a
endif

ifdef USE_MCT
       MCT_INCDIR ?= /usr/local/mct/include
       MCT_LIBDIR ?= /usr/local/mct/lib
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

#--------------------------------------------------------------------------
# ROMS specific rules.
#--------------------------------------------------------------------------

# Turn off bounds checking for function def_var, as "dimension(*)"
# declarations confuse Gnu Fortran 95 bounds-checking code.

ifdef USE_ROMS
 $(SCRATCH_DIR)/def_var.o: FFLAGS += -fno-bounds-check
endif

# Allow integer overflow in ran_state.F.  This is not allowed
# during -O3 optimization. This option should be applied only for
# Gfortran versions >= 4.2.

ifdef USE_ROMS
 FC_TEST := $(findstring $(shell ${FC} --version | head -1 | \
                              awk '{ sub("Fortran 95", "Fortran"); print }' | \
                              cut -d " " -f 4 | \
                              cut -d "." -f 1-2), \
             4.0 4.1)

 ifeq "${FC_TEST}" ""
  $(SCRATCH_DIR)/ran_state.o: FFLAGS += -fno-strict-overflow
 endif
endif

# Set free form format in some ROMS source files to allow long string for
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

#
# Supress free format in SWAN source files since there are comments
# beyond column 72.

ifdef USE_SWAN

$(SCRATCH_DIR)/ocpcre.o:   FFLAGS += -fixed-form
$(SCRATCH_DIR)/ocpids.o:   FFLAGS += -fixed-form
$(SCRATCH_DIR)/ocpmix.o:   FFLAGS += -fixed-form
$(SCRATCH_DIR)/swancom1.o: FFLAGS += -fixed-form
$(SCRATCH_DIR)/swancom2.o: FFLAGS += -fixed-form
$(SCRATCH_DIR)/swancom3.o: FFLAGS += -fixed-form
$(SCRATCH_DIR)/swancom4.o: FFLAGS += -fixed-form
$(SCRATCH_DIR)/swancom5.o: FFLAGS += -fixed-form
$(SCRATCH_DIR)/swanmain.o: FFLAGS += -fixed-form
$(SCRATCH_DIR)/swanout1.o: FFLAGS += -fixed-form
$(SCRATCH_DIR)/swanout2.o: FFLAGS += -fixed-form
$(SCRATCH_DIR)/swanparll.o:FFLAGS += -fixed-form
$(SCRATCH_DIR)/swanpre1.o: FFLAGS += -fixed-form
$(SCRATCH_DIR)/swanpre2.o: FFLAGS += -fixed-form
$(SCRATCH_DIR)/swanser.o:  FFLAGS += -fixed-form
$(SCRATCH_DIR)/swmod1.o:   FFLAGS += -fixed-form
$(SCRATCH_DIR)/swmod2.o:       FFLAGS += -fixed-form
$(SCRATCH_DIR)/SwanSpectPart.o:FFLAGS += -fixed-form
$(SCRATCH_DIR)/m_constants.o:  FFLAGS += -free-form
$(SCRATCH_DIR)/m_fileio.o:     FFLAGS += -free-form
$(SCRATCH_DIR)/mod_xnl4v5.o:   FFLAGS += -free-form
$(SCRATCH_DIR)/serv_xnl4v5.o:  FFLAGS += -free-form
$(SCRATCH_DIR)/nctablemd.o:    FFLAGS += -free-form
$(SCRATCH_DIR)/agioncmd.o:     FFLAGS += -free-form
$(SCRATCH_DIR)/swn_outnc.o:    FFLAGS += -free-form
$(SCRATCH_DIR)/SdsBabanin.o:   FFLAGS += -free-form
$(SCRATCH_DIR)/SwanBpntlist.o: FFLAGS += -free-form
$(SCRATCH_DIR)/SwanCheckGrid.o:FFLAGS += -free-form
$(SCRATCH_DIR)/SwanCompdata.o: FFLAGS += -free-form
$(SCRATCH_DIR)/SwanCompUnstruc.o:       FFLAGS += -free-form
$(SCRATCH_DIR)/SwanComputeForce.o:      FFLAGS += -free-form
$(SCRATCH_DIR)/SwanConvAccur.o:         FFLAGS += -free-form
$(SCRATCH_DIR)/SwanConvStopc.o:         FFLAGS += -free-form
$(SCRATCH_DIR)/SwanCreateEdges.o:       FFLAGS += -free-form
$(SCRATCH_DIR)/SwanCrossObstacle.o:     FFLAGS += -free-form
$(SCRATCH_DIR)/SwanDiffPar.o:           FFLAGS += -free-form
$(SCRATCH_DIR)/SwanDispParm.o:          FFLAGS += -free-form
$(SCRATCH_DIR)/SwanFindObstacles.o:     FFLAGS += -free-form
$(SCRATCH_DIR)/SwanFindPoint.o:         FFLAGS += -free-form
$(SCRATCH_DIR)/SwanGridCell.o:          FFLAGS += -free-form
$(SCRATCH_DIR)/SwanGriddata.o:          FFLAGS += -free-form
$(SCRATCH_DIR)/SwanGridFace.o:          FFLAGS += -free-form
$(SCRATCH_DIR)/SwanGridobjects.o:       FFLAGS += -free-form
$(SCRATCH_DIR)/SwanGridTopology.o:      FFLAGS += -free-form
$(SCRATCH_DIR)/SwanGridVert.o:          FFLAGS += -free-form
$(SCRATCH_DIR)/SwanGSECorr.o:           FFLAGS += -free-form
$(SCRATCH_DIR)/SwanInitCompGrid.o:      FFLAGS += -free-form
$(SCRATCH_DIR)/SwanInterpolateAc.o:     FFLAGS += -free-form
$(SCRATCH_DIR)/SwanInterpolateOutput.o: FFLAGS += -free-form
$(SCRATCH_DIR)/SwanInterpolatePoint.o:  FFLAGS += -free-form
$(SCRATCH_DIR)/SwanIntgratSpc.o:        FFLAGS += -free-form
$(SCRATCH_DIR)/SwanPointinMesh.o:       FFLAGS += -free-form
$(SCRATCH_DIR)/SwanPrepComp.o:          FFLAGS += -free-form
$(SCRATCH_DIR)/SwanPrintGridInfo.o:     FFLAGS += -free-form
$(SCRATCH_DIR)/SwanPropvelS.o:          FFLAGS += -free-form
$(SCRATCH_DIR)/SwanPropvelX.o:          FFLAGS += -free-form
$(SCRATCH_DIR)/SwanReadADCGrid.o:       FFLAGS += -free-form
$(SCRATCH_DIR)/SwanReadEasymeshGrid.o:  FFLAGS += -free-form
$(SCRATCH_DIR)/SwanReadGrid.o:          FFLAGS += -free-form
$(SCRATCH_DIR)/SwanReadTriangleGrid.o:  FFLAGS += -free-form
$(SCRATCH_DIR)/SwanSweepSel.o:          FFLAGS += -free-form
$(SCRATCH_DIR)/SwanThreadBounds.o:      FFLAGS += -free-form
$(SCRATCH_DIR)/SwanTranspAc.o:          FFLAGS += -free-form
$(SCRATCH_DIR)/SwanTranspX.o:           FFLAGS += -free-form
$(SCRATCH_DIR)/SwanVertlist.o:  FFLAGS += -freeform
$(SCRATCH_DIR)/waves_control.o: FFLAGS += -freeform
$(SCRATCH_DIR)/waves_coupler.o: FFLAGS += -freeform

endif
