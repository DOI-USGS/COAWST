# git $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2026 The ROMS Group                                :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Include file for Intel IFORT compiler on Linux
# -------------------------------------------------------------------------
#
# ARPACK_LIBDIR  ARPACK library directory
# FC             Name of the fortran compiler to use
# FFLAGS         Flags to the fortran compiler
# CPP            Name of the C-preprocessor
# CPPFLAGS       Flags to the C-preprocessor
# HDF5_INCDIR    HDF5 include directory
# HDF5_LIBDIR    HDF5 library directory
# HDF5_LIBS      HDF5 library switches
# LIBS           Required libraries during linking
# ROMS_LIB       Directory and name for ROMS library
# NF_CONFIG      NetCDF Fortran configuration script
# NETCDF_INCDIR  NetCDF include directory
# NETCDF_LIBDIR  NetCDF library directory
# NETCDF_LIBS    NetCDF library switches
# PIO_INCDIR     Parallel-IO (PIO) library include directory
# PIO_LIBDIR     Parallel-IO (PIO) library directory
# PIO_LIBS       Parallel-IO (PIO) library switches
# PNETCDF_INCDIR PNetCDF include directory
# PNETCDF_LIBDIR PNetCDF library directory
# PNETCDF_LIBS   PNetCDF library switches

# LD             Program to load the objects into an executable or shared library
# LDFLAGS        Flags to the loader
# RANLIB         Name of ranlib command
# MDEPFLAGS      Flags for sfmakedepend  (-s if you keep .f files)
#
# First the defaults
#
               FC := ifort
           FFLAGS := -fp-model precise
           FFLAGS += -fc=ifort
           FFLAGS += -heap-arrays
       FIXEDFLAGS := -nofree
        FREEFLAGS := -free
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -traditional-cpp -w          # -w turns off warnings
           INCDIR := /usr/include /usr/local/bin
            SLIBS := -L/usr/local/lib -L/usr/lib
            ULIBS :=
             LIBS :=
         ROMS_LIB := -L$(BUILD_DIR) -lROMS
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
      ST_LIB_NAME := libROMS.a
      SH_LIB_NAME := libROMS.so

#--------------------------------------------------------------------------
# Compiling flags for ROMS Applications.
#--------------------------------------------------------------------------

ifdef USE_ROMS
 ifdef USE_DEBUG
           FFLAGS += -g
           FFLAGS += -check all
           FFLAGS += -check bounds
           FFLAGS += -traceback
           FFLAGS += -check uninit
#          FFLAGS += -warn interfaces,nouncalled -gen-interfaces
           FFLAGS += -gen-interfaces
 else
           FFLAGS += -ip -O3
           FFLAGS += -traceback
#          FFLAGS += -check uninit
 endif
 ifdef SHARED
          LDFLAGS += -Wl,-rpath,$(BUILD_DIR)

           FFLAGS += -fPIC
       SH_LDFLAGS += -shared
 endif
endif
        MDEPFLAGS := --cpp --fext=f90 --file=- --objdir=$(BUILD_DIR)

#--------------------------------------------------------------------------
# Compiling flags for CICE Applications.
#--------------------------------------------------------------------------

ifdef CICE_APPLICATION
          CPPDEFS := -DLINUS $(MY_CPP_FLAGS)
 ifdef USE_DEBUG
           FFLAGS := -g
#          FFLAGS += -O2
#          FFLAGS += -r8 -i4 -align all -w
           FFLAGS += -check bounds
           FFLAGS += -traceback
           FFLAGS += -check uninit
           FFLAGS += -ftz -convert big_endian -assume byterecl
           FFLAGS += -warn interfaces,nouncalled
           FFLAGS += -gen-interfaces
 else
           FFLAGS := -r8 -i4 -O2 -align all -w
           FFLAGS += -ftz -convert big_endian -assume byterecl
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


ifdef USE_PIO
       PIO_INCDIR ?= /opt/intelsoft/openmpi/pio/include
       PIO_LIBDIR ?= /opt/intelsoft/openmpi/pio/lib
           FFLAGS += -I$(PIO_INCDIR)
             LIBS += -L$(PIO_LIBDIR) -lpiof -lpioc

   PNETCDF_INCDIR ?= /opt/intelsoft/openmpi/pnetcdf/include
   PNETCDF_LIBDIR ?= /opt/intelsoft/openmpi/pnetcdf/lib
           FFLAGS += -I$(PNETCDF_INCDIR)
             LIBS += -L$(PNETCDF_LIBDIR) -lpnetcdf
endif

ifdef USE_NETCDF4
        NC_CONFIG ?= nc-config
   TEST_NC_CONFIG := $(shell which $(NC_CONFIG))
  ifneq ($(TEST_NC_CONFIG),)
             LIBS += $(shell $(NC_CONFIG) --libs)
  endif
        NF_CONFIG ?= nf-config
    NETCDF_INCDIR ?= $(shell $(NF_CONFIG) --prefix)/include
             LIBS += $(shell $(NF_CONFIG) --flibs)
           INCDIR += $(NETCDF_INCDIR) $(INCDIR)
else
    NETCDF_INCDIR ?= /opt/intelsoft/serial/netcdf3/include
    NETCDF_LIBDIR ?= /opt/intelsoft/serial/netcdf3/lib
      NETCDF_LIBS ?= -lnetcdf
             LIBS += -L$(NETCDF_LIBDIR) $(NETCDF_LIBS)
           INCDIR += $(NETCDF_INCDIR) $(INCDIR)
endif

ifdef USE_HDF5
      HDF5_INCDIR ?= /opt/intelsoft/serial/hdf5/include
      HDF5_LIBDIR ?= /opt/intelsoft/serial/hdf5/lib
        HDF5_LIBS ?= -lhdf5_fortran -lhdf5hl_fortran -lhdf5 -lz
             LIBS += -L$(HDF5_LIBDIR) $(HDF5_LIBS)
           INCDIR += $(HDF5_INCDIR)
endif

ifdef USE_ARPACK
 ifdef USE_MPI
   PARPACK_LIBDIR ?= /opt/intelsoft/PARPACK
             LIBS += -L$(PARPACK_LIBDIR) -lparpack
 endif
    ARPACK_LIBDIR ?= /opt/intelsoft/ARPACK
             LIBS += -L$(ARPACK_LIBDIR) -larpack
endif

ifdef USE_MPI
         CPPFLAGS += -DMPI
 ifdef USE_MPIF90
  ifeq ($(which_MPI), intel)
               FC := mpiifort
  else
               FC := mpif90
  endif
 else
             LIBS += -lfmpi -lmpi
 endif
endif

ifdef USE_OpenMP
         CPPFLAGS += -D_OPENMP
           FFLAGS += -qopenmp -fpp
             LIBS += -liomp5
endif

ifndef USE_SCRIP
             LIBS += $(MCT_PARAMS_DIR)/mct_coupler_params.o
             LIBS += $(MCT_PARAMS_DIR)/mod_coupler_iounits.o
             LIBS += $(MCT_PARAMS_DIR)/get_sparse_matrix.o
endif

ifdef USE_SWAN
           FFLAGS += -assume byterecl
           FFLAGS += -I$(MY_ROOT_DIR)/SWAN/build/mod
           LIBS += $(MY_ROOT_DIR)/SWAN/build/lib/libswan41.51.a
endif

ifdef USE_WW3
             LIBS += WW3/build/model/src/CMakeFiles/ww3_shel.dir/ww3_shel.F90.o
             LIBS += WW3/build/model/src/CMakeFiles/ww3_multi.dir/ww3_multi.F90.o
             LIBS += WW3/build/lib/libww3.a
endif

ifdef USE_MCT
       MCT_INCDIR ?= /opt/intelsoft/mct/include
       MCT_LIBDIR ?= /opt/intelsoft/mct/lib
           FFLAGS += -I$(MCT_INCDIR)
             LIBS += -L$(MCT_LIBDIR) -lmct -lmpeu
           INCDIR += $(MCT_INCDIR) $(INCDIR)
endif

ifdef USE_ESMF
                     include $(ESMFMKFILE)
          ESMF_OS ?= $(OS)
      ESMF_SUBDIR := $(ESMF_OS).$(ESMF_COMPILER).$(ESMF_ABI).$(ESMF_COMM).$(ESMF_SITE)
      ESMF_MK_DIR ?= $(ESMF_DIR)/lib/lib$(ESMF_BOPT)/$(ESMF_SUBDIR)
           FFLAGS += $(ESMF_F90COMPILEPATHS)
             LIBS += $(ESMF_F90LINKPATHS) $(ESMF_F90ESMFLINKLIBS)
             LIBS += -liomp5 -lstdc++
endif

ifdef USE_WRF
 ifeq "$(strip $(WRF_LIB_DIR))" "$(WRF_SRC_DIR)"
             FFLAGS += -I$(WRF_DIR)/main -I$(WRF_DIR)/external/esmf_time_f90 -I$(WRF_DIR)/frame -I$(WRF_DIR)/share
             LIBS += $(WRF_LIB_DIR)/main/module_wrf_top.o
             LIBS += $(WRF_LIB_DIR)/main/libwrflib.a
             LIBS += $(WRF_LIB_DIR)/external/fftpack/fftpack5/libfftpack.a
             LIBS += $(WRF_LIB_DIR)/external/io_grib1/libio_grib1.a
             LIBS += $(WRF_LIB_DIR)/external/io_grib_share/libio_grib_share.a
             LIBS += $(WRF_LIB_DIR)/external/io_int/libwrfio_int.a
             LIBS += $(WRF_LIB_DIR)/external/esmf_time_f90/libesmf_time.a
             LIBS += $(WRF_LIB_DIR)/external/RSL_LITE/librsl_lite.a
             LIBS += $(WRF_LIB_DIR)/frame/module_internal_header_util.o
             LIBS += $(WRF_LIB_DIR)/frame/pack_utils.o
             LIBS += $(WRF_LIB_DIR)/external/io_netcdf/libwrfio_nf.a
     WRF_MOD_DIRS  = main frame phys share external/esmf_time_f90
 else
             LIBS += $(WRF_LIB_DIR)/module_wrf_top.o
             LIBS += $(WRF_LIB_DIR)/libwrflib.a
             LIBS += $(WRF_LIB_DIR)/libfftpack.a
             LIBS += $(WRF_LIB_DIR)/libio_grib1.a
             LIBS += $(WRF_LIB_DIR)/libio_grib_share.a
             LIBS += $(WRF_LIB_DIR)/libwrfio_int.a
             LIBS += $(WRF_LIB_DIR)/libesmf_time.a
             LIBS += $(WRF_LIB_DIR)/librsl_lite.a
             LIBS += $(WRF_LIB_DIR)/module_internal_header_util.o
             LIBS += $(WRF_LIB_DIR)/pack_utils.o
  ifneq ($(NETCDFPAR),)
             LIBS += $(WRF_LIB_DIR)/libwrfio_nfpar.a
  else
             LIBS += $(WRF_LIB_DIR)/libwrfio_nf.a
  endif
 endif
endif

ifdef USE_WRFHYDRO
             FFLAGS += -I $(WRFHYDRO_DIR)/Land_models/NoahMP/IO_code
             LIBS += $(WRFHYDRO_DIR)/Land_models/NoahMP/IO_code/main_hrldas_driver.o
             LIBS += $(WRFHYDRO_DIR)/Land_models/NoahMP/IO_code/module_hrldas_netcdf_io.o 
             LIBS += $(WRFHYDRO_DIR)/Land_models/NoahMP/phys/module_sf_noahmpdrv.o
             LIBS += $(WRFHYDRO_DIR)/Land_models/NoahMP/phys/module_sf_noahmplsm.o
             LIBS += $(WRFHYDRO_DIR)/Land_models/NoahMP/phys/module_sf_noahmp_glacier.o
             LIBS += $(WRFHYDRO_DIR)/Land_models/NoahMP/phys/module_sf_noahmp_groundwater.o
             LIBS += $(WRFHYDRO_DIR)/Land_models/NoahMP/Utility_routines/module_wrf_utilities.o
             LIBS += $(WRFHYDRO_DIR)/Land_models/NoahMP/Utility_routines/module_model_constants.o
             LIBS += $(WRFHYDRO_DIR)/Land_models/NoahMP/Utility_routines/module_date_utilities.o
             LIBS += $(WRFHYDRO_DIR)/Land_models/NoahMP/Utility_routines/kwm_string_utilities.o
             LIBS += $(WRFHYDRO_DIR)/CPL/COAWST_cpl/hydro_coupler.o
             LIBS += $(WRFHYDRO_DIR)/Land_models/NoahMP/IO_code/module_NoahMP_hrldas_driver.o
             LIBS += $(WRFHYDRO_DIR)/lib/libHYDRO.a
#            LIBS +=  $(WRFHYDRO_DIR)/Land_models/NoahMP/Noah/module_sf_myjsfc.o
#            LIBS +=  $(WRFHYDRO_DIR)/Land_models/NoahMP/Noah/module_sf_sfclay.o
endif

# Use full path of compiler.

               FC := $(shell which ${FC})
               LD := $(FC)

#--------------------------------------------------------------------------
# ROMS specific rules.
#--------------------------------------------------------------------------

# Set free form format in some ROMS source files to allow long string for
# local directory and compilation flags inside the code.

ifdef USE_ROMS
 $(BUILD_DIR)/mod_ncparam.o: FFLAGS += $(FREEFLAGS)
 $(BUILD_DIR)/mod_strings.o: FFLAGS += $(FREEFLAGS)
 $(BUILD_DIR)/analytical.o: FFLAGS += $(FREEFLAGS)
 $(BUILD_DIR)/biology.o: FFLAGS += $(FREEFLAGS)

 ifdef USE_ADJOINT
  $(BUILD_DIR)/ad_biology.o: FFLAGS += $(FREEFLAGS)
 endif
 ifdef USE_REPRESENTER
  $(BUILD_DIR)/rp_biology.o: FFLAGS += $(FREEFLAGS)
 endif
 ifdef USE_TANGENT
  $(BUILD_DIR)/tl_biology.o: FFLAGS += $(FREEFLAGS)
 endif
endif

#--------------------------------------------------------------------------
# Model coupling specific rules.
#--------------------------------------------------------------------------

# Add COAMPS library directory to include path of ESMF coupling files.

ifdef USE_COAMPS
 $(BUILD_DIR)/esmf_atm.o: FFLAGS += -I$(COAMPS_LIB_DIR)
 $(BUILD_DIR)/esmf_esm.o: FFLAGS += -I$(COAMPS_LIB_DIR)
endif

# Add WRF library directory to include path of ESMF coupling files.

ifdef USE_WRF
 $(BUILD_DIR)/esmf_atm.o: FFLAGS += -I$(WRF_LIB_DIR)
endif

# Supress free format in SWAN source files since there are comments
# beyond column 72.

ifdef USE_SWAN
$(BUILD_DIR)/ocpcre.o:   FFLAGS += -nofree
$(BUILD_DIR)/ocpids.o:   FFLAGS += -nofree
$(BUILD_DIR)/ocpmix.o:   FFLAGS += -nofree
$(BUILD_DIR)/swancom1.o: FFLAGS += -nofree
$(BUILD_DIR)/swancom2.o: FFLAGS += -nofree
$(BUILD_DIR)/swancom3.o: FFLAGS += -nofree
$(BUILD_DIR)/swancom4.o: FFLAGS += -nofree
$(BUILD_DIR)/swancom5.o: FFLAGS += -nofree
$(BUILD_DIR)/swanmain.o: FFLAGS += -nofree
$(BUILD_DIR)/swanout1.o: FFLAGS += -nofree
$(BUILD_DIR)/swanout2.o: FFLAGS += -nofree
$(BUILD_DIR)/swanparll.o:FFLAGS += -nofree
$(BUILD_DIR)/swanpre1.o: FFLAGS += -nofree
$(BUILD_DIR)/swanpre2.o: FFLAGS += -nofree
$(BUILD_DIR)/swanser.o:  FFLAGS += -nofree
$(BUILD_DIR)/swmod1.o:   FFLAGS += -nofree
$(BUILD_DIR)/swmod2.o:       FFLAGS += -nofree
$(BUILD_DIR)/SwanSpectPart.o:FFLAGS += -nofree
$(BUILD_DIR)/m_constants.o:  FFLAGS += -free
$(BUILD_DIR)/m_fileio.o:     FFLAGS += -free
$(BUILD_DIR)/mod_xnl4v5.o:   FFLAGS += -free
$(BUILD_DIR)/serv_xnl4v5.o:  FFLAGS += -free
$(BUILD_DIR)/nctablemd.o:    FFLAGS += -free
$(BUILD_DIR)/agioncmd.o:     FFLAGS += -free
$(BUILD_DIR)/swn_outnc.o:    FFLAGS += -free
$(BUILD_DIR)/SdsBabanin.o:   FFLAGS += -free
$(BUILD_DIR)/SwanBpntlist.o: FFLAGS += -free
$(BUILD_DIR)/SwanCheckGrid.o:FFLAGS += -free
$(BUILD_DIR)/SwanCompdata.o: FFLAGS += -free
$(BUILD_DIR)/SwanCompUnstruc.o:       FFLAGS += -free
$(BUILD_DIR)/SwanComputeForce.o:      FFLAGS += -free
$(BUILD_DIR)/SwanConvAccur.o:         FFLAGS += -free
$(BUILD_DIR)/SwanConvStopc.o:         FFLAGS += -free
$(BUILD_DIR)/SwanCreateEdges.o:       FFLAGS += -free
$(BUILD_DIR)/SwanCrossObstacle.o:     FFLAGS += -free
$(BUILD_DIR)/SwanDiffPar.o:           FFLAGS += -free
$(BUILD_DIR)/SwanDispParm.o:          FFLAGS += -free
$(BUILD_DIR)/SwanFindObstacles.o:     FFLAGS += -free
$(BUILD_DIR)/SwanFindPoint.o:         FFLAGS += -free
$(BUILD_DIR)/SwanGridCell.o:          FFLAGS += -free
$(BUILD_DIR)/SwanGriddata.o:          FFLAGS += -free
$(BUILD_DIR)/SwanGridFace.o:          FFLAGS += -free
$(BUILD_DIR)/SwanGridobjects.o:       FFLAGS += -free
$(BUILD_DIR)/SwanGridTopology.o:      FFLAGS += -free
$(BUILD_DIR)/SwanGridVert.o:          FFLAGS += -free
$(BUILD_DIR)/SwanGSECorr.o:           FFLAGS += -free
$(BUILD_DIR)/SwanInitCompGrid.o:      FFLAGS += -free
$(BUILD_DIR)/SwanInterpolateAc.o:     FFLAGS += -free
$(BUILD_DIR)/SwanInterpolateOutput.o: FFLAGS += -free
$(BUILD_DIR)/SwanInterpolatePoint.o:  FFLAGS += -free
$(BUILD_DIR)/SwanIntgratSpc.o:        FFLAGS += -free
$(BUILD_DIR)/SwanPointinMesh.o:       FFLAGS += -free
$(BUILD_DIR)/SwanPrepComp.o:          FFLAGS += -free
$(BUILD_DIR)/SwanPrintGridInfo.o:     FFLAGS += -free
$(BUILD_DIR)/SwanPropvelS.o:          FFLAGS += -free
$(BUILD_DIR)/SwanPropvelX.o:          FFLAGS += -free
$(BUILD_DIR)/SwanReadADCGrid.o:       FFLAGS += -free
$(BUILD_DIR)/SwanReadEasymeshGrid.o:  FFLAGS += -free
$(BUILD_DIR)/SwanReadGrid.o:          FFLAGS += -free
$(BUILD_DIR)/SwanReadTriangleGrid.o:  FFLAGS += -free
$(BUILD_DIR)/SwanSweepSel.o:          FFLAGS += -free
$(BUILD_DIR)/SwanThreadBounds.o:      FFLAGS += -free
$(BUILD_DIR)/SwanTranspAc.o:          FFLAGS += -free
$(BUILD_DIR)/SwanTranspX.o:           FFLAGS += -free
$(BUILD_DIR)/SwanVertlist.o:  FFLAGS += -free
$(BUILD_DIR)/waves_control.o: FFLAGS += -free
$(BUILD_DIR)/waves_coupler.o: FFLAGS += -free

endif
