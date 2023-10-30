# git $Id$
# svn $Id: MINGW-gfortran.mk 1151 2023-02-09 03:08:53Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Include file for GNU on MSYS/MinGW
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
# ROMS_LIB       Directory and name for ROMS library
# NF_CONFIG      NetCDF Fortran configuration script
# NETCDF_INCDIR  NetCDF include directory
# NETCDF_LIBDIR  NetCDF library directory
# NETCDF_LIBS    NetCDF library switches

# LD             Program to load the objects into an executable or shared library
# LDFLAGS        Flags to the loader
# RANLIB         Name of ranlib command
# MDEPFLAGS      Flags for sfmakedepend  (-s if you keep .f files)
#
# First the defaults
#

# This needs an if block to prevent building the executable
# when none is requested.
ifdef BIN
              BIN := $(BIN).exe
endif

               FC := gfortran
           FFLAGS := -frepack-arrays
       FIXEDFLAGS := -ffixed-form
        FREEFLAGS := -ffree-form -ffree-line-length-none
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -traditional -w              # -w turns off warnings
           INCDIR := /usr/include /usr/local/bin
            SLIBS := -L/usr/local/lib -L/usr/lib -L/usr/lib64
            ULIBS :=
             LIBS :=
         ROMS_LIB := -L$(SCRATCH_DIR) -lROMS
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
      ST_LIB_NAME := libROMS.a
      SH_LIB_NAME := cygROMS.dll

#--------------------------------------------------------------------------
# Checking for minimum version and versions requiring extra compile flags
#--------------------------------------------------------------------------

# Because of the recursive derived types in Utility/yaml_parser.F, gfortran
# version 7 or greater is required.

          MIN_VER := $(shell expr `$(FC) -dumpversion | cut -f1 -d.` \<= 6)
ifeq "$(MIN_VER)" "1"
      $(error gfortran version 7 or greater is required)
endif

# Starting with gfortran version 10, GNU takes a much stricter approach to argument
# type and size matching. This conflicts with ROMS high level MPI calls so this flag
# is added for compilation to complete successfully.

       STRICT_VER := $(shell expr `$(FC) -dumpversion | cut -f1 -d.` \>= 10)
ifeq "$(STRICT_VER)" "1"
           FFLAGS += -fallow-argument-mismatch       # needed for gfortran 10 and higher
endif

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
           FFLAGS += -ffast-math
 endif
 ifdef SHARED
           FFLAGS += -fPIC
       SH_LDFLAGS += -shared
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
#          FFLAGS += -fcheck=all
           FFLAGS += -fsanitize=address -fsanitize=undefined
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


ifdef USE_NETCDF4
        NF_CONFIG ?= nf-config
    NETCDF_INCDIR ?= $(shell $(NF_CONFIG) --prefix)/include
             LIBS += $(shell $(NF_CONFIG) --flibs)
           INCDIR += $(NETCDF_INCDIR) $(INCDIR)
else
    NETCDF_INCDIR ?= /usr/include
    NETCDF_LIBDIR ?= /usr/lib
      NETCDF_LIBS ?= -lnetcdf
             LIBS += -L$(NETCDF_LIBDIR) $(NETCDF_LIBS)
           INCDIR += $(NETCDF_INCDIR) $(INCDIR)
endif

ifdef USE_HDF5
      HDF5_INCDIR ?= /usr/include
      HDF5_LIBDIR ?= /usr/lib
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
 ifdef USE_MPIF90
               FC := mpif90
 else
             LIBS += -lfmpi -lmpi
 endif
endif

ifdef USE_OpenMP
         CPPFLAGS += -D_OPENMP
           FFLAGS += -fopenmp -static-libgcc
endif

ifndef USE_SCRIP
             LIBS += $(MCT_PARAMS_DIR)/mct_coupler_params.o
             LIBS += $(MCT_PARAMS_DIR)/mod_coupler_iounits.o
endif

ifdef USE_SWAN
           FFLAGS += -assume byterecl
           FFLAGS += -I$(MY_ROOT_DIR)/SWAN/build/mod
           LIBS += $(MY_ROOT_DIR)/SWAN/build/lib/libswan41.45.a
endif

ifdef USE_WW3
             FFLAGS += -frecord-marker=4 -fconvert=big-endian
             LIBS += WW3/build/model/src/CMakeFiles/ww3_shel.dir/ww3_shel.F90.o
             LIBS += WW3/build/lib/libww3.a
endif

ifdef USE_MCT
       MCT_INCDIR ?= /usr/local/mct/include
       MCT_LIBDIR ?= /usr/local/mct/lib
           FFLAGS += -I$(MCT_INCDIR)
             LIBS += -L$(MCT_LIBDIR) -lmct -lmpeu
           INCDIR += $(MCT_INCDIR) $(INCDIR)
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
             LIBS += $(WRF_LIB_DIR)/libwrfio_nf.a
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

#
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

# Add WRF library directory to include path of ESMF coupling files.

ifdef USE_WRF
 ifeq "$(strip $(WRF_LIB_DIR))" "$(WRF_SRC_DIR)"
  $(SCRATCH_DIR)/esmf_atm.o: FFLAGS += $(addprefix -I$(WRF_LIB_DIR)/,$(WRF_MOD_DIRS))
 else
  $(SCRATCH_DIR)/esmf_atm.o: FFLAGS += -I$(WRF_LIB_DIR)
 endif
endif

# Supress free format in SWAN source files since there are comments
# beyond column 72.

ifdef USE_SWAN

$(SCRATCH_DIR)/ocpcre.o: FFLAGS += -ffixed-form
$(SCRATCH_DIR)/ocpids.o: FFLAGS += -ffixed-form
$(SCRATCH_DIR)/ocpmix.o: FFLAGS += -ffixed-form
$(SCRATCH_DIR)/swancom1.o: FFLAGS += -ffixed-form
$(SCRATCH_DIR)/swancom2.o: FFLAGS += -ffixed-form
$(SCRATCH_DIR)/swancom3.o: FFLAGS += -ffixed-form
$(SCRATCH_DIR)/swancom4.o: FFLAGS += -ffixed-form
$(SCRATCH_DIR)/swancom5.o: FFLAGS += -ffixed-form
$(SCRATCH_DIR)/swanmain.o: FFLAGS += -ffixed-form
$(SCRATCH_DIR)/swanout1.o: FFLAGS += -ffixed-form
$(SCRATCH_DIR)/swanout2.o: FFLAGS += -ffixed-form
$(SCRATCH_DIR)/swanparll.o: FFLAGS += -ffixed-form
$(SCRATCH_DIR)/swanpre1.o: FFLAGS += -ffixed-form
$(SCRATCH_DIR)/swanpre2.o: FFLAGS += -ffixed-form
$(SCRATCH_DIR)/swanser.o: FFLAGS += -ffixed-form
$(SCRATCH_DIR)/swmod1.o: FFLAGS += -ffixed-form
$(SCRATCH_DIR)/swmod2.o: FFLAGS += -ffixed-form
$(SCRATCH_DIR)/SwanSpectPart.o:  FFLAGS += -ffixed-form
$(SCRATCH_DIR)/m_constants.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/m_fileio.o:    FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/mod_xnl4v5.o:  FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/serv_xnl4v5.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/nctablemd.o:   FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/agioncmd.o:    FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/swn_outnc.o:   FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SdsBabanin.o:  FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanBpntlist.o:  FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanCheckGrid.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanCompdata.o:  FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanCompUnstruc.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanComputeForce.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanConvAccur.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanConvStopc.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanCreateEdges.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanCrossObstacle.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanDiffPar.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanDispParm.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanFindObstacles.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanFindPoint.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanGridCell.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanGriddata.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanGridFace.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanGridobjects.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanGridTopology.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanGridVert.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanGSECorr.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanInitCompGrid.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanInterpolateAc.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanInterpolateOutput.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanInterpolatePoint.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanIntgratSpc.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanPointinMesh.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanPrepComp.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanPrintGridInfo.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanPropvelS.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanPropvelX.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanReadADCGrid.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanReadEasymeshGrid.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanReadGrid.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanReadTriangleGrid.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanSweepSel.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanThreadBounds.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanTranspAc.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanTranspX.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/SwanVertlist.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/waves_control.o: FFLAGS += -ffree-form -ffree-line-length-none
$(SCRATCH_DIR)/waves_coupler.o: FFLAGS += -ffree-form -ffree-line-length-none

endif
