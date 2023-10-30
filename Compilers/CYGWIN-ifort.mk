# git $Id$
# svn $Id: CYGWIN-ifort.mk 1151 2023-02-09 03:08:53Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#
# Include file for Intel compiler on Cygwin
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

               FC := ifort
           FFLAGS := /align /G7 /MD
       FIXEDFLAGS := /nofree
        FREEFLAGS := /free
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -DCYGWIN -DCYGWIN_ifort -traditional
           INCDIR := /usr/include /usr/local/bin
            SLIBS := -L/usr/local/lib -L/usr/lib -L/usr/lib64
            ULIBS :=
             LIBS :=
         ROMS_LIB := -L$(SCRATCH_DIR) -lROMS
       MOD_SUFFIX := mod
               LD := $(FC)
          LDFLAGS := /link /stack:67108864
               AR := ar
          ARFLAGS := r
            MKDIR := mkdir -p
               CP := cp -p -v
               RM := rm -f
           RANLIB := ranlib
             PERL := perl
             TEST := test
      ST_LIB_NAME := libROMS.a
      SH_LIB_NAME := cygROMS.dll

#--------------------------------------------------------------------------
# Compiling flags for ROMS Applications.
#--------------------------------------------------------------------------

ifdef USE_ROMS
 ifdef USE_DEBUG
           FFLAGS += /debug:full
           FFLAGS += /traceback
           FFLAGS += /Od
           FFLAGS += /Zi
 else
           FFLAGS += /Ox
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
# These are be specified in Unix form and will be converted as
# necessary to Windows form for Windows-native commands. The default
# values below assume that Cygwin mounts have been defined pointing to
# the NETCDF and ARPACK library locations.
#

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
       MCT_LIBDIR ?= c:\\work\\models\\MCT_v2.2\\mct
      MPEU_LIBDIR ?= c:\\work\\models\\MCT_v2.2\\mpeu
       MPI_INCDIR ?= c:\\work\\models\\MPICH2\\include
       MPI_LIBDIR ?= c:\\work\\models\\MPICH2\\lib
       LIBS_WIN32 += "$(MPI_LIBDIR)\fmpich2.lib "
         CPPFLAGS += -DMPI -I$(MPI_INCDIR)
           FFLAGS += -I$(MPI_INCDIR)
endif

ifdef USE_OpenMP
         CPPFLAGS += -D_OPENMP
           FFLAGS += /Qopenmp /Qopenmp_report1
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
             LIBS += WW3/build/model/src/CMakeFiles/ww3_shel.dir/ww3_shel.F90.o
             LIBS += WW3/build/lib/libww3.a
endif

ifdef USE_MCT
         CPPFLAGS += -traditional-cpp
           FFLAGS += -I$(MCT_LIBDIR) -I$(MPEU_LIBDIR)
           FFLAGS += /noextend_source -assume:byterecl
       LIBS_WIN32 += "$(MCT_LIBDIR)\libmct.a" "$(MPEU_LIBDIR)\libmpeu.a"
           INCDIR += $(MCT_INCDIR) $(INCDIR)
endif

ifdef USE_ESMF
          ESMF_OS ?= $(OS)
      ESMF_SUBDIR := $(ESMF_OS).$(ESMF_COMPILER).$(ESMF_ABI).$(ESMF_COMM).$(ESMF_SITE)
      ESMF_MK_DIR ?= $(ESMF_DIR)/lib/lib$(ESMF_BOPT)/$(ESMF_SUBDIR)
                     include $(ESMF_MK_DIR)/esmf.mk
           FFLAGS += $(ESMF_F90COMPILEPATHS)
       LIBS_WIN32 += $(ESMF_F90LINKPATHS) $(ESMF_F90ESMFLINKLIBS)
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

#
# For a Windows compiler, create variables pointing to the Windows
# file names needed when linking. Use of the "=" sign means that
# variables will be evaluated only when needed.
#
         BIN_WIN32 = "$$(cygpath --windows $(BIN))"
        LIBS_WIN32 += "$$(cygpath --windows $(NETCDF_LIB))"
        LIBS_WIN32 += "c:\cygwin\lib\gcc\i686-pc-cygwin\4.9.2\libgcc.a"
ifdef USE_ARPACK
        LIBS_WIN32 += "$$(cygpath --windows $(ARPACK_LIB))"
endif

        LD_WINDOWS := on

#
# Use full path of compiler. Notice that a very special editing is
# done for the FC defintion to allow blank spaces and parenthesis in the
# compiler path. For example, we can have:
#
#      c:\\Software\\My Compilers (64 bit)
# or
#      /cygdrive/c/Program Files (x86)/Intel/Compiler/11.1/051/bin/ia32/ifort
# or
#      /cygdrive/c/Program Files/Intel/Compiler/11.1/051/bin/ia32/ifort
#
               FC := $(shell which ${FC} | sed 's|\([ |(|)]\)|\\\1|g')
               LD := $(FC)

#--------------------------------------------------------------------------
# ROMS specific rules.
#--------------------------------------------------------------------------

# CYGWIN can only load user compiled .dll files located in the same
# directory as the executable. This rule will copy the cygROMS.dll
# to the $(BINDIR) so the executable can run. This is only needed when
# EXEC and SHARED are set and STATIC is NOT. If STATIC is set then
# CYGWIN will automatically link with the static library.

# Blank it out to be sure
       CYG_DLL_CP :=

ifdef SHARED
 ifdef EXEC
  ifndef STATIC
	CYG_DLL_CP := cyg_dll_cp

.PHONY: cyg_dll_cp
cyg_dll_cp: $(BIN)
	$(CP) $(SCRATCH_DIR)/$(SH_LIB_NAME) $(BINDIR)

  endif
 endif
endif

#
# For a Windows compiler, override the compilation rule
#

%.o: %.f90
	cd $(SCRATCH_DIR); $(FC) -c $(FFLAGS) $(notdir $<) /object:$(notdir $@)

#
# Override rule for object files. Use .o instead of .obj
#

define one-compile-rule
  $1: $2 $3
	cd $$(SCRATCH_DIR); $$(FC) -c $$(FFLAGS) $(notdir $2) /object:$(notdir $1)

  $2: $3
	$$(CPP) $$(CPPFLAGS) $$(MY_CPP_FLAGS) $$< > $$@
	$$(CLEAN) $$@

endef

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

$(SCRATCH_DIR)/ocpcre.o:   FFLAGS += -nofree
$(SCRATCH_DIR)/ocpids.o:   FFLAGS += -nofree
$(SCRATCH_DIR)/ocpmix.o:   FFLAGS += -nofree
$(SCRATCH_DIR)/swancom1.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swancom2.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swancom3.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swancom4.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swancom5.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swanmain.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swanout1.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swanout2.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swanparll.o:FFLAGS += -nofree
$(SCRATCH_DIR)/swanpre1.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swanpre2.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swanser.o:  FFLAGS += -nofree
$(SCRATCH_DIR)/swmod1.o:   FFLAGS += -nofree
$(SCRATCH_DIR)/swmod2.o:       FFLAGS += -nofree
$(SCRATCH_DIR)/SwanSpectPart.o:FFLAGS += -nofree
$(SCRATCH_DIR)/m_constants.o:  FFLAGS += -free
$(SCRATCH_DIR)/m_fileio.o:     FFLAGS += -free
$(SCRATCH_DIR)/mod_xnl4v5.o:   FFLAGS += -free
$(SCRATCH_DIR)/serv_xnl4v5.o:  FFLAGS += -free
$(SCRATCH_DIR)/nctablemd.o:    FFLAGS += -free
$(SCRATCH_DIR)/agioncmd.o:     FFLAGS += -free
$(SCRATCH_DIR)/swn_outnc.o:    FFLAGS += -free
$(SCRATCH_DIR)/SdsBabanin.o:   FFLAGS += -free
$(SCRATCH_DIR)/SwanBpntlist.o: FFLAGS += -free
$(SCRATCH_DIR)/SwanCheckGrid.o:FFLAGS += -free
$(SCRATCH_DIR)/SwanCompdata.o: FFLAGS += -free
$(SCRATCH_DIR)/SwanCompUnstruc.o:       FFLAGS += -free
$(SCRATCH_DIR)/SwanComputeForce.o:      FFLAGS += -free
$(SCRATCH_DIR)/SwanConvAccur.o:         FFLAGS += -free
$(SCRATCH_DIR)/SwanConvStopc.o:         FFLAGS += -free
$(SCRATCH_DIR)/SwanCreateEdges.o:       FFLAGS += -free
$(SCRATCH_DIR)/SwanCrossObstacle.o:     FFLAGS += -free
$(SCRATCH_DIR)/SwanDiffPar.o:           FFLAGS += -free
$(SCRATCH_DIR)/SwanDispParm.o:          FFLAGS += -free
$(SCRATCH_DIR)/SwanFindObstacles.o:     FFLAGS += -free
$(SCRATCH_DIR)/SwanFindPoint.o:         FFLAGS += -free
$(SCRATCH_DIR)/SwanGridCell.o:          FFLAGS += -free
$(SCRATCH_DIR)/SwanGriddata.o:          FFLAGS += -free
$(SCRATCH_DIR)/SwanGridFace.o:          FFLAGS += -free
$(SCRATCH_DIR)/SwanGridobjects.o:       FFLAGS += -free
$(SCRATCH_DIR)/SwanGridTopology.o:      FFLAGS += -free
$(SCRATCH_DIR)/SwanGridVert.o:          FFLAGS += -free
$(SCRATCH_DIR)/SwanGSECorr.o:           FFLAGS += -free
$(SCRATCH_DIR)/SwanInitCompGrid.o:      FFLAGS += -free
$(SCRATCH_DIR)/SwanInterpolateAc.o:     FFLAGS += -free
$(SCRATCH_DIR)/SwanInterpolateOutput.o: FFLAGS += -free
$(SCRATCH_DIR)/SwanInterpolatePoint.o:  FFLAGS += -free
$(SCRATCH_DIR)/SwanIntgratSpc.o:        FFLAGS += -free
$(SCRATCH_DIR)/SwanPointinMesh.o:       FFLAGS += -free
$(SCRATCH_DIR)/SwanPrepComp.o:          FFLAGS += -free
$(SCRATCH_DIR)/SwanPrintGridInfo.o:     FFLAGS += -free
$(SCRATCH_DIR)/SwanPropvelS.o:          FFLAGS += -free
$(SCRATCH_DIR)/SwanPropvelX.o:          FFLAGS += -free
$(SCRATCH_DIR)/SwanReadADCGrid.o:       FFLAGS += -free
$(SCRATCH_DIR)/SwanReadEasymeshGrid.o:  FFLAGS += -free
$(SCRATCH_DIR)/SwanReadGrid.o:          FFLAGS += -free
$(SCRATCH_DIR)/SwanReadTriangleGrid.o:  FFLAGS += -free
$(SCRATCH_DIR)/SwanSweepSel.o:          FFLAGS += -free
$(SCRATCH_DIR)/SwanThreadBounds.o:      FFLAGS += -free
$(SCRATCH_DIR)/SwanTranspAc.o:          FFLAGS += -free
$(SCRATCH_DIR)/SwanTranspX.o:           FFLAGS += -free
$(SCRATCH_DIR)/SwanVertlist.o:  FFLAGS += -free
$(SCRATCH_DIR)/waves_control.o: FFLAGS += -free
$(SCRATCH_DIR)/waves_coupler.o: FFLAGS += -free

endif
