# svn $Id: CYGWIN-ifort.mk 734 2008-09-07 01:58:06Z jcwarner $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2016 The ROMS/TOMS Group                           :::
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

              BIN := $(BIN).exe

               FC := ifort
           FFLAGS := /align /G7 /MD
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -DCYGWIN -DCYGWIN_ifort -traditional
               CC := gcc   # ifc?
              CXX := g++
           CFLAGS :=
         CXXFLAGS :=
          LDFLAGS := /link /stack:67108864
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
    NETCDF_INCDIR ?= /usr/local/netcdf4/include
    NETCDF_LIBDIR ?= /usr/local/netcdf4/lib
      HDF5_LIBDIR ?= /usr/local/hdf5/lib
else
    NETCDF_INCDIR ?= /netcdf-win32/include
    NETCDF_LIBDIR ?= /netcdf-win32/lib
endif
       NETCDF_LIB := $(NETCDF_LIBDIR)/libnetcdf.a
ifdef USE_NETCDF4
       NETCDF_LIB += -L$(HDF5_LIBDIR) -lhdf5_hl -lhdf5 -lz
endif

ifdef USE_ARPACK
    ARPACK_LIBDIR ?= /arpack-win32/lib
       ARPACK_LIB := $(ARPACK_LIBDIR)/arpack.lib
endif

#
# Compiler flags
#

ifdef USE_OpenMP
         CPPFLAGS += -D_OPENMP
           FFLAGS += /Qopenmp /Qopenmp_report1
endif

ifdef USE_DEBUG
           FFLAGS += /debug:full /traceback /Od /Zi /check:bounds
else
           FFLAGS += /fp:precise /fp:source /O3
endif

ifdef USE_SWAN
           FFLAGS += /noextend_source -assume:byterecl
endif

ifdef USE_MPI
       MPI_INCDIR ?= c:\\work\\models\\MPICH2\\include
       MPI_LIBDIR ?= c:\\work\\models\\MPICH2\\lib
       LIBS_WIN32 += "$(MPI_LIBDIR)\fmpich2.lib "
         CPPFLAGS += -DMPI -I$(MPI_INCDIR)
           FFLAGS += -I$(MPI_INCDIR)
endif

ifdef USE_MCT
       MCT_LIBDIR ?= c:\\work\\models\\MCT_v2.2\\lib
       MCT_INCDIR ?= c:\\work\\models\\MCT_v2.2\\include
         CPPFLAGS += -traditional-cpp
           FFLAGS += -I$(MCT_INCDIR)
       LIBS_WIN32 += "$(MCT_LIBDIR)\libmct.a" "$(MCT_LIBDIR)\libmpeu.a"
endif

ifdef USE_ESMF
      ESMF_SUBDIR := $(ESMF_OS).$(ESMF_COMPILER).$(ESMF_ABI).$(ESMF_COMM).$(ESMF_SITE)
      ESMF_MK_DIR ?= $(ESMF_DIR)/lib/lib$(ESMF_BOPT)/$(ESMF_SUBDIR)
                     include $(ESMF_MK_DIR)/esmf.mk
           FFLAGS += $(ESMF_F90COMPILEPATHS)
       LIBS_WIN32 += $(ESMF_F90LINKPATHS) -lesmf -lC
endif

ifdef USE_REFDIF
#         CPPFLAGS += -traditional-cpp
           FFLAGS += -I$(MCT_LIBDIR) -I$(MPEU_LIBDIR) 
#           FFLAGS += /noextend_source -assume:byterecl
       LIBS_WIN32 += "$(MCT_LIBDIR)\libmct.a" "$(MPEU_LIBDIR)\libmpeu.a"
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

$(SCRATCH_DIR)/ocpcre.o:    FFLAGS += -nofree
$(SCRATCH_DIR)/ocpids.o:    FFLAGS += -nofree
$(SCRATCH_DIR)/ocpmix.o:    FFLAGS += -nofree
$(SCRATCH_DIR)/swancom1.o:  FFLAGS += -nofree
$(SCRATCH_DIR)/swancom2.o:  FFLAGS += -nofree
$(SCRATCH_DIR)/swancom3.o:  FFLAGS += -nofree
$(SCRATCH_DIR)/swancom4.o:  FFLAGS += -nofree
$(SCRATCH_DIR)/swancom5.o:  FFLAGS += -nofree
$(SCRATCH_DIR)/swanmain.o:  FFLAGS += -nofree
$(SCRATCH_DIR)/swanout1.o:  FFLAGS += -nofree
$(SCRATCH_DIR)/swanout2.o:  FFLAGS += -nofree
$(SCRATCH_DIR)/swanparll.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swanpre1.o:  FFLAGS += -nofree
$(SCRATCH_DIR)/swanpre2.o:  FFLAGS += -nofree
$(SCRATCH_DIR)/swanser.o:   FFLAGS += -nofree
$(SCRATCH_DIR)/swmod1.o:    FFLAGS += -nofree
$(SCRATCH_DIR)/swmod2.o:    FFLAGS += -nofree
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
