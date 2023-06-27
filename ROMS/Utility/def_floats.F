#include "cppdefs.h"
      MODULE def_floats_mod
#ifdef FLOATS
!
!git $Id$
!svn $Id: def_floats.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module creates output FLOAT data file using either the         !
!  standard NetCDF library or the Parallel-IO (PIO) library.  It       !
!  defines its dimensions, attributes, and variables.                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef BIOLOGY
      USE mod_biology
# endif
      USE mod_floats
# ifdef FOUR_DVAR
      USE mod_fourdvar
# endif
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
# ifdef SEDIMENT
      USE mod_sediment
# endif
!
      USE def_dim_mod,  ONLY : def_dim
      USE def_info_mod, ONLY : def_info
      USE def_var_mod,  ONLY : def_var
      USE strings_mod,  ONLY : FoundError
      USE wrt_info_mod, ONLY : wrt_info
!
      implicit none
!
      PUBLIC  :: def_floats
      PRIVATE :: def_floats_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: def_floats_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE def_floats (ng, ldef)
!***********************************************************************
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Create a new history file according to IO type.
!-----------------------------------------------------------------------
!
      SELECT CASE (FLT(ng)%IOtype)
        CASE (io_nf90)
          CALL def_floats_nf90 (ng, ldef)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL def_floats_pio (ng, ldef)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) FLT(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' DEF_FLOATS - Illegal output file type, io_type = ',i0,  &
     &        /,14x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE def_floats
!
!***********************************************************************
      SUBROUTINE def_floats_nf90 (ng, ldef)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
      logical, intent(in) :: ldef
!
!  Local variable declarations.
!
      logical :: got_var(-6:NV)
!
      integer, parameter :: Natt = 25

      integer :: fltdim, i, itrc, j, l
      integer :: recdim, status

      integer :: DimIDs(nDimID)
      integer :: fgrd(2), start(2), total(2)
!
      real(r8) :: Aval(6), Tinp(Nfloats(ng))
!
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", def_floats_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=FLT(ng)%name
!
      IF (Master) THEN
        IF (ldef) THEN
          WRITE (stdout,10) ng, TRIM(ncname)
        ELSE
          WRITE (stdout,20) ng, TRIM(ncname)
        END IF
      END IF
!
!=======================================================================
!  Create a new floats data file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL netcdf_create (ng, iNLM, TRIM(ncname), FLT(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,30) TRIM(ncname)
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define file dimensions.
!-----------------------------------------------------------------------
!
        DimIDs=0
!
# ifdef SOLVE3D
        status=def_dim(ng, iNLM, FLT(ng)%ncid, ncname, 's_rho',         &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, FLT(ng)%ncid, ncname, 's_w',           &
     &                 N(ng)+1, DimIDs(10))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, FLT(ng)%ncid, ncname, 'tracer',        &
     &                 NT(ng), DimIDs(11))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SEDIMENT
        status=def_dim(ng, iNLM, FLT(ng)%ncid, ncname, 'NST',           &
     &                 NST, DimIDs(32))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, FLT(ng)%ncid, ncname, 'Nbed',          &
     &                 Nbed, DimIDs(16))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef ECOSIM
        status=def_dim(ng, iNLM, FLT(ng)%ncid, ncname, 'Nbands',        &
     &                 NBands, DimIDs(33))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, FLT(ng)%ncid, ncname, 'Nphy',          &
     &                 Nphy, DimIDs(25))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, FLT(ng)%ncid, ncname, 'Nbac',          &
     &                 Nbac, DimIDs(26))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, FLT(ng)%ncid, ncname, 'Ndom',          &
     &                 Ndom, DimIDs(27))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, FLT(ng)%ncid, ncname, 'Nfec',          &
     &                 Nfec, DimIDs(28))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

        status=def_dim(ng, iNLM, FLT(ng)%ncid, ncname, 'drifter' ,      &
     &                 Nfloats(ng), DimIDs(15))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, FLT(ng)%ncid, ncname, 'boundary',      &
     &                 4, DimIDs(14))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef FOUR_DVAR
        status=def_dim(ng, iNLM, FLT(ng)%ncid, ncname, 'Nstate',        &
     &                 NstateVar(ng), DimIDs(29))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

        status=def_dim(ng, iNLM, FLT(ng)%ncid, ncname,                  &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        recdim=DimIDs(12)
        fltdim=DimIDs(15)
!
!  Define dimension vectors for point variables.
!
        fgrd(1)=DimIDs(15)
        fgrd(2)=DimIDs(12)
!
!  Initialize unlimited time record dimension.
!
        FLT(ng)%Rindex=0
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,LEN(Vinfo(1))
            Vinfo(i)(j:j)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
!-----------------------------------------------------------------------
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, iNLM, FLT(ng)%ncid, ncname, DimIDs)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define variables and their attributes.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        Vinfo( 2)=Vname(2,idtime)
        WRITE (Vinfo( 3),'(a,a)') 'seconds since ', TRIM(Rclock%string)
        Vinfo( 4)=TRIM(Rclock%calendar)
        Vinfo(14)=Vname(4,idtime)
        status=def_var(ng, iNLM, FLT(ng)%ncid, FLT(ng)%Vid(idtime),     &
     &                 NF_TOUT, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define floats X-grid locations.
!
        Vinfo( 1)='Xgrid'
        Vinfo( 2)='x-grid floats locations'
        Vinfo( 5)='valid_min'
        Vinfo( 6)='valid_max'
        Aval(2)=0.0_r8
        Aval(3)=REAL(Lm(ng)+1,r8)
        Vinfo(14)='Xgrid, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FLT(ng)%ncid, FLT(ng)%Vid(idXgrd),     &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define floats Y-grid locations.
!
        Vinfo( 1)='Ygrid'
        Vinfo( 2)='Y-grid floats locations'
        Vinfo( 5)='valid_min'
        Vinfo( 6)='valid_max'
        Aval(2)=0.0_r8
        Aval(3)=REAL(Mm(ng)+1,r8)
        Vinfo(14)='Ygrid, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FLT(ng)%ncid, FLT(ng)%Vid(idYgrd),     &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef SOLVE3D
!
!  Define floats Z-grid locations.
!
        Vinfo( 1)='Zgrid'
        Vinfo( 2)='Z-grid floats locations'
        Vinfo( 5)='valid_min'
        Vinfo( 6)='valid_max'
        Aval(2)=0.0_r8
        Aval(3)=REAL(N(ng),r8)
        Vinfo(14)='Zgrid, scalar, series'
        Vinfo(16)=Vname(1,idtime)
#  ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  endif
        status=def_var(ng, iNLM, FLT(ng)%ncid, FLT(ng)%Vid(idZgrd),     &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
!  Define floats (lon,lat) or (x,y) locations.
!
        IF (spherical) THEN
          Vinfo( 1)='lon'
          Vinfo( 2)='longitude of floats trajectories'
          Vinfo( 3)='degree_east'
          Vinfo( 5)='valid_min'
          Vinfo( 6)='valid_max'
          Vinfo(14)='lon, scalar, series'
          Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
          Vinfo(24)='_FillValue'
          Aval(6)=spval
# endif
          Aval(2)=-180.0_r8
          Aval(3)=180.0_r8
          status=def_var(ng, iNLM, FLT(ng)%ncid, FLT(ng)%Vid(idglon),   &
     &                   NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

          Vinfo( 1)='lat'
          Vinfo( 2)='latitude of floats trajectories'
          Vinfo( 3)='degree_north'
          Vinfo( 5)='valid_min'
          Vinfo( 6)='valid_max'
          Vinfo(14)='lat, scalar, series'
          Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
          Vinfo(24)='_FillValue'
          Aval(6)=spval
# endif
          Aval(2)=-90.0_r8
          Aval(3)=90.0_r8
          status=def_var(ng, iNLM, FLT(ng)%ncid, FLT(ng)%Vid(idglat),   &
     &                   NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
       ELSE
          Vinfo( 1)='x'
          Vinfo( 2)='x-location of floats trajectories'
          Vinfo( 3)='meter'
          Vinfo(14)='x, scalar, series'
          Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
          Vinfo(24)='_FillValue'
          Aval(6)=spval
# endif
          status=def_var(ng, iNLM, FLT(ng)%ncid, FLT(ng)%Vid(idglon),   &
     &                   NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

          Vinfo( 1)='y'
          Vinfo( 2)='y-location of floats trajectories'
          Vinfo( 3)='meter'
          Vinfo(14)='y, scalar, series'
          Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
          Vinfo(24)='_FillValue'
          Aval(6)=spval
# endif
          status=def_var(ng, iNLM, FLT(ng)%ncid, FLT(ng)%Vid(idglat),   &
     &                   NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# ifdef SOLVE3D
!
!  Define floats depths.
!
        Vinfo( 1)='depth'
        Vinfo( 2)='depth of floats trajectories'
        Vinfo( 3)='meter'
        Vinfo(14)='depth, scalar, series'
        Vinfo(16)=Vname(1,idtime)
#  ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  endif
        status=def_var(ng, iNLM, FLT(ng)%ncid, FLT(ng)%Vid(iddpth),     &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define density anomaly.
!
        Vinfo( 1)=Vname(1,idDano)
        Vinfo( 2)=Vname(2,idDano)
        Vinfo( 3)=Vname(3,idDano)
        Vinfo(14)=Vname(4,idDano)
        Vinfo(16)=Vname(1,idtime)
#  ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  endif
        status=def_var(ng, iNLM, FLT(ng)%ncid, FLT(ng)%Vid(idDano),     &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          Vinfo( 1)=Vname(1,idTvar(itrc))
          Vinfo( 2)=Vname(2,idTvar(itrc))
          Vinfo( 3)=Vname(3,idTvar(itrc))
          Vinfo(14)=Vname(4,idTvar(itrc))
          Vinfo(16)=Vname(1,idtime)
#  ifndef NO_4BYTE_REALS
          Vinfo(24)='_FillValue'
          Aval(6)=spval
#  endif
#  ifdef SEDIMENT
          DO i=1,NST
            IF (itrc.eq.idsed(i)) THEN
              WRITE (Vinfo(19),40) 1000.0_r8*Sd50(i,ng)
            END IF
          END DO
#  endif
          status=def_var(ng, iNLM, FLT(ng)%ncid, FLT(ng)%Tid(itrc),     &
     &                   NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
# endif

# ifdef FLOAT_OYSTER
!
!  Define biological float swimming time.
!
        Vinfo( 1)='swim_time'
        Vinfo( 2)='biological float swimming time'
        Vinfo( 3)='s'
        Vinfo(14)='swim_time, scalar, series'
        Vinfo(16)=Vname(1,idtime)
#  ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  endif
        status=def_var(ng, iNLM, FLT(ng)%ncid, FLT(ng)%Vid(idswim),     &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define biological float vertical velocity.
!
        Vinfo( 1)='w_bio'
        Vinfo( 2)='biological float vertical velocity'
        Vinfo( 3)='m/s'
        Vinfo(14)='w_bio, scalar, series'
        Vinfo(16)=Vname(1,idtime)
#  ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  endif
        status=def_var(ng, iNLM, FLT(ng)%ncid, FLT(ng)%Vid(idwbio),     &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define biological float size (length).
!
        Vinfo( 1)='bio_size'
        Vinfo( 2)='biological float size'
        Vinfo( 3)='um'
        Vinfo(14)='bio_size, scalar, series'
        Vinfo(16)=Vname(1,idtime)
#  ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  endif
        status=def_var(ng, iNLM, FLT(ng)%ncid, FLT(ng)%Vid(idsize),     &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define biological float sinking velocity.
!
        Vinfo( 1)='bio_sink'
        Vinfo( 2)='biological float sinking velocity'
        Vinfo( 3)='m/s'
        Vinfo(14)='bio_sink, scalar, series'
        Vinfo(16)=Vname(1,idtime)
#  ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  endif
        status=def_var(ng, iNLM, FLT(ng)%ncid, FLT(ng)%Vid(idwsin),     &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
!  Initialize unlimited time record dimension.
!
        FLT(ng)%Rindex=0
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, iNLM, ncname, FLT(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, iNLM, FLT(ng)%ncid, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing floats file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=FLT(ng)%name
!
!  Open floats file for read/write.
!
        CALL netcdf_open (ng, iNLM, ncname, 1, FLT(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,50) TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, iNLM, ncname,                        &
     &                         ncid = FLT(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Get the size of the drifter dimension.
!
        DO i=1,n_dim
          IF (TRIM(dim_name(i)).eq.'drifter') THEN
            Nfloats(ng)=dim_size(i)
            EXIT
          END IF
        END DO
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, iNLM, ncname,                          &
     &                       ncid = FLT(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  float variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            FLT(ng)%Vid(idtime)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Xgrid') THEN
            got_var(idXgrd)=.TRUE.
            FLT(ng)%Vid(idXgrd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Ygrid') THEN
            got_var(idYgrd)=.TRUE.
            FLT(ng)%Vid(idYgrd)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.'Zgrid') THEN
            got_var(idZgrd)=.TRUE.
            FLT(ng)%Vid(idZgrd)=var_id(i)
# endif
          ELSE IF (spherical.and.TRIM(var_name(i)).eq.'lon') THEN
            got_var(idglon)=.TRUE.
            FLT(ng)%Vid(idglon)=var_id(i)
          ELSE IF (spherical.and.TRIM(var_name(i)).eq.'lat') THEN
            got_var(idglat)=.TRUE.
            FLT(ng)%Vid(idglat)=var_id(i)
          ELSE IF (.not.spherical.and.TRIM(var_name(i)).eq.'x') THEN
            got_var(idglon)=.TRUE.
            FLT(ng)%Vid(idglon)=var_id(i)
          ELSE IF (.not.spherical.and.TRIM(var_name(i)).eq.'y') THEN
            got_var(idglat)=.TRUE.
            FLT(ng)%Vid(idglat)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.'depth') THEN
            got_var(iddpth)=.TRUE.
            FLT(ng)%Vid(iddpth)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            FLT(ng)%Vid(idDano)=var_id(i)
# endif
# ifdef FLOAT_OYSTER
          ELSE IF (TRIM(var_name(i)).eq.'swim_time') THEN
            got_var(idswim)=.TRUE.
            FLT(ng)%Vid(idswim)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'w_bio') THEN
            got_var(idwbio)=.TRUE.
            FLT(ng)%Vid(idwbio)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'bio_size') THEN
            got_var(idsize)=.TRUE.
            FLT(ng)%Vid(idsize)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'bio_sink') THEN
            got_var(idwsin)=.TRUE.
            FLT(ng)%Vid(idwsin)=var_id(i)
# endif
          END IF
# ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              FLT(ng)%Tid(itrc)=var_id(i)
            END IF
          END DO
# endif
        END DO
!
!  Check if floats variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idXgrd)) THEN
          IF (Master) WRITE (stdout,60) 'Xgrid', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idYgrd)) THEN
          IF (Master) WRITE (stdout,60) 'Ygrid', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idZgrd)) THEN
          IF (Master) WRITE (stdout,60) 'Zgrid', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
        IF (.not.got_var(idglon)) THEN
          IF (spherical) THEN
            IF (Master) WRITE (stdout,60) 'lon', TRIM(ncname)
          ELSE
            IF (Master) WRITE (stdout,60) 'x', TRIM(ncname)
          END IF
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idglat)) THEN
          IF (spherical) THEN
            IF (Master) WRITE (stdout,60) 'lat', TRIM(ncname)
          ELSE
            IF (Master) WRITE (stdout,60) 'y', TRIM(ncname)
          END IF
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(iddpth)) THEN
          IF (Master) WRITE (stdout,60) 'depth', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDano)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef FLOAT_OYSTER
        IF (.not.got_var(idswim)) THEN
          IF (Master) WRITE (stdout,60) 'swim_time', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idwbio)) THEN
          IF (Master) WRITE (stdout,60) 'w_bio', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idsize)) THEN
          IF (Master) WRITE (stdout,60) 'bio_size', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idwsin)) THEN
          IF (Master) WRITE (stdout,60) 'bio_sink', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc))) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Initialize floats positions to the appropriate values.
!-----------------------------------------------------------------------
!
!  Set-up floats time record.
!
        IF (frrec(ng).lt.0) THEN
          FLT(ng)%Rindex=rec_size
        ELSE
          FLT(ng)%Rindex=ABS(frrec(ng))
        END IF
!
!  Read in floats nondimentional horizontal positions.  If the floats
!  have not been released yet at restart time, the values of Xgrid,
!  Ygrid, and Zgrid will be _FillValue (1.0E+37) in the FLOATS NetCDF
!  file. The calls to 'netcdf_get_fvar' will replace such values with
!  zero.  Therefore, we need to read Zgrid first so the bounded switch
!  is false in such cases tp trigger release.  Then, the bounded switch
!  is set correctly when reading Xgrid and/or Ygrid since the lower
!  bound is 0.5 in fractional coordinates.
!
# ifdef SOLVE3D
        CALL netcdf_get_fvar (ng, iNLM, ncname, 'Zgrid',                &
     &                        Tinp,                                     &
     &                        ncid = FLT(ng)%ncid,                      &
     &                        start = (/1,FLT(ng)%Rindex/),             &
     &                        total = (/Nfloats(ng),1/))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        DO l=1,Nfloats(ng)
          IF ((Tinp(l).gt.REAL(N(ng),r8)).or.                           &
     &        (Tinp(l).lt.0.0_r8)) THEN
            DRIFTER(ng)%bounded(l)=.FALSE.
          ELSE
            DRIFTER(ng)%bounded(l)=.TRUE.
            DO i=0,NFT
              DRIFTER(ng)%track(izgrd,i,l)=Tinp(l)
              DRIFTER(ng)%track(izrhs,i,l)=0.0_r8
            END DO
          END IF
        END DO
!
# endif
        CALL netcdf_get_fvar (ng, iNLM, ncname, 'Xgrid',                &
     &                        Tinp,                                     &
     &                        ncid = FLT(ng)%ncid,                      &
     &                        start = (/1,FLT(ng)%Rindex/),             &
     &                        total = (/Nfloats(ng),1/))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        DO l=1,Nfloats(ng)
          IF ((Tinp(l).gt.REAL(Lm(ng)+1,r8)-0.5_r8).or.                 &
     &        (Tinp(l).lt.0.5_r8)) THEN
            DRIFTER(ng)%bounded(l)=.FALSE.
          ELSE
            DRIFTER(ng)%bounded(l)=.TRUE.
            DO i=0,NFT
              DRIFTER(ng)%track(ixgrd,i,l)=Tinp(l)
              DRIFTER(ng)%track(ixrhs,i,l)=0.0_r8
            END DO
          END IF
        END DO
!
        CALL netcdf_get_fvar (ng, iNLM, ncname, 'Ygrid',                &
     &                        Tinp,                                     &
     &                        ncid = FLT(ng)%ncid,                      &
     &                        start = (/1,FLT(ng)%Rindex/),             &
     &                        total = (/Nfloats(ng),1/))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        DO l=1,Nfloats(ng)
          IF ((Tinp(l).gt.REAL(Mm(ng)+1,r8)-0.5_r8).or.                 &
     &        (Tinp(l).lt.0.5_r8)) THEN
            DRIFTER(ng)%bounded(l)=.FALSE.
          ELSE
            DRIFTER(ng)%bounded(l)=.TRUE.
            DO i=0,NFT
              DRIFTER(ng)%track(iygrd,i,l)=Tinp(l)
              DRIFTER(ng)%track(iyrhs,i,l)=0.0_r8
            END DO
          END IF
        END DO

      END IF QUERY
!
  10  FORMAT (2x,'DEF_FLOATS_NF90  - creating floats file,',t56,        &
     &        'Grid ',i2.2,': ',a)
  20  FORMAT (2x,'DEF_FLOATS_NF90  - inquiring floats file,',t56,       &
     &        'Grid ',i2.2,': ',a)
  30  FORMAT (/,' DEF_FLOATS_NF90 - unable to create floats NetCDF',    &
     &        ' file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_FLOATS_NF90 - unable to open floats NetCDF',      &
     &        ' file: ',a)
  60  FORMAT (/,' DEF_FLOATS_NF90 - unable to find variable: ',a,2x,    &
     &        ' in floats NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE def_floats_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE def_floats_pio (ng, ldef)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
      logical, intent(in) :: ldef
!
!  Local variable declarations.
!
      logical :: got_var(-6:NV)
!
      integer, parameter :: Natt = 25

      integer :: fltdim, i, itrc, j, l
      integer :: recdim, status

      integer :: DimIDs(nDimID)
      integer :: fgrd(2), start(2), total(2)
!
      real(r8) :: Aval(6), Tinp(Nfloats(ng))
!
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", def_floats_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=FLT(ng)%name
!
      IF (Master) THEN
        IF (ldef) THEN
          WRITE (stdout,10) ng, TRIM(ncname)
        ELSE
          WRITE (stdout,20) ng, TRIM(ncname)
        END IF
      END IF
!
!=======================================================================
!  Create a new floats data file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL pio_netcdf_create (ng, iNLM, TRIM(ncname), FLT(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,30) TRIM(ncname)
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define file dimensions.
!-----------------------------------------------------------------------
!
        DimIDs=0
!
#  ifdef SOLVE3D
        status=def_dim(ng, iNLM, FLT(ng)%pioFile, ncname, 's_rho',      &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, FLT(ng)%pioFile, ncname, 's_w',        &
     &                 N(ng)+1, DimIDs(10))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, FLT(ng)%pioFile, ncname, 'tracer',     &
     &                 NT(ng), DimIDs(11))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   ifdef SEDIMENT
        status=def_dim(ng, iNLM, FLT(ng)%pioFile, ncname, 'NST',        &
     &                 NST, DimIDs(32))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, FLT(ng)%pioFile, ncname, 'Nbed',       &
     &                 Nbed, DimIDs(16))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif

#   ifdef ECOSIM
        status=def_dim(ng, iNLM, FLT(ng)%pioFile, ncname, 'Nbands',     &
     &                 NBands, DimIDs(33))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, FLT(ng)%pioFile, ncname, 'Nphy',       &
     &                 Nphy, DimIDs(25))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, FLT(ng)%pioFile, ncname, 'Nbac',       &
     &                 Nbac, DimIDs(26))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, FLT(ng)%pioFile, ncname, 'Ndom',       &
     &                 Ndom, DimIDs(27))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, FLT(ng)%pioFile, ncname, 'Nfec',       &
     &                 Nfec, DimIDs(28))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif

        status=def_dim(ng, iNLM, FLT(ng)%pioFile, ncname, 'drifter' ,   &
     &                 Nfloats(ng), DimIDs(15))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, FLT(ng)%pioFile, ncname, 'boundary',   &
     &                 4, DimIDs(14))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef FOUR_DVAR
        status=def_dim(ng, iNLM, FLT(ng)%pioFile, ncname, 'Nstate',     &
     &                 NstateVar(ng), DimIDs(29))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

        status=def_dim(ng, iNLM, FLT(ng)%pioFile, ncname,               &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 PIO_unlimited, DimIDs(12))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        recdim=DimIDs(12)
        fltdim=DimIDs(15)
!
!  Define dimension vectors for point variables.
!
        fgrd(1)=DimIDs(15)
        fgrd(2)=DimIDs(12)
!
!  Initialize unlimited time record dimension.
!
        FLT(ng)%Rindex=0
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,LEN(Vinfo(1))
            Vinfo(i)(j:j)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
!-----------------------------------------------------------------------
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, iNLM, FLT(ng)%pioFile, ncname, DimIDs)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define variables and their attributes.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        Vinfo( 2)=Vname(2,idtime)
        WRITE (Vinfo( 3),'(a,a)') 'seconds since ', TRIM(Rclock%string)
        Vinfo( 4)=TRIM(Rclock%calendar)
        Vinfo(14)=Vname(4,idtime)
        FLT(ng)%pioVar(idtime)%dkind=PIO_TOUT
        FLT(ng)%pioVar(idtime)%gtype=0
!
        status=def_var(ng, iNLM, FLT(ng)%pioFile,                       &
     &                 FLT(ng)%pioVar(idtime)%vd,                       &
     &                 PIO_TOUT, 1, (/recdim/), Aval, Vinfo, ncname,    &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define floats X-grid locations.
!
        Vinfo( 1)='Xgrid'
        Vinfo( 2)='x-grid floats locations'
        Vinfo( 5)='valid_min'
        Vinfo( 6)='valid_max'
        Aval(2)=0.0_r8
        Aval(3)=REAL(Lm(ng)+1,r8)
        Vinfo(14)='Xgrid, scalar, series'
        Vinfo(16)=Vname(1,idtime)
#  ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  endif
        FLT(ng)%pioVar(idXgrd)%dkind=PIO_FOUT
        FLT(ng)%pioVar(idXgrd)%gtype=0
!
        status=def_var(ng, iNLM, FLT(ng)%pioFile,                       &
     &                 FLT(ng)%pioVar(idXgrd)%vd,                       &
     &                 PIO_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define floats Y-grid locations.
!
        Vinfo( 1)='Ygrid'
        Vinfo( 2)='Y-grid floats locations'
        Vinfo( 5)='valid_min'
        Vinfo( 6)='valid_max'
        Aval(2)=0.0_r8
        Aval(3)=REAL(Mm(ng)+1,r8)
        Vinfo(14)='Ygrid, scalar, series'
        Vinfo(16)=Vname(1,idtime)
#  ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  endif
        FLT(ng)%pioVar(idYgrd)%dkind=PIO_FOUT
        FLT(ng)%pioVar(idYgrd)%gtype=0
!
        status=def_var(ng, iNLM, FLT(ng)%pioFile,                       &
     &                 FLT(ng)%pioVar(idYgrd)%vd,                       &
     &                 PIO_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SOLVE3D
!
!  Define floats Z-grid locations.
!
        Vinfo( 1)='Zgrid'
        Vinfo( 2)='Z-grid floats locations'
        Vinfo( 5)='valid_min'
        Vinfo( 6)='valid_max'
        Aval(2)=0.0_r8
        Aval(3)=REAL(N(ng),r8)
        Vinfo(14)='Zgrid, scalar, series'
        Vinfo(16)=Vname(1,idtime)
#   ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#   endif
        FLT(ng)%pioVar(idZgrd)%dkind=PIO_FOUT
        FLT(ng)%pioVar(idZgrd)%gtype=0
!
        status=def_var(ng, iNLM, FLT(ng)%pioFile,                       &
     &                 FLT(ng)%pioVar(idZgrd)%vd,                       &
     &                 PIO_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
!
!  Define floats (lon,lat) or (x,y) locations.
!
        IF (spherical) THEN
          Vinfo( 1)='lon'
          Vinfo( 2)='longitude of floats trajectories'
          Vinfo( 3)='degree_east'
          Vinfo( 5)='valid_min'
          Vinfo( 6)='valid_max'
          Vinfo(14)='lon, scalar, series'
          Vinfo(16)=Vname(1,idtime)
#  ifndef NO_4BYTE_REALS
          Vinfo(24)='_FillValue'
          Aval(6)=spval
#  endif
          Aval(2)=-180.0_r8
          Aval(3)=180.0_r8
          FLT(ng)%pioVar(idglon)%dkind=PIO_FOUT
          FLT(ng)%pioVar(idglon)%gtype=0
!
          status=def_var(ng, iNLM, FLT(ng)%pioFile,                     &
     &                   FLT(ng)%pioVar(idglon)%vd,                     &
     &                   PIO_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          Vinfo( 1)='lat'
          Vinfo( 2)='latitude of floats trajectories'
          Vinfo( 3)='degree_north'
          Vinfo( 5)='valid_min'
          Vinfo( 6)='valid_max'
          Vinfo(14)='lat, scalar, series'
          Vinfo(16)=Vname(1,idtime)
#  ifndef NO_4BYTE_REALS
          Vinfo(24)='_FillValue'
          Aval(6)=spval
#  endif
          Aval(2)=-90.0_r8
          Aval(3)=90.0_r8
          FLT(ng)%pioVar(idglat)%dkind=PIO_FOUT
          FLT(ng)%pioVar(idglat)%gtype=0
!
          status=def_var(ng, iNLM, FLT(ng)%pioFile,                     &
     &                   FLT(ng)%pioVar(idglat)%vd,                     &
     &                   PIO_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
       ELSE
          Vinfo( 1)='x'
          Vinfo( 2)='x-location of floats trajectories'
          Vinfo( 3)='meter'
          Vinfo(14)='x, scalar, series'
          Vinfo(16)=Vname(1,idtime)
#  ifndef NO_4BYTE_REALS
          Vinfo(24)='_FillValue'
          Aval(6)=spval
#  endif
          FLT(ng)%pioVar(idglon)%dkind=PIO_FOUT
          FLT(ng)%pioVar(idglon)%gtype=0
!
          status=def_var(ng, iNLM, FLT(ng)%pioFile,                     &
     &                   FLT(ng)%pioVar(idglon)%vd,                     &
     &                   PIO_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          Vinfo( 1)='y'
          Vinfo( 2)='y-location of floats trajectories'
          Vinfo( 3)='meter'
          Vinfo(14)='y, scalar, series'
          Vinfo(16)=Vname(1,idtime)
#  ifndef NO_4BYTE_REALS
          Vinfo(24)='_FillValue'
          Aval(6)=spval
#  endif
          FLT(ng)%pioVar(idglat)%dkind=PIO_FOUT
          FLT(ng)%pioVar(idglat)%gtype=0
!
          status=def_var(ng, iNLM, FLT(ng)%pioFile,                     &
     &                   FLT(ng)%pioVar(idglat)%vd,                     &
     &                   PIO_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  ifdef SOLVE3D
!
!  Define floats depths.
!
        Vinfo( 1)='depth'
        Vinfo( 2)='depth of floats trajectories'
        Vinfo( 3)='meter'
        Vinfo(14)='depth, scalar, series'
        Vinfo(16)=Vname(1,idtime)
#   ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#   endif
        FLT(ng)%pioVar(iddpth)%dkind=PIO_FOUT
        FLT(ng)%pioVar(iddpth)%gtype=0
!
        status=def_var(ng, iNLM, FLT(ng)%pioFile,                       &
     &                 FLT(ng)%pioVar(iddpth)%vd,                       &
     &                 PIO_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define density anomaly.
!
        Vinfo( 1)=Vname(1,idDano)
        Vinfo( 2)=Vname(2,idDano)
        Vinfo( 3)=Vname(3,idDano)
        Vinfo(14)=Vname(4,idDano)
        Vinfo(16)=Vname(1,idtime)
#   ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#   endif
        FLT(ng)%pioVar(idDano)%dkind=PIO_FOUT
        FLT(ng)%pioVar(idDano)%gtype=0
!
        status=def_var(ng, iNLM, FLT(ng)%pioFile,                       &
     &                 FLT(ng)%pioVar(idDano)%vd,                       &
     &                 PIO_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          Vinfo( 1)=Vname(1,idTvar(itrc))
          Vinfo( 2)=Vname(2,idTvar(itrc))
          Vinfo( 3)=Vname(3,idTvar(itrc))
          Vinfo(14)=Vname(4,idTvar(itrc))
          Vinfo(16)=Vname(1,idtime)
#   ifndef NO_4BYTE_REALS
          Vinfo(24)='_FillValue'
          Aval(6)=spval
#   endif
#   ifdef SEDIMENT
          DO i=1,NST
            IF (itrc.eq.idsed(i)) THEN
              WRITE (Vinfo(19),40) 1000.0_r8*Sd50(i,ng)
            END IF
          END DO
#   endif
          FLT(ng)%pioTrc(itrc)%dkind=PIO_FOUT
          FLT(ng)%pioTrc(itrc)%gtype=0
!
          status=def_var(ng, iNLM, FLT(ng)%pioFile,                     &
     &                   FLT(ng)%pioTrc(itrc)%vd,                       &
     &                   PIO_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
#  endif

#  ifdef FLOAT_OYSTER
!
!  Define biological float swimming time.
!
        Vinfo( 1)='swim_time'
        Vinfo( 2)='biological float swimming time'
        Vinfo( 3)='s'
        Vinfo(14)='swim_time, scalar, series'
        Vinfo(16)=Vname(1,idtime)
#   ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#   endif
        FLT(ng)%pioVar(idswim)%dkind=PIO_FOUT
        FLT(ng)%pioVar(idswim)%gtype=0
!
        status=def_var(ng, iNLM, FLT(ng)%pioFile,                       &
     &                 FLT(ng)%pioVar(idswim)%vd,                       &
     &                 PIO_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define biological float vertical velocity.
!
        Vinfo( 1)='w_bio'
        Vinfo( 2)='biological float vertical velocity'
        Vinfo( 3)='m/s'
        Vinfo(14)='w_bio, scalar, series'
        Vinfo(16)=Vname(1,idtime)
#   ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#   endif
        FLT(ng)%pioVar(idwbio)%dkind=PIO_FOUT
        FLT(ng)%pioVar(idwbio)%gtype=0
!
        status=def_var(ng, iNLM, FLT(ng)%pioFile,                       &
     &                 FLT(ng)%pioVar(idwbio)%vd,                       &
     &                 PIO_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define biological float size (length).
!
        Vinfo( 1)='bio_size'
        Vinfo( 2)='biological float size'
        Vinfo( 3)='um'
        Vinfo(14)='bio_size, scalar, series'
        Vinfo(16)=Vname(1,idtime)
#   ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#   endif
        FLT(ng)%pioVar(idsize)%dkind=PIO_FOUT
        FLT(ng)%pioVar(idsize)%gtype=0
!
        status=def_var(ng, iNLM, FLT(ng)%pioFile,                       &
     &                 FLT(ng)%pioVar(idsize)%vd,                       &
     &                 PIO_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define biological float sinking velocity.
!
        Vinfo( 1)='bio_sink'
        Vinfo( 2)='biological float sinking velocity'
        Vinfo( 3)='m/s'
        Vinfo(14)='bio_sink, scalar, series'
        Vinfo(16)=Vname(1,idtime)
#   ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#   endif
        FLT(ng)%pioVar(idwsin)%dkind=PIO_FOUT
        FLT(ng)%pioVar(idwsin)%gtype=0
!
        status=def_var(ng, iNLM, FLT(ng)%pioFile,                       &
     &                 FLT(ng)%pioVar(idwsin)%vd,                       &
     &                 PIO_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
!
!  Initialize unlimited time record dimension.
!
        FLT(ng)%Rindex=0
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL pio_netcdf_enddef (ng, iNLM, ncname, FLT(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, iNLM, FLT(ng)%pioFile, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing floats file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=FLT(ng)%name
!
!  Open floats file for read/write.
!
        CALL pio_netcdf_open (ng, iNLM, ncname, 1, FLT(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,50) TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL pio_netcdf_check_dim (ng, iNLM, ncname,                    &
     &                             pioFile = FLT(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Get the size of the drifter dimension.
!
        DO i=1,n_dim
          IF (TRIM(dim_name(i)).eq.'drifter') THEN
            Nfloats(ng)=dim_size(i)
            EXIT
          END IF
        END DO
!
!  Inquire about the variables.
!
        CALL pio_netcdf_inq_var (ng, iNLM, ncname,                      &
     &                           pioFile = FLT(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  float variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            FLT(ng)%pioVar(idtime)%vd=var_desc(i)
            FLT(ng)%pioVar(idtime)%dkind=PIO_TOUT
            FLT(ng)%pioVar(idtime)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.'Xgrid') THEN
            got_var(idXgrd)=.TRUE.
            FLT(ng)%pioVar(idXgrd)%vd=var_desc(i)
            FLT(ng)%pioVar(idXgrd)%dkind=PIO_FOUT
            FLT(ng)%pioVar(idXgrd)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.'Ygrid') THEN
            got_var(idYgrd)=.TRUE.
            FLT(ng)%pioVar(idYgrd)%vd=var_desc(i)
            FLT(ng)%pioVar(idYgrd)%dkind=PIO_FOUT
            FLT(ng)%pioVar(idYgrd)%gtype=0
#  ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.'Zgrid') THEN
            got_var(idZgrd)=.TRUE.
            FLT(ng)%pioVar(idZgrd)%vd=var_desc(i)
            FLT(ng)%pioVar(idZgrd)%dkind=PIO_FOUT
            FLT(ng)%pioVar(idZgrd)%gtype=0
#  endif
          ELSE IF (spherical.and.TRIM(var_name(i)).eq.'lon') THEN
            got_var(idglon)=.TRUE.
            FLT(ng)%pioVar(idglon)%vd=var_desc(i)
            FLT(ng)%pioVar(idglon)%dkind=PIO_FOUT
            FLT(ng)%pioVar(idglon)%gtype=0
          ELSE IF (spherical.and.TRIM(var_name(i)).eq.'lat') THEN
            got_var(idglat)=.TRUE.
            FLT(ng)%pioVar(idglat)%vd=var_desc(i)
            FLT(ng)%pioVar(idglat)%dkind=PIO_FOUT
            FLT(ng)%pioVar(idglat)%gtype=0
          ELSE IF (.not.spherical.and.TRIM(var_name(i)).eq.'x') THEN
            got_var(idglon)=.TRUE.
            FLT(ng)%pioVar(idglon)%vd=var_desc(i)
            FLT(ng)%pioVar(idglon)%dkind=PIO_FOUT
            FLT(ng)%pioVar(idglon)%gtype=0
          ELSE IF (.not.spherical.and.TRIM(var_name(i)).eq.'y') THEN
            got_var(idglat)=.TRUE.
            FLT(ng)%pioVar(idglat)%vd=var_desc(i)
            FLT(ng)%pioVar(idglat)%dkind=PIO_FOUT
            FLT(ng)%pioVar(idglat)%gtype=0
#  ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.'depth') THEN
            got_var(iddpth)=.TRUE.
            FLT(ng)%pioVar(iddpth)%vd=var_desc(i)
            FLT(ng)%pioVar(iddpth)%dkind=PIO_FOUT
            FLT(ng)%pioVar(iddpth)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            FLT(ng)%pioVar(idDano)%vd=var_desc(i)
            FLT(ng)%pioVar(idDano)%dkind=PIO_FOUT
            FLT(ng)%pioVar(idDano)%gtype=0
#  endif
#  ifdef FLOAT_OYSTER
          ELSE IF (TRIM(var_name(i)).eq.'swim_time') THEN
            got_var(idswim)=.TRUE.
            FLT(ng)%pioVar(idswim)%vd=var_desc(i)
            FLT(ng)%pioVar(idswim)%dkind=PIO_FOUT
            FLT(ng)%pioVar(idswim)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.'w_bio') THEN
            got_var(idwbio)=.TRUE.
            FLT(ng)%pioVar(idwbio)%vd=var_desc(i)
            FLT(ng)%pioVar(idwbio)%dkind=PIO_FOUT
            FLT(ng)%pioVar(idwbio)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.'bio_size') THEN
            got_var(idsize)=.TRUE.
            FLT(ng)%pioVar(idsize)%vd=var_desc(i)
            FLT(ng)%pioVar(idsize)%dkind=PIO_FOUT
            FLT(ng)%pioVar(idsize)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.'bio_sink') THEN
            got_var(idwsin)=.TRUE.
            FLT(ng)%pioVar(idwsin)%vd=var_desc(i)
            FLT(ng)%pioVar(idwsin)%dkind=PIO_FOUT
            FLT(ng)%pioVar(idwsin)%gtype=0
#  endif
          END IF
#  ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              FLT(ng)%pioTrc(itrc)%vd=var_desc(i)
              FLT(ng)%pioTrc(itrc)%dkind=PIO_FOUT
              FLT(ng)%pioTrc(itrc)%gtype=0
            END IF
          END DO
#  endif
        END DO
!
!  Check if floats variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idXgrd)) THEN
          IF (Master) WRITE (stdout,60) 'Xgrid', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idYgrd)) THEN
          IF (Master) WRITE (stdout,60) 'Ygrid', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef SOLVE3D
        IF (.not.got_var(idZgrd)) THEN
          IF (Master) WRITE (stdout,60) 'Zgrid', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idglon)) THEN
          IF (spherical) THEN
            IF (Master) WRITE (stdout,60) 'lon', TRIM(ncname)
          ELSE
            IF (Master) WRITE (stdout,60) 'x', TRIM(ncname)
          END IF
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idglat)) THEN
          IF (spherical) THEN
            IF (Master) WRITE (stdout,60) 'lat', TRIM(ncname)
          ELSE
            IF (Master) WRITE (stdout,60) 'y', TRIM(ncname)
          END IF
          exit_flag=3
          RETURN
        END IF
#  ifdef SOLVE3D
        IF (.not.got_var(iddpth)) THEN
          IF (Master) WRITE (stdout,60) 'depth', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDano)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef FLOAT_OYSTER
        IF (.not.got_var(idswim)) THEN
          IF (Master) WRITE (stdout,60) 'swim_time', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idwbio)) THEN
          IF (Master) WRITE (stdout,60) 'w_bio', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idsize)) THEN
          IF (Master) WRITE (stdout,60) 'bio_size', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idwsin)) THEN
          IF (Master) WRITE (stdout,60) 'bio_sink', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc))) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  endif
!
!-----------------------------------------------------------------------
!  Initialize floats positions to the appropriate values.
!-----------------------------------------------------------------------
!
!  Set-up floats time record.
!
        IF (frrec(ng).lt.0) THEN
          FLT(ng)%Rindex=rec_size
        ELSE
          FLT(ng)%Rindex=ABS(frrec(ng))
        END IF
!
!  Read in floats nondimentional horizontal positions.  If the floats
!  have not been released yet at restart time, the values of Xgrid,
!  Ygrid, and Zgrid will be _FillValue (1.0E+37) in the FLOATS NetCDF
!  file. The calls to 'netcdf_get_fvar' will replace such values with
!  zero.  Therefore, we need to read Zgrid first so the bounded switch
!  is false in such cases tp trigger release.  Then, the bounded switch
!  is set correctly when reading Xgrid and/or Ygrid since the lower
!  bound is 0.5 in fractional coordinates.
!
#  ifdef SOLVE3D
        CALL pio_netcdf_get_fvar (ng, iNLM, ncname, 'Zgrid',            &
     &                            Tinp,                                 &
     &                            pioFile = FLT(ng)%pioFile,            &
     &                            start = (/1,FLT(ng)%Rindex/),         &
     &                            total = (/Nfloats(ng),1/))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        DO l=1,Nfloats(ng)
          IF ((Tinp(l).gt.REAL(N(ng),r8)).or.                           &
     &        (Tinp(l).lt.0.0_r8)) THEN
            DRIFTER(ng)%bounded(l)=.FALSE.
          ELSE
            DRIFTER(ng)%bounded(l)=.TRUE.
            DO i=0,NFT
              DRIFTER(ng)%track(izgrd,i,l)=Tinp(l)
              DRIFTER(ng)%track(izrhs,i,l)=0.0_r8
            END DO
          END IF
        END DO
!
#  endif
        CALL pio_netcdf_get_fvar (ng, iNLM, ncname, 'Xgrid',            &
     &                            Tinp,                                 &
     &                            pioFile = FLT(ng)%pioFile,            &
     &                            start = (/1,FLT(ng)%Rindex/),         &
     &                            total = (/Nfloats(ng),1/))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        DO l=1,Nfloats(ng)
          IF ((Tinp(l).gt.REAL(Lm(ng)+1,r8)-0.5_r8).or.                 &
     &        (Tinp(l).lt.0.5_r8)) THEN
            DRIFTER(ng)%bounded(l)=.FALSE.
          ELSE
            DRIFTER(ng)%bounded(l)=.TRUE.
            DO i=0,NFT
              DRIFTER(ng)%track(ixgrd,i,l)=Tinp(l)
              DRIFTER(ng)%track(ixrhs,i,l)=0.0_r8
            END DO
          END IF
        END DO
!
        CALL pio_netcdf_get_fvar (ng, iNLM, ncname, 'Ygrid',            &
     &                            Tinp,                                 &
     &                            pioFile = FLT(ng)%pioFile,            &
     &                            start = (/1,FLT(ng)%Rindex/),         &
     &                            total = (/Nfloats(ng),1/))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        DO l=1,Nfloats(ng)
          IF ((Tinp(l).gt.REAL(Mm(ng)+1,r8)-0.5_r8).or.                 &
     &        (Tinp(l).lt.0.5_r8)) THEN
            DRIFTER(ng)%bounded(l)=.FALSE.
          ELSE
            DRIFTER(ng)%bounded(l)=.TRUE.
            DO i=0,NFT
              DRIFTER(ng)%track(iygrd,i,l)=Tinp(l)
              DRIFTER(ng)%track(iyrhs,i,l)=0.0_r8
            END DO
          END IF
        END DO

      END IF QUERY
!
  10  FORMAT (2x,'DEF_FLOATS_PIO   - creating floats file,',t56,        &
     &        'Grid ',i2.2,': ',a)
  20  FORMAT (2x,'DEF_FLOATS_PIO   - inquiring floats file,',t56,       &
     &        'Grid ',i2.2,': ',a)
  30  FORMAT (/,' DEF_FLOATS_PIO - unable to create floats NetCDF',     &
     &        ' file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_FLOATS_PIO - unable to open floats NetCDF',       &
     &        ' file: ',a)
  60  FORMAT (/,' DEF_FLOATS_PIO - unable to find variable: ',a,2x,     &
     &        ' in floats NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE def_floats_pio
# endif
#endif
      END MODULE def_floats_mod
