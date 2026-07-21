#include "cppdefs.h"
      MODULE sediment_output_mod

#if (defined SEDIMENT || defined BBL_MODEL) && defined SOLVE3D
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module defines/writes Waves Effect on Currents variables into  !
!  output NetCDF files.                                                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef AVERAGES
      USE mod_average
# endif
# ifdef BBL_MODEL
      USE mod_bbl
# endif
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
      USE mod_sedbed
# if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
# endif
      USE mod_stepping
!
      USE def_var_mod,     ONLY : def_var
# ifdef STATIONS
      USE extract_sta_mod, ONLY : extract_sta2d
#  ifdef SOLVE3D
      USE extract_sta_mod, ONLY : extract_sta3d
#  endif
# endif
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
# ifdef SOLVE3D
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
      USE omega_mod,       ONLY : scale_omega
# endif
      USE strings_mod,     ONLY : FoundError
!
      implicit none
!
      PUBLIC :: sediment_def_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: sediment_def_pio
# endif
# ifdef STATIONS
      PUBLIC :: sediment_def_station_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: sediment_def_station_pio
#  endif
# endif
      PUBLIC :: sediment_wrt_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: sediment_wrt_pio
# endif
# ifdef STATIONS
      PUBLIC :: sediment_wrt_station_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: sediment_wrt_station_pio
#  endif
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE sediment_def_nf90 (ng, model, ldef, VarOut, S,         &
     &                              t2dgrd, u2dgrd, v2dgrd,             &
     &                              b3dgrd)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef, VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: t2dgrd(:), u2dgrd(:), v2dgrd(:)
      integer, intent(in), optional :: b3dgrd(:)
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer :: i, itrc, j, nvd3, nvd4, status
!
      real(r8) :: Aval(6)
!
# ifdef ADJOINT
      character (len=21)  :: Prefix
# else
      character (len=13)  :: Prefix
# endif
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", sediment_def_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define sediment output variables.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
!
      DEFINE : IF (ldef) THEN
!
!  Set number of dimensions for output variables.
!
# if defined WRITE_WATER && defined MASKING
        nvd3=2
        nvd4=2
# else
        nvd3=3
        nvd4=4
# endif
!
!  Set long name prefix string.
!
# ifdef ADJOINT
!!      Prefix='time-averaged adjoint'
        Prefix='adjoint'
# else
!!      Prefix='time-averaged'
        Prefix=CHAR(32)                                   ! blank space
# endif
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

# if defined SEDIMENT && defined SED_MORPH
!
!  Define time-varying bathymetry.
!
        IF (VarOut(idbath,ng)) THEN
          Vinfo( 1)=Vname(1,idbath)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idbath))
          ELSE
            Vinfo( 3)=Vname(3,idbath)
          END IF
          Vinfo(14)=Vname(4,idbath)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idbath)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idbath,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idbath),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif

# if defined SEDIMENT && defined BEDLOAD
!
!  Define Bedload transport U-direction.
!
        DO i=1,NST
          IF (VarOut(idUbld(i),ng)) THEN
            Vinfo( 1)=Vname(1,idUbld(i))
            IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                      &
     &                                     TRIM(Vname(2,idUbld(i)))
            ELSE
              Vinfo( 2)=Vname(2,idUbld(i))
            END IF
            Vinfo( 3)=Vname(3,idUbld(i))
            Vinfo(14)=Vname(4,idUbld(i))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_u'
#  endif
            Vinfo(21)=Vname(6,idUbld(i))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idUbld(i),ng),r8)
            status=def_var(ng, model, S(ng)%ncid,                       &
     &                     S(ng)%Vid(idUbld(i)), NF_FOUT,               &
     &                     nvd3, u2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!  Define Bedload transport V-direction.
!
          IF (VarOut(idVbld(i),ng)) THEN
            Vinfo( 1)=Vname(1,idVbld(i))
            IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                      &
     &                                     TRIM(Vname(2,idVbld(i)))
            ELSE
              Vinfo( 2)=Vname(2,idVbld(i))
            END IF
            Vinfo( 3)=Vname(3,idVbld(i))
            Vinfo(14)=Vname(4,idVbld(i))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_v'
#  endif
            Vinfo(21)=Vname(6,idVbld(i))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idVbld(i),ng),r8)
            status=def_var(ng, model, S(ng)%ncid,                       &
     &                     S(ng)%Vid(idVbld(i)), NF_FOUT,               &
     &                     nvd3, v2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
# endif

# ifdef SEDIMENT
!
!  Define sediment fraction of each size class in each bed layer.
!
        DO i=1,NST
          IF (VarOut(idfrac(i),ng)) THEN
            Vinfo( 1)=Vname(1,idfrac(i))
            IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                      &
     &                                     TRIM(Vname(2,idfrac(i)))
            ELSE
              Vinfo( 2)=Vname(2,idfrac(i))
            END IF
            Vinfo( 3)=Vname(3,idfrac(i))
            Vinfo(14)=Vname(4,idfrac(i))
            Vinfo(16)=Vname(1,idtime)
            WRITE (Vinfo(19),10) 1000.0_r8*Sd50(i,ng)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(21)=Vname(6,idfrac(i))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idfrac(i),ng))
            status=def_var(ng, model, S(ng)%ncid,                       &
     &                     S(ng)%Vid(idfrac(i)), NF_FOUT,               &
     &                     nvd4, b3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define sediment mass of each size class in each bed layer.
!
        DO i=1,NST
          IF (VarOut(idBmas(i),ng)) THEN
            Vinfo( 1)=Vname(1,idBmas(i))
            IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                      &
     &                                     TRIM(Vname(2,idBmas(i)))
            ELSE
              Vinfo( 2)=Vname(2,idBmas(i))
            END IF
            Vinfo( 3)=Vname(3,idBmas(i))
            Vinfo(14)=Vname(4,idBmas(i))
            Vinfo(16)=Vname(1,idtime)
            WRITE (Vinfo(19),10) 1000.0_r8*Sd50(i,ng)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(21)=Vname(6,idBmas(i))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idBmas(i),ng),r8)
            status=def_var(ng, model, S(ng)%ncid,                       &
     &                     S(ng)%Vid(idBmas(i)), NF_FOUT,               &
     &                     nvd4, b3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define sediment properties in each bed layer.
!
        DO i=1,MBEDP
          IF (VarOut(idSbed(i),ng)) THEN
            Vinfo( 1)=Vname(1,idSbed(i))
            IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                      &
     &                                     TRIM(Vname(2,idSbed(i)))
            ELSE
              Vinfo( 2)=Vname(2,idSbed(i))
            END IF
            Vinfo( 3)=Vname(3,idSbed(i))
            Vinfo(14)=Vname(4,idSbed(i))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(21)=Vname(6,idSbed(i))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idSbed(i),ng),r8)
            status=def_var(ng, model, S(ng)%ncid,                       &
     &                     S(ng)%Vid(idSbed(i)), NF_FOUT,               &
     &                     nvd4, b3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
# endif

# if defined SEDIMENT || defined BBL_MODEL
!
!  Define exposed sediment layer properties.
!
        DO i=1,MBOTP
          IF (VarOut(idBott(i),ng)) THEN
            Vinfo( 1)=Vname(1,idBott(i))
            IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                      &
     &                                     TRIM(Vname(2,idBott(i)))
            ELSE
              Vinfo( 2)=Vname(2,idBott(i))
            END IF
            Vinfo( 3)=Vname(3,idBott(i))
            Vinfo(14)=Vname(4,idBott(i))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(21)=Vname(6,idBott(i))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idBott(i),ng),r8)
            status=def_var(ng, model, S(ng)%ncid,                       &
     &                     S(ng)%Vid(idBott(i)), NF_FOUT,               &
     &                     nvd3, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
# endif

# if defined SEDIMENT && defined BEDLOAD && defined BEDLOAD_VANDERA
!
!  Define Ursell number from skewed waveform parameters.
!
        IF (VarOut(idsurs,ng)) THEN
          Vinfo( 1)=Vname(1,idsurs)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idsurs))
          ELSE
            Vinfo( 2)=Vname(2,idsurs)
          END IF
          Vinfo( 3)=Vname(3,idsurs)
          Vinfo(14)=Vname(4,idsurs)
          Vinfo(16)=Vname(1,idsurs)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idsurs)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idsurs,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idsurs),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wave velocity skewness coefficient for asymmetric bedload.
!
        IF (VarOut(idsrrw,ng)) THEN
          Vinfo( 1)=Vname(1,idsrrw)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idsrrw))
          ELSE
            Vinfo( 2)=Vname(2,idsrrw)
          END IF
          Vinfo( 3)=Vname(3,idsrrw)
          Vinfo(14)=Vname(4,idsrrw)
          Vinfo(16)=Vname(1,idsrrw)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idsrrw)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idsrrw,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idsrrw),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wave acceleration skewness coefficient for asymmetric bedload.
!
        IF (VarOut(idsbtw,ng)) THEN
          Vinfo( 1)=Vname(1,idsbtw)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idsbtw))
          ELSE
            Vinfo( 2)=Vname(2,idsbtw)
          END IF
          Vinfo( 3)=Vname(3,idsbtw)
          Vinfo(14)=Vname(4,idsbtw)
          Vinfo(16)=Vname(1,idsbtw)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idsbtw)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idsbtw,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idsbtw),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define peak wave crest cycle velocity for asymmetric bedload.
!
        IF (VarOut(idsucr,ng)) THEN
          Vinfo( 1)=Vname(1,idsucr)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idsucr))
          ELSE
            Vinfo( 2)=Vname(2,idsucr)
          END IF
          Vinfo( 3)=Vname(3,idsucr)
          Vinfo(14)=Vname(4,idsucr)
          Vinfo(16)=Vname(1,idsucr)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idsucr)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idsucr,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idsucr),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define peak wave trough cycle velocity for asymmetric bedload.
!
        IF (VarOut(idsutr,ng)) THEN
          Vinfo( 1)=Vname(1,idsutr)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idsutr))
          ELSE
            Vinfo( 2)=Vname(2,idsutr)
          END IF
          Vinfo( 3)=Vname(3,idsutr)
          Vinfo(14)=Vname(4,idsutr)
          Vinfo(16)=Vname(1,idsutr)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idsutr)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idsutr,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idsutr),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wave crest cycle period for asymmetric bedload.
!
        IF (VarOut(idstcr,ng)) THEN
          Vinfo( 1)=Vname(1,idstcr)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idstcr))
          ELSE
            Vinfo( 2)=Vname(2,idstcr)
          END IF
          Vinfo( 3)=Vname(3,idstcr)
          Vinfo(14)=Vname(4,idstcr)
          Vinfo(16)=Vname(1,idstcr)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idstcr)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idstcr,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idstcr),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wave trough cycle period for asymmetric bedload.
!
        IF (VarOut(idsttr,ng)) THEN
          Vinfo( 1)=Vname(1,idsttr)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idsttr))
          ELSE
            Vinfo( 2)=Vname(2,idsttr)
          END IF
          Vinfo( 3)=Vname(3,idsttr)
          Vinfo(14)=Vname(4,idsttr)
          Vinfo(16)=Vname(1,idsttr)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idsttr)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idsttr,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idsttr),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif

      END IF DEFINE
!
!-----------------------------------------------------------------------
!  Otherwise, check existing output file and prepare for appending
!  data.
!-----------------------------------------------------------------------
!
      QUERY : IF (.not.ldef) THEN
!
!  Initialize local logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  Waves Effect on Currents variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            S(ng)%Vid(idtime)=var_id(i)
# if defined SEDIMENT && defined SED_MORPH
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idbath))) THEN
            got_var(idbath)=.TRUE.
            S(ng)%Vid(idbath)=var_id(i)
# endif
# if defined SEDIMENT && defined BEDLOAD && defined BEDLOAD_VANDERA
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idsurs))) THEN
            got_var(idsurs)=.TRUE.
            S(ng)%Vid(idsurs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idsrrw))) THEN
            got_var(idsrrw)=.TRUE.
            S(ng)%Vid(idsrrw)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idsbtw))) THEN
            got_var(idsbtw)=.TRUE.
            S(ng)%Vid(idsbtw)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idsucr))) THEN
            got_var(idsucr)=.TRUE.
            S(ng)%Vid(idsucr)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idsutr))) THEN
            got_var(idsutr)=.TRUE.
            S(ng)%Vid(idsutr)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idstcr))) THEN
            got_var(idstcr)=.TRUE.
            S(ng)%Vid(idstcr)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idsttr))) THEN
            got_var(idsttr)=.TRUE.
            S(ng)%Vid(idsttr)=var_id(i)
# endif
          END IF
# ifdef SEDIMENT
          DO itrc=1,NST
            IF (TRIM(var_name(i)).eq.                                   &
     &               TRIM(Vname(1,idfrac(itrc)))) THEN
              got_var(idfrac(itrc))=.TRUE.
              S(ng)%Vid(idfrac(itrc))=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idBmas(itrc)))) THEN
              got_var(idBmas(itrc))=.TRUE.
              S(ng)%Vid(idBmas(itrc))=var_id(i)
#  ifdef BEDLOAD
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idUbld(itrc)))) THEN
              got_var(idUbld(itrc))=.true.
              S(ng)%Vid(idUbld(itrc))=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idVbld(itrc)))) THEN
              got_var(idVbld(itrc))=.true.
              S(ng)%Vid(idVbld(itrc))=var_id(i)
#  endif
            END IF
          END DO
          DO itrc=1,MBEDP
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSbed(itrc)))) THEN
              got_var(idSbed(itrc))=.TRUE.
              S(ng)%Vid(idSbed(itrc))=var_id(i)
            END IF
          END DO
# endif
# if defined SEDIMENT || defined BBL_MODEL
          DO itrc=1,MBOTP
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idBott(itrc)))) THEN
              got_var(idBott(itrc))=.TRUE.
              S(ng)%Vid(idBott(itrc))=var_id(i)
            END IF
          END DO
# endif
        END DO
!
!  Check if output variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# if defined SEDIMENT && defined SED_MORPH
        IF (.not.got_var(idbath).and.VarOut(idbath,ng)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idbath)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined SEDIMENT && defined BEDLOAD && defined BEDLOAD_VANDERA
        IF (.not.got_var(idsurs).and.VarOut(idsurs,ng)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idsurs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idsrrw).and.VarOut(idsrrw,ng)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idsrrw)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idsbtw).and.VarOut(idsbtw,ng)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idsbtw)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idsucr).and.VarOut(idsucr,ng)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idsucr)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idsutr).and.VarOut(idsutr,ng)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idsutr)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idstcr).and.VarOut(idstcr,ng)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idstcr)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idsttr).and.VarOut(idsttr,ng)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idsttr)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef SEDIMENT
        DO i=1,NST
          IF (.not.got_var(idfrac(i)).and.VarOut(idfrac(i),ng)) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idfrac(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF

          IF(.not.got_var(idBmas(i)).and.VarOut(idBmas(i),ng)) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idBmas(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#  ifdef BEDLOAD
          IF (.not.got_var(idUbld(i)).and.VarOut(idUbld(i),ng)) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idUbld(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idVbld(i)).and.VarOut(idVbld(i),ng)) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idVbld(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#  endif
        END DO
        DO i=1,MBEDP
          IF (.not.got_var(idSbed(i)).and.VarOut(idSbed(i),ng)) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idSbed(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
# if defined SEDIMENT || defined BBL_MODEL
        DO i=1,MBOTP
          IF (.not.got_var(idBott(i)).and.VarOut(idBott(i),ng)) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idBott(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
      END IF QUERY
!
  10  FORMAT (1pe11.4,1x,'millimeter')
  20  FORMAT (/,' SEDIMENT_DEF_NF90 - unable to find variable: ',       &
     &        a,2x,' in output NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE sediment_def_nf90

# ifdef STATIONS
!
!***********************************************************************
      SUBROUTINE sediment_def_station_nf90 (ng, model, ldef, VarOut, S, &
     &                                      bgrd, pgrd, rgrd)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef, VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: bgrd(:), pgrd(:), rgrd(:)
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer :: i, itrc, j, status
!
      real(r8) :: Aval(6)
!
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", sediment_def_station_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define sediment output stations variables.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
!
      DEFINE : IF (ldef) THEN
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

#  if defined SEDIMENT && defined SED_MORPH
!
!  Define time-varying bathymetry.
!
        IF (VarOut(idbath,ng)) THEN
          Vinfo( 1)=Vname(1,idbath)
          Vinfo( 2)=Vname(2,idbath)
          Vinfo( 3)=Vname(3,idbath)
          Vinfo(14)=Vname(4,idbath)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idbath),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .FALSE.,                          &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#  ifdef SEDIMENT
!
!  Define sediment fraction of each size class in each bed layer.
!
        DO i=1,NST
          IF (VarOut(idfrac(i),ng)) THEN
            Vinfo( 1)=Vname(1,idfrac(i))
            Vinfo( 2)=Vname(2,idfrac(i))
            Vinfo( 3)=Vname(3,idfrac(i))
            Vinfo(14)=Vname(4,idfrac(i))
            Vinfo(16)=Vname(1,idtime)
            WRITE (Vinfo(19),10) 1000.0_r8*Sd50(i,ng)
            status=def_var(ng, model, S(ng)%ncid,                       &
     &                     S(ng)%Vid(idfrac(i)), NF_FOUT,               &
     &                     3, bgrd, Aval, Vinfo, ncname,                &
     &                     SetFillVal = .TRUE.,                         &
     &                     SetParAccess = .TRUE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!  Define sediment mass of each size class in each bed layer.
!
          IF (VarOut(idBmas(i),ng)) THEN
            Vinfo( 1)=Vname(1,idBmas(i))
            Vinfo( 2)=Vname(2,idBmas(i))
            Vinfo( 3)=Vname(3,idBmas(i))
            Vinfo(14)=Vname(4,idBmas(i))
            Vinfo(16)=Vname(1,idtime)
            WRITE (Vinfo(19),10) 1000.0_r8*Sd50(i,ng)
            status=def_var(ng, model, S(ng)%ncid,                       &
     &                     S(ng)%Vid(idBmas(i)), NF_FOUT,               &
     &                     3, bgrd, Aval, Vinfo, ncname,                &
     &                     SetFillVal = .TRUE.,                         &
     &                     SetParAccess = .TRUE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define sediment properties in each bed layer.
!
        DO i=1,MBEDP
          IF (VarOut(idSbed(i),ng)) THEN
            Vinfo( 1)=Vname(1,idSbed(i))
            Vinfo( 2)=Vname(2,idSbed(i))
            Vinfo( 3)=Vname(3,idSbed(i))
            Vinfo(14)=Vname(4,idSbed(i))
            Vinfo(16)=Vname(1,idtime)
            status=def_var(ng, model, S(ng)%ncid,                       &
     &                     S(ng)%Vid(idSbed(i)), NF_FOUT,               &
     &                     3, bgrd, Aval, Vinfo, ncname,                &
     &                     SetFillVal = .TRUE.,                         &
     &                     SetParAccess = .TRUE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#  endif

#  if defined SEDIMENT || defined BBL_MODEL
!
!  Define exposed sediment layer properties.
!
        DO i=1,MBOTP
          IF (VarOut(idBott(i),ng)) THEN
            Vinfo( 1)=Vname(1,idBott(i))
            Vinfo( 2)=Vname(2,idBott(i))
            Vinfo( 3)=Vname(3,idBott(i))
            Vinfo(14)=Vname(4,idBott(i))
            Vinfo(16)=Vname(1,idtime)
            status=def_var(ng, model, S(ng)%ncid,                       &
     &                     S(ng)%Vid(idBott(i)), NF_FOUT,               &
     &                     2, pgrd, Aval, Vinfo, ncname,                &
     &                     SetFillVal = .TRUE.,                         &
     &                     SetParAccess = .TRUE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#  endif

      END IF DEFINE
!
!-----------------------------------------------------------------------
!  Otherwise, check existing output file and prepare for appending
!  data.
!-----------------------------------------------------------------------
!
      QUERY : IF (.not.ldef) THEN
!
!  Initialize locallogical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  Waves Effect on Currents variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            S(ng)%Vid(idtime)=var_id(i)
#  if defined SEDIMENT && defined SED_MORPH
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idbath))) THEN
            got_var(idbath)=.TRUE.
            S(ng)%Vid(idbath)=var_id(i)
#  endif
          END IF
#  ifdef SEDIMENT
          DO itrc=1,NST
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idfrac(itrc)))) THEN
              got_var(idfrac(itrc))=.TRUE.
              S(ng)%Vid(idfrac(itrc))=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idBmas(itrc)))) THEN
              got_var(idBmas(itrc))=.TRUE.
              S(ng)%Vid(idBmas(itrc))=var_id(i)
            END IF
          END DO
          DO itrc=1,MBEDP
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSbed(itrc)))) THEN
              got_var(idSbed(itrc))=.TRUE.
              S(ng)%Vid(idSbed(itrc))=var_id(i)
            END IF
          END DO
#  endif
#  if defined SEDIMENT || defined BBL_MODEL
          DO itrc=1,MBOTP
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idBott(itrc)))) THEN
              got_var(idBott(itrc))=.TRUE.
              S(ng)%Vid(idBott(itrc))=var_id(i)
            END IF
          END DO
#  endif
        END DO
!
!  Check if output variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  if defined SEDIMENT && defined SED_MORPH
        IF (.not.got_var(idbath).and.VarOut(idbath,ng)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idbath)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef SEDIMENT
        DO i=1,NST
          IF (.not.got_var(idfrac(i)).and.VarOut(idfrac(i),ng)) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idfrac(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idBmas(i)).and.VarOut(idBmas(i),ng)) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idBmas(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO i=1,MBEDP
          IF (.not.got_var(idSbed(i)).and.VarOut(idSbed(i),ng)) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idSbed(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  endif
#  if defined SEDIMENT || defined BBL_MODEL
        DO i=1,MBOTP
          IF (.not.got_var(idBott(i)).and.VarOut(idBott(i),ng)) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idBott(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  endif

      END IF QUERY
!
  10  FORMAT (1pe11.4,1x,'millimeter')
  20  FORMAT (/,' SEDIMENT_DEF_STATION_NF90 - unable to find variable:',&
     &        1x,a,2x,' in output NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE sediment_def_station_nf90
# endif
!
!***********************************************************************
      SUBROUTINE sediment_wrt_nf90 (ng, model, tile,                    &
    &                               LBi, UBi, LBj, UBj,                 &
    &                               VarOut, S)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: Linstataneous
!
      integer :: gfactor, gtype, i, status
!
      real(dp) :: scale
!
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", sediment_wrt_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out Waves Effect on Currents output variables into specified
!  output NetCDF file.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
# if defined WRITE_WATER && defined MASKING
      gfactor=-1
# else
      gfactor=1
# endif
!
!  Set instantaneous fields.
!
      IF ((S(ng)%ncid.eq.HIS(ng)%ncid).or.                              &
     &    (S(ng)%ncid.eq.QCK(ng)%ncid)) THEN
        Linstataneous=.TRUE.
      ELSE
        Linstataneous=.FALSE.                  ! time-averged fiels
      END IF

# if defined SEDIMENT && defined SED_MORPH
!
!  Write out time-dependent bathymetry (m)
!
      IF (VarOut(idBath,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idbath,             &
     &                       S(ng)%Vid(idbath),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       GRID(ng) % h,                              &
     &                       SetFillVal = .FALSE.)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idbath)), S(ng)%Rindex,    &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
# endif

# if defined SEDIMENT && defined BEDLOAD
!
!  Write out bed load transport in U-direction.
!
      DO i=1,NST
        IF (VarOut(idUbld(i),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*u2dvar
          IF (Linstataneous) THEN
            status=nf_fwrite2d(ng, model, S(ng)%ncid, idUbld(i),        &
     &                         S(ng)%Vid(idUbld(i)),                    &
     &                         S(ng)%Rindex, gtype,                     &
     &                         LBi, UBi, LBj, UBj, scale,               &
#  ifdef MASKING
     &                         GRID(ng) % umask,                        &
#  endif
     &                         SEDBED(ng) % bedldu(:,:,i))
#  ifdef AVERAGES
          ELSE
            status=nf_fwrite2d(ng, model, S(ng)%ncid, idUbld(i),        &
     &                         S(ng)%Vid(idUbld(i)),                    &
     &                         S(ng)%Rindex, gtype,                     &
     &                         LBi, UBi, LBj, UBj, scale,               &
#   ifdef MASKING
     &                         GRID(ng) % umask,                        &
#   endif
     &                         SEDBED(ng) % avgbedldu(:,:,i))
#  endif
          END IF
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUbld(i))), S(ng)%Rindex, &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
!  Write out bed load transport in V-direction.
!
        IF (VarOut(idVbld(i),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*v2dvar
          IF (Linstataneous) THEN
            status=nf_fwrite2d(ng, model, S(ng)%ncid, idVbld(i),        &
     &                         S(ng)%Vid(idVbld(i)),                    &
     &                         S(ng)%Rindex, gtype,                     &
     &                         LBi, UBi, LBj, UBj, scale,               &
#  ifdef MASKING
     &                         GRID(ng) % vmask,                        &
#  endif
     &                         SEDBED(ng) % bedldv(:,:,i))
#  ifdef AVERAGES
          ELSE
            status=nf_fwrite2d(ng, model, S(ng)%ncid, idVbld(i),        &
     &                         S(ng)%Vid(idVbld(i)),                    &
     &                         S(ng)%Rindex, gtype,                     &
     &                         LBi, UBi, LBj, UBj, scale,               &
#   ifdef MASKING
     &                         GRID(ng) % vmask,                        &
#   endif
     &                         SEDBED(ng) % avgbedldv(:,:,i))
#  endif
          END IF
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVbld(i))), S(ng)%Rindex, &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif

# ifdef SEDIMENT
!
!  Write out sediment fraction of each size class in each bed layer.
!
      DO i=1,NST
        IF (VarOut(idfrac(i),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*b3dvar
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idfrac(i),          &
     &                       S(ng)%Vid(idfrac(i)),                      &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, Nbed, scale,        &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       SEDBED(ng) % bed_frac(:,:,:,i))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idfrac(i))), S(ng)%Rindex, &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out sediment mass of each size class in each bed layer.
!
      DO i=1,NST
        IF (VarOut(idBmas(i),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*b3dvar
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idBmas(i),          &
     &                       S(ng)%Vid(idBmas(i)),                      &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, Nbed, scale,        &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       SEDBED(ng) % bed_mass(:,:,:,NOUT,i))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idBmas(i))), S(ng)%Rindex, &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out sediment properties in each bed layer.
!
      DO i=1,MBEDP
        IF (VarOut(idSbed(i),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*b3dvar
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idSbed(i),          &
     &                       S(ng)%Vid(idSbed(i)),                      &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, Nbed, scale,        &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       SEDBED(ng) % bed(:,:,:,i))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idSbed(i))), S(ng)%Rindex, &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif

# if defined SEDIMENT || defined BBL_MODEL
!
!  Write out exposed sediment layer properties.
!
      DO i=1,MBOTP
        IF (VarOut(idBott(i),ng)) THEN
          IF (i.eq.itauc) THEN
            scale=rho0
          ELSE
            scale=1.0_dp
          END IF
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idBott(i),          &
     &                       S(ng)%Vid(idBott(i)),                      &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       SEDBED(ng) % bottom(:,:,i))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idBott(i))), S(ng)%Rindex, &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif

# if defined SEDIMENT && defined BEDLOAD && defined BEDLOAD_VANDERA
!
!  Write out Ursell number from skewed waveform parameters..
!
      IF (VarOut(idsurs,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, S(ng)%ncid, idsurs,               &
     &                     S(ng)%Vid(idsurs),                           &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     SEDBED(ng) % ursell_no)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idsurs)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wave velocity skewness coefficient for asymmetric bedload.
!
      IF (VarOut(idsrrw,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, S(ng)%ncid, idsrrw,               &
     &                     S(ng)%Vid(idsrrw),                           &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     SEDBED(ng) % RR_asymwave)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idsrrw)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wave acceleration skewness coefficient for asymmetric
!  bedload.
!
      IF (VarOut(idsbtw,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, S(ng)%ncid, idsbtw,               &
     &                     S(ng)%Vid(idsbtw),                           &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     SEDBED(ng) % beta_asymwave)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idsbtw)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out peak wave crest cycle velocity for asymmetric bedload.
!
      IF (VarOut(idsucr,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, S(ng)%ncid, idsucr,               &
     &                     S(ng)%Vid(idsucr),                           &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     SEDBED(ng) % ucrest_r)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idsucr)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out peak wave trough cycle velocity for asymmetric bedload.
!
      IF (VarOut(idsutr,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, S(ng)%ncid, idsutr,               &
     &                     S(ng)%Vid(idsutr),                           &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     SEDBED(ng) % utrough_r)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idsutr)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wave crest cycle period for asymmetric bedload.
!
      IF (VarOut(idstcr,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, S(ng)%ncid, idstcr,               &
     &                     S(ng)%Vid(idstcr),                           &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     SEDBED(ng) % T_crest)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idstcr)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wave trough cycle period for asymmetric bedload.
!
      IF (VarOut(idsttr,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, S(ng)%ncid, idsttr,               &
     &                     S(ng)%Vid(idsttr),                           &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     SEDBED(ng) % T_trough)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idsttr)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
  10  FORMAT (/," SEDIMENT_WRT_NF90 - error while writing variable '",  &
     &        a,"', time record = ",i0,/,11x,'into file: ',a)
!
      RETURN
      END SUBROUTINE sediment_wrt_nf90

# ifdef STATIONS
!
!***********************************************************************
      SUBROUTINE sediment_wrt_station_nf90 (ng, model, tile,            &
     &                                      LBi, UBi, LBj, UBj,         &
     &                                      VarOut, S)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: Cgrid
!
      integer :: NposB
      integer :: i, k, np, status
!
      real(dp) :: scale
!
      real(r8), dimension(Nstation(ng)) :: Xpos, Ypos, Zpos, psta
#  ifdef SEDIMENT
      real(r8), dimension(Nstation(ng)*Nbed) :: XposB, YposB, ZposB
      real(r8), dimension(Nstation(ng)*Nbed) :: bsta
#  endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", sediment_wrt_station_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out sediment output variables into specified stations NetCDF
!  file.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set switch to extract station data at native C-grid position (TRUE)
!  or at RHO-points (FALSE).
!
#  ifdef STATIONS_CGRID
      Cgrid=.TRUE.
#  else
      Cgrid=.FALSE.
#  endif
!
!  Set positions for generic extraction routine.
!
      NposB=Nstation(ng)*Nbed
      DO i=1,Nstation(ng)
        Xpos(i)=SCALARS(ng)%SposX(i)
        Ypos(i)=SCALARS(ng)%SposY(i)
        Zpos(i)=1.0_r8
#  ifdef SEDIMENT
        DO k=1,Nbed
          np=k+(i-1)*Nbed
          XposB(np)=SCALARS(ng)%SposX(i)
          YposB(np)=SCALARS(ng)%SposY(i)
          ZposB(np)=REAL(k,r8)
        END DO
#  endif
      END DO

#  if defined SEDIMENT && defined SED_MORPH
!
!  Write out time-varying bathymetry.
!
      IF (VarOut(idbath,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idbath, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%h,                          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idbath)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idbath))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif

#  ifdef SEDIMENT
!
!  Write out sediment fraction of each size class in each bed layer.
!
      DO i=1,NST
        IF (VarOut(idfrac(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, model, Cgrid, idfrac(i), b3dvar,      &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        scale, SEDBED(ng)%bed_frac(:,:,:,i),      &
     &                        NposB, XposB, YposB, ZposB, bsta)
          CALL netcdf_put_fvar (ng, model, S(ng)%name,                  &
     &                          TRIM(Vname(1,idfrac(i))), bsta,         &
     &                          (/1,1,S(ng)%Rindex/),                   &
     &                          (/Nbed,Nstation(ng),1/),                &
     &                          ncid = S(ng)%ncid,                      &
     &                          varid = S(ng)%Vid(idfrac(i)))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Write out sediment mass of each size class in each bed layer.
!
        IF (VarOut(idBmas(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, model, Cgrid, idBmas(i), b3dvar,      &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        scale,                                    &
     &                        SEDBED(ng)%bed_mass(:,:,:,NOUT,i),        &
     &                        NposB, XposB, YposB, ZposB, bsta)
          CALL netcdf_put_fvar (ng, model, S(ng)%name,                  &
     &                          TRIM(Vname(1,idBmas(i))), bsta,         &
     &                          (/1,1,S(ng)%Rindex/),                   &
     &                          (/Nbed,Nstation(ng),1/),                &
     &                          ncid = S(ng)%ncid,                      &
     &                          varid = S(ng)%Vid(idBmas(i)))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
!
!  Write out sediment properties in each bed layer.
!
      DO i=1,MBEDP
        IF (VarOut(idSbed(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, model, Cgrid, idSbed(i), b3dvar,      &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        scale, SEDBED(ng)%bed(:,:,:,i),           &
     &                        NposB, XposB, YposB, ZposB, bsta)
          CALL netcdf_put_fvar (ng, model, S(ng)%name,                  &
     &                          TRIM(Vname(1,idSbed(i))), bsta,         &
     &                          (/1,1,S(ng)%Rindex/),                   &
     &                          (/Nbed,Nstation(ng),1/),                &
     &                          ncid = S(ng)%ncid,                      &
     &                          varid = S(ng)%Vid(idSbed(i)))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
#  endif

#  if defined SEDIMENT || defined BBL_MODEL
!
!  Write out exposed sediment layer properties.
!
      DO i=1,MBEDP
        IF (VarOut(idBott(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta2d (ng, model, Cgrid, idBott(i), r2dvar,      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        scale, SEDBED(ng)%bottom(:,:,i),          &
     &                        Nstation(ng), Xpos, Ypos, psta)
          CALL netcdf_put_fvar (ng, model, S(ng)%name,                  &
     &                          TRIM(Vname(1,idBott(i))), psta,         &
     &                          (/1,S(ng)%Rindex/),                     &
     &                          (/Nstation(ng),1/),                     &
     &                          ncid = S(ng)%ncid,                      &
     &                          varid = S(ng)%Vid(idBott(i)))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
#  endif
!
      RETURN
      END SUBROUTINE sediment_wrt_station_nf90
# endif

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE sediment_def_pio (ng, model, ldef, VarOut, S,          &
     &                             t2dgrd, u2dgrd, v2dgrd,              &
     &                             b3dgrd)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef, VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: t2dgrd(:), u2dgrd(:), v2dgrd(:)
      integer, intent(in), optional :: b3dgrd(:)
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer :: i, itrc, j, nvd3, nvd4, status
!
      real(r8) :: Aval(6)
!
#  ifdef ADJOINT
      character (len=21)  :: Prefix
#  else
      character (len=13)  :: Prefix
#  endif
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", sediment_def_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define Waves Effect on Currents output variables.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
!
      DEFINE : IF (ldef) THEN
!
!  Set number of dimensions for output variables.
!
#  if defined WRITE_WATER && defined MASKING
        nvd3=2
        nvd4=2
#  else
        nvd3=3
        nvd4=4
#  endif
!
!  Set long name prefix string.
!
#  ifdef ADJOINT
!!      Prefix='time-averaged adjoint'
        Prefix='adjoint'
#  else
!!      Prefix='time-averaged'
        Prefix=CHAR(32)                                   ! blank space
#  endif
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

#  if defined SEDIMENT && defined SED_MORPH
!
!  Define time-varying bathymetry.
!
        IF (VarOut(idbath,ng)) THEN
          Vinfo( 1)=Vname(1,idbath)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idbath))
          ELSE
            Vinfo( 2)=Vname(2,idbath)
          END IF
          Vinfo( 3)=Vname(3,idbath)
          Vinfo(14)=Vname(4,idbath)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idbath)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idbath,ng),r8)
          S(ng)%pioVar(idbath)%dkind=PIO_FOUT
          S(ng)%pioVar(idbath)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idbath)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#  if defined SEDIMENT && defined BEDLOAD
!
!  Define Bedload transbort U-direction.
!
        DO i=1,NST
          IF (VarOut(idUbld(i),ng)) THEN
            Vinfo( 1)=Vname(1,idUbld(i))
            IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                      &
     &                                     TRIM(Vname(2,idUbld(i)))
            ELSE
              Vinfo( 2)=Vname(2,idUbld(i))
            END IF
            Vinfo( 3)=Vname(3,idUbld(i))
            Vinfo(14)=Vname(4,idUbld(i))
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_u'
#   endif
            Vinfo(21)=Vname(6,idUbld(i))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idUbld(i),ng),r8)
            S(ng)%pioVar(idUbld(i))%dkind=PIO_FOUT
            S(ng)%pioVar(idUbld(i))%gtype=u2dvar
!
            status=def_var(ng, model, S(ng)%pioFile,                    &
     &                     S(ng)%pioVar(idUbld(i))%vd,                  &
     &                     PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!  Define Bedload transport V-direction.
!
          IF (VarOut(idVbld(i),ng)) THEN
            Vinfo( 1)=Vname(1,idVbld(i))
            IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                      &
     &                                     TRIM(Vname(2,idVbld(i)))
            ELSE
              Vinfo( 2)=Vname(2,idVbld(i))
            END IF
            Vinfo( 3)=Vname(3,idVbld(i))
            Vinfo(14)=Vname(4,idVbld(i))
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_v'
#   endif
            Vinfo(21)=Vname(6,idVbld(i))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idVbld(i),ng),r8)
            S(ng)%pioVar(idVbld(i))%dkind=PIO_FOUT
            S(ng)%pioVar(idVbld(i))%gtype=v2dvar
!
            status=def_var(ng, model, S(ng)%pioFile,                    &
     &                     S(ng)%pioVar(idVbld(i))%vd,                  &
     &                     PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#  endif

#  ifdef SEDIMENT
!
!  Define sediment fraction of each size class in each bed layer.
!
        DO i=1,NST
          IF (VarOut(idfrac(i),ng)) THEN
            Vinfo( 1)=Vname(1,idfrac(i))
            IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                      &
     &                                     TRIM(Vname(2,idfrac(i)))
            ELSE
              Vinfo( 2)=Vname(2,idfrac(i))
            END IF
            Vinfo( 3)=Vname(3,idfrac(i))
            Vinfo(14)=Vname(4,idfrac(i))
            Vinfo(16)=Vname(1,idtime)
            WRITE (Vinfo(19),10) 1000.0_r8*Sd50(i,ng)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#   endif
            Vinfo(21)=Vname(6,idfrac(i))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idfrac(i),ng))
            S(ng)%pioVar(idfrac(i))%dkind=PIO_FOUT
            S(ng)%pioVar(idfrac(i))%gtype=b3dvar
!
            status=def_var(ng, model, S(ng)%pioFile,                    &
     &                     S(ng)%pioVar(idfrac(i))%vd,                  &
     &                     PIO_FOUT, nvd4, b3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define sediment mass of each size class in each bed layer.
!
        DO i=1,NST
          IF (VarOut(idBmas(i),ng)) THEN
            Vinfo( 1)=Vname(1,idBmas(i))
            IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                      &
     &                                     TRIM(Vname(2,idBmas(i)))
            ELSE
              Vinfo( 2)=Vname(2,idBmas(i))
            END IF
            Vinfo( 3)=Vname(3,idBmas(i))
            Vinfo(14)=Vname(4,idBmas(i))
            Vinfo(16)=Vname(1,idtime)
            WRITE (Vinfo(19),10) 1000.0_r8*Sd50(i,ng)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#   endif
            Vinfo(21)=Vname(6,idBmas(i))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idBmas(i),ng),r8)
            S(ng)%pioVar(idBmas(i))%dkind=PIO_FOUT
            S(ng)%pioVar(idBmas(i))%gtype=b3dvar
!
            status=def_var(ng, model, S(ng)%pioFile,                    &
     &                     S(ng)%pioVar(idBmas(i))%vd,                  &
     &                     PIO_FOUT, nvd4, b3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define sediment properties in each bed layer.
!
        DO i=1,MBEDP
          IF (VarOut(idSbed(i),ng)) THEN
            Vinfo( 1)=Vname(1,idSbed(i))
            IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                      &
     &                                     TRIM(Vname(2,idSbed(i)))
            ELSE
              Vinfo( 2)=Vname(2,idSbed(i))
            END IF
            Vinfo( 3)=Vname(3,idSbed(i))
            Vinfo(14)=Vname(4,idSbed(i))
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#   endif
            Vinfo(21)=Vname(6,idSbed(i))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idSbed(i),ng),r8)
            S(ng)%pioVar(idSbed(i))%dkind=PIO_FOUT
            S(ng)%pioVar(idSbed(i))%gtype=b3dvar
!
            status=def_var(ng, model, S(ng)%pioFile,                    &
     &                     S(ng)%pioVar(idSbed(i))%vd,                  &
                           PIO_FOUT, nvd4, b3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#  endif

#  if defined SEDIMENT || defined BBL_MODEL
!
!  Define exposed sediment layer properties.
!
        DO i=1,MBOTP
          IF (VarOut(idBott(i),ng)) THEN
            Vinfo( 1)=Vname(1,idBott(i))
            IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                      &
     &                                     TRIM(Vname(2,idBott(i)))
            ELSE
              Vinfo( 2)=Vname(2,idBott(i))
            END IF
            Vinfo( 3)=Vname(3,idBott(i))
            Vinfo(14)=Vname(4,idBott(i))
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#   endif
            Vinfo(21)=Vname(6,idBott(i))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idBott(i),ng),r8)
            S(ng)%pioVar(idBott(i))%dkind=PIO_FOUT
            S(ng)%pioVar(idBott(i))%gtype=r2dvar
!
            status=def_var(ng, model, S(ng)%pioFile,                    &
     &                     S(ng)%pioVar(idBott(i))%vd,                  &
     &                     PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#  endif

#  if defined SEDIMENT && defined BEDLOAD && defined BEDLOAD_VANDERA
!
!  Define Ursell number from skewed waveform parameters.
!
        IF (VarOut(idsurs,ng)) THEN
          Vinfo( 1)=Vname(1,idsurs)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idsurs))
          ELSE
            Vinfo( 2)=Vname(2,idsurs)
          END IF
          Vinfo( 3)=Vname(3,idsurs)
          Vinfo(14)=Vname(4,idsurs)
          Vinfo(16)=Vname(1,idsurs)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idsurs)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idsurs,ng),r8)
          S(ng)%pioVar(idsurs)%dkind=PIO_FOUT
          S(ng)%pioVar(idsurs)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idsurs)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wave velocity skewness coefficient for asymmetric bedload.
!
        IF (VarOut(idsrrw,ng)) THEN
          Vinfo( 1)=Vname(1,idsrrw)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idsrrw))
          ELSE
            Vinfo( 2)=Vname(2,idsrrw)
          END IF
          Vinfo( 3)=Vname(3,idsrrw)
          Vinfo(14)=Vname(4,idsrrw)
          Vinfo(16)=Vname(1,idsrrw)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idsrrw)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idsrrw,ng),r8)
          S(ng)%pioVar(idsrrw)%dkind=PIO_FOUT
          S(ng)%pioVar(idsrrw)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idsrrw)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wave acceleration skewness coefficient for asymmetric bedload.
!
        IF (VarOut(idsbtw,ng)) THEN
          Vinfo( 1)=Vname(1,idsbtw)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idsbtw))
          ELSE
            Vinfo( 2)=Vname(2,idsbtw)
          END IF
          Vinfo( 3)=Vname(3,idsbtw)
          Vinfo(14)=Vname(4,idsbtw)
          Vinfo(16)=Vname(1,idsbtw)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idsbtw)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idsbtw,ng),r8)
          S(ng)%pioVar(idsbtw)%dkind=PIO_FOUT
          S(ng)%pioVar(idsbtw)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idsbtw)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define peak wave crest cycle velocity for asymmetric bedload.
!
        IF (VarOut(idsucr,ng)) THEN
          Vinfo( 1)=Vname(1,idsucr)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idsucr))
          ELSE
            Vinfo( 2)=Vname(2,idsucr)
          END IF
          Vinfo( 3)=Vname(3,idsucr)
          Vinfo(14)=Vname(4,idsucr)
          Vinfo(16)=Vname(1,idsucr)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idsucr)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idsucr,ng),r8)
          S(ng)%pioVar(idsucr)%dkind=PIO_FOUT
          S(ng)%pioVar(idsucr)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idsucr)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define peak wave trough cycle velocity for asymmetric bedload.
!
        IF (VarOut(idsutr,ng)) THEN
          Vinfo( 1)=Vname(1,idsutr)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idsutr))
          ELSE
            Vinfo( 2)=Vname(2,idsutr)
          END IF
          Vinfo( 3)=Vname(3,idsutr)
          Vinfo(14)=Vname(4,idsutr)
          Vinfo(16)=Vname(1,idsutr)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idsutr)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idsutr,ng),r8)
          S(ng)%pioVar(idsutr)%dkind=PIO_FOUT
          S(ng)%pioVar(idsutr)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idsutr)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wave crest cycle period for asymmetric bedload.
!
        IF (VarOut(idstcr,ng)) THEN
          Vinfo( 1)=Vname(1,idstcr)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idstcr))
          ELSE
            Vinfo( 2)=Vname(2,idstcr)
          END IF
          Vinfo( 3)=Vname(3,idstcr)
          Vinfo(14)=Vname(4,idstcr)
          Vinfo(16)=Vname(1,idstcr)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idstcr)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idstcr,ng),r8)
          S(ng)%pioVar(idstcr)%dkind=PIO_FOUT
          S(ng)%pioVar(idstcr)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idstcr)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wave trough cycle period for asymmetric bedload.
!
        IF (VarOut(idsttr,ng)) THEN
          Vinfo( 1)=Vname(1,idsttr)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idsttr))
          ELSE
            Vinfo( 2)=Vname(2,idsttr)
          END IF
          Vinfo( 3)=Vname(3,idsttr)
          Vinfo(14)=Vname(4,idsttr)
          Vinfo(16)=Vname(1,idsttr)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idsttr)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idsttr,ng),r8)
          S(ng)%pioVar(idsttr)%dkind=PIO_FOUT
          S(ng)%pioVar(idsttr)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idsttr)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

      END IF DEFINE
!
!-----------------------------------------------------------------------
!  Otherwise, check existing output file and prepare for appending
!  data.
!-----------------------------------------------------------------------
!
      QUERY : IF (.not.ldef) THEN
!
!  Initialize locallogical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  Waves Effect on Currents variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            S(ng)%pioVar(idtime)%vd=var_desc(i)
            S(ng)%pioVar(idtime)%dkind=PIO_TOUT
            S(ng)%pioVar(idtime)%gtype=0
#  if defined SEDIMENT && defined SED_MORPH
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idbath))) THEN
            got_var(idbath)=.TRUE.
            S(ng)%pioVar(idbath)%vd=var_desc(i)
            S(ng)%pioVar(idbath)%dkind=PIO_FOUT
            S(ng)%pioVar(idbath)%gtype=r2dvar
#  endif
#  if defined SEDIMENT && defined BEDLOAD && defined BEDLOAD_VANDERA
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idsurs))) THEN
            got_var(idsurs)=.TRUE.
            S(ng)%pioVar(idsurs)%vd=var_desc(i)
            S(ng)%pioVar(idsurs)%dkind=PIO_FOUT
            S(ng)%pioVar(idsurs)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idsrrw))) THEN
            got_var(idsrrw)=.TRUE.
            S(ng)%pioVar(idsrrw)%vd=var_desc(i)
            S(ng)%pioVar(idsrrw)%dkind=PIO_FOUT
            S(ng)%pioVar(idsrrw)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idsbtw))) THEN
            got_var(idsbtw)=.TRUE.
            S(ng)%pioVar(idsbtw)%vd=var_desc(i)
            S(ng)%pioVar(idsbtw)%dkind=PIO_FOUT
            S(ng)%pioVar(idsbtw)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idsucr))) THEN
            got_var(idsucr)=.TRUE.
            S(ng)%pioVar(idsucr)%vd=var_desc(i)
            S(ng)%pioVar(idsucr)%dkind=PIO_FOUT
            S(ng)%pioVar(idsucr)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idsutr))) THEN
            got_var(idsutr)=.TRUE.
            S(ng)%pioVar(idsutr)%vd=var_desc(i)
            S(ng)%pioVar(idsutr)%dkind=PIO_FOUT
            S(ng)%pioVar(idsutr)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idstcr))) THEN
            got_var(idstcr)=.TRUE.
            S(ng)%pioVar(idstcr)%vd=var_desc(i)
            S(ng)%pioVar(idstcr)%dkind=PIO_FOUT
            S(ng)%pioVar(idstcr)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idsttr))) THEN
            got_var(idsttr)=.TRUE.
            S(ng)%pioVar(idsttr)%vd=var_desc(i)
            S(ng)%pioVar(idsttr)%dkind=PIO_FOUT
            S(ng)%pioVar(idsttr)%gtype=r2dvar
#  endif
          END IF
#  ifdef SEDIMENT
          DO itrc=1,NST
            IF (TRIM(var_name(i)).eq.                                   &
     &               TRIM(Vname(1,idfrac(itrc)))) THEN
              got_var(idfrac(itrc))=.TRUE.
              S(ng)%pioVar(idfrac(itrc))%vd=var_desc(i)
              S(ng)%pioVar(idfrac(itrc))%dkind=PIO_FOUT
              S(ng)%pioVar(idfrac(itrc))%gtype=b3dvar
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idBmas(itrc)))) THEN
              got_var(idBmas(itrc))=.TRUE.
              S(ng)%pioVar(idBmas(itrc))%vd=var_desc(i)
              S(ng)%pioVar(idBmas(itrc))%dkind=PIO_FOUT
              S(ng)%pioVar(idBmas(itrc))%gtype=b3dvar
#   ifdef BEDLOAD
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idUbld(itrc)))) THEN
              got_var(idUbld(itrc))=.true.
              S(ng)%pioVar(idUbld(itrc))%vd=var_desc(i)
              S(ng)%pioVar(idUbld(itrc))%dkind=PIO_FOUT
              S(ng)%pioVar(idUbld(itrc))%gtype=u2dvar
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idVbld(itrc)))) THEN
              got_var(idVbld(itrc))=.true.
              S(ng)%pioVar(idVbld(itrc))%vd=var_desc(i)
              S(ng)%pioVar(idVbld(itrc))%dkind=PIO_FOUT
              S(ng)%pioVar(idVbld(itrc))%gtype=v2dvar
#   endif
            END IF
          END DO
          DO itrc=1,MBEDP
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSbed(itrc)))) THEN
              got_var(idSbed(itrc))=.TRUE.
              S(ng)%pioVar(idSbed(itrc))%vd=var_desc(i)
              S(ng)%pioVar(idSbed(itrc))%dkind=PIO_FOUT
              S(ng)%pioVar(idSbed(itrc))%gtype=b3dvar
            END IF
          END DO
#  endif
#  if defined SEDIMENT || defined BBL_MODEL
          DO itrc=1,MBOTP
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idBott(itrc)))) THEN
              got_var(idBott(itrc))=.TRUE.
              S(ng)%pioVar(idBott(itrc))%vd=var_desc(i)
              S(ng)%pioVar(idBott(itrc))%dkind=PIO_FOUT
              S(ng)%pioVar(idBott(itrc))%gtype=r2dvar
            END IF
          END DO
#  endif
        END DO
!
!  Check if output variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  if defined SEDIMENT && defined SED_MORPH
        IF (.not.got_var(idbath).and.VarOut(idbath,ng)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idbath)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined SEDIMENT && defined BEDLOAD && defined BEDLOAD_VANDERA
        IF (.not.got_var(idsurs).and.VarOut(idsurs,ng)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idsurs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idsrrw).and.VarOut(idsrrw,ng)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idsrrw)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idsbtw).and.VarOut(idsbtw,ng)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idsbtw)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idsucr).and.VarOut(idsucr,ng)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idsucr)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idsutr).and.VarOut(idsutr,ng)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idsutr)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idstcr).and.VarOut(idstcr,ng)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idstcr)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idsttr).and.VarOut(idsttr,ng)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idsttr)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef SEDIMENT
        DO i=1,NST
          IF (.not.got_var(idfrac(i)).and.VarOut(idfrac(i),ng)) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idfrac(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF

          IF(.not.got_var(idBmas(i)).and.VarOut(idBmas(i),ng)) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idBmas(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#   ifdef BEDLOAD
          IF (.not.got_var(idUbld(i)).and.VarOut(idUbld(i),ng)) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idUbld(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idVbld(i)).and.VarOut(idVbld(i),ng)) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idVbld(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#   endif
        END DO
        DO i=1,MBEDP
          IF (.not.got_var(idSbed(i)).and.VarOut(idSbed(i),ng)) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idSbed(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  endif
#  if defined SEDIMENT || defined BBL_MODEL
        DO i=1,MBOTP
          IF (.not.got_var(idBott(i)).and.VarOut(idBott(i),ng)) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idBott(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  endif
      END IF QUERY
!
  10  FORMAT (1pe11.4,1x,'millimeter')
  20  FORMAT (/,' SEDIMENT_DEF_PIO  - unable to find variable: ',       &
     &        a,2x,' in output NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE sediment_def_pio
!
!***********************************************************************
      SUBROUTINE sediment_wrt_pio (ng, model, tile,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             VarOut, S)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: Linstataneous
!
      integer :: i, status
!
      real(dp) :: scale
!
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", sediment_wrt_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out sediment output variables into specified NetCDF file.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
!
!  Set instantaneous fields.
!
      IF ((S(ng)%ncid.eq.S(ng)%ncid).or.                              &
     &    (S(ng)%ncid.eq.QCK(ng)%ncid)) THEN
        Linstataneous=.TRUE.
      ELSE
        Linstataneous=.FALSE.                  ! time-averged fiels
      END IF

#  if defined SEDIMENT && defined SED_MORPH
!
!  Write out time-dependent bathymetry (m)
!
      IF (VarOut(idBath,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idbath)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, idbath,            &
     &                     S(ng)%pioVar(idbath),                        &
     &                     S(ng)%Rindex,                                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     GRID(ng) % h,                                &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idbath)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif

#  if defined SEDIMENT && defined BEDLOAD
!
!  Write out bed load transport in U-direction.
!
      DO i=1,NST
        IF (VarOut(idUbld(i),ng)) THEN
          scale=1.0_dp
          IF (S(ng)%pioVar(idUbld(i))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_u2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_u2dvar(ng)
          END IF
          IF (Linstataneous) THEN
            status=nf_fwrite2d(ng, model, S(ng)%pioFile, idUbld(i),     &
     &                         S(ng)%pioVar(idUbld(i)),                 &
     &                         S(ng)%Rindex,                            &
     &                         ioDesc,                                  &
     &                         LBi, UBi, LBj, UBj, scale,               &
#   ifdef MASKING
     &                         GRID(ng) % umask,                        &
#   endif
     &                         SEDBED(ng) % bedldu(:,:,i))
#   ifdef AVERAGES
          ELSE
            status=nf_fwrite2d(ng, model, S(ng)%pioFile, idUbld(i),     &
     &                         S(ng)%pioVar(idUbld(i)),                 &
     &                         S(ng)%Rindex,                            &
     &                         ioDesc,                                  &
     &                         LBi, UBi, LBj, UBj, scale,               &
#    ifdef MASKING
     &                         GRID(ng) % umask,                        &
#    endif
     &                         SEDBED(ng) % avgbedldu(:,:,i))
#   endif
          END IF
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUbld(i))), S(ng)%Rindex, &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
!  Write out bed load transport in V-direction.
!
        IF (VarOut(idVbld(i),ng)) THEN
          scale=1.0_dp
          IF (S(ng)%pioVar(idVbld(i))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_v2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_v2dvar(ng)
          END IF
          IF (Linstataneous) THEN
            status=nf_fwrite2d(ng, model, S(ng)%pioFile, idVbld(i),     &
     &                         S(ng)%pioVar(idVbld(i)),                 &
     &                         S(ng)%Rindex,                            &
     &                         ioDesc,                                  &
     &                         LBi, UBi, LBj, UBj, scale,               &
#   ifdef MASKING
     &                         GRID(ng) % vmask,                        &
#   endif
     &                         SEDBED(ng) % bedldv(:,:,i))
#   ifdef AVERAGES
          ELSE
            status=nf_fwrite2d(ng, model, S(ng)%pioFile, idVbld(i),     &
     &                         S(ng)%pioVar(idVbld(i)),                 &
     &                         S(ng)%Rindex,                            &
     &                         ioDesc,                                  &
     &                         LBi, UBi, LBj, UBj, scale,               &
#    ifdef MASKING
     &                         GRID(ng) % vmask,                        &
#    endif
     &                         SEDBED(ng) % avgbedldv(:,:,i))
#   endif
          END IF
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVbld(i))), S(ng)%Rindex, &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif

#  ifdef SEDIMENT
!
!  Write out sediment fraction of each size class in each bed layer.
!
      DO i=1,NST
        IF (VarOut(idfrac(i),ng)) THEN
          scale=1.0_dp
          IF (S(ng)%pioVar(idfrac(i))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_b3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_b3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, S(ng)%pioFile, idfrac(i),       &
     &                       S(ng)%pioVar(idfrac(i)),                   &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, Nbed, scale,        &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       SEDBED(ng) % bed_frac(:,:,:,i))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idfrac(i))), S(ng)%Rindex, &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out sediment mass of each size class in each bed layer.
!
      DO i=1,NST
        IF (VarOut(idBmas(i),ng)) THEN
          scale=1.0_dp
          IF (S(ng)%pioVar(idBmas(i))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_b3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_b3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, S(ng)%pioFile, idBmas(i),       &
     &                       S(ng)%pioVar(idBmas(i)),                   &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, Nbed, scale,        &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       SEDBED(ng) % bed_mass(:,:,:,NOUT,i))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idBmas(i))), S(ng)%Rindex, &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out sediment properties in each bed layer.
!
      DO i=1,MBEDP
        IF (VarOut(idSbed(i),ng)) THEN
          scale=1.0_dp
          IF (S(ng)%pioVar(idSbed(i))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_b3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_b3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, S(ng)%pioFile, idSbed(i),       &
     &                       S(ng)%pioVar(idSbed(i)),                   &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, Nbed, scale,        &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       SEDBED(ng) % bed(:,:,:,i))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idSbed(i))), S(ng)%Rindex, &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif

#  if defined SEDIMENT || defined BBL_MODEL
!
!  Write out exposed sediment layer properties.
!
      DO i=1,MBOTP
        IF (VarOut(idBott(i),ng)) THEN
          IF (i.eq.itauc) THEN
            scale=rho0
          ELSE
            scale=1.0_dp
          END IF
          IF (S(ng)%pioVar(idBott(i))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dvar(ng)
          END IF
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idBott(i),       &
     &                       S(ng)%pioVar(idBott(i)),                   &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       SEDBED(ng) % bottom(:,:,i))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idBott(i))), S(ng)%Rindex, &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif

#  if defined SEDIMENT && defined BEDLOAD && defined BEDLOAD_VANDERA
!
!  Write out Ursell number from skewed waveform parameters.
!
      IF (VarOut(idsurs,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idsurs)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, idsurs,            &
     &                     S(ng)%pioVar(idsurs),                        &
     &                     S(ng)%Rindex,                                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     SEDBED(ng) % ursell_no)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idsurs)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wave velocity skewness coefficient for asymmetric bedload.
!
      IF (VarOut(idsrrw,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idsrrw)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, idsrrw,            &
     &                     S(ng)%pioVar(idsrrw),                        &
     &                     S(ng)%Rindex,                                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     SEDBED(ng) % RR_asymwave)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idsrrw)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wave acceleration skewness coefficient for asymmetric
!  bedload.
!
      IF (VarOut(idsbtw,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idsbtw)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, idsbtw,            &
     &                     S(ng)%pioVar(idsbtw),                        &
     &                     S(ng)%Rindex,                                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     SEDBED(ng) % beta_asymwave)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idsbtw)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out peak wave crest cycle velocity for asymmetric bedload.
!
      IF (VarOut(idsucr,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idsucr)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, idsucr,            &
     &                     S(ng)%pioVar(idsucr),                        &
     &                     S(ng)%Rindex,                                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     SEDBED(ng) % ucrest_r)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idsucr)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out peak wave trough cycle velocity for asymmetric bedload.
!
      IF (VarOut(idsutr,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idsutr)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, idsutr,            &
     &                     S(ng)%pioVar(idsutr),                        &
     &                     S(ng)%Rindex,                                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     SEDBED(ng) % utrough_r)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idsutr)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wave crest cycle period for asymmetric bedload.
!
      IF (VarOut(idstcr,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idstcr)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, idstcr,            &
     &                     S(ng)%pioVar(idstcr),                        &
     &                     S(ng)%Rindex,                                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     SEDBED(ng) % T_crest)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idstcr)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wave trough cycle period for asymmetric bedload.
!
      IF (VarOut(idsttr,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idsttr)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, idsttr,            &
     &                     S(ng)%pioVar(idsttr),                        &
     &                     S(ng)%Rindex,                                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     SEDBED(ng) % T_trough)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idsttr)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
  10  FORMAT (/," SEDIMENT_WRT_PIO  - error while writing variable '",  &
     &        a,"', time record = ",i0,/,11x,'into file: ',a)
!
      RETURN
      END SUBROUTINE sediment_wrt_pio

#  ifdef STATIONS
!
!***********************************************************************
      SUBROUTINE sediment_def_station_pio (ng, model, ldef, VarOut, S,  &
     &                                     bgrd, pgrd, rgrd)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef, VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: bgrd(:), pgrd(:), rgrd(:)
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer :: i, itrc, j, status
!
      real(r8) :: Aval(6)
!
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", sediment_def_station_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define sediment output stations variables.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
!
      DEFINE : IF (ldef) THEN
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

#   if defined SEDIMENT && defined SED_MORPH
!
!  Define time-varying bathymetry.
!
        IF (VarOut(idbath,ng)) THEN
          Vinfo( 1)=Vname(1,idbath)
          Vinfo( 2)=Vname(2,idbath)
          Vinfo( 3)=Vname(3,idbath)
          Vinfo(14)=Vname(4,idbath)
          Vinfo(16)=Vname(1,idtime)
          S(ng)%pioVar(idbath)%dkind=PIO_FOUT
          S(ng)%pioVar(idbath)%gtype=0
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idbath)%vd,                       &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .FALSE.,                          &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif

#   ifdef SEDIMENT
!
!  Define sediment fraction of each size class in each bed layer.
!
        DO i=1,NST
          IF (VarOut(idfrac(i),ng)) THEN
            Vinfo( 1)=Vname(1,idfrac(i))
            Vinfo( 2)=Vname(2,idfrac(i))
            Vinfo( 3)=Vname(3,idfrac(i))
            Vinfo(14)=Vname(4,idfrac(i))
            Vinfo(16)=Vname(1,idtime)
            WRITE (Vinfo(19),10) 1000.0_r8*Sd50(i,ng)
            S(ng)%pioVar(idfrac(i))%dkind=PIO_FOUT
            S(ng)%pioVar(idfrac(i))%gtype=0
!
            status=def_var(ng, model, S(ng)%pioFile,                    &
     &                     S(ng)%pioVar(idfrac(i))%vd,                  &
     &                     PIO_FOUT, 3, bgrd, Aval, Vinfo, ncname,      &
     &                     SetFillVal = .TRUE.,                         &
     &                     SetParAccess = .TRUE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!  Define sediment mass of each size class in each bed layer.
!
          IF (VarOut(idBmas(i),ng)) THEN
            Vinfo( 1)=Vname(1,idBmas(i))
            Vinfo( 2)=Vname(2,idBmas(i))
            Vinfo( 3)=Vname(3,idBmas(i))
            Vinfo(14)=Vname(4,idBmas(i))
            Vinfo(16)=Vname(1,idtime)
            WRITE (Vinfo(19),10) 1000.0_r8*Sd50(i,ng)
            S(ng)%pioVar(idBmas(i))%dkind=PIO_FOUT
            S(ng)%pioVar(idBmas(i))%gtype=0
!
            status=def_var(ng, model, S(ng)%pioFile,                    &
     &                     S(ng)%pioVar(idBmas(i))%vd,                  &
     &                     PIO_FOUT, 3, bgrd, Aval, Vinfo, ncname,      &
     &                     SetFillVal = .TRUE.,                         &
     &                     SetParAccess = .TRUE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define sediment properties in each bed layer.
!
        DO i=1,MBEDP
          IF (VarOut(idSbed(i),ng)) THEN
            Vinfo( 1)=Vname(1,idSbed(i))
            Vinfo( 2)=Vname(2,idSbed(i))
            Vinfo( 3)=Vname(3,idSbed(i))
            Vinfo(14)=Vname(4,idSbed(i))
            Vinfo(16)=Vname(1,idtime)
            S(ng)%pioVar(idSbed(i))%dkind=PIO_FOUT
            S(ng)%pioVar(idSbed(i))%gtype=0
!
            status=def_var(ng, model, S(ng)%pioFile,                    &
     &                     S(ng)%pioVar(idSbed(i))%vd,                  &
     &                     PIO_FOUT, 3, bgrd, Aval, Vinfo, ncname,      &
     &                     SetFillVal = .TRUE.,                         &
     &                     SetParAccess = .TRUE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#   endif

#   if defined SEDIMENT || defined BBL_MODEL
!
!  Define exposed sediment layer properties.
!
        DO i=1,MBOTP
          IF (VarOut(idBott(i),ng)) THEN
            Vinfo( 1)=Vname(1,idBott(i))
            Vinfo( 2)=Vname(2,idBott(i))
            Vinfo( 3)=Vname(3,idBott(i))
            Vinfo(14)=Vname(4,idBott(i))
            Vinfo(16)=Vname(1,idtime)
            S(ng)%pioVar(idBott(i))%dkind=PIO_FOUT
            S(ng)%pioVar(idBott(i))%gtype=0
!
            status=def_var(ng, model, S(ng)%pioFile,                    &
     &                     S(ng)%pioVar(idBott(i))%vd,                  &
     &                     PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,      &
     &                     SetFillVal = .TRUE.,                         &
     &                     SetParAccess = .TRUE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#   endif

      END IF DEFINE
!
!-----------------------------------------------------------------------
!  Open an existing stations file, check its contents, and prepare for
!  appending data.
!-----------------------------------------------------------------------
!
      QUERY : IF (.not.ldef) THEN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  stations variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            S(ng)%pioVar(idtime)%vd=var_desc(i)
            S(ng)%pioVar(idtime)%dkind=PIO_TOUT
            S(ng)%pioVar(idtime)%gtype=0
#   if defined SEDIMENT && defined SED_MORPH
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idbath))) THEN
            got_var(idbath)=.TRUE.
            S(ng)%pioVar(idbath)%vd=var_desc(i)
            S(ng)%pioVar(idbath)%dkind=PIO_FOUT
            S(ng)%pioVar(idbath)%gtype=0
#   endif
          END IF
#   ifdef SEDIMENT
          DO itrc=1,NST
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idfrac(itrc)))) THEN
              got_var(idfrac(itrc))=.TRUE.
              S(ng)%pioVar(idfrac(itrc))%vd=var_desc(i)
              S(ng)%pioVar(idfrac(itrc))%dkind=PIO_FOUT
              S(ng)%pioVar(idfrac(itrc))%gtype=0
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idBmas(itrc)))) THEN
              got_var(idBmas(itrc))=.TRUE.
              S(ng)%pioVar(idBmas(itrc))%vd=var_desc(i)
              S(ng)%pioVar(idBmas(itrc))%dkind=PIO_FOUT
              S(ng)%pioVar(idBmas(itrc))%gtype=0
            END IF
          END DO
          DO itrc=1,MBEDP
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSbed(itrc)))) THEN
              got_var(idSbed(itrc))=.TRUE.
              S(ng)%pioVar(idSbed(itrc))%vd=var_desc(i)
              S(ng)%pioVar(idSbed(itrc))%dkind=PIO_FOUT
              S(ng)%pioVar(idSbed(itrc))%gtype=0
            END IF
          END DO
#   endif
#   if defined SEDIMENT || defined BBL_MODEL
          DO itrc=1,MBOTP
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idBott(itrc)))) THEN
              got_var(idBott(itrc))=.TRUE.
              S(ng)%pioVar(idBott(itrc))%vd=var_desc(i)
              S(ng)%pioVar(idBott(itrc))%dkind=PIO_FOUT
              S(ng)%pioVar(idBott(itrc))%gtype=0
            END IF
          END DO
#   endif
        END DO
!
!  Check if station variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   if defined SEDIMENT && defined SED_MORPH
        IF (.not.got_var(idbath).and.VarOut(idbath,ng)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idbath)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   ifdef SEDIMENT
        DO i=1,NST
          IF (.not.got_var(idfrac(i)).and.VarOut(idfrac(i),ng)) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idfrac(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idBmas(i)).and.VarOut(idBmas(i),ng)) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idBmas(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO i=1,MBEDP
          IF (.not.got_var(idSbed(i)).and.VarOut(idSbed(i),ng)) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idSbed(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#   endif
#   if defined SEDIMENT || defined BBL_MODEL
        DO i=1,MBOTP
          IF (.not.got_var(idBott(i)).and.VarOut(idBott(i),ng)) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idBott(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#   endif

      END IF QUERY
!
  10  FORMAT (1pe11.4,1x,'millimeter')
  20  FORMAT (/,' SEDIMENT_DEF_STATION_PIO - unable to find variable:', &
     &        1x,a,2x,' in stations NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE sediment_def_station_pio
!
!***********************************************************************
      SUBROUTINE sediment_wrt_station_pio (ng, model, tile,             &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     VarOut, S)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: Cgrid
!
      integer :: NposB
      integer :: i, k, np, status
!
      real(dp) :: scale
!
      real(r8), dimension(Nstation(ng)) :: Xpos, Ypos, Zpos, psta

#   ifdef SEDIMENT
      real(r8), dimension(Nstation(ng)*Nbed) :: XposB, YposB, ZposB
      real(r8), dimension(Nstation(ng)*Nbed) :: bsta
#   endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", sediment_wrt_station_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out sediment output variables into specified stations NetCDF
!  file.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set switch to extract station data at native C-grid position (TRUE)
!  or at RHO-points (FALSE).
!
#   ifdef STATIONS_CGRID
      Cgrid=.TRUE.
#   else
      Cgrid=.FALSE.
#   endif
!
!  Set positions for generic extraction routine.
!
      NposB=Nstation(ng)*Nbed
      DO i=1,Nstation(ng)
        Xpos(i)=SCALARS(ng)%SposX(i)
        Ypos(i)=SCALARS(ng)%SposY(i)
        Zpos(i)=1.0_r8
#   ifdef SEDIMENT
        DO k=1,Nbed
          np=k+(i-1)*Nbed
          XposB(np)=SCALARS(ng)%SposX(i)
          YposB(np)=SCALARS(ng)%SposY(i)
          ZposB(np)=REAL(k,r8)
        END DO
#   endif
      END DO

#   if defined SEDIMENT && defined SED_MORPH
!
!  Write out time-varying bathymetry.
!
      IF (VarOut(idbath,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idbath, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%h,                          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idbath)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idbath)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif

#   ifdef SEDIMENT
!
!  Write out sediment fraction of each size class in each bed layer.
!
      DO i=1,NST
        IF (VarOut(idfrac(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, model, Cgrid, idfrac(i), b3dvar,      &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        scale, SEDBED(ng)%bed_frac(:,:,:,i),      &
     &                        NposB, XposB, YposB, ZposB, bsta)
          CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,              &
     &                              TRIM(Vname(1,idfrac(i))), bsta,     &
     &                              (/1,1,S(ng)%Rindex/),               &
     &                              (/Nbed,Nstation(ng),1/),            &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idfrac(i))%vd)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Write out sediment mass of each size class in each bed layer.
!
        IF (VarOut(idBmas(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, model, Cgrid, idBmas(i), b3dvar,      &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        scale,                                    &
     &                        SEDBED(ng)%bed_mass(:,:,:,NOUT,i),        &
     &                        NposB, XposB, YposB, ZposB, bsta)
          CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,              &
     &                              TRIM(Vname(1,idBmas(i))), bsta,     &
     &                              (/1,1,S(ng)%Rindex/),               &
     &                              (/Nbed,Nstation(ng),1/),            &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idBmas(i))%vd)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
!
!  Write out sediment properties in each bed layer.
!
      DO i=1,MBEDP
        IF (VarOut(idSbed(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, model, Cgrid, idSbed(i), b3dvar,      &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        scale, SEDBED(ng)%bed(:,:,:,i),           &
     &                        NposB, XposB, YposB, ZposB, bsta)
          CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,              &
     &                              TRIM(Vname(1,idSbed(i))), bsta,     &
     &                              (/1,1,S(ng)%Rindex/),               &
     &                              (/Nbed,Nstation(ng),1/),            &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idSbed(i))%vd)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
#   endif

#   if defined SEDIMENT || defined BBL_MODEL
!
!  Write out exposed sediment layer properties.
!
      DO i=1,MBEDP
        IF (VarOut(idBott(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta2d (ng, model, Cgrid, idBott(i), r2dvar,      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        scale, SEDBED(ng)%bottom(:,:,i),          &
     &                        Nstation(ng), Xpos, Ypos, psta)
          CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,              &
     &                              TRIM(Vname(1,idBott(i))), psta,     &
     &                              (/1,S(ng)%Rindex/),                 &
     &                              (/Nstation(ng),1/),                 &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idBott(i))%vd)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
#   endif
!
      RETURN
      END SUBROUTINE sediment_wrt_station_pio
#  endif
# endif
!
#endif
      END MODULE sediment_output_mod
