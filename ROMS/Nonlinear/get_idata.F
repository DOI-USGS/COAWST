#include "cppdefs.h"
#ifdef TIMELESS_DATA
      SUBROUTINE get_idata (ng)
!
!svn $Id: get_idata.F 1597 2008-07-25 02:31:35Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads input data that needs to be obtained only once.  !
!                                                                      !
!  Currently,  this routine is only executed in serial mode by the     !
!  main thread.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      USE mod_netcdf
      USE mod_parallel
# endif
      USE mod_scalars
# if defined UV_PSOURCE || defined TS_PSOURCE || defined Q_PSOURCE
      USE mod_sources
# endif
      USE mod_stepping
# if defined SSH_TIDES || defined UV_TIDES
      USE mod_tides
# endif
!
      USE nf_fread3d_mod, ONLY : nf_fread3d
# ifdef SOLVE3D
      USE nf_fread4d_mod, ONLY : nf_fread4d
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical, dimension(3) :: update =                                 &
     &         (/ .FALSE., .FALSE., .FALSE. /)

      integer :: LBi, UBi, LBj, UBj
      integer :: itrc, is
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      integer :: gtype, status, varid, Vsize(4)

      real(r8), parameter :: Fscl = 1.0_r8

      real(r8) :: Fmin, Fmax, Htime
# endif
      real(r8) :: time_save = 0.0_r8
!
      SourceFile='get_idata.F'
!
!  Lower and upper bounds for tiled arrays.
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 3)
# endif

# ifdef REFINED_GRID
      IF (ng.eq.1) THEN
# endif

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Set Vsize to zero to deactivate interpolation of input data to model
!  grid in "nf_fread2d" and "nf_fread3d".
!
      DO is=1,4
        Vsize(is)=0
      END DO
# endif

# if defined SSH_TIDES || defined UV_TIDES
!
!-----------------------------------------------------------------------
!  Tide period, amplitude, phase, and currents.
!-----------------------------------------------------------------------
!
!  Tidal Period.
!
      IF (iic(ng).eq.0) THEN
        CALL get_ngfld (ng, iNLM, idTper, ncFRCid(idTper,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  1, MTC, 1, 1, 1, NTC(ng), 1,                    &
     &                  TIDES(ng) % Tperiod(1))
        IF (exit_flag.ne.NoError) RETURN
      END IF
# endif

# ifdef SSH_TIDES
!
!  Tidal elevation amplitude and phase. In order to read data as a
!  function of tidal period, we need to reset the model time variables
!  temporarily.
!
      IF (iic(ng).eq.0) THEN
        time_save=time(ng)
        time(ng)=8640000.0_r8
        tdays(ng)=time(ng)*sec2day

        CALL get_2dfld (ng, iNLM, idTzam, ncFRCid(idTzam,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  LBi, UBi, LBj, UBj, MTC, NTC(ng),               &
#  ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#  endif
     &                  TIDES(ng) % SSH_Tamp(LBi,LBj,1))
        IF (exit_flag.ne.NoError) RETURN

        CALL get_2dfld (ng, iNLM, idTzph, ncFRCid(idTzph,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  LBi, UBi, LBj, UBj, MTC, NTC(ng),               &
#  ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#  endif
     &                  TIDES(ng) % SSH_Tphase(LBi,LBj,1))
        IF (exit_flag.ne.NoError) RETURN

        time(ng)=time_save
        tdays(ng)=time(ng)*sec2day
      END IF
# endif

# ifdef UV_TIDES
!
!  Tidal currents angle, phase, major and minor ellipse axis.
!
      IF (iic(ng).eq.0) THEN
        time_save=time(ng)
        time(ng)=8640000.0_r8
        tdays(ng)=time(ng)*sec2day

        CALL get_2dfld (ng, iNLM, idTvan, ncFRCid(idTvan,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  LBi, UBi, LBj, UBj, MTC, NTC(ng),               &
#  ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#  endif
     &                  TIDES(ng) % UV_Tangle(LBi,LBj,1))
        IF (exit_flag.ne.NoError) RETURN

        CALL get_2dfld (ng, iNLM, idTvph, ncFRCid(idTvph,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  LBi, UBi, LBj, UBj, MTC, NTC(ng),               &
#  ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#  endif
     &                  TIDES(ng) % UV_Tphase(LBi,LBj,1))
        IF (exit_flag.ne.NoError) RETURN

        CALL get_2dfld (ng, iNLM, idTvma, ncFRCid(idTvma,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  LBi, UBi, LBj, UBj, MTC, NTC(ng),               &
#  ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#  endif
     &                  TIDES(ng) % UV_Tmajor(LBi,LBj,1))
        IF (exit_flag.ne.NoError) RETURN

        CALL get_2dfld (ng, iNLM, idTvmi, ncFRCid(idTvmi,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  LBi, UBi, LBj, UBj, MTC, NTC(ng),               &
#  ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#  endif
     &                  TIDES(ng) % UV_Tminor(LBi,LBj,1))
        IF (exit_flag.ne.NoError) RETURN

        time(ng)=time_save
        tdays(ng)=time(ng)*sec2day
      END IF
# endif

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!-----------------------------------------------------------------------
!  If detiding and applicable, define additional variable to store
!  time-accumulated tide harmonics variables.  This variable are
!  defined and written into input tide forcing NetCDF file.
!-----------------------------------------------------------------------
!
      CALL def_tides (ng, LdefTIDE(ng))
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  If restarting, read in time-accumulated tide harmonics variables.
!-----------------------------------------------------------------------
!
      IF (.not.LdefTIDE(ng).and.(nrrec(ng).ne.0)) THEN
!
!  For consistency, check time of written accumulate harmonics and
!  compare to current time.
!
        CALL netcdf_get_fvar (ng, iNLM, TIDEname(ng), Vname(1,idtime),  &
     &                        Htime,                                    &
     &                        ncid = ncTIDEid(ng))
        IF (exit_flag.ne.NoError) RETURN

        IF (time(ng).ne.Htime) THEN
          IF (Master) THEN
            WRITE (stdout,20) tdays(ng), Htime*sec2day
          END IF
          exit_flag=2
          ioerror=0
          RETURN
        END IF
!
!  Number of time-acummulate tide harmonics.
!
        CALL netcdf_get_ivar (ng, iNLM, TIDEname(ng), 'Hcount',         &
     &                        Hcount(ng),                               &
     &                        ncid = ncTIDEid(ng))
        IF (exit_flag.ne.NoError) RETURN
!
!  Time-accumulated COS(omega(k)*t) harmonics.
!
        CALL get_ngfld (ng, iNLM, idCosW, ncFRCid(idCosW,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  1, MTC, 1, 1, 1, NTC(ng), 1,                    &
     &                  TIDES(ng) % CosW_sum(1))
        IF (exit_flag.ne.NoError) RETURN
!
!  Time-accumulated SIN(omega(k)*t) harmonics.
!
        CALL get_ngfld (ng, iNLM, idSinW, ncFRCid(idSinW,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  1, MTC, 1, 1, 1, NTC(ng), 1,                    &
     &                  TIDES(ng) % SinW_sum(1))
        IF (exit_flag.ne.NoError) RETURN
!
!  Time-accumulated COS(omega(k)*t)*COS(omega(k)*t) harmonics.
!
        CALL get_ngfld (ng, iNLM, idCos2, ncFRCid(idCos2,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  1, MTC, MTC, 1, 1, NTC(ng), NTC(ng),            &
     &                  TIDES(ng) % CosWCosW(1,1))
        IF (exit_flag.ne.NoError) RETURN
!
!  Time-accumulated SIN(omega(k)*t)*SIN(omega(k)*t) harmonics.
!
        CALL get_ngfld (ng, iNLM, idSin2, ncFRCid(idSin2,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  1, MTC, MTC, 1, 1, NTC(ng), NTC(ng),            &
     &                  TIDES(ng) %  SinWSinW(1,1))
        IF (exit_flag.ne.NoError) RETURN
!
!  Time-accumulated SIN(omega(k)*t)*COS(omega(k)*t) harmonics.
!
        CALL get_ngfld (ng, iNLM, idSWCW, ncFRCid(idSWCW,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  1, MTC, MTC, 1, 1, NTC(ng), NTC(ng),            &
     &                  TIDES(ng) %  SinWCosW(1,1))
        IF (exit_flag.ne.NoError) RETURN
!
!  Time-accumulated free-surface tide harmonics.
!
        gtype=r3dvar
        status=nf_fread3d(ng, iNLM, TIDEname(ng), ncTIDEid(ng),         &
     &                    Vname(1,idFsuH), tideVid(idFsuH,ng),          &
     &                    0, gtype, Vsize,                              &
     &                    LBi, UBi, LBj, UBj, 0, 2*NTC(ng),             &
     &                    Fscl, Fmin, Fmax,                             &
#  ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#  endif
     &                    TIDES(ng) % zeta_tide)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsuH)), TRIM(TIDEname(ng))
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(2,idFsuH)), Fmin, Fmax
          END IF
        END IF
!
!  Time-accumulated 2D u-momentum tide harmonics.
!
        gtype=u3dvar
        status=nf_fread3d(ng, iNLM, TIDEname(ng), ncTIDEid(ng),         &
     &                    Vname(1,idu2dH), tideVid(idu2dH,ng),          &
     &                    0, gtype, Vsize,                              &
     &                    LBi, UBi, LBj, UBj, 0, 2*NTC(ng),             &
     &                    Fscl, Fmin, Fmax,                             &
#  ifdef MASKING
     &                    GRID(ng) % umask,                             &
#  endif
     &                    TIDES(ng) % ubar_tide)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu2dH)), TRIM(TIDEname(ng))
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(2,idu2dH)), Fmin, Fmax
          END IF
        END IF
!
!  Time-accumulated 2D v-momentum tide harmonics.
!
        gtype=v3dvar
        status=nf_fread3d(ng, iNLM, TIDEname(ng), ncTIDEid(ng),         &
     &                    Vname(1,idv2dH), tideVid(idv2dH,ng),          &
     &                    0, gtype, Vsize,                              &
     &                    LBi, UBi, LBj, UBj, 0, 2*NTC(ng),             &
     &                    Fscl, Fmin, Fmax,                             &
#  ifdef MASKING
     &                    GRID(ng) % vmask,                             &
#  endif
     &                    TIDES(ng) % vbar_tide)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv2dH)), TRIM(TIDEname(ng))
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(2,idv2dH)), Fmin, Fmax
          END IF
        END IF

#  ifdef SOLVE3D
!
!  Time-accumulated 3D u-momentum tide harmonics.
!
        gtype=u3dvar
        status=nf_fread4d(ng, iNLM, TIDEname(ng), ncTIDEid(ng),         &
     &                    Vname(1,idu3dH), tideVid(idu3dH,ng),          &
     &                    0, gtype, Vsize,                              &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 0, 2*NTC(ng),   &
     &                    Fscl, Fmin, Fmax,                             &
#   ifdef MASKING
     &                    GRID(ng) % umask,                             &
#   endif
     &                    TIDES(ng) % u_tide)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu3dH)), TRIM(TIDEname(ng))
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(2,idu3dH)), Fmin, Fmax
          END IF
        END IF
!
!  Time-accumulated 3D v-momentum tide harmonics.
!
        gtype=v3dvar
        status=nf_fread4d(ng, iNLM, TIDEname(ng), ncTIDEid(ng),         &
     &                    Vname(1,idv3dH), tideVid(idv3dH,ng),          &
     &                    0, gtype, Vsize,                              &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 0, 2*NTC(ng),   &
     &                    Fscl, Fmin, Fmax,                             &
#   ifdef MASKING
     &                    GRID(ng) % vmask,                             &
#   endif
     &                    TIDES(ng) % v_tide)
#  endif
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv3dH)), TRIM(TIDEname(ng))
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(2,idv3dH)), Fmin, Fmax
          END IF
        END IF
      END IF
# endif

# ifdef REFINED_GRID
      END IF
# endif

# if !defined ANA_PSOURCE && (defined UV_PSOURCE || \
                              defined TS_PSOURCE || defined Q_PSOURCE)
!
!-----------------------------------------------------------------------
!  Point Sources/Sinks position, direction, special flag, and mass
!  transport nondimensional shape profile.  Point sources are at U-
!  and V-points.
!-----------------------------------------------------------------------
!
      IF (iic(ng).eq.0) THEN
        CALL get_ngfld (ng, iNLM, idRxpo, ncFRCid(idRxpo,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  1, Nsrc(ng), 1, 1, 1, Nsrc(ng), 1,              &
     &                  SOURCES(ng) % Xsrc(1))
        IF (exit_flag.ne.NoError) RETURN

        CALL get_ngfld (ng, iNLM, idRepo, ncFRCid(idRepo,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  1, Nsrc(ng), 1, 1, 1, Nsrc(ng), 1,              &
     &                  SOURCES(ng) % Ysrc(1))
        IF (exit_flag.ne.NoError) RETURN

        CALL get_ngfld (ng, iNLM, idRdir, ncFRCid(idRdir,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  1, Nsrc(ng), 1, 1, 1, Nsrc(ng), 1,              &
     &                  SOURCES(ng) % Dsrc(1))
        IF (exit_flag.ne.NoError) RETURN

        CALL get_ngfld (ng, iNLM, idRvsh, ncFRCid(idRvsh,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  1, Nsrc(ng), N(ng), 1, 1, Nsrc(ng), N(ng),      &
     &                  SOURCES(ng) % Qshape(1,1))
        IF (exit_flag.ne.NoError) RETURN

        DO is=1,Nsrc(ng)
          SOURCES(ng)%Isrc(is)=                                         &
     &                MAX(1,MIN(NINT(SOURCES(ng)%Xsrc(is)),Lm(ng)+1))
          SOURCES(ng)%Jsrc(is)=                                         &
     &                MAX(1,MIN(NINT(SOURCES(ng)%Ysrc(is)),Mm(ng)+1))
        END DO
      END IF
# endif

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iNLM, 3)
# endif

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
  10  FORMAT (/,' GET_IDATA - error while reading variable: ',a,        &
     &        /,13x,'in input NetCDF file: ',a)
  20  FORMAT (/,' GET_IDATA - incosistent restart and harmonics time:', &
     &        /,13x,f15.4,2x,f15.4)
  30  FORMAT (16x,'- ',a,/,19x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,')')
# endif

      RETURN
      END SUBROUTINE get_idata
#else
      SUBROUTINE get_idata
      RETURN
      END SUBROUTINE get_idata
#endif
