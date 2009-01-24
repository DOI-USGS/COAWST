#include "cppdefs.h"
#if defined TL_IOMS && defined TIMELESS_DATA
      SUBROUTINE rp_get_idata (ng)
!
!svn $Id: rp_get_idata.F 503 2008-01-10 00:11:51Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
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
      USE mod_scalars
# if defined UV_PSOURCE || defined TS_PSOURCE || defined Q_PSOURCE
      USE mod_sources
# endif
      USE mod_stepping
# if defined SSH_TIDES_NOT_YET || defined UV_TIDES_NOT_YET
      USE mod_tides
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

      real(r8) :: time_save = 0.0_r8
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
      CALL wclock_on (ng, iRPM, 3)
# endif
# if defined SSH_TIDES_NOT_YET || defined UV_TIDES_NOT_YET
!
!-----------------------------------------------------------------------
!  Tide period, amplitude, phase, and currents.
!-----------------------------------------------------------------------
!
!  Tidal Period.
!
      IF (iic(ng).eq.0) THEN
        NTC(ng)=0
        CALL get_ngfld (ng, iRPM, idTper, ncFRCid(idTper,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  1, MTC, 1, 1, 1, NTC(ng), 1,                    &
     &                  TIDES(ng) % Tperiod(1))
      END IF
# endif
# ifdef SSH_TIDES_NOT_YET
!
!  Tidal elevation amplitude and phase. In order to read data as a
!  function of tidal period, we need to reset the model time variables
!  temporarily.
!
      IF (iic(ng).eq.0) THEN
        time_save=time(ng)
        time(ng)=8640000.0_r8
        tdays(ng)=time(ng)*sec2day
        CALL get_2dfld (ng, iRPM, idTzam, ncFRCid(idTzam,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  LBi, UBi, LBj, UBj, MTC, NTC(ng),               &
#  ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#  endif
     &                  TIDES(ng) % SSH_Tamp(LBi,LBj,1))
        CALL get_2dfld (ng, iRPM, idTzph, ncFRCid(idTzph,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  LBi, UBi, LBj, UBj, MTC, NTC(ng),               &
#  ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#  endif
     &                  TIDES(ng) % SSH_Tphase(LBi,LBj,1))
        time(ng)=time_save
        tdays(ng)=time(ng)*sec2day
      END IF
# endif
# ifdef UV_TIDES_NOT_YET
!
!  Tidal currents angle, phase, major and minor ellipse axis.
!
      IF (iic(ng).eq.0) THEN
        time_save=time(ng)
        time(ng)=8640000.0_r8
        tdays(ng)=time(ng)*sec2day
        CALL get_2dfld (ng, iRPM, idTvan, ncFRCid(idTvan,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  LBi, UBi, LBj, UBj, MTC, NTC(ng),               &
#  ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#  endif
     &                  TIDES(ng) % UV_Tangle(LBi,LBj,1))
        CALL get_2dfld (ng, iRPM, idTvph, ncFRCid(idTvph,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  LBi, UBi, LBj, UBj, MTC, NTC(ng),               &
#  ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#  endif
     &                  TIDES(ng) % UV_Tphase(LBi,LBj,1))
        CALL get_2dfld (ng, iRPM, idTvma, ncFRCid(idTvma,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  LBi, UBi, LBj, UBj, MTC, NTC(ng),               &
#  ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#  endif
     &                  TIDES(ng) % UV_Tmajor(LBi,LBj,1))
        CALL get_2dfld (ng, iRPM, idTvmi, ncFRCid(idTvmi,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  LBi, UBi, LBj, UBj, MTC, NTC(ng),               &
#  ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#  endif
     &                  TIDES(ng) % UV_Tminor(LBi,LBj,1))
        time(ng)=time_save
        tdays(ng)=time(ng)*sec2day
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
        CALL get_ngfld (ng, iRPM, idRxpo, ncFRCid(idRxpo,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  1, Nsrc(ng), 1, 1, 1, Nsrc(ng), 1,              &
     &                  SOURCES(ng) % Xsrc(1))
        CALL get_ngfld (ng, iRPM, idRepo, ncFRCid(idRepo,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  1, Nsrc(ng), 1, 1, 1, Nsrc(ng), 1,              &
     &                  SOURCES(ng) % Ysrc(1))
        CALL get_ngfld (ng, iRPM, idRdir, ncFRCid(idRdir,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  1, Nsrc(ng), 1, 1, 1, Nsrc(ng), 1,              &
     &                  SOURCES(ng) % Dsrc(1))
        CALL get_ngfld (ng, iRPM, idRvsh, ncFRCid(idRvsh,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  1, Nsrc(ng), N(ng), 1, 1, Nsrc(ng), N(ng),      &
     &                  SOURCES(ng) % Qshape(1,1))
#  ifdef TS_PSOURCE
        CALL get_ngfld (ng, iRPM, idRflg, ncFRCid(idRflg,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  1, Nsrc(ng), 1, 1, 1, Nsrc(ng), 1,              &
     &                  SOURCES(ng) % Fsrc(1))
        IF (update(1)) THEN
          DO itrc=1,NT(ng)
            SOURCES(ng)%Ltracer(itrc)=.FALSE.
            DO is=1,Nsrc(ng)
              SOURCES(ng)%Lsrc(is,itrc)=.FALSE.
            END DO
          END DO
          DO is=1,Nsrc(ng)
            IF (SOURCES(ng)%Fsrc(is).eq.1.0_r8) THEN
              SOURCES(ng)%Lsrc(is,itemp)=.TRUE.
              SOURCES(ng)%Ltracer(itemp)=.TRUE.
            END IF
            IF (SOURCES(ng)%Fsrc(is).eq.2.0_r8) THEN
              SOURCES(ng)%Lsrc(is,isalt)=.TRUE.
              SOURCES(ng)%Ltracer(isalt)=.TRUE.
            END IF
            IF (SOURCES(ng)%Fsrc(is).ge.3.0_r8) THEN
              SOURCES(ng)%Lsrc(is,itemp)=.TRUE.
              SOURCES(ng)%Lsrc(is,isalt)=.TRUE.
              SOURCES(ng)%Ltracer(itemp)=.TRUE.
              SOURCES(ng)%Ltracer(isalt)=.TRUE.
            END IF
#   if defined RIVER_SEDIMENT && defined SEDIMENT
            IF (SOURCES(ng)%Fsrc(is).ge.4.0_r8) THEN
              DO itrc=1,NST
                SOURCES(ng)%Lsrc(is,idsed(itrc))=.TRUE.
                SOURCES(ng)%Ltracer(idsed(itrc))=.TRUE.
              END DO
            END IF
#   endif
#   if defined RIVER_BIOLOGY && defined BIOLOGY
            IF (SOURCES(ng)%Fsrc(is).ge.5.0_r8) THEN
              DO itrc=1,NBT
                SOURCES(ng)%Lsrc(is,idbio(itrc))=.TRUE.
                SOURCES(ng)%Ltracer(idbio(itrc))=.TRUE.
              END DO
            END IF
#   endif
          END DO
        END IF
#  endif
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
      CALL wclock_off (ng, iRPM, 3)
# endif
      RETURN
      END SUBROUTINE rp_get_idata
#else 
      SUBROUTINE rp_get_idata
      RETURN
      END SUBROUTINE rp_get_idata
#endif
