#include "cppdefs.h"
#if defined TL_IOMS && defined TIMELESS_DATA
      SUBROUTINE rp_get_idata (ng)
!
!svn $Id: rp_get_idata.F 429 2009-12-20 17:30:26Z arango $
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
        IF (exit_flag.ne.NoError) RETURN
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
        IF (exit_flag.ne.NoError) RETURN

        CALL get_2dfld (ng, iRPM, idTzph, ncFRCid(idTzph,ng),           &
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
        IF (exit_flag.ne.NoError) RETURN

        CALL get_2dfld (ng, iRPM, idTvph, ncFRCid(idTvph,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  LBi, UBi, LBj, UBj, MTC, NTC(ng),               &
#  ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#  endif
     &                  TIDES(ng) % UV_Tphase(LBi,LBj,1))
        IF (exit_flag.ne.NoError) RETURN

        CALL get_2dfld (ng, iRPM, idTvma, ncFRCid(idTvma,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  LBi, UBi, LBj, UBj, MTC, NTC(ng),               &
#  ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#  endif
     &                  TIDES(ng) % UV_Tmajor(LBi,LBj,1))
        IF (exit_flag.ne.NoError) RETURN

        CALL get_2dfld (ng, iRPM, idTvmi, ncFRCid(idTvmi,ng),           &
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
        IF (exit_flag.ne.NoError) RETURN

        CALL get_ngfld (ng, iRPM, idRepo, ncFRCid(idRepo,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  1, Nsrc(ng), 1, 1, 1, Nsrc(ng), 1,              &
     &                  SOURCES(ng) % Ysrc(1))
        IF (exit_flag.ne.NoError) RETURN

        CALL get_ngfld (ng, iRPM, idRdir, ncFRCid(idRdir,ng),           &
     &                  nFfiles(ng), FRCname(1,ng), update(1),          &
     &                  1, Nsrc(ng), 1, 1, 1, Nsrc(ng), 1,              &
     &                  SOURCES(ng) % Dsrc(1))
        IF (exit_flag.ne.NoError) RETURN

        CALL get_ngfld (ng, iRPM, idRvsh, ncFRCid(idRvsh,ng),           &
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
      CALL wclock_off (ng, iRPM, 3)
# endif

      RETURN
      END SUBROUTINE rp_get_idata
#else
      SUBROUTINE rp_get_idata
      RETURN
      END SUBROUTINE rp_get_idata
#endif
