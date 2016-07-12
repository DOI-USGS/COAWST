#include "cppdefs.h"
      SUBROUTINE get_idata (ng)
!
!svn $Id: get_idata.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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
#ifdef RED_TIDE
      USE mod_biology
#endif
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
#if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      USE mod_netcdf
#endif
#ifdef RED_TIDE
      USE mod_ocean
#endif
      USE mod_parallel
      USE mod_scalars
      USE mod_sources
      USE mod_stepping
#if defined SSH_TIDES || defined UV_TIDES || defined POT_TIDES
      USE mod_tides
#endif
#if defined M2TIDE_DIFF
      USE mod_mixing
#endif
!
      USE nf_fread3d_mod, ONLY : nf_fread3d
#ifdef SOLVE3D
      USE nf_fread4d_mod, ONLY : nf_fread4d
#endif
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

#if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      integer :: gtype, status, varid, Vsize(4)

      real(r8), parameter :: Fscl = 1.0_r8

      real(r8) :: Fmin, Fmax, Htime
#endif
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

#if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Set Vsize to zero to deactivate interpolation of input data to model
!  grid in "nf_fread2d" and "nf_fread3d".
!
      DO is=1,4
        Vsize(is)=0
      END DO
#endif

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 3)
#endif
#if defined SSH_TIDES || defined UV_TIDES || defined POT_TIDES
!
!-----------------------------------------------------------------------
!  Tide period, amplitude, phase, and currents.  In grid refinement,
!  only the coarser grid (RefineScale(ng)=0) tidal forcing data is
!  processed and needed.
!-----------------------------------------------------------------------
!
!  Tidal Period.
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (iic(ng).eq.0) THEN
          CALL get_ngfld (ng, iNLM, idTper, ncFRCid(idTper,ng),         &
     &                    nFfiles(ng), FRC(1,ng), update(1),            &
     &                    1, MTC, 1, 1, 1, NTC(ng), 1,                  &
     &                    TIDES(ng) % Tperiod)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
#endif

#ifdef SSH_TIDES
!
!  Tidal elevation amplitude and phase. In order to read data as a
!  function of tidal period, we need to reset the model time variables
!  temporarily.
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (iic(ng).eq.0) THEN
          time_save=time(ng)
          time(ng)=8640000.0_r8
          tdays(ng)=time(ng)*sec2day

          CALL get_2dfld (ng, iNLM, idTzam, ncFRCid(idTzam,ng),         &
     &                    nFfiles(ng), FRC(1,ng), update(1),            &
     &                    LBi, UBi, LBj, UBj, MTC, NTC(ng),             &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    TIDES(ng) % SSH_Tamp)
          IF (exit_flag.ne.NoError) RETURN

          CALL get_2dfld (ng, iNLM, idTzph, ncFRCid(idTzph,ng),         &
     &                    nFfiles(ng), FRC(1,ng), update(1),            &
     &                    LBi, UBi, LBj, UBj, MTC, NTC(ng),             &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    TIDES(ng) % SSH_Tphase)
          IF (exit_flag.ne.NoError) RETURN

          time(ng)=time_save
          tdays(ng)=time(ng)*sec2day
        END IF
      END IF
#endif

#ifdef UV_TIDES
!
!  Tidal currents angle, phase, major and minor ellipse axis.
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (iic(ng).eq.0) THEN
          time_save=time(ng)
          time(ng)=8640000.0_r8
          tdays(ng)=time(ng)*sec2day

          CALL get_2dfld (ng, iNLM, idTvan, ncFRCid(idTvan,ng),         &
     &                    nFfiles(ng), FRC(1,ng), update(1),            &
     &                    LBi, UBi, LBj, UBj, MTC, NTC(ng),             &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    TIDES(ng) % UV_Tangle)
          IF (exit_flag.ne.NoError) RETURN

          CALL get_2dfld (ng, iNLM, idTvph, ncFRCid(idTvph,ng),         &
     &                    nFfiles(ng), FRC(1,ng), update(1),            &
     &                    LBi, UBi, LBj, UBj, MTC, NTC(ng),             &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    TIDES(ng) % UV_Tphase)
          IF (exit_flag.ne.NoError) RETURN

          CALL get_2dfld (ng, iNLM, idTvma, ncFRCid(idTvma,ng),         &
     &                    nFfiles(ng), FRC(1,ng), update(1),            &
     &                    LBi, UBi, LBj, UBj, MTC, NTC(ng),             &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    TIDES(ng) % UV_Tmajor)
          IF (exit_flag.ne.NoError) RETURN

          CALL get_2dfld (ng, iNLM, idTvmi, ncFRCid(idTvmi,ng),         &
     &                    nFfiles(ng), FRC(1,ng), update(1),            &
     &                    LBi, UBi, LBj, UBj, MTC, NTC(ng),             &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    TIDES(ng) % UV_Tminor)
          IF (exit_flag.ne.NoError) RETURN

          time(ng)=time_save
          tdays(ng)=time(ng)*sec2day
        END IF
      END IF
#endif

#ifdef M2TIDE_DIFF
!
!  Tidal velocity amplitude, M2 only.
!
      IF (iic(ng).eq.0) THEN
        CALL get_2dfld (ng, iNLM, idM2am, ncFRCid(idM2am,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#  endif
     &                  MIXING(ng) % m2_amp)
        IF (exit_flag.ne.NoError) RETURN
      END IF
#endif

#if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
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
        CALL netcdf_get_fvar (ng, iNLM, TIDE(ng)%name, Vname(1,idtime), &
     &                        Htime,                                    &
     &                        ncid = TIDE(ng)%ncid)
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
        CALL netcdf_get_ivar (ng, iNLM, TIDE(ng)%name, 'Hcount',        &
     &                        Hcount(ng),                               &
     &                        ncid = TIDE(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
!  Time-accumulated COS(omega(k)*t) harmonics.
!
        CALL get_ngfld (ng, iNLM, idCosW, ncFRCid(idCosW,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  1, MTC, 1, 1, 1, NTC(ng), 1,                    &
     &                  TIDES(ng) % CosW_sum)
        IF (exit_flag.ne.NoError) RETURN
!
!  Time-accumulated SIN(omega(k)*t) harmonics.
!
        CALL get_ngfld (ng, iNLM, idSinW, ncFRCid(idSinW,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  1, MTC, 1, 1, 1, NTC(ng), 1,                    &
     &                  TIDES(ng) % SinW_sum)
        IF (exit_flag.ne.NoError) RETURN
!
!  Time-accumulated COS(omega(k)*t)*COS(omega(k)*t) harmonics.
!
        CALL get_ngfld (ng, iNLM, idCos2, ncFRCid(idCos2,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  1, MTC, MTC, 1, 1, NTC(ng), NTC(ng),            &
     &                  TIDES(ng) % CosWCosW)
        IF (exit_flag.ne.NoError) RETURN
!
!  Time-accumulated SIN(omega(k)*t)*SIN(omega(k)*t) harmonics.
!
        CALL get_ngfld (ng, iNLM, idSin2, ncFRCid(idSin2,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  1, MTC, MTC, 1, 1, NTC(ng), NTC(ng),            &
     &                  TIDES(ng) %  SinWSinW)
        IF (exit_flag.ne.NoError) RETURN
!
!  Time-accumulated SIN(omega(k)*t)*COS(omega(k)*t) harmonics.
!
        CALL get_ngfld (ng, iNLM, idSWCW, ncFRCid(idSWCW,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  1, MTC, MTC, 1, 1, NTC(ng), NTC(ng),            &
     &                  TIDES(ng) %  SinWCosW)
        IF (exit_flag.ne.NoError) RETURN
!
!  Time-accumulated free-surface tide harmonics.
!
        IF (Aout(idFsuD,ng)) THEN
          gtype=r3dvar
          status=nf_fread3d(ng, iNLM, TIDE(ng)%name, TIDE(ng)%ncid,     &
     &                      Vname(1,idFsuH), TIDE(ng)%Vid(idFsuH),      &
     &                      0, gtype, Vsize,                            &
     &                      LBi, UBi, LBj, UBj, 0, 2*NTC(ng),           &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      TIDES(ng) % zeta_tide)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idFsuH)),                  &
     &                          TRIM(TIDE(ng)%name)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(2,idFsuH)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Time-accumulated 2D u-momentum tide harmonics.
!
        IF (Aout(idu2dD,ng)) THEN
          gtype=u3dvar
          status=nf_fread3d(ng, iNLM, TIDE(ng)%name, TIDE(ng)%ncid,     &
     &                      Vname(1,idu2dH), TIDE(ng)%Vid(idu2dH),      &
     &                      0, gtype, Vsize,                            &
     &                      LBi, UBi, LBj, UBj, 0, 2*NTC(ng),           &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % umask,                           &
# endif
     &                      TIDES(ng) % ubar_tide)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idu2dH)),                  &
     &                          TRIM(TIDE(ng)%name)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(2,idu2dH)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Time-accumulated 2D v-momentum tide harmonics.
!
        IF (Aout(idv2dD,ng)) THEN
          gtype=v3dvar
          status=nf_fread3d(ng, iNLM, TIDE(ng)%name, TIDE(ng)%ncid,     &
     &                      Vname(1,idv2dH), TIDE(ng)%Vid(idv2dH),      &
     &                      0, gtype, Vsize,                            &
     &                      LBi, UBi, LBj, UBj, 0, 2*NTC(ng),           &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % vmask,                           &
# endif
     &                      TIDES(ng) % vbar_tide)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idv2dH)),                  &
     &                          TRIM(TIDE(ng)%name)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(2,idv2dH)), Fmin, Fmax
            END IF
          END IF
        END IF

# ifdef SOLVE3D
!
!  Time-accumulated 3D u-momentum tide harmonics.
!
        IF (Aout(idu3dD,ng)) THEN
          gtype=u3dvar
          status=nf_fread4d(ng, iNLM, TIDE(ng)%name, TIDE(ng)%ncid,     &
     &                      Vname(1,idu3dH), TIDE(ng)%Vid(idu3dH),      &
     &                      0, gtype, Vsize,                            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng), 0, 2*NTC(ng), &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % umask,                           &
#  endif
     &                      TIDES(ng) % u_tide)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idu3dH)),                  &
     &                          TRIM(TIDE(ng)%name)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(2,idu3dH)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Time-accumulated 3D v-momentum tide harmonics.
!
        IF (Aout(idv3dD,ng)) THEN
          gtype=v3dvar
          status=nf_fread4d(ng, iNLM, TIDE(ng)%name, TIDE(ng)%ncid,     &
     &                      Vname(1,idv3dH), TIDE(ng)%Vid(idv3dH),      &
     &                      0, gtype, Vsize,                            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng), 0, 2*NTC(ng), &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % vmask,                           &
#  endif
     &                      TIDES(ng) % v_tide)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idv3dH)),                  &
     &                          TRIM(TIDE(ng)%name)
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
!
!  Time-accumulated temperature and salinity tide harmonics.
!
        DO itrc=1,NAT
          IF (Aout(idTrcD(itrc),ng)) THEN
            gtype=r3dvar
            status=nf_fread4d(ng, iNLM, TIDE(ng)%name, TIDE(ng)%ncid,   &
     &                        Vname(1,idTrcH(itrc)),                    &
     &                        TIDE(ng)%Vid(idTrcH(itrc)),               &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        0, 2*NTC(ng),                             &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
     &                        TIDES(ng) % t_tide(:,:,:,:,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,10) TRIM(Vname(1,idTrcH(itrc))),          &
     &                            TRIM(TIDE(ng)%name)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) TRIM(Vname(2,idTrcH(itrc))),          &
     &                            Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
# endif
      END IF
#endif
#ifdef POT_TIDES
!
!  Tidal potential amplitude and phase. In order to read data as a
!  function of tidal period, we need to reset the model time variables
!  temporarily.
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (iic(ng).eq.0) THEN
          time_save=time(ng)
          time(ng)=8640000.0_r8
          tdays(ng)=time(ng)*sec2day
          CALL get_2dfld (ng, iNLM, idTpam, ncFRCid(idTpam,ng),         &
     &                    nFfiles(ng), FRC(1,ng), update(1),            &
     &                    LBi, UBi, LBj, UBj, MTC, NTC(ng),             &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    TIDES(ng) % POT_Tamp)
          IF (exit_flag.ne.NoError) RETURN

          CALL get_2dfld (ng, iNLM, idTpph, ncFRCid(idTpph,ng),         &
     &                    nFfiles(ng), FRC(1,ng), update(1),            &
     &                    LBi, UBi, LBj, UBj, MTC, NTC(ng),             &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    TIDES(ng) % POT_Tphase)
          IF (exit_flag.ne.NoError) RETURN

          time(ng)=time_save
          tdays(ng)=time(ng)*sec2day
        END IF
      END IF
#endif

#ifndef ANA_PSOURCE
!
!-----------------------------------------------------------------------
!  Read in point Sources/Sinks position, direction, special flag, and
!  mass transport nondimensional shape profile.  Point sources are at
!  U- and V-points.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).eq.0).and.                                           &
     &    (LuvSrc(ng).or.LwSrc(ng).or.ANY(LtracerSrc(:,ng)))) THEN
        CALL get_ngfld (ng, iNLM, idRxpo, SSF(ng)%ncid,                 &
     &                  1, SSF(ng), update(1),                          &
     &                  1, Nsrc(ng), 1, 1, 1, Nsrc(ng), 1,              &
     &                  SOURCES(ng) % Xsrc)
        IF (exit_flag.ne.NoError) RETURN

        CALL get_ngfld (ng, iNLM, idRepo, SSF(ng)%ncid,                 &
     &                  1, SSF(ng), update(1),                          &
     &                  1, Nsrc(ng), 1, 1, 1, Nsrc(ng), 1,              &
     &                  SOURCES(ng) % Ysrc)
        IF (exit_flag.ne.NoError) RETURN

        CALL get_ngfld (ng, iNLM, idRdir, SSF(ng)%ncid,                 &
     &                  1, SSF(ng), update(1),                          &
     &                  1, Nsrc(ng), 1, 1, 1, Nsrc(ng), 1,              &
     &                  SOURCES(ng) % Dsrc)
        IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
        CALL get_ngfld (ng, iNLM, idRvsh, SSF(ng)%ncid,                 &
     &                  1, SSF(ng), update(1),                          &
     &                  1, Nsrc(ng), N(ng), 1, 1, Nsrc(ng), N(ng),      &
     &                  SOURCES(ng) % Qshape)
        IF (exit_flag.ne.NoError) RETURN
# endif

        DO is=1,Nsrc(ng)
          SOURCES(ng)%Isrc(is)=                                         &
     &                MAX(1,MIN(NINT(SOURCES(ng)%Xsrc(is)),Lm(ng)+1))
          SOURCES(ng)%Jsrc(is)=                                         &
     &                MAX(1,MIN(NINT(SOURCES(ng)%Ysrc(is)),Mm(ng)+1))
        END DO
      END IF
#endif

#ifdef RED_TIDE
!
!-----------------------------------------------------------------------
!  Red tides.
!-----------------------------------------------------------------------
!
!  Read initial dynoflagellate bottom cyst concentration.
!
      CALL get_2dfld (ng, iNLM, idCyst, ncFRCid(idCyst,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
     &                OCEAN(ng) % CystIni)
      IF (exit_flag.ne.NoError) RETURN
      IF (.not.Linfo(1,idCyst,ng)) THEN             ! if not gridded,
        OCEAN(ng) % CystIni = Fpoint(1,idCyst,ng)   ! load point data
      END IF
#endif

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iNLM, 3)
#endif

#if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
  10  FORMAT (/,' GET_IDATA - error while reading variable: ',a,        &
     &        /,13x,'in input NetCDF file: ',a)
  20  FORMAT (/,' GET_IDATA - incosistent restart and harmonics time:', &
     &        /,13x,f15.4,2x,f15.4)
  30  FORMAT (16x,'- ',a,/,19x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,')')
#endif

      RETURN
      END SUBROUTINE get_idata
