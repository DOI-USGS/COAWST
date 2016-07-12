#include "cppdefs.h"
#if (defined FOUR_DVAR || defined VERIFICATION) && defined OBSERVATIONS
      SUBROUTINE obs_read (ng, model, backward)
!
!svn $Id: obs_read.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine reads in observations data when appropriate from    !
!  observations input NetCDF file.  The observations data is stored    !
!  for use elsewhere.                                                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: backward

      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      logical :: readNLmod, readTLmod

      integer :: Mstr, Mend
      integer :: i, iobs, itrc, status
!
      SourceFile='obs_read.F'
!
!---------------------------------------------------------------------
!  Read observation variables needed for interpolating the model
!  state at the observation locations.
!---------------------------------------------------------------------
!
      IF (ProcessObs(ng)) THEN
# if defined TLM_OBS
#  if defined IS4DVAR_SENSITIVITY
        readNLmod=.FALSE.
#  else
        readNLmod=.TRUE.
#  endif
        readTLmod=.TRUE.
# else
        readNLmod=.FALSE.
        readTLmod=.FALSE.
# endif
!
!  Initialize observations processing counters.
!
        DO i=1,NstateVar(ng)
          FOURDVAR(ng)%ObsCount(i)=0
          FOURDVAR(ng)%ObsReject(i)=0
        END DO
        IF (backward) THEN
          ObsSurvey(ng)=ObsSurvey(ng)-1
        ELSE
          ObsSurvey(ng)=ObsSurvey(ng)+1
        END IF
!
!  Set number of observations to process.
!
        Nobs(ng)=FOURDVAR(ng)%NobsSurvey(ObsSurvey(ng))
!
!  Set number of datum to process at current time-step.
!
        IF (backward) THEN
          NendObs(ng)=NstrObs(ng)-1
          NstrObs(ng)=NstrObs(ng)-Nobs(ng)
        ELSE
          NstrObs(ng)=NendObs(ng)+1
          NendObs(ng)=NstrObs(ng)+Nobs(ng)-1
        END IF
!
!  Set starting index of obervation vectors for reading.  In weak
!  constraint, the entire observation data is loaded. Otherwise,
!  only the observartion for the current time window are loaded
!  and started from vector index one.
!
# ifdef WEAK_CONSTRAINT
        Mstr=NstrObs(ng)
        Mend=NendObs(ng)
# else
        Mstr=1
        Mend=Nobs(ng)
# endif
!
!  Read in observation type identifier.
!
        CALL netcdf_get_ivar (ng, model, OBS(ng)%name, Vname(1,idOtyp), &
     &                        ObsType(Mstr:),                           &
     &                        ncid = OBS(ng)%ncid,                      &
     &                        start = (/NstrObs(ng)/),                  &
     &                        total = (/Nobs(ng)/))
        IF (exit_flag.ne.NoError) RETURN
!
!  Read in observation provenance identifier.
!
        CALL netcdf_get_ivar (ng, model, OBS(ng)%name, Vname(1,idOpro), &
     &                        ObsProv(Mstr:),                           &
     &                        ncid = OBS(ng)%ncid,                      &
     &                        start = (/NstrObs(ng)/),                  &
     &                        total = (/Nobs(ng)/))
        IF (exit_flag.ne.NoError) RETURN
!
!  Read in observation time (days).
!
        CALL netcdf_get_fvar (ng, model, OBS(ng)%name, Vname(1,idObsT), &
     &                        Tobs(Mstr:),                              &
     &                        ncid = OBS(ng)%ncid,                      &
     &                        start = (/NstrObs(ng)/),                  &
     &                        total = (/Nobs(ng)/))
        IF (exit_flag.ne.NoError) RETURN
!
!  Read in observation X-location (grid units).
!
        CALL netcdf_get_fvar (ng, model, OBS(ng)%name, Vname(1,idObsX), &
     &                        Xobs(Mstr:),                              &
     &                        ncid = OBS(ng)%ncid,                      &
     &                        start = (/NstrObs(ng)/),                  &
     &                        total = (/Nobs(ng)/))
        IF (exit_flag.ne.NoError) RETURN
!
!  Read in observation Y-location (grid units).
!
        CALL netcdf_get_fvar (ng, model, OBS(ng)%name, Vname(1,idObsY), &
     &                        Yobs(Mstr:),                              &
     &                        ncid = OBS(ng)%ncid,                      &
     &                        start = (/NstrObs(ng)/),                  &
     &                        total = (/Nobs(ng)/))
        IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
!
!  Read in observation depth, Zobs.  If negative, depth is meter. If
!  greater than zero, Zobs is in model fractional vertical levels
!  (1 <= Zobs <= N). If Zobs < 0, its fractional level value is
!  computed in routine "extract_obs3d" and over-written so it can
!  be written into the observation NetCDF file for latter use.
!
        IF (wrote_Zobs(ng)) THEN
          CALL netcdf_get_fvar (ng, model, OBS(ng)%name,                &
     &                          Vname(1,idObsZ),                        &
     &                          Zobs(Mstr:),                            &
     &                          ncid = OBS(ng)%ncid,                    &
     &                          start = (/NstrObs(ng)/),                &
     &                          total = (/Nobs(ng)/))
          IF (exit_flag.ne.NoError) RETURN
        ELSE
          CALL netcdf_get_fvar (ng, model, OBS(ng)%name,                &
     &                          Vname(1,idObsD),                        &
     &                          Zobs(Mstr:),                            &
     &                          ncid = OBS(ng)%ncid,                    &
     &                          start = (/NstrObs(ng)/),                &
     &                          total = (/Nobs(ng)/))
          IF (exit_flag.ne.NoError) RETURN
        END IF
        Load_Zobs(ng)=.FALSE.
        IF ((MINVAL(Zobs).lt.0.0_r8).or.                                &
     &      (MAXVAL(Zobs).lt.0.0_r8)) THEN
          Load_Zobs(ng)=.TRUE.
        END IF

#  ifdef DISTRIBUTE
!
!  If distributed-memory and Zobs in meters (Zobs < 0),  zero-out
!  Zobs values in all nodes but itself to facilitate exchages between
!  tiles latter before writting into observation NetCDF file.
!
        IF (.not.wrote_Zobs(ng)) THEN
          CALL obs_depth (ng, MyRank, model)
        END IF
#  endif
# endif
!
!  Read in observation values.
!
        CALL netcdf_get_fvar (ng, model, OBS(ng)%name, Vname(1,idOval), &
     &                        ObsVal(Mstr:),                            &
     &                        ncid = OBS(ng)%ncid,                      &
     &                        start = (/NstrObs(ng)/),                  &
     &                        total = (/Nobs(ng)/))
        IF (exit_flag.ne.NoError) RETURN
!
# ifdef WEAK_CONSTRAINT
!  Read in observation error covariance.
# else
!  Read in observation error covariance. To avoid successive divisions,
!  convert to inverse observation error covariance.
# endif
!
        CALL netcdf_get_fvar (ng, model, OBS(ng)%name, Vname(1,idOerr), &
     &                        ObsErr(Mstr:),                            &
     &                        ncid = OBS(ng)%ncid,                      &
     &                        start = (/NstrObs(ng)/),                  &
     &                        total = (/Nobs(ng)/))
        IF (exit_flag.ne.NoError) RETURN

# ifndef WEAK_CONSTRAINT
        DO iobs=1,Nobs(ng)
          ObsErr(iobs)=1.0_r8/ObsErr(iobs)
        END DO
# endif
!
!  Read in nonlinear model values at observation locations.
!
        IF (readNLmod.and.haveNLmod(ng)) THEN
# ifdef VERIFICATION
          CALL netcdf_get_fvar (ng, model, DAV(ng)%name,                &
     &                          Vname(1,idNLmo),                        &
     &                          NLmodVal(Mstr:),                        &
     &                          ncid = DAV(ng)%ncid,                    &
     &                          start = (/NstrObs(ng)/),                &
     &                          total = (/Nobs(ng)/))
          IF (exit_flag.ne.NoError) RETURN
# else
          CALL netcdf_get_fvar (ng, model, DAV(ng)%name,                &
     &                          Vname(1,idNLmo),                        &
     &                          NLmodVal(Mstr:),                        &
     &                          ncid = DAV(ng)%ncid,                    &
     &                          start = (/NstrObs(ng),outer/),          &
     &                          total = (/Nobs(ng),1/))
          IF (exit_flag.ne.NoError) RETURN
# endif
# ifdef BGQC
          CALL netcdf_get_fvar (ng, model, DAV(ng)%name,                &
     &                          Vname(1,idBgEr),                        &
     &                          BgErr(Mstr:),                           &
     &                          ncid = DAV(ng)%ncid,                    &
     &                          start = (/NstrObs(ng)/),                &
     &                          total = (/Nobs(ng)/))
          IF (exit_flag.ne.NoError) RETURN
          CALL netcdf_get_fvar (ng, model, DAV(ng)%name,                &
     &                          Vname(1,idObsS),                        &
     &                          ObsScale(Mstr:),                        &
     &                          ncid = DAV(ng)%ncid,                    &
     &                          start = (/NstrObs(ng)/),                &
     &                          total = (/Nobs(ng)/))
          IF (exit_flag.ne.NoError) RETURN
# endif
        END IF

# if defined TLM_OBS && !defined WEAK_CONSTRAINT
!
!  If adjoint pass and incremental 4DVar, read in tangent linear model
!  values at observation locations.
!
        IF (readTLmod.and.haveTLmod(ng)) THEN
          CALL netcdf_get_fvar (ng, model, DAV(ng)%name,                &
     &                          Vname(1,idTLmo),                        &
     &                          TLmodVal(Mstr:),                        &
     &                          ncid = DAV(ng)%ncid,                    &
     &                          start = (/NstrObs(ng)/),                &
     &                          total = (/Nobs(ng)/))
          IF (exit_flag.ne.NoError) RETURN

#  if defined IS4DVAR
!
!  Reset TLM values to zero at the first pass of the inner loop.
!
          IF (inner.eq.0) THEN
            DO iobs=1,Mobs
              TLmodVal(iobs)=0.0_r8
            END DO
          END IF
#  endif
        END IF
# endif
!
!-----------------------------------------------------------------------
!  Set counters for number of observations to processed for each state
!  variable.
!-----------------------------------------------------------------------
!
        DO iobs=Mstr,Mend
          IF  (ObsType(iobs).eq.isFsur) THEN
            FOURDVAR(ng)%ObsCount(isFsur)=                              &
     &                           FOURDVAR(ng)%ObsCount(isFsur)+1
          ELSE IF (ObsType(iobs).eq.isUbar) THEN
            FOURDVAR(ng)%ObsCount(isUbar)=                              &
     &                           FOURDVAR(ng)%ObsCount(isUbar)+1
          ELSE IF (ObsType(iobs).eq.isVbar) THEN
            FOURDVAR(ng)%ObsCount(isVbar)=                              &
     &                           FOURDVAR(ng)%ObsCount(isVbar)+1
# ifdef SOLVE3D
          ELSE IF (ObsType(iobs).eq.isUvel) THEN
            FOURDVAR(ng)%ObsCount(isUvel)=                              &
     &                           FOURDVAR(ng)%ObsCount(isUvel)+1
          ELSE IF (ObsType(iobs).eq.isVvel) THEN
            FOURDVAR(ng)%ObsCount(isVvel)=                              &
     &                           FOURDVAR(ng)%ObsCount(isVvel)+1
          ELSE
            DO itrc=1,NT(ng)
              IF (ObsType(iobs).eq.isTvar(itrc)) THEN
                i=isTvar(itrc)
                FOURDVAR(ng)%ObsCount(i)=FOURDVAR(ng)%ObsCount(i)+1
              END IF
            END DO
# endif
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  If applicable, set next observation survey time to process.
!-----------------------------------------------------------------------
!
        IF (Master) THEN
          WRITE (stdout,10) ObsTime(ng)*sec2day
        END IF
        IF (backward) THEN
          IF ((ObsSurvey(ng)-1).ge.1) THEN
            ObsTime(ng)=FOURDVAR(ng)%SurveyTime(ObsSurvey(ng)-1)*day2sec
          END IF
        ELSE
          IF ((ObsSurvey(ng)+1).le.Nsurvey(ng)) THEN
            ObsTime(ng)=FOURDVAR(ng)%SurveyTime(ObsSurvey(ng)+1)*day2sec
          END IF
        END IF
      END IF
!
  10  FORMAT (/,' Number of State Observations Processed:',             &
     &        t58,'ObsTime = ',f12.4,/,/,                               &
     &        10x,'Variable',10x,'IstrObs',4x,'IendObs',6x,'Count',     &
     &        3x,'Rejected',/)

      RETURN
      END SUBROUTINE obs_read
#else
      SUBROUTINE obs_read
      RETURN
      END SUBROUTINE obs_read
#endif
