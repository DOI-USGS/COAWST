#include "cppdefs.h"
#if (defined FOUR_DVAR || defined VERIFICATION) && defined OBSERVATIONS
      SUBROUTINE obs_read (ng, model, backward)
!
!svn $Id: obs_read.F 584 2008-03-18 20:01:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
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
# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti, mp_bcastl
# endif
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

      integer, dimension(2) :: start, total
!
!---------------------------------------------------------------------
!  Read observation variables needed for interpolating the model
!  state at the observation locations.
!---------------------------------------------------------------------
!
      IF (ProcessObs(ng)) THEN
# if defined S4DVAR
        readNLmod=backward
# elif defined TLM_OBS
        readNLmod=.TRUE.
        readTLmod=.TRUE.
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
        start(1)=NstrObs(ng)
        total(1)=Nobs(ng)
!
!  Set starting index of obervation vectors for reading.  In weak
!  constraint, the entire observation data is loaded. Otherwise,
!  only the observartion for the current time window are loaded
!  and started from vector index one.
!
# if defined WEAK_CONSTRAINT || defined IOM
        Mstr=NstrObs(ng)
        Mend=NendObs(ng)
# else
        Mstr=1
        Mend=Nobs(ng)
# endif
!
!  Read in observation type identifier.
!
        IF (InpThread) THEN
          status=nf90_get_var(ncOBSid(ng), obsVid(idOtyp,ng),           &
     &                        ObsType(Mstr:), start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOtyp)), ObsSurvey(ng),     &
     &                        TRIM(OBSname(ng))
            exit_flag=2
            ioerror=status
            RETURN
          END IF
        END IF
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, model, ObsType, Mobs)
# endif
!
!  Read in observation time (days).
!
        IF (InpThread) THEN
          status=nf90_get_var(ncOBSid(ng), obsVid(idObsT,ng),           &
     &                        Tobs(Mstr:), start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idObsT)), ObsSurvey(ng),     &
     &                        TRIM(OBSname(ng))
            exit_flag=2
            ioerror=status
            RETURN
          END IF
        END IF
# ifdef DISTRIBUTE
        CALL mp_bcastf (ng, model, Tobs, Mobs)
# endif
!
!  Read in observation X-location (grid units).
!
        IF (InpThread) THEN
          status=nf90_get_var(ncOBSid(ng), obsVid(idObsX,ng),           &
     &                        Xobs(Mstr:), start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idObsX)), ObsSurvey(ng),     &
     &                        TRIM(OBSname(ng))
            exit_flag=2
            ioerror=status
            RETURN
          END IF
        END IF
# ifdef DISTRIBUTE
        CALL mp_bcastf (ng, model, Xobs, Mobs)
# endif
!
!  Read in observation Y-location (grid units).
!
        IF (InpThread) THEN
          status=nf90_get_var(ncOBSid(ng), obsVid(idObsY,ng),           &
     &                        Yobs(Mstr:), start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idObsY)), ObsSurvey(ng),     &
     &                        TRIM(OBSname(ng))
            exit_flag=2
            ioerror=status
            RETURN
          END IF
        END IF
# ifdef DISTRIBUTE
        CALL mp_bcastf (ng, model, Yobs, Mobs)
# endif

# ifdef SOLVE3D
!
!  Read in observation depth, Zobs.  If negative, depth is meter. If
!  greater than zero, Zobs is in model fractional vertical levels
!  (1 <= Zobs <= N). If Zobs < 0, its fractional level value is
!  computed in routine "extract_obs3d" and over-written so it can
!  be written into the observation NetCDF file for latter use.
!
        IF (InpThread) THEN
          IF (wrote_Zobs(ng)) THEN
            status=nf90_get_var(ncOBSid(ng), obsVid(idObsZ,ng),         &
     &                          Zobs(Mstr:), start, total)
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,10) TRIM(Vname(1,idObsZ)), ObsSurvey(ng),   &
     &                          TRIM(OBSname(ng))
              exit_flag=2
              ioerror=status
              RETURN
            END IF
          ELSE
            status=nf90_get_var(ncOBSid(ng), obsVid(idObsD,ng),         &
     &                          Zobs(Mstr:), start, total)
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,10) TRIM(Vname(1,idObsD)), ObsSurvey(ng),   &
     &                          TRIM(OBSname(ng))
              exit_flag=2
              ioerror=status
              RETURN
            END IF
          END IF
          Load_Zobs(ng)=.FALSE.
          IF ((MINVAL(Zobs).lt.0.0_r8).or.                              &
     &        (MAXVAL(Zobs).lt.0.0_r8)) THEN
            Load_Zobs(ng)=.TRUE.
          END IF
        END IF
#  ifdef DISTRIBUTE
        CALL mp_bcastl (ng, model, Load_Zobs(ng), 1)
        CALL mp_bcastf (ng, model, Zobs, Mobs)
!
!  If distributed-memory and Zobs in meters (Zobs < 0),  zero-out
!  Zobs values in all nodes by itself to facilitate exchages between
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
        IF (InpThread) THEN
          status=nf90_get_var(ncOBSid(ng), obsVid(idOval,ng),           &
     &                        ObsVal(Mstr:), start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOval)), ObsSurvey(ng),     &
     &                        TRIM(OBSname(ng))
            exit_flag=2
            ioerror=status
            RETURN
          END IF
        END IF
# ifdef DISTRIBUTE
        CALL mp_bcastf (ng, model, ObsVal, Mobs)
# endif
!
# if defined WEAK_CONSTRAINT || defined IOM
!  Read in observation error covariance.
# else
!  Read in observation error covariance. To avoid successive divisions,
!  convert to inverse observation error covariance.
# endif
!
        IF (InpThread) THEN
          status=nf90_get_var(ncOBSid(ng), obsVid(idOerr,ng),           &
     &                        ObsErr(Mstr:), start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOerr)), ObsSurvey(ng),     &
     &                        TRIM(OBSname(ng))
            exit_flag=2
            ioerror=status
            RETURN
          END IF
# if !(defined WEAK_CONSTRAINT || defined IOM)
          DO iobs=1,Nobs(ng)
            ObsErr(iobs)=1.0_r8/ObsErr(iobs)
          END DO
# endif
        END IF
# ifdef DISTRIBUTE
        CALL mp_bcastf (ng, model, ObsErr, Mobs)
# endif
!
!  Read in nonlinear model values at observation locations.
!
        IF (readNLmod.and.haveNLmod(ng)) THEN
          IF (InpThread) THEN
# ifdef S4DVAR
            start(1)=NstrObs(ng)
            total(1)=Nobs(ng)
            start(2)=Ipass
            total(2)=1
            status=nf90_get_var(ncMODid(ng), modVid(idNLmo,ng),         &
     &                          NLmodVal(:,Ipass), start, total)
# else
            start(1)=NstrObs(ng)
            total(1)=Nobs(ng)
            status=nf90_get_var(ncMODid(ng), modVid(idNLmo,ng),         &
     &                          NLmodVal(Mstr:), start, total)
# endif
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,10) TRIM(Vname(1,idNLmo)), ObsSurvey(ng),   &
     &                          TRIM(MODname(ng))
              exit_flag=2
              ioerror=status
              RETURN
            END IF
          END IF
# ifdef DISTRIBUTE
#  ifdef S4DVAR
          CALL mp_bcastf (ng, model, NLmodVal(:,Ipass), Mobs)
#  else
          CALL mp_bcastf (ng, model, NLmodVal, Mobs)
#  endif
# endif
        END IF

# if defined TLM_OBS && !(defined WEAK_CONSTRAINT || defined IOM)
!
!  If adjoint pass and incremental S4DVAR, read in tangent linear model
!  values at observation locations.
!
#  ifdef IS4DVAR_OLD
        IF (readTLmod.and.haveTLmod(ng).and.(Ipass.ne.0)) THEN
#  else
        IF (readTLmod.and.haveTLmod(ng)) THEN
#  endif
          IF (InpThread) THEN
#  ifdef IS4DVAR_OLD
            start(1)=NstrObs(ng)
            total(1)=Nobs(ng)
            start(2)=Ipass
            total(2)=1
            status=nf90_get_var(ncMODid(ng), modVid(idTLmo,ng),         &
     &                          TLmodVal(:,Ipass), start, total)
#  else
            start(1)=NstrObs(ng)
            total(1)=Nobs(ng)
            status=nf90_get_var(ncMODid(ng), modVid(idTLmo,ng),         &
     &                          TLmodVal(Mstr:), start, total)
#  endif
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTLmo)), ObsSurvey(ng),   &
     &                          TRIM(MODname(ng))
              exit_flag=2
              ioerror=status
              RETURN
            END IF
          END IF
!
!  Reset TLM values to zero at the first pass of the inner loop.
!
#  if defined IS4DVAR_OLD
          IF ((inner.eq.1).and.(Ipass.eq.1)) THEN
            DO iobs=1,Mobs
              TLmodVal(iobs,1)=0.0_r8
              TLmodVal(iobs,2)=0.0_r8
            END DO
          END IF
#  elif defined IS4DVAR
          IF (inner.eq.0) THEN
            DO iobs=1,Mobs
              TLmodVal(iobs)=0.0_r8
            END DO
          END IF
#  endif
#  ifdef DISTRIBUTE
#   ifdef IS4DVAR_OLD
          CALL mp_bcastf (ng, model, TLmodVal(:,1), Mobs)
          CALL mp_bcastf (ng, model, TLmodVal(:,2), Mobs)
#   else
          CALL mp_bcastf (ng, model, TLmodVal, Mobs)
#   endif
#  endif
        END IF
# endif
# if defined IOM && defined ADJOINT
!
!  If multiple executables, read in representer coefficients (or
!  its approximation) and load it in ADmodVal.
!
        IF (backward) THEN
          IF (InpThread) THEN
            start(1)=1
            total(1)=Ndatum(ng)
            status=nf90_get_var(ncMODid(ng), modVid(idRepC,ng),         &
     &                          ADmodVal, start, total)
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,10) TRIM(Vname(1,idRepC)), ObsSurvey(ng),   &
     &                          TRIM(MODname(ng))
              exit_flag=2
              ioerror=status
              RETURN
            END IF
          END IF
#  ifdef DISTRIBUTE
          CALL mp_bcastf (ng, model, ADmodVal, Mobs)
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
          WRITE (stdout,20) ObsTime(ng)*sec2day
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
  10  FORMAT (/,' OBS_READ - error while reading variable: ',a,         &
     &        ' at survey record = ',i4,/,                              &
     &        11x,'in input NetCDF file: ',a)
  20  FORMAT (/,' Number of State Observations Processed:',             &
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
