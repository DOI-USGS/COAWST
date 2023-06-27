#include "cppdefs.h"
      MODULE obs_initial_mod
#if (defined FOUR_DVAR || defined VERIFICATION) && defined OBSERVATIONS
!
!git $Id$
!svn $Id: obs_initial.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module opens and inquires observations file using the standard !
!  NetCDF library or the Parallel-IO (PIO) library.                    !
!                                                                      !
!  It sets various variables needed for the processing of the state    !
!  solution at observations locations in variational data assimilation.!
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE strings_mod, ONLY : FoundError
!
      implicit none
!
      PUBLIC  :: obs_initial
      PRIVATE :: obs_initial_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: obs_initial_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE obs_initial (ng, model, backward)
!***********************************************************************
!
!  Imported variable declarations.
!
      logical, intent(in) :: backward
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Write out time-averaged fields according to IO type.
!-----------------------------------------------------------------------
!
      SELECT CASE (OBS(ng)%IOtype)
        CASE (io_nf90)
          CALL obs_initial_nf90 (ng, model, backward)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL obs_initial_pio (ng, model, backward)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) OBS(ng)%IOtype
          exit_flag=2
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' OBS_INITIAL - Illegal input file type, io_type = ',i0,  &
     &        /,15x,'Check KeyWord ''INP_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE obs_initial
!
!***********************************************************************
      SUBROUTINE obs_initial_nf90 (ng, model, backward)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: backward
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      logical, dimension(NV) :: got_var(NV)
!
      integer :: Ifirst, i, nvd, recdim, status
      integer :: Vsize(4)
!
# ifdef WEAK_CONSTRAINT
      real(r8), parameter :: IniVal = 0.0_r8
# endif
      real(r8) :: tend
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//"obs_initial_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Inquire about the contents of observation NetCDF file:  Inquire about
!  the dimensions and variables.
!-----------------------------------------------------------------------
!
      QUERY : IF (OBS(ng)%ncid.eq.-1) THEN
!
!  Open observations NetCDF file.
!
        CALL netcdf_open (ng, model, OBS(ng)%name, 1, OBS(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,10) TRIM(OBS(ng)%name)
          RETURN
        END IF
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, model, OBS(ng)%name,                   &
     &                       OBS(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from observation NetCDF and activate switches for
!  required variables.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOday))) THEN
            got_var(idOday)=.TRUE.
            OBS(ng)%Vid(idOday)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idNobs))) THEN
            got_var(idNobs)=.TRUE.
            OBS(ng)%Vid(idNobs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOtyp))) THEN
            got_var(idOtyp)=.TRUE.
            OBS(ng)%Vid(idOtyp)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idObsT))) THEN
            got_var(idObsT)=.TRUE.
            OBS(ng)%Vid(idObsT)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idObsD))) THEN
            got_var(idObsD)=.TRUE.
            OBS(ng)%Vid(idObsD)=var_id(i)
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idObsX))) THEN
            got_var(idObsX)=.TRUE.
            OBS(ng)%Vid(idObsX)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idObsY))) THEN
            got_var(idObsY)=.TRUE.
            OBS(ng)%Vid(idObsY)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idObsZ))) THEN
            got_var(idObsZ)=.TRUE.
            OBS(ng)%Vid(idObsZ)=var_id(i)
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOerr))) THEN
            got_var(idOerr)=.TRUE.
            OBS(ng)%Vid(idOerr)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOval))) THEN
            got_var(idOval)=.TRUE.
            OBS(ng)%Vid(idOval)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOmet))) THEN
            got_var(idOmet)=.TRUE.
            haveObsMeta(ng)=.TRUE.
            OBS(ng)%Vid(idOmet)=var_id(i)
          END IF
        END DO
!
!  Check if needed obsrvation variables are available.
!
        IF (.not.got_var(idOday)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idOday)),          &
     &                                  TRIM(OBS(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idNobs)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idNobs)),          &
     &                                  TRIM(OBS(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idOtyp)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idOtyp)),          &
     &                                  TRIM(OBS(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idObsT)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idObsT)),          &
     &                                  TRIM(OBS(ng)%name)
          exit_flag=2
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idObsD)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idObsD)),          &
     &                                  TRIM(OBS(ng)%name)
          exit_flag=2
          RETURN
        END IF
# endif
        IF (.not.got_var(idObsX)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idObsX)),          &
     &                                  TRIM(OBS(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idObsY)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idObsY)),          &
     &                                  TRIM(OBS(ng)%name)
          exit_flag=2
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idObsZ)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idObsZ)),          &
     &                                  TRIM(OBS(ng)%name)
          exit_flag=2
          RETURN
        END IF
# endif
        IF (.not.got_var(idOerr)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idOerr)),          &
     &                                  TRIM(OBS(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idOval)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idOval)),          &
     &                                  TRIM(OBS(ng)%name)
          exit_flag=2
          RETURN
        END IF
      END IF QUERY
!
!-----------------------------------------------------------------------
!  Set observation processing variables.
!-----------------------------------------------------------------------
!
!  Determine if there is any data available in the model time
!  window.  Set first survey record to process.
!
      IF (backward) THEN
        Ifirst=0
# ifdef GENERIC_DSTART
        tend=(time(ng)-ntimes(ng)*dt(ng))*sec2day
# else
        tend=(time(ng)-(ntstart(ng)-1)*dt(ng))*sec2day
# endif
        DO i=1,Nsurvey(ng)
          IF ((tend.le.FOURDVAR(ng)%SurveyTime(i)).and.                 &
     &        (FOURDVAR(ng)%SurveyTime(i).le.tdays(ng))) THEN
            Ifirst=MAX(Ifirst,i)
          END IF
        END DO
# ifndef SP4DVAR
        IF (Ifirst.eq.0) THEN
          IF (Master) WRITE (stdout,30) tend, tdays(ng)
          exit_flag=2
          RETURN
        END IF
# endif
      ELSE
        Ifirst=Nsurvey(ng)
        tend=(time(ng)+ntimes(ng)*dt(ng))*sec2day
        DO i=1,Nsurvey(ng)
          IF ((tdays(ng).le.FOURDVAR(ng)%SurveyTime(i)).and.            &
     &        (FOURDVAR(ng)%SurveyTime(i).le.tend)) THEN
            Ifirst=MIN(Ifirst,i)
          END IF
        END DO
# ifndef SP4DVAR
        IF (Ifirst.eq.0) THEN
          IF (Master) WRITE (stdout,30) tdays(ng), tend
          exit_flag=2
          RETURN
        END IF
# endif
      END IF
      ObsTime(ng)=FOURDVAR(ng)%SurveyTime(Ifirst)*day2sec
!
!  Initialize observation survey counter.  This is the counter of data
!  assimilation cycles or observations survey times on which the model
!  state is extracted (interpolated) at the observation locations.
!
      IF (backward) THEN
        ObsSurvey(ng)=Ifirst+1
      ELSE
        ObsSurvey(ng)=Ifirst-1
      END IF
!
!  Initialize time switch to process model state at observation
!  locations.
!
      ProcessObs(ng)=.FALSE.

# ifdef I4DVAR
!
!  Initialize cost function misfit between model and observations.
!  The IF statement is to avoid rewritting its value before it is
!  written into the initial NetCDF file.
!
      IF (.not.backward) THEN
        DO i=0,NobsVar(ng)
          FOURDVAR(ng)%ObsCost(i)=0.0_r8
        END DO
      END IF
# endif
!
!  Set staring and ending observation indices.
!
      IF (backward) THEN
        NstrObs(ng)=0
        NendObs(ng)=0
        DO i=1,Ifirst
          NstrObs(ng)=NstrObs(ng)+FOURDVAR(ng)%NobsSurvey(i)
        END DO
        NstrObs(ng)=NstrObs(ng)+1
      ELSE
        IF (Ifirst.eq.1) THEN
          NstrObs(ng)=0
          NendObs(ng)=0
        ELSE
          NstrObs(ng)=0
          NendObs(ng)=0
          DO i=1,Ifirst-1
            NendObs(ng)=NendObs(ng)+FOURDVAR(ng)%NobsSurvey(i)
          END DO
        END IF
      END IF
!
!  Initialize total observation counter in structure array.  Notice
!  that the zero index carries the total summation.
!
      FOURDVAR(ng)%ObsCount(0)=0
      FOURDVAR(ng)%ObsReject(0)=0

# ifdef WEAK_CONSTRAINT
!
!  Initialize model values at observation locations.
!
      DO i=1,Mobs
        NLmodVal(i)=IniVal
        ObsVetting(i)=IniVal
#  ifndef SP4DVAR
        TLmodVal(i)=IniVal
#  endif
#  if defined SP4DVAR           && defined DISJOINTED && \
      defined CONCURRENT_KERNEL
        TLmodVal(i)=IniVal
#  endif
      END DO
#  if defined SP4DVAR           && defined DISJOINTED && \
     !defined CONCURRENT_KERNEL
      IF (model.ne.iADM) THEN
        DO i=1,Mobs
          TLmodVal(i)=IniVal
        END DO
      END IF
#  endif
# endif
!
  10  FORMAT (/,' OBS_INITIAL_NF90 - unable to open input NetCDF',      &
     &        ' file: ',a)
  20  FORMAT (/,' OBS_INITIAL_NF90 - unable to find model variable:',   &
     &        1x,a,/,20x,'in input NetCDF file: ',a)
  30  FORMAT (/,' OBS_INITIAL_NF90 - No are observations available',    &
     &        ' for time window (days): ',/,12x,f12.4,' - ',f12.4,/)
!
      RETURN
      END SUBROUTINE obs_initial_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE obs_initial_pio (ng, model, backward)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: backward
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      logical, dimension(NV) :: got_var(NV)
!
      integer :: Ifirst, i, nvd, recdim, status
      integer :: Vsize(4)
!
#  ifdef WEAK_CONSTRAINT
      real(r8), parameter :: IniVal = 0.0_r8
#  endif
      real(r8) :: tend
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//"obs_initial_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Inquire about the contents of observation NetCDF file:  Inquire about
!  the dimensions and variables.
!-----------------------------------------------------------------------
!
      QUERY : IF (OBS(ng)%ncid.eq.-1) THEN
!
!  Open observations NetCDF file.
!
        CALL pio_netcdf_open (ng, model, OBS(ng)%name, 1,               &
     &                        OBS(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,10) TRIM(OBS(ng)%name)
          RETURN
        END IF
!
!  Inquire about the variables.
!
        CALL pio_netcdf_inq_var (ng, model, OBS(ng)%name,               &
     &                           OBS(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from observation NetCDF and activate switches for
!  required variables.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOday))) THEN
            got_var(idOday)=.TRUE.
            OBS(ng)%pioVar(idOday)%vd=var_desc(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idNobs))) THEN
            got_var(idNobs)=.TRUE.
            OBS(ng)%pioVar(idNobs)%vd=var_desc(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOtyp))) THEN
            got_var(idOtyp)=.TRUE.
            OBS(ng)%pioVar(idOtyp)%vd=var_desc(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idObsT))) THEN
            got_var(idObsT)=.TRUE.
            OBS(ng)%pioVar(idObsT)%vd=var_desc(i)
#  ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idObsD))) THEN
            got_var(idObsD)=.TRUE.
            OBS(ng)%pioVar(idObsD)%vd=var_desc(i)
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idObsX))) THEN
            got_var(idObsX)=.TRUE.
            OBS(ng)%pioVar(idObsX)%vd=var_desc(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idObsY))) THEN
            got_var(idObsY)=.TRUE.
            OBS(ng)%pioVar(idObsY)%vd=var_desc(i)
#  ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idObsZ))) THEN
            got_var(idObsZ)=.TRUE.
            OBS(ng)%pioVar(idObsZ)%vd=var_desc(i)
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOerr))) THEN
            got_var(idOerr)=.TRUE.
            OBS(ng)%pioVar(idOerr)%vd=var_desc(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOval))) THEN
            got_var(idOval)=.TRUE.
            OBS(ng)%pioVar(idOval)%vd=var_desc(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOmet))) THEN
            got_var(idOmet)=.TRUE.
            haveObsMeta(ng)=.TRUE.
            OBS(ng)%pioVar(idOmet)%vd=var_desc(i)
          END IF
        END DO
!
!  Check if needed obsrvation variables are available.
!
        IF (.not.got_var(idOday)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idOday)),          &
     &                                  TRIM(OBS(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idNobs)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idNobs)),          &
     &                                  TRIM(OBS(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idOtyp)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idOtyp)),          &
     &                                  TRIM(OBS(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idObsT)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idObsT)),          &
     &                                  TRIM(OBS(ng)%name)
          exit_flag=2
          RETURN
        END IF
#  ifdef SOLVE3D
        IF (.not.got_var(idObsD)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idObsD)),          &
     &                                  TRIM(OBS(ng)%name)
          exit_flag=2
          RETURN
        END IF
#  endif
        IF (.not.got_var(idObsX)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idObsX)),          &
     &                                  TRIM(OBS(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idObsY)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idObsY)),          &
     &                                  TRIM(OBS(ng)%name)
          exit_flag=2
          RETURN
        END IF
#  ifdef SOLVE3D
        IF (.not.got_var(idObsZ)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idObsZ)),          &
     &                                  TRIM(OBS(ng)%name)
          exit_flag=2
          RETURN
        END IF
#  endif
        IF (.not.got_var(idOerr)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idOerr)),          &
     &                                  TRIM(OBS(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idOval)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idOval)),          &
     &                                  TRIM(OBS(ng)%name)
          exit_flag=2
          RETURN
        END IF
      END IF QUERY
!
!-----------------------------------------------------------------------
!  Set observation processing variables.
!-----------------------------------------------------------------------
!
!  Determine if there is any data available in the model time
!  window.  Set first survey record to process.
!
      IF (backward) THEN
        Ifirst=0
#  ifdef GENERIC_DSTART
        tend=(time(ng)-ntimes(ng)*dt(ng))*sec2day
#  else
        tend=(time(ng)-(ntstart(ng)-1)*dt(ng))*sec2day
#  endif
        DO i=1,Nsurvey(ng)
          IF ((tend.le.FOURDVAR(ng)%SurveyTime(i)).and.                 &
     &        (FOURDVAR(ng)%SurveyTime(i).le.tdays(ng))) THEN
            Ifirst=MAX(Ifirst,i)
          END IF
        END DO
#  ifndef SP4DVAR
        IF (Ifirst.eq.0) THEN
          IF (Master) WRITE (stdout,30) tend, tdays(ng)
          exit_flag=2
          RETURN
        END IF
#  endif
      ELSE
        Ifirst=Nsurvey(ng)
        tend=(time(ng)+ntimes(ng)*dt(ng))*sec2day
        DO i=1,Nsurvey(ng)
          IF ((tdays(ng).le.FOURDVAR(ng)%SurveyTime(i)).and.            &
     &        (FOURDVAR(ng)%SurveyTime(i).le.tend)) THEN
            Ifirst=MIN(Ifirst,i)
          END IF
        END DO
#  ifndef SP4DVAR
        IF (Ifirst.eq.0) THEN
          IF (Master) WRITE (stdout,30) tdays(ng), tend
          exit_flag=2
          RETURN
        END IF
#  endif
      END IF
      ObsTime(ng)=FOURDVAR(ng)%SurveyTime(Ifirst)*day2sec
!
!  Initialize observation survey counter.  This is the counter of data
!  assimilation cycles or observations survey times on which the model
!  state is extracted (interpolated) at the observation locations.
!
      IF (backward) THEN
        ObsSurvey(ng)=Ifirst+1
      ELSE
        ObsSurvey(ng)=Ifirst-1
      END IF
!
!  Initialize time switch to process model state at observation
!  locations.
!
      ProcessObs(ng)=.FALSE.

#  ifdef I4DVAR
!
!  Initialize cost function misfit between model and observations.
!  The IF statement is to avoid rewritting its value before it is
!  written into the initial NetCDF file.
!
      IF (.not.backward) THEN
        DO i=0,NobsVar(ng)
          FOURDVAR(ng)%ObsCost(i)=0.0_r8
        END DO
      END IF
#  endif
!
!  Set staring and ending observation indices.
!
      IF (backward) THEN
        NstrObs(ng)=0
        NendObs(ng)=0
        DO i=1,Ifirst
          NstrObs(ng)=NstrObs(ng)+FOURDVAR(ng)%NobsSurvey(i)
        END DO
        NstrObs(ng)=NstrObs(ng)+1
      ELSE
        IF (Ifirst.eq.1) THEN
          NstrObs(ng)=0
          NendObs(ng)=0
        ELSE
          NstrObs(ng)=0
          NendObs(ng)=0
          DO i=1,Ifirst-1
            NendObs(ng)=NendObs(ng)+FOURDVAR(ng)%NobsSurvey(i)
          END DO
        END IF
      END IF
!
!  Initialize total observation counter in structure array.  Notice
!  that the zero index carries the total summation.
!
      FOURDVAR(ng)%ObsCount(0)=0
      FOURDVAR(ng)%ObsReject(0)=0

#  ifdef WEAK_CONSTRAINT
!
!  Initialize model values at observation locations.
!
      DO i=1,Mobs
        NLmodVal(i)=IniVal
        ObsVetting(i)=IniVal
#   ifndef SP4DVAR
        TLmodVal(i)=IniVal
#   endif
#   if defined SP4DVAR           && defined DISJOINTED && \
      defined CONCURRENT_KERNEL
        TLmodVal(i)=IniVal
#   endif
      END DO
#   if defined SP4DVAR           && defined DISJOINTED && \
     !defined CONCURRENT_KERNEL
      IF (model.ne.iADM) THEN
        DO i=1,Mobs
          TLmodVal(i)=IniVal
        END DO
      END IF
#   endif
#  endif
!
  10  FORMAT (/,' OBS_INITIAL_PIO - unable to open input NetCDF',       &
     &        ' file: ',a)
  20  FORMAT (/,' OBS_INITIAL_PIO - unable to find model variable:',    &
     &        1x,a,/,20x,'in input NetCDF file: ',a)
  30  FORMAT (/,' OBS_INITIAL_PIO - No are observations available',     &
     &        ' for time window (days): ',/,12x,f12.4,' - ',f12.4,/)
!
      RETURN
      END SUBROUTINE obs_initial_pio
# endif
#endif
      END MODULE obs_initial_mod
