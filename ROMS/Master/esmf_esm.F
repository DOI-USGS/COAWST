#include "cppdefs.h"
      MODULE esmf_esm_mod

#if defined MODEL_COUPLING && defined ESMF_LIB
!
!git $Id$
!svn $Id: esmf_esm.F 1151 2023-02-09 03:08:53Z arango $
!=======================================================================
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license         Hernan G. Arango     !
!    See License_ROMS.txt                         Ufuk Utku Turuncoglu !
!=======================================================================
!                                                                      !
!  This module sets the Earth System Model (ESM) coupled gridded       !
!  components interface routines using the ESMF/NUOPC layer:           !
!                                                                      !
!    ESM_SetServices        Sets components shared-object entry        !
!                           points using NUOPC generic methods         !
!                                                                      !
!    ESM_SetModelServices   Sets shared-object entry points for each   !
!                           coupled component, its connectors, and     !
!                           coupled system internal clock              !
!                                                                      !
!    ESM_SetRunSequence     Sets connectors (explicit or implicit)     !
!                           dependencies between ESM components        !
!                                                                      !
!  The coupled components are attached to the driver via connectors.   !
!                                                                      !
!  ESMF:   Earth System Modeling Framework (Version 7 or higher)       !
!            https://www.earthsystemcog.org/projects/esmf              !
!                                                                      !
!  NUOPC:  National Unified Operational Prediction Capability          !
!            https://www.earthsystemcog.org/projects/nuopc             !
!                                                                      !
!=======================================================================
!
      USE ESMF
      USE NUOPC
      USE NUOPC_Driver,                                                 &
     &    NUOPC_SetServices            => SetServices,                  &
     &    NUOPC_Label_SetModelServices => label_SetModelServices,       &
     &    NUOPC_Label_SetRunSequence   => label_SetRunSequence
!
      USE mod_esmf_esm          ! ESM coupling structures and variables
!
# ifdef ATM_COUPLING
#  if defined COAMPS_COUPLING
      USE esmf_coamps_mod,  ONLY: ATM_SetServices
#  elif defined REGCM_COUPLING
      USE esmf_regcm_mod,   ONLY: ATM_SetServices
#  elif defined WRF_COUPLING
      USE esmf_wrf_mod,     ONLY: ATM_SetServices
#  else
      USE esmf_atm_mod,     ONLY: ATM_SetServices
#  endif
# endif
      USE esmf_coupler_mod, ONLY: Coupler_SetServices
# ifdef DATA_COUPLING
      USE esmf_data_mod,    ONLY: DATA_SetServices
# endif
# ifdef ICE_COUPLING
#  if defined CICE_COUPLING
      USE esmf_cice_mod,    ONLY: ICE_SetServices
#  else
      USE esmf_ice_mod,     ONLY: ICE_SetServices
#  endif
# endif
      USE esmf_roms_mod,    ONLY: ROMS_SetServices
# ifdef WAV_COUPLING
#  if defined WAM_COUPLING
      USE esmf_wam_mod,     ONLY: WAV_SetServices
#  else
      USE esmf_wav_mod,     ONLY: WAV_SetServices
#  endif
# endif
!
      implicit none
!
      PUBLIC  :: ESM_SetServices
      PRIVATE :: ESM_SetModelServices
      PRIVATE :: ESM_SetRunSequence
!
      CONTAINS
!
      SUBROUTINE ESM_SetServices (driver, rc)
!
!=======================================================================
!                                                                      !
!  Sets the gridded component shared-object entry points using NUOPC   !
!  generic methods followed by the appropriate specializations for     !
!  "component services" and "run sequence".                            !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      TYPE(ESMF_GridComp) :: driver
      TYPE(ESMF_Config)   :: config
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ESM_SetServices"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering ESM_SetServices',         &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Register generic methods.
!-----------------------------------------------------------------------
!
      CALL NUOPC_CompDerive (driver,                                    &
     &                       NUOPC_SetServices,                         &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Attach ESM specializing methods.
!-----------------------------------------------------------------------
!
!  Set gridded component services shared-object entry point.
!
      CALL NUOPC_CompSpecialize (driver,                                &
     &                           specLabel=NUOPC_Label_SetModelServices,&
     &                           specRoutine=ESM_SetModelServices,      &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set component run sequence shared-object entry point.
!
      CALL NUOPC_CompSpecialize (driver,                                &
     &                           specLabel=NUOPC_Label_SetRunSequence,  &
     &                           specRoutine=ESM_SetRunSequence,        &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Create, open, and set ESM configuration.  The ESM application run
!  sequence is read from the input configuration file.
!-----------------------------------------------------------------------
!
!  Create configuration object.
!
      config = ESMF_ConfigCreate(rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Read in ESN application configuration file.
!
      CALL ESMF_ConfigLoadFile(config,                                  &
     &                         TRIM(CONFname),                          &
     &                         rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set configuration object in ESM driver.
!
      CALL ESMF_GridCompSet(driver,                                     &
     &                      config=config,                              &
     &                      rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  ESM_SetServices',         &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      END SUBROUTINE ESM_SetServices
!
      SUBROUTINE ESM_SetModelServices (driver, rc)
!
!=======================================================================
!                                                                      !
!  Sets shared-object entry point for each active coupled component    !
!  services. Then, set the "connectors" between active coupled         !
!  component and internal clocks.                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      TYPE (ESMF_GridComp) :: driver
!
!  Local variable declarations.
!
      integer :: i, j
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ESM_SetModelServices"
!
      TYPE (ESMF_GridComp) :: model
      TYPE (ESMF_Clock)    :: clock
      TYPE (ESMF_CplComp)  :: connector
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering ESM_SetModelServices',    &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Set services for ESM active coupled components.
!-----------------------------------------------------------------------
!
      DO i=1,Nmodels
        IF (MODELS(i)%IsActive) THEN
          SELECT CASE (TRIM(Clabel(i)))
            CASE ('OCN')
              CALL NUOPC_DriverAddComp (driver,                         &
     &                                  TRIM(MODELS(i)%name),           &
     &                                  ROMS_SetServices,               &
     &                                  petList=MODELS(i)%PETlist(:),   &
     &                                  comp=model,                     &
     &                                  rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
# ifdef ATM_COUPLING
            CASE ('ATM')
              CALL NUOPC_DriverAddComp (driver,                         &
     &                                  TRIM(MODELS(i)%name),           &
     &                                  ATM_SetServices,                &
     &                                  petList=MODELS(i)%PETlist(:),   &
     &                                  comp=model,                     &
     &                                  rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
# endif
# ifdef ICE_COUPLING
            CASE ('ICE')
              CALL NUOPC_DriverAddComp (driver,                         &
     &                                  TRIM(MODELS(i)%name),           &
     &                                  ICE_SetServices,                &
     &                                  petList=MODELS(i)%PETlist(:),   &
     &                                  comp=model,                     &
     &                                  rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
# endif
# ifdef WAV_COUPLING
            CASE ('WAV')
              CALL NUOPC_DriverAddComp (driver,                         &
     &                                  TRIM(MODELS(i)%name),           &
     &                                  WAV_SetServices,                &
     &                                  petList=MODELS(i)%PETlist(:),   &
     &                                  comp=model,                     &
     &                                  rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
# endif
# ifdef DATA_COUPLING
            CASE ('DAT')
              CALL NUOPC_DriverAddComp (driver,                         &
     &                                  TRIM(MODELS(i)%name),           &
     &                                  DATA_SetServices,               &
     &                                  petList=MODELS(i)%PETlist(:),   &
     &                                  comp=model,                     &
     &                                  rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
# endif
          END SELECT
!
!  Set debugging flag.
!
          IF (DebugLevel.gt.0) THEN
            CALL ESMF_AttributeSet (model,                              &
     &                              name="Verbosity",                   &
     &                              value="high",                       &
     &                              rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Set services for ESM active connector components.
!-----------------------------------------------------------------------
!
      DO i=1,Nmodels
        DO j=1,Nmodels
          IF (CONNECTORS(i,j)%IsActive) THEN
            CALL NUOPC_DriverAddComp (driver,                           &
     &                      srcCompLabel=TRIM(MODELS(i)%name),          &
     &                      dstCompLabel=TRIM(MODELS(j)%name),          &
     &                      compSetServicesRoutine=Coupler_SetServices, &
     &                                comp=connector,                   &
     &                                rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
            IF (DebugLevel.gt.0) THEN
              CALL ESMF_AttributeSet (connector,                        &
     &                                name="Verbosity",                 &
     &                                value="high",                     &
     &                                rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
            END IF
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Set internal clock for application: coupling starting and stopping
!  times.
!-----------------------------------------------------------------------
!
      IF (ClockInfo(Idriver)%StartTime.ne.                              &
     &    ClockInfo(Idriver)%RestartTime) THEN
        ClockInfo(Idriver)%Restarted=.TRUE.
        Clock=ESMF_ClockCreate(ClockInfo(Idriver)%TimeStep,             &
     &                         ClockInfo(Idriver)%RestartTime,          &
     &                         stopTime=ClockInfo(Idriver)%StopTime,    &
!!   &                         refTime=ClockInfo(Idriver)%ReferenceTime,&
     &                         name='ESM_clock',                        &
     &                         rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
      ELSE
        ClockInfo(Idriver)%Restarted=.FALSE.
        Clock=ESMF_ClockCreate(ClockInfo(Idriver)%TimeStep,             &
     &                         ClockInfo(Idriver)%StartTime,            &
     &                         stopTime=ClockInfo(Idriver)%StopTime,    &
!!   &                         refTime=ClockInfo(Idriver)%ReferenceTime,&
     &                         name='ESM_clock',                        &
     &                         rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
      END IF
      ClockInfo(Idriver)%Clock=clock
!
      CALL ESMF_GridCompSet (driver,                                    &
     &                       clock=ClockInfo(Idriver)%Clock,            &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  ESM_SetModelServices',    &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE ESM_SetModelServices
!
      SUBROUTINE ESM_SetRunSequence (driver, rc)
!
!=======================================================================
!                                                                      !
!  Sets the gridded component shared-object entry points using NUOPC   !
!  generic methods followed by the appropriate specializations for     !
!  "component services" and "run sequence".                            !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      TYPE (ESMF_GridComp)    :: driver
!
!  Local variable declarations.
!
      integer :: localPET
!
      real(r8) :: Time_Step
!
      character (len=100)  :: name

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ESM_SetRunSequence"
!
      TYPE (ESMF_Clock)        :: Clock
      TYPE (ESMF_Config)       :: Config
      TYPE (ESMF_TimeInterval) :: TimeStep
      TYPE (NUOPC_FreeFormat)  :: runSeqFF
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering ESM_SetRunSequence',      &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Set run sequence from input configuration file.
!-----------------------------------------------------------------------
!
!  Query driver for information.
!
      CALL ESMF_GridCompGet(driver,                                     &
     &                      name=name,                                  &
     &                      rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Read in free-format run sequence from configuration object.
!
      CALL ESMF_GridCompGet(driver,                                     &
     &                      config=Config,                              &
     &                      rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      runSeqFF = NUOPC_FreeFormatCreate(Config,                         &
     &                                  label='runSeq::',               &
     &                                  rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Ingest free-format run sequence.
!
      CALL NUOPC_DriverIngestRunSequence (driver,                       &
     &                                    runSeqFF,                     &
     &                                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        IF (localPET.eq.0) THEN
          WRITE (cplout,10) TRIM(CONFname)
        END IF
        RETURN
      END IF
!
!  Get driver clock object.
!
      CALL ESMF_GridCompGet (driver,                                    &
     &                       clock=Clock,                               &
     &                       localPet=localPET,                         &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Inquire clock for the coupling time step set in ingested run
!  sequence.
!
      CALL ESMF_ClockGet (Clock,                                        &
     &                    timeStep=TimeStep,                            &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL ESMF_TimeIntervalGet (TimeStep,                              &
     &                           s_r8=Time_Step,                        &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Check ingested coupling time step value to the one provided in the
!  standard input script. Both values need to be the same.
!
      IF (Time_Step.ne.ClockInfo(Idriver)%Time_Step) THEN
        IF (localPET.eq.0) THEN
          WRITE (cplout,20) ClockInfo(Idriver)%Time_Step,               &
     &                      TRIM(CinpName),                             &
     &                      Time_Step,                                  &
     &                      TRIM(CONFname)
        END IF
        rc=ESMF_RC_VAL_WRONG
        RETURN
      END IF
!
!  Report internal driver information.
!
      IF (DebugLevel.ge.2) THEN
        CALL NUOPC_DriverPrint (driver,                                 &
     &                          orderflag=.TRUE.,                       &
     &                          rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
        CALL my_flush (6)                   ! flush standard output unit
      END IF
!
!  Destroy free format object. All internal memory is deallocated.
!
      CALL NUOPC_FreeFormatDestroy(runSeqFF,                            &
     &                             rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  ESM_SetRunSequence',      &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
 10   FORMAT (/,' ESM_SetRunSequence - Error while ingesting',          &
     &        ' RunSequence configuration file:',/,22x,a,               &
     &        /,22x,'Check if connections (->) between components '     &
     &        ' are needed or not.')
 20   FORMAT (/,' ESM_SetRunSequence - Inconsistent coupling time',     &
     &        ' step (seconds) from ingested RunSequence:',             &
     &        /,22x,'TimeStep = ',f15.8,2x,'(in ',a,')',                &
     &        /,22x,'TimeStep = ',f15.8,2x,'(in ',a,')',                &
     &        /,22x,'Correct either input file to the desired value.',  &
     &        /,22x,'The value needs to be the same in both files!')
!
      RETURN
      END SUBROUTINE ESM_SetRunSequence
#endif
      END MODULE esmf_esm_mod
