      PROGRAM esmf_coupler
!
!svn $Id: esmf_coupler.h 584 2008-03-18 20:01:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Master program to couple ROMS/TOMS to other models using the ESMF   !
!  (Earth System Modeling Framework) library.                          !
!                                                                      !
!  The following models are coupled to ROMS/TOMS:                      !
!                                                                      !
#ifdef WRF_COUPLING
!  WRF, Weather Research and Forecasting model:                        !
!       http://www.wrf-model.org                                       !
!                                                                      !
#endif
#ifdef SWAN_COUPLING
!  SWAN, Simulating WAves Nearshore model:                             !
!        http://vlm089.citg.tudelft.nl/swan/index.htm                  !
!                                                                      !
#endif
!=======================================================================
!
      USE ESMF_mod
      USE mod_parallel
      USE mod_coupler
!
      USE ESMF_ROMS, ONLY : OCEAN_SetServices => ROMS_SetServices
#ifdef SWAN_COUPLING
# define WAVES_COUPLING
      USE ESMF_SWAN, ONLY : WAVES_SetServices => SWAN_SetServices
#endif
#ifdef WRF_COUPLING
# define ATMOS_COUPLING
      USE ESMF_WRF,  ONLY : ATMOS_SetServices => WRF_SetServices
#endif
!
      implicit none
!
!  Local variable declarations.
!
      logical :: CheckError

      integer :: model, status
      integer :: MyMPIcomm, MyPetRank, Nnodes

      character (len=80) :: ExpName, GridName, ImpName, name

      TYPE (ESMF_VM) :: MyVM
!
!-----------------------------------------------------------------------
!  Initialize distributed-memory (MPI) configuration
!-----------------------------------------------------------------------
!
!  Initialize error flag.
!
      status=ESMF_SUCCESS
!
!  Initialize ESMF and Virtual Machine (VM) parallel environment.
!
      CALL ESMF_Initialize (defaultLogType=ESMF_LOG_SINGLE,             &
     &                      vm=MyVM,                                    &
     &                      rc=status)
      IF (CheckError(status, 'ESMF', 'esmf_coupler.h',                  &
     &               'initializing ESMF and VM environment')) THEN
        CALL abort (status)
      END IF
!
!  Inquire VM for the MPI communicator handle.
!
      CALL ESMF_VMGet (MyVM,                                            &
     &                 localPet=MyPetRank,                              &
     &                 petCount=Nnodes,                                 &
     &                 mpiCommunicator=MyMPIcomm,                       &
     &                 rc=status)
      IF (CheckError(status, 'ESMF', 'esmf_coupler.h',                  &
     &               'getting mpiCommunicator handle')) THEN
        CALL abort (status)
      END IF
!
!  Assign ocean global communicator to ocean communicator.
!
      OCN_COMM_WORLD=MyMPIcomm
!
!  Read in coupling models parameters from standard input.
!
      CALL read_CouplePar (iNLM)
!
!  Allocate several coupling variables.
!
      CALL allocate_coupler (Nnodes)
!
!  Store VM object. Notice that all the coupling arrays have been
!  allocated/associated.
!
      VM(0)=MyVM
!
!-----------------------------------------------------------------------
!  Create gridded components. Assign parallel nodes for each component.
!-----------------------------------------------------------------------
!
      DO model=1,Nmodels
        IF (model.eq.Iocean) THEN
          GridName='Ocean Model Component, ROMS'
        ELSE IF (model.eq.Iwaves) THEN
          GridName='Wave Model Component'
        ELSE IF (model.eq.Iatmos) THEN
          GridName='Atmosphere Model Component'
        END IF
        GridComp(model)=ESMF_GridCompCreate (name=TRIM(GridName),       &
     &                                       petList=pets(model)%val,   &
     &                                       rc=status)
        IF (CheckError(status, 'ESMF', 'esmf_coupler.h',                &
     &                 'creating '//TRIM(GridName))) THEN
          CALL abort (status)
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Set gridded components services. Register "initialize", "run", and
!  "finalize" routines for each model components. Notice that the actual
!  subroutine name needs to be passed as argument. This is not a
!  character string, so we cannot use a do loop over all coupled
!  components.
!-----------------------------------------------------------------------
!
!  Ocean model registration.
!
      CALL ESMF_GridCompSetServices (GridComp(Iocean),                  &
     &                               OCEAN_SetServices,                 &
     &                               status)
      IF (CheckError(status, 'ESMF', 'esmf_coupler.h',                  &
     &               'setting ocean model services')) THEN
        CALL abort (status)
      END IF

#ifdef WAVES_COUPLING
!
!  Wave model registration
!
      CALL ESMF_GridCompSetServices (GridComp(Iwaves),                  &
     &                               WAVES_SetServices,                 &
     &                               status)
      IF (CheckError(status, 'ESMF', 'esmf_coupler.h',                  &
     &               'setting wave model services')) THEN
        CALL abort (status)
      END IF
#endif
#ifdef ATMOS_COUPLING
!
!  Atmosphere model registration
!
      CALL ESMF_GridCompSetServices (GridComp(Iatmos),                  &
     &                               ATMOS_SetServices,                 &
     &                               status)
      IF (CheckError(status, 'ESMF', 'esmf_coupler.h',                  &
     &               'setting atmosphere model services')) THEN
        CALL abort (status)
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Create gridded components export/import state objects.
!-----------------------------------------------------------------------
!
      DO model=1,Nmodels
        IF (model.eq.Iocean) THEN
          ExpName='Ocean Model Component Export State'
          ImpName='Ocean Model Component Import State'
        ELSE IF (model.eq.Iwaves) THEN
          ExpName='Wave Model Component Export State'
          ImpName='Wave Model Component Import State'
        ELSE IF (model.eq.Iatmos) THEN
          ExpName='Atmosphere Model Component Export State'
          ImpName='Atmosphere Model Component Import State'
        END IF
        StateExport(model)=ESMF_StateCreate (TRIM(ExpName),             &
     &                                       ESMF_STATE_EXPORT,         &
     &                                       rc=status)
        IF (CheckError(status, 'ESMF', 'esmf_coupler.h',                &
     &                 'creating '//TRIM(ExpName))) THEN
          CALL abort (status)
        END IF
!
        StateImport(model)=ESMF_StateCreate (TRIM(ImpName),             &
     &                                       ESMF_STATE_IMPORT,         &
     &                                       rc=status)
        IF (CheckError(status, 'ESMF', 'esmf_coupler.h',                &
     &                 'creating '//TRIM(ImpName))) THEN
          CALL abort (status)
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Initialize gridded components.
!-----------------------------------------------------------------------
!
      DO model=1,Nmodels
        CALL ESMF_GridCompInitialize (GridComp(model),                  &
     &                                importState=StateImport(model),   &
     &                                exportState=StateExport(model),   &
     &                                rc=status)
        IF (CheckError(status, 'ESMF', 'esmf_coupler.h',                &
     &                 'initializing gridded component')) THEN
          CALL abort (status)
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Check all coupled models internal time clocks and determine minimum
!  start time, maximum stop time, and minimum time step. Then, set
!  these values to vector index zero.  It is assumed that all clocks
!  have the same calendar.
!-----------------------------------------------------------------------
!
      StartTime(0)=StartTime(1)
      StopTime(0)=StopTime(1)
      TimeStep(0)=TimeStep(1)
      DO model=2,Nmodels
        IF (StartTime(model) < StartTime(0)) THEN
          StartTime(0)=StartTime(model)
        END IF
        IF (StopTime(model) < StopTime(0)) THEN
          StopTime(0)=StopTime(model)
        END IF
        IF (TimeStep(model) < TimeStep(0)) THEN
          TimeStep(0)=TimeStep(model)
        END IF
      END DO
!
!  Create external clock.
!
      name='External Time Clock'
      TimeClock(0)=ESMF_ClockCreate(TRIM(name),                         &
     &                              timeStep=TimeStep(0),               &
     &                              startTime=StartTime(0),             &
     &                              stopTime=StopTime(0),               &
     &                              rc=status)
      IF (CheckError(status, 'ESMF', 'esmf_coupler.h',                  &
     &               'creating external time clock')) RETURN
!
!  Validate external time clock.
!
      CALL ESMF_ClockValidate (TimeClock(0),                            &
     &                         rc=status)
      IF (CheckError(status, 'ESMF', 'ROMS_SetClock',                   &
     &               'validating external time clock')) RETURN
!
!-----------------------------------------------------------------------
!  Run gridded components.
!-----------------------------------------------------------------------
!
      DO WHILE (.not.ESMF_ClockIsStopTime(TimeClock(0)))
!
!  Advance external time clock.
!
        CALL ESMF_ClockAdvance (TimeClock(0),                           &
     &                          rc=status)
        IF (CheckError(status, 'ESMF', 'esmf_coupler.h',                &
     &                 'advancing external time clock')) EXIT
!
!  Get external clock current time.
!
        CALL ESMF_ClockGet (TimeClock(0),                               &
     &                      currTime=CurrTime(0),                       &
     &                      rc=status)
        IF (CheckError(status, 'ESMF', 'esmf_coupler.h',                &
     &                 'getting external current time object')) EXIT
        CALL ESMF_TimeGet (CurrTime(0),                                 &
     &                     yy=CurrentTime(0)%year,                      &
     &                     mm=CurrentTime(0)%month,                     &
     &                     dd=CurrentTime(0)%day,                       &
     &                     h=CurrentTime(0)%hour,                       &
     &                     m=CurrentTime(0)%minute,                     &
     &                     s=CurrentTime(0)%second,                     &
     &                     timeZone=CurrentTime(0)%TimeZone,            &
     &                     timeStringISOFrac=CurrentTime(0)%string,     &
     &                     dayOfYear=CurrentTime(0)%YearDay,            &
     &                     rc=status)
        IF (CheckError(status, 'ESMF', 'esmf_coupler.h',                &
     &                 'getting external current time')) EXIT
!
!  Run coupled models.
!
        DO model=1,Nmodels
          CALL ESMF_GridCompRun (GridComp(model),                       &
     &                           exportState=StateExport(model),        &
     &                           importState=StateImport(model),        &
     &                           rc=status)
          IF (CheckError(status, 'ESMF', 'esmf_coupler.h',              &
     &                   'running gridded')) THEN
            CALL abort (status)
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Finalize gridded components.
!-----------------------------------------------------------------------
!
      DO model=1,Nmodels
        CALL ESMF_GridCompFinalize (GridComp(model),                    &
     &                              exportState=StateExport(model),     &
     &                              importState=StateImport(model),     &
     &                              rc=status)
        IF (CheckError(status, 'ESMF', 'esmf_coupler.h',                &
     &                 'finalizing gridded component')) THEN
            CALL abort (status)
        END IF
      END DO
!
!  Terminates all the ESMF/MPI processing.
!
      CALL esmf_finalize (rc=status)

      END PROGRAM esmf_coupler
