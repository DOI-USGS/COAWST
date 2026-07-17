      PROGRAM esmf_driver
!
!git $Id$
!svn $Id: esmf_driver.h 1151 2023-02-09 03:08:53Z arango $
!=======================================================================
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license         Hernan G. Arango     !
!    See License_ROMS.txt                         Ufuk Utku Turuncoglu !
!=======================================================================
!                                                                      !
!  Master program to couple ROMS/TOMS to other Earth System Models     !
!  (ESM) using the ESMF library with the NUOPC layer.                  !
!                                                                      !
!  ESMF:   Earth System Modeling Framework (Version 7 or higher)       !
!  NUOPC:  National Unified Operational Prediction Capability          !
!                                                                      !
!=======================================================================
!
      USE ESMF
!
      USE mod_esmf_esm          ! ESM coupling structures and variables
!
      USE esmf_esm_mod, ONLY : ESM_SetServices
!
      implicit none
!
!  Local variable declarations.
!
      integer :: rc, urc
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
      TYPE (ESMF_GridComp)     :: driver
      TYPE (ESMF_LogKind_flag) :: LogKindFlag
!
!-----------------------------------------------------------------------
!  Initialize ESMF coupling.
!-----------------------------------------------------------------------
!
!  It is optimal to have a single log file using ESMF_LOGKIND_SINGLE
!  combining messages from all PETs. However, it is not supported on
!  some platforms.  If this is the case, use ESMF_LOGKIND_MULTI to
!  create one log file per PET. This is kind of annoying. In
!  applications using s large number of processors, opening a large
!  number of log files and writing messages can be a bottleneck.
!
      LogKindFlag=ESMF_LOGKIND_SINGLE
!!    LogKindFlag=ESMF_LOGKIND_MULTI
!
!  Initialize ESMF.
!
      CALL ESMF_Initialize (defaultLogFileName=TRIM(ESMnameLog),        &
     &                      logappendflag=.FALSE.,                      &
     &                      logkindflag=LogKindFlag,                    &
     &                      vm=VMdriver,                                &
     &                      rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        CALL ESMF_Finalize (endflag=ESMF_END_ABORT)
      END IF
!
!-----------------------------------------------------------------------
!  Create ESM component.
!-----------------------------------------------------------------------
!
      driver=ESMF_GridCompCreate(name="roms_esm",                       &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        CALL ESMF_Finalize (endflag=ESMF_END_ABORT)
      END IF
!
!-----------------------------------------------------------------------
!  Coupling ESM configuration.
!-----------------------------------------------------------------------
!
!  Read in coupling ESM configuration parameters from standard input.
!  (See "ROMS/External/coupling_esmf.in").
!
      CALL read_ESMconfig (VMdriver, rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        CALL ESMF_Finalize (endflag=ESMF_END_ABORT)
      END IF
!
!  Set user defined ESM log.
!
      IF (TraceLevel.ge.2) THEN
        CALL ESMF_LogSet (flush=.TRUE.,                                 &
     &                    logmsgList=(/ ESMF_LOGMSG_ALL /),             &
     &                    trace=.TRUE.,                                 &
     &                    rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          CALL ESMF_Finalize (endflag=ESMF_END_ABORT)
        END IF
      ELSE
        CALL ESMF_LogSet (flush=.TRUE.,                                 &
     &                    logmsgList=(/ ESMF_LOGMSG_NOTRACE /),         &
     &                    trace=.FALSE.,                                &
     &                    rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          CALL ESMF_Finalize (endflag=ESMF_END_ABORT)
        END IF
      END IF
!
!  Read in and set import and export fields metadata. Add additional
!  fields to NOUPC dictionary.
!
      CALL set_metadata (VMdriver, rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        CALL ESMF_Finalize (endflag=ESMF_END_ABORT)
      END IF
!
!-----------------------------------------------------------------------
!  Register ESM components.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompSetServices (driver,                            &
     &                               ESM_SetServices,                   &
     &                               userRc=urc,                        &
     &                               rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        CALL ESMF_Finalize (endflag=ESMF_END_ABORT)
      END IF
      IF (ESMF_LogFoundError(rcToCheck=urc,                             &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        CALL ESMF_Finalize (endflag=ESMF_END_ABORT)
      END IF
!
!-----------------------------------------------------------------------
!  Initialize ESM component.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompInitialize (driver,                             &
     &                              userRc=urc,                         &
     &                              rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        CALL ESMF_Finalize (endflag=ESMF_END_ABORT)
      END IF
      IF (ESMF_LogFoundError(rcToCheck=urc,                             &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        CALL ESMF_Finalize (endflag=ESMF_END_ABORT)
      END IF
!
!  Wait for the initialization phase to finish.
!
      CALL ESMF_VMBarrier (VMdriver, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        CALL ESMF_Finalize (endflag=ESMF_END_ABORT)
      END IF
!
!-----------------------------------------------------------------------
!  Run ESM components.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompRun (driver,                                    &
     &                       userRc=urc,                                &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        CALL ESMF_Finalize (endflag=ESMF_END_ABORT)
      END IF
      IF (ESMF_LogFoundError(rcToCheck=urc,                             &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        CALL ESMF_Finalize (endflag=ESMF_END_ABORT)
      END IF
!
!-----------------------------------------------------------------------
!  Finalize ESM components.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompFinalize (driver,                               &
     &                            userRc=urc,                           &
     &                            rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        CALL ESMF_Finalize (endflag=ESMF_END_ABORT)
      END IF
      IF (ESMF_LogFoundError(rcToCheck=urc,                             &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        CALL ESMF_Finalize (endflag=ESMF_END_ABORT)
      END IF
!
!-----------------------------------------------------------------------
!  Destroy ESM components.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompDestroy (driver,                                &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        CALL ESMF_Finalize (endflag=ESMF_END_ABORT)
      END IF
!
!-----------------------------------------------------------------------
!  Finalize ESMF coupling.
!-----------------------------------------------------------------------
!
      CALL ESMF_Finalize (rc=rc)
!
!  Flush and lose coupling standard out unit.
!
      CALL my_flush (cplout)
      CLOSE (cplout)
!
      STOP
      END PROGRAM esmf_driver
