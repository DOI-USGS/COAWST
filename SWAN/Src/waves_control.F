#include "cppdefs.h"
      MODULE waves_control_mod
!
!svn $Id: waves_control.F 814 2008-10-29 01:42:17Z jcwarner $
!================================================== John C. Warner  ====
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  SWAN model:                                                         !
!                                                                      !
!  This driver executes SWAN by controlling initialization,            !
!  time-stepping, and finalization for nested grids.                   !
!                                                                      !
!     SWAN_initialize                                                  !
!     SWAN_run                                                         !
!     SWAN_finalize                                                    !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: SWAN_driver
      PUBLIC  :: SWAN_driver_run
      PUBLIC  :: SWAN_driver_finalize

      CONTAINS

      SUBROUTINE SWAN_driver (MyCOMM, CouplingTime, SWAN_inplist)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes ROMS/TOMS state variables    !
!  and internal and external parameters.                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
!     USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE mod_coupler
      USE M_GENARR
      USE M_PARALL
      USE M_PARALL2
      USE SWCOMM3
      USE M_MPI
!
!  Imported variable declarations.
!
      integer, intent(in) :: MyCOMM
      real(r4), intent(in)  :: CouplingTime             ! single precision
      character*80, intent(in) :: SWAN_inplist
!
!  Local variable declarations.
!
      integer :: ng, inp, ModelId
      integer :: MyColor, MyKey
      INTEGER :: MySize, MyError, MyTag, MyStatus, MyHandle
      integer :: MyCOMM_wave
      integer :: count1, MyRank
#if defined REFINED_GRID
      integer :: count2, count3, count4
#endif
!      ALLOCATE(WAV_COMM_WORLDS(Ngrids))
      MyKey=0
      MyTag=52
      inp=22
      Numspec=17
!      Numspecl=0
      NBGRPT_last=0
!
!-----------------------------------------------------------------------
!  Set distribute-memory (MPI) world communicator.
!-----------------------------------------------------------------------
!
      CALL mpi_barrier(MyCOMM, MyError)
      CALL mpi_comm_rank (MyCOMM, MyRank, MyError)
      CALL mpi_comm_size (MyCOMM, MySize, MyError)
!
!     IF (MyRank.eq.MyMaster) THEN
      OPEN (inp, FILE=TRIM(SWAN_inplist), FORM='formatted',             &
     &      STATUS='old', ERR=20)
      DO ng=1,Ngrids
        READ(inp,'(a80)') SWAN_infile(ng)
      END DO
      CLOSE(inp)
      GO TO 40
 20   WRITE (stdout,30)
      STOP
 30   FORMAT (/,' swan_inp - Unable to open SWAN file list')
 40   CONTINUE

!     END IF
!
!  Initialize coupling streams between model(s).
!
!  NEED to be consistent with mysize/ngridsand mpi_send / recv below.
!      DO ng=1,Ngrids
!        WAV_COMM_WORLDS(ng)=0
!      END DO
      IF (MyRank.le.(MySize/Ngrids-1)) THEN
        ng=1
        Modelid=Ngrids+ng
        MyColor=Modelid
      ELSE
        ng=2
        Modelid=Ngrids+ng
        MyColor=Modelid
      END IF
      CALL mpi_comm_split (MyCOMM, MyColor, MyKey, MyCOMM_wave,         &
     &                     MyError)
! send all wave worlds to all processors
!      CALL mpi_allgather(MyCOMM_WAVE,1,MPI_INTEGER,WAV_COMM_WORLDS,1,   &
!     &                   MPI_INTEGER,MyCOMM,MyError)

!       if (ng.eq.1) then 
!         CALL MPI_Intercomm_create(WAV_COMM_WORLDS(1),0,       &
!     &            WAV_COMM_WORLDS(2),0,MyTag,WAVCOM,MyError)
!         CALL mpi_comm_rank (WAVCOM, MyRank, MyError)
!         CALL mpi_comm_size (WAVCOM, MySize, MyError)
!       endif

      IF (.not.ALLOCATED(BGRIDPG)) ALLOCATE(BGRIDPG(100000))
      BGRIDPG=9999

      CALL SWAN_INITIALIZE (ng, MyCOMM_wave, SWAN_infile(ng),           &
     &                      NumSWANsteps)
      CALL SWOUTP ( AC2   , SPCSIG, SPCDIR, COMPDA, XYTST ,             &
     &              KGRPNT, XCGRID, YCGRID, KGRBND, OURQT ,             &
     &              0, 0, ng, -1)

! Broadcast the number of boundary arrays to be passed

      ac2size=MSC*MDC*Numspec
!
! Start big time loop here.
!
      NL_LOOP : DO count1=0, NumSWANsteps
        MyTag=0
!
! Do a barrier so we all start together.
!
        CALL mpi_barrier(MyCOMM, MyError)
!        write(*,*) 'mycomm barrier ', ng
!        CALL mpi_comm_rank (MyCOMM, MyRank, MyError)
!        CALL mpi_comm_size (MyCOMM, MySize, MyError)
!
! First grid.
!
        sc(1)=count1
        IF (ng.eq.1) THEN
!          write(*,*) 'call swan run ng = ', ng
          CALL SWAN_RUN (CouplingTime, sc(1), ng, count1)
!         IF (count1.gt.1) THEN
!            CALL MPI_SEND(AC2_parent2,ac2size,MPI_REAL,              &
!     &                         1,MyTag,MyCOMM, MyError)
!          END IF
        END IF
#if defined REFINED_GRID
!        CALL SWBROADC ( AC2_parent2, ac2size, SWREAL )

        CALL MPI_BCAST(AC2_parent2,ac2size,SWREAL,                 &
     &                0,MyCOMM,MyError)
!
! Second grid.
!
!        CALL mpi_barrier(MyCOMM, MyError)
        IF (ng.eq.2) THEN
          DO count2=1,5           !nrefined(2)
!            write(*,*) 'call swan run ng = ', ng
            IF (count2.eq.1) THEN
!              CALL MPI_RECV(AC2_parent2,ac2size,MPI_REAL,               &
!     &                      0,MyTag,MyCOMM, MyStatus,MyError)
            END IF
!           sc(2)=sc(1)*nrefined(1)*nrefined(2)+count2
            sc(2)=sc(1)*1*5+count2-1
            CALL SWAN_RUN (CouplingTime, sc(2), ng, count2-1)
          END DO
        END IF
#endif
        IF (exit_flag.ne.NoError) THEN
          IF (MyRank.eq.0) THEN
            WRITE (stdout,'(/,a,i3,/)') Rerror(exit_flag), exit_flag
          END IF
          RETURN
        END IF
      END DO NL_LOOP
!      DEALLOCATE(WAV_COMM_WORLDS)

      RETURN
      END SUBROUTINE SWAN_driver

      SUBROUTINE SWAN_driver_run (CouplingTime)
!
!=======================================================================
!                                                                      !
!  This routine runs SWAN.                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#ifdef VERIFICATION
      USE mod_fourdvar
#endif
      USE mod_iounits
      USE mod_scalars
      USE mod_coupler
      USE M_PARALL
!
!  Imported variable declarations.
!
      real(r4), intent(in)  :: CouplingTime             ! single precision
!
!  Local variable declarations.
!
      integer :: ng, my_iic, MyError
      integer :: count1, MyCOMM_wave
#if defined REFINED_GRID
      integer :: count2, count3, count4
#endif
!
!-----------------------------------------------------------------------
!  Run model for all nested grids, if any.
!-----------------------------------------------------------------------

      NL_LOOP : DO count1=0, NumSWANsteps

        IF (ng.eq.1) THEN
          sc(1)=count1
          CALL SWAN_RUN (CouplingTime, sc(1))
        END IF

#if defined REFINED_GRID
        IF (ng.eq.2) THEN
          DO count2=1,nrefined(2)
            CALL mpi_barrier(MyCOMM_wave, MyError)
            sc(2)=sc(1)*nrefined(1)*nrefined(2)+count2
            CALL SWAN_RUN (CouplingTime, sc(ng))
          END DO
        END IF


!        IF (Ngrids.ge.2) THEN
!          DO count2=1,nrefined(2)
!            CALL mpi_barrier(MyCOMM_wave, MyError)
!            ng=2
!            sc(2)=sc(1)*nrefined(1)*nrefined(2)+count2

!            CALL SWAN_RUN (CouplingTime, sc(ng))

!            IF (Ngrids.ge.3) THEN
!              ng=3
!              DO count3=1,nrefined(ng)
!                CALL mpi_barrier(MyCOMM_wave, MyError)
!                sc(ng)=sc(1)*nrefined(1)*nrefined(2)*                  &
!     &                nrefined(ng)+(count2-1)*nrefined(ng)+count3

!                CALL SWAN_RUN (CouplingTime, sc(ng))

!                IF (Ngrids.ge.4) THEN
!                  ng=4
!                  DO count4=1,nrefined(ng)
!                    CALL mpi_barrier(MyCOMM_wave, MyError)
!                    sc(ng)=sc(1)*nrefined(1)*nrefined(2)*              &
!     &                      nrefined(ng-1)*nrefined(ng)+                  &
!     &                      (count2-1)*nrefined(2)*nrefined(ng-1)+       &
!     &                      (count3-1)*nrefined(ng-1)+count4

!                    CALL SWAN_RUN (CouplingTime, sc(ng))

!                  END DO
!                END IF
!              END DO
!            END IF
!          END DO
!        END IF
#endif
        IF (exit_flag.ne.NoError) THEN
          IF (MyRank.eq.0) THEN
            WRITE (stdout,'(/,a,i3,/)') Rerror(exit_flag), exit_flag
          END IF
          RETURN
        END IF

      END DO NL_LOOP


      RETURN
      END SUBROUTINE SWAN_driver_run

      SUBROUTINE SWAN_driver_finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates SWAN.                                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
!  Local variable declarations.
!
      integer :: ng

      DO ng=1,Ngrids
        CALL SWAN_FINALIZE
      END DO

      RETURN
      END SUBROUTINE SWAN_driver_finalize

      END MODULE waves_control_mod
