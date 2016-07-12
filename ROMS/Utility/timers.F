#include "cppdefs.h"
      SUBROUTINE wclock_on (ng, model, region)
!
!svn $Id: timers.F 795 2016-05-11 01:42:43Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine turns on wall clock to meassure the elapsed time in    !
!  seconds spend by each parallel thread in requested model region.    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_strings
#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_barrier
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) ::  ng, model, region
!
!  Local variable declarations.
!
      integer :: iregion, MyModel, NSUB

      integer :: my_getpid
#ifndef DISTRIBUTE
      integer :: my_threadnum
#endif

      real(r8), dimension(2) :: wtime

      real(r8) :: my_wtime
!
!-----------------------------------------------------------------------
! Initialize timing for all threads.
!-----------------------------------------------------------------------
!
!  Set number of subdivisions, same as for global reductions.
!
#ifdef MPI
      NSUB=1
#else
      NSUB=numthreads
#endif
!
!  Ensure that MyModel is not zero.  Notice that zero value is used to
!  indicate restart of the nonlinear model.
!
      MyModel=MAX(1,model)
      Cstr(region,MyModel,ng)=my_wtime(wtime)
      IF ((region.eq.0).and.(proc(1,MyModel,ng).eq.0)) THEN
        DO iregion=1,Nregion
          Cend(iregion,MyModel,ng)=0.0_r8
          Csum(iregion,MyModel,ng)=0.0_r8
        END DO
        proc(1,MyModel,ng)=1
        proc(0,MyModel,ng)=my_getpid()

!$OMP CRITICAL (START_WCLOCK)
#ifdef DISTRIBUTE
        CALL mp_barrier (ng)
        WRITE (stdout,10) ' Node #', MyRank,                            &
     &                    ' (pid=',proc(0,MyModel,ng),') is active.'
        CALL my_flush (stdout)
#else
        WRITE (stdout,10) ' Thread #', MyThread,                        &
     &                    ' (pid=',proc(0,MyModel,ng),') is active.'
#endif
 10     FORMAT (a,i3,a,i8,a)
        thread_count=thread_count+1
        IF (thread_count.eq.NSUB) thread_count=0
!$OMP END CRITICAL (START_WCLOCK)
      END IF
      RETURN
      END SUBROUTINE wclock_on
      SUBROUTINE wclock_off (ng, model, region)
!
!=======================================================================
!                                                                      !
!  This routine turns off wall clock to meassure the elapsed time in   !
!  seconds spend by each parallel thread in requested model region.    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_strings
#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_barrier, mp_reduce
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) ::  ng, model, region
!
!  Local variable declarations.
!
      integer :: imodel, iregion, MyModel, NSUB

      integer :: my_threadnum

      real(r8) :: percent, sumcpu, sumper, total

      real(r8), dimension(2) :: wtime

      real(r8) :: my_wtime

#ifdef DISTRIBUTE
      real(r8), dimension(0:Nregion) :: buffer

      character (len= 3), dimension(0:Nregion) :: op_handle
#endif
      character (len=14), dimension(4) :: label
!
!-----------------------------------------------------------------------
!  Compute elapsed wall time for all threads.
!-----------------------------------------------------------------------
!
!  Set number of subdivisions, same as for global reductions.
!
#ifdef MPI
      NSUB=1
#else
      NSUB=numthreads
#endif
!
!  Insure that MyModel is not zero.  Notice that zero value is used to
!  indicate restart of the nonlinear model.
!
      MyModel=MAX(1,model)
      IF (region.ne.0) THEN
        Cend(region,MyModel,ng)=Cend(region,MyModel,ng)+                &
     &                          (my_wtime(wtime)-                       &
     &                           Cstr(region,MyModel,ng))
      END IF
!
!  Report elapsed wall time.
!
      IF ((region.eq.0).and.(proc(1,MyModel,ng).eq.1)) THEN
        Cend(region,MyModel,ng)=Cend(region,MyModel,ng)+                &
     &                          (my_wtime(wtime)-                       &
     &                           Cstr(region,MyModel,ng))
        DO imodel=1,4
          proc(1,imodel,ng)=0
        END DO

!$OMP CRITICAL (FINALIZE_WCLOCK)
!
!  Report total elapsed time (seconds) for each CPU.  We get the same
!  time for all nested grids.
!
        IF (ng.eq.1) THEN
#ifdef DISTRIBUTE
          CALL mp_barrier (ng)
          WRITE (stdout,10) ' Node   #', MyRank, ' CPU:',               &
     &                      Cend(region,MyModel,ng)
         CALL my_flush (stdout)
#else
         WRITE (stdout,10) ' Thread #', MyThread, ' CPU:',              &
     &                     Cend(region,MyModel,ng)
#endif
 10      FORMAT (a,i3,a,f12.3)
       END IF
!
! Report elapsed time profile for each region of the code.
!
        thread_count=thread_count+1
        DO imodel=1,4
          Csum(region,imodel,ng)=Csum(region,imodel,ng)+                &
     &                           Cend(region,imodel,ng)
          DO iregion=1,Nregion
            Csum(iregion,imodel,ng)=Csum(iregion,imodel,ng)+            &
     &                              Cend(iregion,imodel,ng)
          END DO
        END DO
        DO imodel=1,4
          IF (imodel.ne.MyModel) THEN
            DO iregion=1,Nregion
              Csum(region,imodel,ng)=Csum( region,imodel,ng)+           &
     &                               Csum(iregion,imodel,ng)
            END DO
          END IF
        END DO
        IF (thread_count.eq.NSUB) THEN
          thread_count=0
#ifdef DISTRIBUTE
          op_handle(0:Nregion)='SUM'
          DO imodel=1,4
            DO iregion=0,Nregion
              buffer(iregion)=Csum(iregion,imodel,ng)
            END DO
            CALL mp_reduce (ng, MyModel, Nregion+1, buffer(0:),         &
     &                      op_handle(0:))
            DO iregion=0,Nregion
              Csum(iregion,imodel,ng)=buffer(iregion)
            END DO
          END DO
#endif
          IF (Master.and.(ng.eq.1)) THEN
            total=0.0_r8
            DO imodel=1,4
              total=total+Csum(region,imodel,ng)
            END DO
            WRITE (stdout,20) ' Total:', total
 20         FORMAT (a,8x,f14.3)
          END IF
#ifdef PROFILE
!
!  Report profiling times.
!
          label(iNLM)='Nonlinear     '
          label(iTLM)='Tangent linear'
          label(iRPM)='Representer   '
          label(iADM)='Adjoint       '
          DO imodel=1,4
            IF (Master.and.(Csum(region,imodel,ng).gt.0.0_r8)) THEN
              WRITE (stdout,30) TRIM(label(imodel)),                    &
     &                    'ocean model elapsed time profile, Grid:',    &
     &                          ng
 30           FORMAT (/,1x,a,1x,a,1x,i2.2/)
            END IF
            sumcpu=0.0_r8
            sumper=0.0_r8
            DO iregion=1,Mregion-1
              IF (Csum(iregion,imodel,ng).gt.0.0_r8) THEN
                percent=100.0_r8*Csum(iregion, imodel,ng)/              &
     &                           Csum( region,MyModel,ng)
                IF (Master) WRITE (stdout,40) Pregion(iregion),         &
     &                                        Csum(iregion,imodel,ng),  &
     &                                        percent
                sumcpu=sumcpu+Csum(iregion,imodel,ng)
                sumper=sumper+percent
              END IF
            END DO
 40         FORMAT (2x,a,t53,f14.3,2x,'(',f7.4,' %)')
            IF (Master.and.(Csum(region,imodel,ng).gt.0.0_r8)) THEN
              WRITE (stdout,50) sumcpu, sumper
 50           FORMAT (t47,'Total:',f14.3,2x,f8.4)
            END IF
          END DO
# ifdef ICE_MODEL
!
!  Report elapsed time for message passage communications.
!
          DO imodel=1,4
            total=0.0_r8
            DO iregion=Iceregion,Nregion
              total=total+Csum(iregion,imodel,ng)
            END DO
            IF (Master.and.(total.gt.0.0_r8)) THEN
              WRITE (stdout,30) TRIM(label(imodel)),                    &
     &                'sea-ice model elapsed time profile, Grid:', ng
            END IF
            sumcpu=0.0_r8
            sumper=0.0_r8
            IF (total.gt.0.0_r8) THEN
              DO iregion=Iceregion,Nregion
                IF (Csum(iregion,imodel,ng).gt.0.0_r8) THEN
                  percent=100.0_r8*Csum(iregion, imodel,ng)/            &
     &                             Csum( region,Mymodel,ng)
                  IF (Master) WRITE (stdout,40) Pregion(iregion),       &
     &                                          Csum(iregion,imodel,ng),&
     &                                          percent
                  sumcpu=sumcpu+Csum(iregion,imodel,ng)
                  sumper=sumper+percent
                END IF
              END DO
              IF (Master.and.(total.gt.0.0_r8)) THEN
                WRITE (stdout,50) sumcpu, sumper
              END IF
            END IF
          END DO
# endif
# ifdef DISTRIBUTE
!
!  Report elapsed time for message passage communications.
!
          DO imodel=1,4
            total=0.0_r8
            DO iregion=Mregion,Iceregion-1
              total=total+Csum(iregion,imodel,ng)
            END DO
            IF (Master.and.(total.gt.0.0_r8)) THEN
              WRITE (stdout,30) TRIM(label(imodel)),                    &
     &                          'model message Passage profile, Grid:', &
     &                          ng
            END IF
            sumcpu=0.0_r8
            sumper=0.0_r8
            IF (total.gt.0.0_r8) THEN
              DO iregion=Mregion,Iceregion-1
                IF (Csum(iregion,imodel,ng).gt.0.0_r8) THEN
                  percent=100.0_r8*Csum(iregion, imodel,ng)/            &
     &                             Csum( region,Mymodel,ng)
                  IF (Master) WRITE (stdout,40) Pregion(iregion),       &
     &                                          Csum(iregion,imodel,ng),&
     &                                          percent
                  sumcpu=sumcpu+Csum(iregion,imodel,ng)
                  sumper=sumper+percent
                END IF
              END DO
              IF (Master.and.(total.gt.0.0_r8)) THEN
                WRITE (stdout,50) sumcpu, sumper
              END IF
            END IF
          END DO
# endif
          IF (Master) WRITE (stdout,60) Csum(region,MyModel,ng)
  60      FORMAT (/,' All percentages are with respect to total time =',&
     &            5x,f12.3,/)
#endif
        END IF
!$OMP END CRITICAL (FINALIZE_WCLOCK)
      END IF
      RETURN
      END SUBROUTINE wclock_off
