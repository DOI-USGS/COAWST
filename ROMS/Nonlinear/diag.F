#include "cppdefs.h"
      MODULE diag_mod
!
!svn $Id: diag.F 790 2008-10-14 19:32:07Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes various diagnostic fields.                    !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: diag

      CONTAINS
!
!***********************************************************************
      SUBROUTINE diag (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
#ifdef ANA_DIAG
      USE analytical_mod, ONLY : ana_diag
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 7)
#endif
      CALL diag_tile (ng, tile,                                         &
     &                LBi, UBi, LBj, UBj,                               &
     &                nstp(ng), krhs(ng),                               &
     &                GRID(ng) % h,                                     &
     &                GRID(ng) % omn,                                   &
#ifdef SOLVE3D
     &                GRID(ng) % Hz,                                    &
     &                GRID(ng) % z_r,                                   &
     &                GRID(ng) % z_w,                                   &
     &                OCEAN(ng) % rho,                                  &
     &                OCEAN(ng) % u,                                    &
     &                OCEAN(ng) % v,                                    &
#endif
     &                OCEAN(ng) % ubar,                                 &
     &                OCEAN(ng) % vbar,                                 &
     &                OCEAN(ng) % zeta)
#ifdef ANA_DIAG
      CALL ana_diag (ng, tile, iNLM)
#endif
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 7)
#endif
      RETURN
      END SUBROUTINE diag
!
!***********************************************************************
      SUBROUTINE diag_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      nstp, krhs,                                 &
     &                      h, omn,                                     &
#ifdef SOLVE3D
     &                      Hz, z_r, z_w,                               &
     &                      rho, u, v,                                  &
#endif
     &                      ubar, vbar, zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nstp, krhs

#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: omn(LBi:,LBj:)
# ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
# endif
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#else
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
# ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
# endif
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: NSUB, i, j, k, trd
#ifdef DISTRIBUTE
# ifdef SOLVE3D
      integer, parameter :: Nreduce = 5
# else
      integer, parameter :: Nreduce = 4
# endif
      real(r8), dimension(Nreduce) :: buffer
      character (len=3), dimension(Nreduce) :: op_handle
#else
      integer :: my_threadnum      
#endif

      real(r8) :: cff, my_avgke, my_avgpe, my_volume
      real(r8) :: my_maxspeed, u2v2
#ifdef SOLVE3D
      real(r8) :: my_maxrho
#endif
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: ke2d
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: pe2d

      character (len=8) :: kechar, pechar

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute and report out volume averaged kinetic, potential
!  total energy, and volume.
!-----------------------------------------------------------------------
!
      IF (MOD(iic(ng)-1,ninfo(ng)).eq.0) THEN
        my_maxspeed=0.0_r8
#ifdef SOLVE3D
        my_maxrho=-1.0E+20_r8
#endif
        DO j=Jstr,Jend
#ifdef SOLVE3D
          DO i=Istr,Iend
            ke2d(i,j)=0.0_r8
            pe2d(i,j)=0.5_r8*g*z_w(i,j,N(ng))*z_w(i,j,N(ng))
          END DO
          cff=g/rho0
          DO k=N(ng),1,-1
            DO i=Istr,Iend
              u2v2=u(i  ,j,k,nstp)*u(i  ,j,k,nstp)+                     &
     &             u(i+1,j,k,nstp)*u(i+1,j,k,nstp)+                     &
     &             v(i,j  ,k,nstp)*v(i,j  ,k,nstp)+                     &
     &             v(i,j+1,k,nstp)*v(i,j+1,k,nstp)
              ke2d(i,j)=ke2d(i,j)+                                      &
     &                  Hz(i,j,k)*0.25_r8*u2v2
              pe2d(i,j)=pe2d(i,j)+                                      &
     &                  cff*Hz(i,j,k)*(rho(i,j,k)+1000.0_r8)*           &
     &                  (z_r(i,j,k)-z_w(i,j,0))
              my_maxspeed=MAX(my_maxspeed,SQRT(0.5_r8*u2v2))
              my_maxrho=MAX(my_maxrho,rho(i,j,k))
            END DO
          END DO
#else
          cff=0.5_r8*g
          DO i=Istr,Iend
            u2v2=ubar(i  ,j,krhs)*ubar(i  ,j,krhs)+                     &
     &           ubar(i+1,j,krhs)*ubar(i+1,j,krhs)+                     &
     &           vbar(i,j  ,krhs)*vbar(i,j  ,krhs)+                     &
     &           vbar(i,j+1,krhs)*vbar(i,j+1,krhs)
            ke2d(i,j)=(zeta(i,j,krhs)+h(i,j))*0.25_r8*u2v2
            pe2d(i,j)=cff*zeta(i,j,krhs)*zeta(i,j,krhs)
            my_maxspeed=MAX(my_maxspeed,SQRT(0.5_r8*u2v2))
          END DO
#endif
        END DO
!
!  Integrate horizontally within one tile. In order to reduce the
!  round-off errors, the summation is performed in two stages. First,
!  the index j is collapsed and then the accumulation is carried out
!  along index i. In this order, the partial sums consist on much
!  fewer number of terms than in a straightforward two-dimensional
!  summation. Thus, adding numbers which are orders of magnitude
!  apart is avoided.
!
        DO i=Istr,Iend
          pe2d(i,Jend+1)=0.0_r8
          pe2d(i,Jstr-1)=0.0_r8
          ke2d(i,Jstr-1)=0.0_r8
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend
#ifdef SOLVE3D
            pe2d(i,Jend+1)=pe2d(i,Jend+1)+                              &
     &                     omn(i,j)*(z_w(i,j,N(ng))-z_w(i,j,0))
#else
            pe2d(i,Jend+1)=pe2d(i,Jend+1)+                              &
     &                     omn(i,j)*(zeta(i,j,krhs)+h(i,j))
#endif
            pe2d(i,Jstr-1)=pe2d(i,Jstr-1)+omn(i,j)*pe2d(i,j)
            ke2d(i,Jstr-1)=ke2d(i,Jstr-1)+omn(i,j)*ke2d(i,j)
          END DO
        END DO
        my_volume=0.0_r8
        my_avgpe=0.0_r8
        my_avgke=0.0_r8
        DO i=Istr,Iend
          my_volume=my_volume+pe2d(i,Jend+1)
          my_avgpe =my_avgpe +pe2d(i,Jstr-1)
          my_avgke =my_avgke +ke2d(i,Jstr-1)
        END DO
!
!  Perform global summation: whoever gets first to the critical region
!  resets global sums before global summation starts; after the global
!  summation is completed, thread, which is the last one to enter the
!  critical region, finalizes the computation of diagnostics and prints
!  them out.
!
        IF (SOUTH_WEST_CORNER.and.                                      &
     &      NORTH_EAST_CORNER) THEN
          NSUB=1                         ! non-tiled application
        ELSE
          NSUB=NtileX(ng)*NtileE(ng)     ! tiled application
        END IF
!$OMP CRITICAL (NL_DIAGNOSTICS)
        IF (tile_count.eq.0) THEN
          volume=0.0_r8
          avgke=0.0_r8
          avgpe=0.0_r8
          maxspeed(ng)=0.0_r8
#ifdef SOLVE3D
          maxrho(ng)=0.0_r8
#endif
        END IF
        volume=volume+my_volume
        avgke=avgke+my_avgke
        avgpe=avgpe+my_avgpe
        maxspeed(ng)=MAX(maxspeed(ng),my_maxspeed)
#ifdef SOLVE3D
        maxrho(ng)=MAX(maxrho(ng),my_maxrho)
#endif
        tile_count=tile_count+1
        IF (tile_count.eq.NSUB) THEN
          tile_count=0
#ifdef DISTRIBUTE
          buffer(1)=volume
          buffer(2)=avgke
          buffer(3)=avgpe
          buffer(4)=maxspeed(ng)
# ifdef SOLVE3D
          buffer(5)=maxrho(ng)
# endif
          op_handle(1)='SUM'
          op_handle(2)='SUM'
          op_handle(3)='SUM'
          op_handle(4)='MAX'
# ifdef SOLVE3D
          op_handle(5)='MAX'
# endif
          CALL mp_reduce (ng, iNLM, Nreduce, buffer, op_handle)
          volume=buffer(1)
          avgke=buffer(2)
          avgpe=buffer(3)
          maxspeed(ng)=buffer(4)
# ifdef SOLVE3D
          maxrho(ng)=buffer(5)
# endif
          trd=MyMaster
#else
          trd=my_threadnum()
#endif
          avgke=avgke/volume
          avgpe=avgpe/volume
          avgkp=avgke+avgpe
          IF (first_time(ng).eq.0) THEN
            first_time(ng)=1
            IF (Master) WRITE (stdout,10) 'GRID','STEP', 'Day HH:MM:SS',&
     &                                    'KINETIC_ENRG', 'POTEN_ENRG', &
     &                                    'TOTAL_ENRG', 'NET_VOLUME'
 10         FORMAT (/,2x,a,2x,a,3x,a,2x,a,3x,a,4x,a,4x,a,/)
          END IF
          IF (Master) WRITE (stdout,20) ng,iic(ng)-1, time_code(ng),    &
     &                                  avgke, avgpe, avgkp, volume
 20       FORMAT (i4,1x,i7,1x,a,4(1pe14.6))
          IF (Master) CALL my_flush (stdout)
!
!  If blowing-up, set exit_flag to stop computations.
!
          WRITE (kechar,'(1pe8.1)') avgke
          WRITE (pechar,'(1pe8.1)') avgpe
          DO i=1,8
            IF ((kechar(i:i).eq.'N').or.(pechar(i:i).eq.'N').or.        &
     &          (kechar(i:i).eq.'n').or.(pechar(i:i).eq.'n').or.        &
     &          (kechar(i:i).eq.'*').or.(pechar(i:i).eq.'*')) THEN
              exit_flag=1
            END IF
          END DO
!
!  Stop computations if exceeding maximum speed allowed.  This will be
!  useful during debugging to avoid NaNs in output NetCDF files.
!
          IF (maxspeed(ng).gt.max_speed) THEN
            exit_flag=1
          END IF
#ifdef SOLVE3D
!
!  Stop computation if exceeding maximum density anomaly allowed.
!  Recall that  density is computed from potential temperature and
!  salinity.  This is a good way to screen for very bad values which
!  indicates that the model is blowing-up.
!
          IF (maxrho(ng).gt.max_rho) THEN
            exit_flag=1
          END IF
#endif

        END IF
!$OMP END CRITICAL (NL_DIAGNOSTICS)
      END IF
      RETURN
      END SUBROUTINE diag_tile
      END MODULE diag_mod
