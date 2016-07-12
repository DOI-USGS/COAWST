#include "cppdefs.h"
      MODULE rp_diag_mod
#ifdef TL_IOMS
!
!svn $Id: rp_diag.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes various representer tangent linear diagnostic !
!  fields.                                                             !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: rp_diag

      CONTAINS
!
!***********************************************************************
      SUBROUTINE rp_diag (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iRPM, 7)
# endif
      CALL rp_diag_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp(ng), kstp(ng),                            &
     &                   GRID(ng) % h,                                  &
     &                   GRID(ng) % omn,                                &
# ifdef SOLVE3D
     &                   GRID(ng) % tl_Hz,                              &
     &                   GRID(ng) % tl_z_r,                             &
     &                   GRID(ng) % tl_z_w,                             &
     &                   OCEAN(ng) % tl_rho,                            &
     &                   OCEAN(ng) % tl_u,                              &
     &                   OCEAN(ng) % tl_v,                              &
# endif
     &                   OCEAN(ng) % tl_ubar,                           &
     &                   OCEAN(ng) % tl_vbar,                           &
     &                   OCEAN(ng) % tl_zeta)
# ifdef PROFILE
      CALL wclock_off (ng, iRPM, 7)
# endif
      RETURN
      END SUBROUTINE rp_diag
!
!***********************************************************************
      SUBROUTINE rp_diag_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nstp, kstp,                              &
     &                         h, omn,                                  &
# ifdef SOLVE3D
     &                         tl_Hz, tl_z_r, tl_z_w,                   &
     &                         tl_rho, tl_u, tl_v,                      &
# endif
     &                         tl_ubar, tl_vbar, tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, kstp
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: omn(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(in) :: tl_Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: tl_rho(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_v(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(in) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_zeta(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(in) :: tl_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: tl_rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  endif
      real(r8), intent(in) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: NSUB, i, j, k, trd

# ifdef DISTRIBUTE
      real(r8), dimension(3) :: buffer
      character (len=3), dimension(3) :: op_handle
# else
      integer :: my_threadnum
# endif

      real(r8) :: cff, my_avgke, my_avgpe, my_volume

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ke2d
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: pe2d

      character (len=8) :: kechar, pechar

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute and report out volume averaged kinetic, potential total
!  energy, and volume of tangent linear fields. Notice that the
!  proper adjoint of these quantities are not coded.
!-----------------------------------------------------------------------
!
      IF (MOD(iic(ng)-1,ninfo(ng)).eq.0) THEN
        DO j=Jstr,Jend
# ifdef SOLVE3D
          DO i=Istr,Iend
            ke2d(i,j)=0.0_r8
            pe2d(i,j)=0.5_r8*g*tl_z_w(i,j,N(ng))*tl_z_w(i,j,N(ng))
          END DO
          cff=g/rho0
          DO k=N(ng),1,-1
            DO i=Istr,Iend
              ke2d(i,j)=ke2d(i,j)+                                      &
     &                  tl_Hz(i,j,k)*                                   &
     &                  0.25_r8*(tl_u(i  ,j,k,nstp)*tl_u(i  ,j,k,nstp)+ &
     &                           tl_u(i+1,j,k,nstp)*tl_u(i+1,j,k,nstp)+ &
     &                           tl_v(i,j  ,k,nstp)*tl_v(i,j  ,k,nstp)+ &
     &                           tl_v(i,j+1,k,nstp)*tl_v(i,j+1,k,nstp))
              pe2d(i,j)=pe2d(i,j)+                                      &
     &                  cff*tl_Hz(i,j,k)*(tl_rho(i,j,k)+1000.0_r8)*     &
     &                  (tl_z_r(i,j,k)-tl_z_w(i,j,0))
            END DO
          END DO
# else
          cff=0.5_r8*g
          DO i=Istr,Iend
            ke2d(i,j)=(tl_zeta(i,j,kstp)+h(i,j))*                       &
     &                0.25_r8*(tl_ubar(i  ,j,kstp)*tl_ubar(i  ,j,kstp)+ &
     &                         tl_ubar(i+1,j,kstp)*tl_ubar(i+1,j,kstp)+ &
     &                         tl_vbar(i,j  ,kstp)*tl_vbar(i,j  ,kstp)+ &
     &                         tl_vbar(i,j+1,kstp)*tl_vbar(i,j+1,kstp))
            pe2d(i,j)=cff*tl_zeta(i,j,kstp)*tl_zeta(i,j,kstp)
          END DO
# endif
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
# ifdef SOLVE3D
!!          pe2d(i,Jend+1)=pe2d(i,Jend+1)+                              &
!!   &                     omn(i,j)*(z_w(i,j,N(ng))-z_w(i,j,0))
# else
!!          pe2d(i,Jend+1)=pe2d(i,Jend+1)+                              &
!!   &                     omn(i,j)*(zeta(i,j,kstp)+h(i,j))
# endif
            pe2d(i,Jend+1)=pe2d(i,Jend+1)+omn(i,j)*h(i,j)
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
# ifdef DISTRIBUTE
        NSUB=1                           ! distributed-memory
# else
        IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                      &
     &      DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          NSUB=1                         ! non-tiled application
        ELSE
          NSUB=NtileX(ng)*NtileE(ng)     ! tiled application
        END IF
# endif
!$OMP CRITICAL (RP_DIAGNOSTICS)
        IF (tile_count.eq.0) THEN
          volume=0.0_r8
          avgke=0.0_r8
          avgpe=0.0_r8
        END IF
        volume=volume+my_volume
        avgke=avgke+my_avgke
        avgpe=avgpe+my_avgpe
        tile_count=tile_count+1
        IF (tile_count.eq.NSUB) THEN
          tile_count=0
# ifdef DISTRIBUTE
          buffer(1)=volume
          buffer(2)=avgke
          buffer(3)=avgpe
          op_handle(1)='SUM'
          op_handle(2)='SUM'
          op_handle(3)='SUM'
          CALL mp_reduce (ng, iRPM, 3, buffer, op_handle)
          volume=buffer(1)
          avgke=buffer(2)
          avgpe=buffer(3)
          trd=MyMaster
# else
          trd=my_threadnum()
# endif
          avgke=avgke/volume
          avgpe=avgpe/volume
          avgkp=avgke+avgpe
          IF (first_time(ng).eq.0) THEN
            first_time(ng)=1
            IF (Master.and.(ng.eq.1)) THEN
              WRITE (stdout,10) 'STEP', 'Day HH:MM:SS',                 &
     &                          'KINETIC_ENRG', 'POTEN_ENRG',           &
# ifdef NESTING
     &                          'TOTAL_ENRG', 'NET_VOLUME', 'Grid'
# else
     &                          'TOTAL_ENRG', 'NET_VOLUME'
# endif
# ifdef NESTING
 10           FORMAT (/,3x,a,3x,a,2x,a,3x,a,4x,a,4x,a,2x,a)
# else
 10           FORMAT (/,3x,a,3x,a,2x,a,3x,a,4x,a,4x,a)
# endif
            END IF
          END IF
          IF (Master) THEN
            WRITE (stdout,20) iic(ng)-1, time_code(ng),                 &
# ifdef NESTING
     &                        avgke, avgpe, avgkp, volume, ng
# else
     &                        avgke, avgpe, avgkp, volume
# endif
# ifdef NESTING
 20         FORMAT (i7,1x,a,4(1pe14.6),2x,i2.2)
# else
 20         FORMAT (i7,1x,a,4(1pe14.6))
# endif
            CALL my_flush (stdout)
          END IF
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
        END IF
!$OMP END CRITICAL (RP_DIAGNOSTICS)
      END IF
      RETURN
      END SUBROUTINE rp_diag_tile
#endif
      END MODULE rp_diag_mod
