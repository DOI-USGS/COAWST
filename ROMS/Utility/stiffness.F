#include "cppdefs.h"
      MODULE stiffness_mod
!
!svn $Id: stiffness.F 795 2016-05-11 01:42:43Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine surveys the 3D grid in order to determine maximum      !
!  grid stiffness ratio:                                               !
!                                                                      !
!             z(i,j,k)-z(i-1,j,k)+z(i,j,k-1)-z(i-1,j,k-1)              !
!      r_x = ---------------------------------------------             !
!             z(i,j,k)+z(i-1,j,k)-z(i,j,k-1)-z(i-1,j,k-1)              !
!                                                                      !
!  This is done for diagnostic purposes and it does not affect the     !
!  computations.                                                       !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: stiffness

      CONTAINS
!
!***********************************************************************
      SUBROUTINE stiffness (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
#include "tile.h"
!
      CALL stiffness_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
#ifdef MASKING
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
#endif
     &                     GRID(ng) % h,                                &
     &                     GRID(ng) % omn,                              &
#ifdef SOLVE3D
     &                     GRID(ng) % Hz,                               &
     &                     GRID(ng) % z_w,                              &
#else
# ifdef ICESHELF
     &                     GRID(ng) % zice,                             &
# endif
#endif 
     &                     OCEAN(ng)% zeta)
      RETURN
      END SUBROUTINE stiffness
!
!***********************************************************************
      SUBROUTINE stiffness_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
#ifdef MASKING
     &                           rmask, umask, vmask,                   &
#endif
     &                           h, omn,                                &
#ifdef SOLVE3D
     &                           Hz, z_w,                               &
#else
# ifdef ICESHELF
     &                           zice,                                  &
# endif
#endif
     &                           zeta)
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
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: omn(LBi:,LBj:)
# ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
# else
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#  endif
# endif
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
# ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
# else
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#  endif
# endif
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: NSUB, i, j, k

      real(r8) :: cff, ratio

#ifdef SOLVE3D
      real(r8) :: my_rx0, my_rx1
#endif
      real(r8) :: my_volume0, my_volume1, my_volume2

#ifdef DISTRIBUTE
# ifdef SOLVE3D
      real(r8), dimension(5) :: buffer
      character (len=3), dimension(5) :: op_handle
# else
      real(r8), dimension(3) :: buffer
      character (len=3), dimension(3) :: op_handle
# endif
#endif

#include "set_bounds.h"

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute grid stiffness.
!-----------------------------------------------------------------------
!
      my_rx0=0.0_r8
      my_rx1=0.0_r8
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
# ifdef MASKING
          IF (umask(i,j).gt.0.0_r8) THEN
# endif
            my_rx0=MAX(my_rx0,ABS((z_w(i,j,0)-z_w(i-1,j,0))/            &
     &                            (z_w(i,j,0)+z_w(i-1,j,0))))
            DO k=1,N(ng)
              my_rx1=MAX(my_rx1,ABS((z_w(i,j,k  )-z_w(i-1,j,k  )+       &
     &                               z_w(i,j,k-1)-z_w(i-1,j,k-1))/      &
     &                              (z_w(i,j,k  )+z_w(i-1,j,k  )-       &
     &                               z_w(i,j,k-1)-z_w(i-1,j,k-1))))
            END DO
# ifdef MASKING
          END IF
# endif
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
# ifdef MASKING
          IF (vmask(i,j).gt.0.0_r8) THEN
# endif
            my_rx0=MAX(my_rx0,ABS((z_w(i,j,0)-z_w(i,j-1,0))/            &
     &                            (z_w(i,j,0)+z_w(i,j-1,0))))
            DO k=1,N(ng)
              my_rx1=MAX(my_rx1,ABS((z_w(i,j,k  )-z_w(i,j-1,k  )+       &
     &                               z_w(i,j,k-1)-z_w(i,j-1,k-1))/      &
     &                              (z_w(i,j,k  )+z_w(i,j-1,k  )-       &
     &                               z_w(i,j,k-1)-z_w(i,j-1,k-1))))
            END DO
# ifdef MASKING
          END IF
# endif
        END DO
      END DO
#endif
!
!-------------------------------------------------------------------------
!  Compute initial basin volume and grid cell minimum and maximum volumes.
!-------------------------------------------------------------------------
!
      my_volume0=0.0_r8
      my_volume1=1.0E+20_r8
      my_volume2=0.0_r8
!
#ifdef SOLVE3D
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
# ifdef MASKING
            IF (rmask(i,j).gt.0.0_r8) THEN
# endif
              cff=omn(i,j)*Hz(i,j,k)
              my_volume0=my_volume0+cff
              my_volume1=MIN(my_volume1,cff)
              my_volume2=MAX(my_volume2,cff)
# ifdef MASKING
            END IF
# endif
          END DO
        END DO
      END DO
#else
      DO j=Jstr,Jend
        DO i=Istr,Iend
# ifdef MASKING
          IF (rmask(i,j).gt.0.0_r8) THEN
# endif
# ifdef ICESHELF
            cff=omn(i,j)*(zeta(i,j,1)+h(i,j)-ABS(zice(i,j)))
# else
            cff=omn(i,j)*(zeta(i,j,1)+h(i,j))
# endif
            my_volume0=my_volume0+cff
            my_volume1=MIN(my_volume1,cff)
            my_volume2=MAX(my_volume2,cff)
# ifdef MASKING
          END IF
# endif
        END DO
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Compute global values.
!-----------------------------------------------------------------------
!
#ifdef DISTRIBUTE
      NSUB=1                             ! distributed-memory
#else
      IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                        &
     &    DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
#endif
!$OMP CRITICAL (R_FACTOR)
      IF (tile_count.eq.0) THEN
        TotVolume=my_volume0
        MinVolume=my_volume1
        MaxVolume=my_volume2
#ifdef SOLVE3D
        rx0=my_rx0
        rx1=my_rx1
#endif
      ELSE
        TotVolume=TotVolume+my_volume0
        MinVolume=MIN(MinVolume,my_volume1)
        MaxVolume=MAX(MaxVolume,my_volume2)
#ifdef SOLVE3D
        rx0=MAX(rx0,my_rx0)
        rx1=MAX(rx1,my_rx1)
#endif
      END IF
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
#ifdef DISTRIBUTE
        buffer(1)=TotVolume
        buffer(2)=MinVolume
        buffer(3)=MaxVolume
# ifdef SOLVE3D
        buffer(4)=rx0
        buffer(5)=rx1
# endif
        op_handle(1)='SUM'
        op_handle(2)='MIN'
        op_handle(3)='MAX'
# ifdef SOLVE3D
        op_handle(4)='MAX'
        op_handle(5)='MAX'
# endif
# ifdef SOLVE3D
        CALL mp_reduce (ng, model, 5, buffer, op_handle)
# else
        CALL mp_reduce (ng, model, 3, buffer, op_handle)
# endif
        TotVolume=buffer(1)
        MinVolume=buffer(2)
        MaxVolume=buffer(3)
# ifdef SOLVE3D
        rx0=buffer(4)
        rx1=buffer(5)
# endif
#endif
        IF (Master) THEN
          WRITE (stdout,10) ng
  10      FORMAT (/,' Basin information for Grid ',i2.2,':',/)
#ifdef SOLVE3D
          WRITE (stdout,20) rx0, rx1
  20      FORMAT (' Maximum grid stiffness ratios:  rx0 = ',1pe14.6,    &
     &            ' (Beckmann and Haidvogel)',/,t34,'rx1 = ',1pe14.6,   &
     &            ' (Haney)')
#endif
          IF (MinVolume.ne.0.0_r8) THEN
            ratio=MaxVolume/MinVolume
          ELSE
            ratio=0.0_r8
          END IF
          WRITE (stdout,30) TotVolume, MinVolume, MaxVolume, ratio
  30      FORMAT (/,' Initial basin volumes: TotVolume = ',1p,e17.10,   &
     &            0p,' m3',/,t25,'MinVolume = ',1p,e17.10,0p,' m3',     &
     &            /,t25,'MaxVolume = ',1p,e17.10,0p,' m3',              &
     &            /,t25,'  Max/Min = ',1p,e17.10,0p)
        END IF
      END IF
!$OMP END CRITICAL (R_FACTOR)
      RETURN
      END SUBROUTINE stiffness_tile

      END MODULE stiffness_mod
