      SUBROUTINE ana_vmix (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets vertical mixing coefficients for momentum "Akv"   !
!  and tracers "Akt" (m2/s) using analytical expressions.              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_vmix_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    knew(ng),                                     &
     &                    GRID(ng) % h,                                 &
     &                    GRID(ng) % z_r,                               &
     &                    GRID(ng) % z_w,                               &
     &                    OCEAN(ng) % zeta,                             &
     &                    MIXING(ng) % Akv,                             &
     &                    MIXING(ng) % Akt)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(35)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_vmix
!
!***********************************************************************
      SUBROUTINE ana_vmix_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          knew,                                   &
     &                          h, z_r, z_w, zeta, Akv, Akt)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_3d_mod, ONLY : exchange_w3d_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: knew
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(out) :: Akv(LBi:,LBj:,0:)
      real(r8), intent(out) :: Akt(LBi:,LBj:,0:,:)
#else
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(out) :: Akv(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(out) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
#endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set vertical viscosity coefficient (m2/s).
!-----------------------------------------------------------------------
!
#if defined CANYON
      DO k=1,N(ng)-1
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Akv(i,j,k)=1.0E-03_r8+95.0E-04_r8*EXP(z_w(i,j,k)/50.0_r8)+  &
     &                 95.0E-04_r8*EXP(-(z_w(i,j,k)+h(i,j))/50.0_r8)
          END DO
        END DO
      END DO
#elif defined CHANNEL_NECK
      DO k=1,N(ng)-1
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Akv(i,j,k)=2.0E-04_r8+8.0E-04_r8*EXP(z_w(i,j,k)/5.0_r8)
          END DO
        END DO
      END DO
#elif defined COUPLING_TEST
      DO k=1,N(ng)-1
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Akv(i,j,k)=2.0E-03_r8+8.0E-03_r8*EXP(z_w(i,j,k)/1500.0_r8)
          END DO
        END DO
      END DO
#elif defined ESTUARY_TEST
      DO k=1,N(ng)-1
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Akv(i,j,k)=0.002_r8
          END DO
        END DO
      END DO
#elif defined LAKE_SIGNELL
      DO k=1,N(ng)-1
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Akv(i,j,k)=0.0005_r8
          END DO
        END DO
      END DO
#elif defined NJ_BIGHT
      DO k=1,N(ng)-1
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Akv(i,j,k)=1.0E-03_r8+2.0E-04_r8*EXP(z_r(i,j,k)/10.0_r8)
          END DO
        END DO
      END DO
#elif defined SED_TEST1
      DO k=1,N(ng)-1                         !  vonkar*ustar*z*(1-z/D)
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Akv(i,j,k)=0.025_r8*(h(i,j)+z_w(i,j,k))*                    &
     &                 (1.0_r8-(h(i,j)+z_w(i,j,k))/                     &
     &                  (h(i,j)+zeta(i,j,knew)))
            Akt(i,j,k,itemp)=Akv(i,j,k)*0.49_r8/0.39_r8
            Akt(i,j,k,isalt)=Akt(i,j,k,itemp)
          END DO
        END DO
      END DO
#elif defined SED_TOY
      DO k=1,N(ng)-1                         !  vonkar*ustar*z*(1-z/D)
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Akv(i,j,k)=0.41_r8*0.01_r8*(h(i,j)+z_w(i,j,k))*             &
     &                 (1.0_r8-(h(i,j)+z_w(i,j,k))/                     &
     &                  (h(i,j)+zeta(i,j,knew)))
          END DO
        END DO
      END DO
#elif defined SHOREFACE
      DO k=1,N(ng)-1
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Akv(i,j,k)=0.025_r8*(h(i,j)+z_w(i,j,k))*                    &
     &                 (1.0_r8-(h(i,j)+z_w(i,j,k))/                     &
     &                  (h(i,j)+zeta(i,j,knew)))
          END DO
        END DO
      END DO
#elif defined TEST_CHAN
      DO k=1,N(ng)-1                         !  vonkar*ustar*z*(1-z/D)
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Akv(i,j,k)=0.41_r8*0.0625_r8*(h(i,j)+z_w(i,j,k))*           &
     &                 (1.0_r8-(h(i,j)+z_w(i,j,k))/                     &
     &                  (h(i,j)+zeta(i,j,knew)))
          END DO
        END DO
      END DO
#elif defined UPWELLING
      DO k=1,N(ng)-1
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Akv(i,j,k)=2.0E-03_r8+8.0E-03_r8*EXP(z_w(i,j,k)/150.0_r8)
          END DO
        END DO
      END DO
#else
      ana_vmix.h: no values provided for Akv.
#endif
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          Akv)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Akv)
#endif
!
!-----------------------------------------------------------------------
!  Set vertical diffusion coefficient (m2/s).
!-----------------------------------------------------------------------
!
#if defined CANYON
      DO k=1,N(ng)-1
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Akt(i,j,k,itemp)=Akt_bak(itemp,ng)
          END DO
        END DO
      END DO
#elif defined CHANNEL_NECK
      DO k=1,N(ng)-1
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Akt(i,j,k,itemp)=2.0E-06_r8+                                &
     &                       8.0E-06_r8*EXP(z_w(i,j,k)/5.0_r8)
          END DO
        END DO
      END DO
#elif defined COUPLING_TEST
      DO k=1,N(ng)-1
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Akt(i,j,k,itemp)=Akt_bak(itemp,ng)
            Akt(i,j,k,isalt)=Akt_bak(isalt,ng)
          END DO
        END DO
      END DO
#elif defined ESTUARY_TEST
      DO k=1,N(ng)-1
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Akt(i,j,k,itemp)=Akv(i,j,k)
            Akt(i,j,k,isalt)=Akv(i,j,k)
          END DO
        END DO
      END DO
#elif defined LAKE_SIGNELL
      DO k=1,N(ng)-1
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Akt(i,j,k,itemp)=Akt_bak(itemp,ng)
            Akt(i,j,k,isalt)=Akt_bak(isalt,ng)
          END DO
        END DO
      END DO
#elif defined NJ_BIGHT
      DO k=1,N(ng)-1
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Akt(i,j,k,itemp)=1.0E-05_r8+                                &
     &                       2.0E-06_r8*EXP(z_r(i,j,k)/10.0_r8)
            Akt(i,j,k,isalt)=Akt(i,j,k,itemp)
          END DO
        END DO
      END DO
#elif defined SED_TOY
      DO k=1,N(ng)-1                         !  vonkar*ustar*z*(1-z/D)
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Akt(i,j,k,itemp)=Akv(i,j,k)
            Akt(i,j,k,isalt)=Akv(i,j,k)
          END DO
        END DO
      END DO
#elif defined SHOREFACE
      DO k=1,N(ng)-1
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Akt(i,j,k,itemp)=Akv(i,j,k)
            Akt(i,j,k,isalt)=Akv(i,j,k)
          END DO
        END DO
      END DO
#elif defined TEST_CHAN
      DO k=1,N(ng)-1
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Akt(i,j,k,itemp)=Akv(i,j,k)*0.49_r8/0.39_r8
            Akt(i,j,k,isalt)=Akt(i,j,k,itemp)
          END DO
        END DO
      END DO
#elif defined UPWELLING
      DO k=1,N(ng)-1
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Akt(i,j,k,itemp)=Akt_bak(itemp,ng)
            Akt(i,j,k,isalt)=Akt_bak(isalt,ng)
          END DO
        END DO
      END DO
#else
      ana_vmix.h: no values provided for Akt.
#endif
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        DO itrc=1,NAT
          CALL exchange_w3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 0, N(ng),         &
     &                            Akt(:,:,:,itrc))
        END DO
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 0, N(ng), 1, NAT,         &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Akt)
#endif

      RETURN
      END SUBROUTINE ana_vmix_tile
