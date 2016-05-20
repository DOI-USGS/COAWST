#include "cppdefs.h"
      MODULE ice_evp_sig_mod
#if defined ICE_MOMENTUM && defined ICE_EVP
!
!================================================ W. Paul Budgell ======
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine time steps the EVP stresses                            !
!  term.                                                               !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC ice_evp_sig

      CONTAINS

      SUBROUTINE ice_evp_sig (ng, tile)

      USE mod_param
      USE mod_grid
      USE mod_ice
      USE mod_stepping

      integer, intent(in) :: ng, tile
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 55)
# endif
!
      CALL ice_evp_sig_tile (ng, tile,                                  &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold(ng), lieol(ng), lienw(ng),            &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
     &                      ICE(ng) % uie,                              &
     &                      ICE(ng) % vie,                              &
     &                      ICE(ng) % hi,                               &
     &                      ICE(ng) % pice,                             &
     &                      ICE(ng) % zetai,                            &
     &                      ICE(ng) % eta,                              &
     &                      ICE(ng) % sig11,                            &
     &                      ICE(ng) % sig22,                            &
     &                      ICE(ng) % sig12                             &
     &                      )
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 55)
# endif
      RETURN
      END SUBROUTINE ice_evp_sig
!
!***********************************************************************
      SUBROUTINE ice_evp_sig_tile (ng, tile,                            &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        liold, lieol, lienw,                      &
     &                        pm, pn,                                   &
     &                        uie, vie,                                 &
     &                        hi, pice, zetai, eta,                     &
     &                        sig11, sig22, sig12                       &
     &                        ) 
!***********************************************************************
!

      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_boundary
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
      USE i2d_bc_mod, ONLY : i2d_bc_tile
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: liold, lieol, lienw

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: uie(LBi:,LBj:,:)
      real(r8), intent(in) :: vie(LBi:,LBj:,:)
      real(r8), intent(in) :: hi(LBi:,LBj:,:)
      real(r8), intent(in) :: pice(LBi:,LBj:)
      real(r8), intent(in) :: zetai(LBi:,LBj:)
      real(r8), intent(in) :: eta(LBi:,LBj:)
      real(r8), intent(inout) :: sig11(LBi:,LBj:,:)
      real(r8), intent(inout) :: sig22(LBi:,LBj:,:)
      real(r8), intent(inout) :: sig12(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: uie(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: vie(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: hi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: pice(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: zetai(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: eta(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: sig11(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: sig22(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: sig12(LBi:UBi,LBj:UBj,2)
# endif

! Local variable definitions
!
      integer :: i, j

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: eps11
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: eps22
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: eps12

      real(r8) :: e
      real(r8) :: e0
      real(r8) :: ep
      real(r8) :: ee
      real(r8) :: ees 
      real(r8) :: gamma
      real(r8) :: f1
      real(r8) :: f2
      real(r8) :: f3 
      real(r8) :: alfa
      real(r8) :: beta
      real(r8) :: epx
      real(r8) :: epy

# include "set_bounds.h"

!----------------
!
!.....initial value for Youngs modulus (between 0 and 1) 
!
      e0=0.25_r8

!
!...........stress tensor 
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          eps11(i,j) = (uie(i+1,j,lieol)-uie(i,j,lieol))*pm(i,j)
          eps22(i,j) = (vie(i,j+1,lieol)-vie(i,j,lieol))*pn(i,j)
          epx = 0.25_r8*( vie(i+1,j+1,lieol)+vie(i+1,j,lieol)           &
     &             - vie(i-1,j+1,lieol)-vie(i-1,j,lieol) )*pm(i,j)
          epy = 0.25_r8*( uie(i+1,j+1,lieol)+uie(i,j+1,lieol)           &
     &             - uie(i+1,j-1,lieol)-uie(i,j-1,lieol) )*pn(i,j)
          eps12(i,j) = 0.5_r8*(epx + epy)
        END DO
      END DO
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
!
          IF (hi(i,j,liold).gt.0.01_r8) THEN 
            e = 2.0_r8*e0*rhoice(ng)*hi(i,j,liold)/(pm(i,j)*dte(ng))**2 
            ep = e*pice(i,j)/(4.0_r8*zetai(i,j)+1.0E-8_r8) 
            ee = e/(2.0_r8*eta(i,j)+1.0E-8_r8) 
            ees= e*(eta(i,j)-zetai(i,j))/                               &
     &             (4.0_r8*eta(i,j)*zetai(i,j)+1.0E-8_r8) 
!
            alfa = 1.0_r8/dte(ng) + ee + ees 
            beta = ees 
            gamma = 1/dte(ng) + ee 
            f1 = e*eps11(i,j) - ep + sig11(i,j,lieol)/dte(ng) 
            f2 = e*eps22(i,j) - ep + sig22(i,j,lieol)/dte(ng) 
            f3 = e*eps12(i,j) + sig12(i,j,lieol)/dte(ng) 
            sig11(i,j,lienw) = (alfa*f1 - beta*f2)/(alfa**2 - beta**2) 
            sig12(i,j,lienw) = f3/gamma 
            sig22(i,j,lienw) = (alfa*f2 - beta*f1)/(alfa**2 - beta**2) 
!
          ELSE
! 
            sig11(i,j,lienw) = 2.0_r8*eta(i,j)*eps11(i,j)+              &
     &        (zetai(i,j)-eta(i,j))*                                    &
     &        (eps11(i,j)+eps22(i,j))  - pice(i,j)/2.0_r8
!     
            sig22(i,j,lienw) = 2.0_r8*eta(i,j)*eps22(i,j)+              &
     &        (zetai(i,j)-eta(i,j))*                                    &
     &        (eps11(i,j)+eps22(i,j))  - pice(i,j)/2.0_r8 
!
            sig12(i,j,lienw) = 2.0_r8*eta(i,j)*eps12(i,j)
!
          END IF 
        END DO 
      END DO

      CALL i2d_bc_tile (ng, tile, iNLM,                                 &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  lieol, lienw,                                   &
     &                  BOUNDARY(ng)%sig11_west(LBj:UBj),               &
     &                  BOUNDARY(ng)%sig11_east(LBj:UBj),               &
     &                  BOUNDARY(ng)%sig11_north(LBi:UBi),              &
     &                  BOUNDARY(ng)%sig11_south(LBi:UBi),              &
     &                  uie, vie, sig11, LBC(:,isSig11,ng))
      CALL i2d_bc_tile (ng, tile, iNLM,                                 &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  lieol, lienw,                                   &
     &                  BOUNDARY(ng)%sig22_west(LBj:UBj),               &
     &                  BOUNDARY(ng)%sig22_east(LBj:UBj),               &
     &                  BOUNDARY(ng)%sig22_north(LBi:UBi),              &
     &                  BOUNDARY(ng)%sig22_south(LBi:UBi),              &
     &                  uie, vie, sig22, LBC(:,isSig22,ng))
      CALL i2d_bc_tile (ng, tile, iNLM,                                 &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  lieol, lienw,                                   &
     &                  BOUNDARY(ng)%sig12_west(LBj:UBj),               &
     &                  BOUNDARY(ng)%sig12_east(LBj:UBj),               &
     &                  BOUNDARY(ng)%sig12_north(LBi:UBi),              &
     &                  BOUNDARY(ng)%sig12_south(LBi:UBi),              &
     &                  uie, vie, sig12, LBC(:,isSig12,ng))
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sig11(:,:,lienw))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sig22(:,:,lienw))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sig12(:,:,lienw))
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    sig11(:,:,lienw), sig22(:,:,lienw),           &
     &                    sig12(:,:,lienw))
# endif
!
      RETURN
      END SUBROUTINE ice_evp_sig_tile
#endif
      END MODULE ice_evp_sig_mod
