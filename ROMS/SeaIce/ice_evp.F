#include "cppdefs.h"
      MODULE ice_evp_mod
#if defined ICE_MOMENTUM && defined ICE_EVP
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes the parameters for elastic-viscous-plastic    !
!  rheology.                                                           !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC ice_evp

      CONTAINS

      SUBROUTINE ice_evp (ng, tile)

      USE mod_param
      USE mod_grid
      USE mod_ice
      USE mod_stepping
!
      implicit none
!
      integer, intent(in) :: ng, tile
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 45)
# endif
!
      CALL ice_evp_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   liold(ng), lieol(ng),                          &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
# ifdef WET_DRY
     &                   GRID(ng) % rmask_wet,                          &
# endif
     &                   GRID(ng) % pm,                                 &
     &                   GRID(ng) % pn,                                 &
     &                   ICE(ng) % uie,                                 &
     &                   ICE(ng) % vie,                                 &
     &                   ICE(ng) % ai,                                  &
     &                   ICE(ng) % hi,                                  &
     &                   ICE(ng) % pice,                                &
     &                   ICE(ng) % zetai,                               &
     &                   ICE(ng) % eta                                  &
     &                   )
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 45)
# endif
      RETURN
      END SUBROUTINE ice_evp
!
!***********************************************************************
      SUBROUTINE ice_evp_tile (ng, tile,                                &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        liold, lieol,                             &
# ifdef MASKING
     &                        rmask,                                    &
# endif
# ifdef WET_DRY
     &                        rmask_wet,                                &
# endif
     &                        pm, pn, uie, vie,                         &
     &                        ai, hi, pice, zetai, eta)
!***********************************************************************
!

      USE mod_param
      USE mod_scalars
!
      USE bc_2d_mod, ONLY : bc_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: liold, lieol
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: uie(LBi:,LBj:,:)
      real(r8), intent(in) :: vie(LBi:,LBj:,:)
      real(r8), intent(in) :: ai(LBi:,LBj:,:)
      real(r8), intent(in) :: hi(LBi:,LBj:,:)
      real(r8), intent(out) :: pice(LBi:,LBj:)
      real(r8), intent(out) :: zetai(LBi:,LBj:)
      real(r8), intent(out) :: eta(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: uie(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: vie(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: ai(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: hi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(out) :: pice(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: zetai(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: eta(LBi:UBi,LBj:UBj)
# endif
! Local variable definitions
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      integer :: i, j

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: eps11
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: eps22
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: eps12

      real(r8) :: eone
      real(r8) :: etwos
      real(r8) :: epx
      real(r8) :: epy
      real(r8) :: e2r
      real(r8) :: delta
      real(r8) :: zmax
      real(r8), parameter :: epso = 1.E-12_r8

# include "set_bounds.h"
! ------------------------------------------------------------
!

      e2r = 1.0_r8/(ellip_sq(ng))
!
! *** Compute strain rates
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
         eps11(i,j) = (uie(i+1,j,lieol)-uie(i,j,lieol))*pm(i,j)
         eps22(i,j) = (vie(i,j+1,lieol)-vie(i,j,lieol))*pn(i,j)
         epx = 0.25_r8*( vie(i+1,j+1,lieol)+vie(i+1,j,lieol)            &
     &             - vie(i-1,j+1,lieol)-vie(i-1,j,lieol) )*pm(i,j)
         epy = 0.25_r8*( uie(i+1,j+1,lieol)+uie(i,j+1,lieol)            &
     &             - uie(i+1,j-1,lieol)-uie(i,j-1,lieol) )*pn(i,j)
         eps12(i,j) = 0.5_r8*(epx + epy)
!
         eone=eps11(i,j)+eps22(i,j) 
         etwos=(eps11(i,j)-eps22(i,j))*(eps11(i,j)-eps22(i,j))+         &
     &         4.0_r8*eps12(i,j)*eps12(i,j)
! 
         delta=abs(eone**2+e2r*etwos)
         delta=max(sqrt(delta),epso)
         pice(i,j)=pstar(ng)*hi(i,j,liold)                              &
     &                *exp(-astren(ng)*(1.0_r8-ai(i,j,liold)))
         zetai(i,j)=pice(i,j)/(2.0_r8*delta)
         zmax = 2.5E+8_r8*pice(i,j)
         zetai(i,j)= min(zetai(i,j),zetamax(ng))
         zetai(i,j)= max(zetai(i,j),zetamin(ng))
         eta(i,j)=e2r*zetai(i,j)
        ENDDO
       ENDDO

        CALL bc_r2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    pice)
        CALL bc_r2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    zetai)
        CALL bc_r2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    eta)
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    pice, zetai, eta)
#endif

      RETURN
      END SUBROUTINE ice_evp_tile

#endif
      END MODULE ice_evp_mod
