#include "cppdefs.h"
      MODULE atm2ocn_flux_mod
#ifdef ATM2OCN_FLUXES
!
!svn $Id: atm2ocn_flux.F 732 2008-09-07 01:55:51Z jcwarner $
!================================================== John C. Warner   ===
!                                                   Hernan G. Arango   =
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the total heat fluxes.                        !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: atm2ocn_flux

      CONTAINS
!
!***********************************************************************
      SUBROUTINE atm2ocn_flux (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
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
      CALL wclock_on (ng, iNLM, 17)
# endif
      CALL atm2ocn_flux_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % umask,                         &
     &                        GRID(ng) % vmask,                         &
# endif
# ifdef EMINUSP
     &                        FORCES(ng) % EminusP,                     &
     &                        FORCES(ng) % evap,                        &
# endif
     &                        FORCES(ng) % rain,                        &
     &                        FORCES(ng) % lhflx,                       &
     &                        FORCES(ng) % lrflx,                       &
     &                        FORCES(ng) % shflx,                       &
     &                        FORCES(ng) % srflx,                       &
     &                        FORCES(ng) % stflx,                       &
     &                        FORCES(ng) % sustr,                       &
     &                        FORCES(ng) % svstr,                       &
     &                        FORCES(ng) % Taux,                        &
     &                        FORCES(ng) % Tauy)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 17)
# endif
      RETURN
      END SUBROUTINE atm2ocn_flux
!
!***********************************************************************
      SUBROUTINE atm2ocn_flux_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
# ifdef MASKING
     &                              rmask, umask, vmask,                &
# endif
# ifdef EMINUSP
     &                              EminusP, evap,                      &
# endif
     &                              rain, lhflx,                        &
     &                              lrflx, shflx,                       &
     &                              srflx, stflx,                       &
     &                              sustr, svstr,                       &
     &                              Taux, Tauy)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef EMINUSP
      real(r8), intent(out) :: EminusP(LBi:,LBj:)
      real(r8), intent(inout) :: evap(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: rain(LBi:,LBj:)
      real(r8), intent(inout) :: lhflx(LBi:,LBj:)
      real(r8), intent(inout) :: lrflx(LBi:,LBj:)
      real(r8), intent(inout) :: shflx(LBi:,LBj:)
      real(r8), intent(inout) :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
      real(r8), intent(inout) :: sustr(LBi:,LBj:)
      real(r8), intent(inout) :: svstr(LBi:,LBj:)
      real(r8), intent(in) :: Taux(LBi:,LBj:)
      real(r8), intent(in) :: Tauy(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef EMINUSP
      real(r8), intent(out) :: EminusP(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: evap(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: rain(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: lhflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: lrflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: shflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(inout) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Taux(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Tauy(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j

      real(r8) :: cff

# include "set_bounds.h"
!
!=======================================================================
!  Atmosphere-Ocean flux computations.
!=======================================================================
!
      cff=1.0_r8/rhow
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          stflx(i,j,itemp)=(srflx(i,j)+lrflx(i,j)+                      &
     &                      lhflx(i,j)+shflx(i,j))
# ifdef EMINUSP
          stflx(i,j,isalt)=cff*(evap(i,j)-rain(i,j))
# endif
# ifdef MASKING
          stflx(i,j,itemp)=stflx(i,j,itemp)*rmask(i,j)
          rain(i,j)=rain(i,j)*rmask(i,j)
#   ifdef EMINUSP
          evap(i,j)=evap(i,j)*rmask(i,j)
          stflx(i,j,isalt)=stflx(i,j,isalt)*rmask(i,j)
#   endif
# endif
# ifdef EMINUSP
          EminusP(i,j)=stflx(i,j,isalt)
# endif
        END DO
      END DO
!
!  Transfer Taux and Tauy from rho points to u and v points.
!
      cff=0.5_r8
      DO j=JstrR,JendR
        DO i=Istr,IendR
          sustr(i,j)=cff*(Taux(i-1,j)+Taux(i,j))
# ifdef MASKING
          sustr(i,j)=sustr(i,j)*umask(i,j)
# endif
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          svstr(i,j)=cff*(Tauy(i,j-1)+Tauy(i,j))
# ifdef MASKING
          svstr(i,j)=svstr(i,j)*vmask(i,j)
# endif
        END DO
      END DO
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          stflx(:,:,itemp))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rain)
# ifdef EMINUSP
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          evap)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          EminusP)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          stflx(:,:,isalt))
# endif
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sustr)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          svstr)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rain, stflx(:,:,itemp))
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    sustr, svstr)
#  ifdef EMINUSP
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    evap, EminusP,                                &
     &                    stflx(:,:,isalt))
#  endif
# endif
      RETURN
      END SUBROUTINE atm2ocn_flux_tile
#endif
      END MODULE atm2ocn_flux_mod
