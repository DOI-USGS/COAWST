#include "cppdefs.h"
      MODULE ini_strengthcoef_mod
!
! svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine initializes pstar_grid based on the grid spacing.      !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ini_strengthcoef

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ini_strengthcoef (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
#ifdef ICE_MODEL
      USE mod_ice
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
#ifdef SOLVE3D
      real(r8) :: diffusion2(MT), diffusion4(MT)
#endif
      real(r8) :: viscosity2, viscosity4
!
#include "tile.h"

      CALL ini_strengthcoef_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            GRID(ng) % grdscl                     &
#ifdef ICE_STRENGTH_QUAD
     &                           ,ICE(ng) % pstar_grid                  &
#endif
     &                            )

      RETURN
      END SUBROUTINE ini_strengthcoef
!
!***********************************************************************
      SUBROUTINE ini_strengthcoef_tile (ng, tile, model,                &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  IminS, ImaxS, JminS, JmaxS,     &
     &                                  grdscl                          &
#ifdef ICE_STRENGTH_QUAD
     &                                  ,pstar_grid                     &
#endif
     &                              )
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: grdscl(LBi:,LBj:)
# ifdef ICE_STRENGTH_QUAD
      real(r8), intent(inout) :: pstar_grid(LBi:,LBj:)
# endif
#else
      real(r8), intent(in) :: grdscl(LBi:UBi,LBj:UBj)
# ifdef ICE_STRENGTH_QUAD
      real(r8), intent(inout) :: pstar_grid(LBi:UBi,LBj:UBj)
# endif
#endif
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j
      real(r8) :: cff

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set horizontal, constant, mixing coefficient according to model flag.
!-----------------------------------------------------------------------
!
#ifdef ICE_STRENGTH_QUAD
!
!-----------------------------------------------------------------------
!  Scale ice strength according to the grid size. Scale by 1 km.
!-----------------------------------------------------------------------
!
      cff=pstar_e(ng)/1000.0_r8
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          pstar_grid(i,j)=cff*grdscl(i,j)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pstar_grid)
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pstar_grid)
# endif
#endif

      RETURN
      END SUBROUTINE ini_strengthcoef_tile

      END MODULE ini_strengthcoef_mod
