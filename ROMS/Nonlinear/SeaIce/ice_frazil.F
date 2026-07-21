#include "cppdefs.h"
      MODULE ice_frazil_mod

#if defined ICE_MODEL && defined ICE_THERMO
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group           Paul Budgell       !
!    Licensed under a MIT/X style license           Katherine Hedstrom !
!    See License_ROMS.md                            Scott M. Durski    !
!=======================================================================
!                                                                      !
!  This routine computes the frazil ice growth in the water when the   !
!  water temperature gets below freezing. It adjusts both the water    !
!  temperature and salinity accordingly.                               !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Steele, M., G.L. Mellor, and M.G. Mcphee, 1989: Role of the       !
!      molecular sublayer in the melting or freezing of sea ice,       !
!      J. Phys. Oceanogr., 19, 139-147.                                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ice
      USE mod_ocean
      USE mod_scalars
!
      USE bc_2d_mod,       ONLY : bc_r2d_tile
# ifdef DISTRIBUTE
      USE distribute_mod,  ONLY : mp_reduce
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
      implicit none
!
      PRIVATE
      PUBLIC  :: freezing_point
      PUBLIC  :: ice_frazil
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_frazil (ng, tile, model)
!***********************************************************************
!
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 42, __LINE__, MyFile)
# endif
      CALL ice_frazil_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nnew(ng),                                   &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
# endif
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % z_r,                             &
     &                      OCEAN(ng) % rho,                            &
     &                      OCEAN(ng) % t,                              &
     &                      ICE(ng) % Fi)
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE ice_frazil
!
!***********************************************************************
      SUBROUTINE ice_frazil_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nnew,                                 &
# ifdef MASKING
     &                            rmask,                                &
# endif
# ifdef WET_DRY
     &                            rmask_wet,                            &
# endif
     &                            Hz, z_r, rho, t,                      &
     &                            Fi)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nnew

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: Fi(LBi:,LBj:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: Fi(LBi:UBi,LBj:UBj,nIceF)
# endif
!
! Local variable definitions
!
      logical :: overLand
!
      integer :: i, j, k, itrc
!
      real(r8), parameter :: Lhat = 79.2_r8
      real(r8), parameter :: r = 0.5_r8

      real(r8) :: delta_wfr, gamma_k, orhoi, pfac, t_fr

# ifdef DISTRIBUTE
      real(r8), allocatable           :: buffer(:)
      character (len=3), allocatable  :: op_handle(:)
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute frazil ice growth.
!-----------------------------------------------------------------------
!
!  Initialize rate of frazil ice growth (m3/s).
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Fi(i,j,icW_fr)=0.0_r8
        END DO
      END DO
!
!  Original formulation.
!
      pfac=1.0_r8                                    ! adjusting factor
      DO j=Jstr,Jend
        DO i=Istr,Iend
          DO k=1,N(ng)
# ifdef MASKING
            overLand=rmask(i,j).lt.1.0_r8
#  ifdef WET_DRY
            overLand=overland.or.(rmask_wet(i,j).lt.1.0_r8)
#  endif
# else
            overLand=.FALSE.
# endif
            IF (.not.overLand) THEN
              orhoi=1.0_r8/IceRho(ng)
              t_fr=freezing_point(t(i,j,k,nnew,isalt), z_r(i,j,k))
              IF (t(i,j,k,nnew,itemp).lt.t_fr) THEN
                gamma_k=pfac*(t_fr-t(i,j,k,nnew,itemp))/                &
     &                       (Lhat+                                     &
                              t(i,j,k,nnew,itemp)*(1.0_r8-r)+           &
     &                        0.0543_r8*t(i,j,k,nnew,isalt))
                IF ((gamma_k.lt.0.0_r8).and.(k.eq.N(ng))) THEN
                  gamma_k=0.0_r8
                END IF
                Fi(i,j,icW_fr)=Fi(i,j,icW_fr)+                          &
     &                         gamma_k*Hz(i,j,k)*                       &
     &                         (1000.0_r8+rho(i,j,k))*orhoi
                t(i,j,k,nnew,itemp)=t(i,j,k,nnew,itemp)+                &
     &                              gamma_k*                            &
     &                              (Lhat+                              &
     &                               t(i,j,k,nnew,itemp)*(1.0_r8-r))
!
!  Use heat at this level to melt some ice from below (gamma_k becomes
!  negative here). The salinity is not adjusted since salt flux from
!  frazil ice is considered a surface flux (SMD).
!
              ELSE IF ((Fi(i,j,icW_fr).gt.0.0_r8).and.                  &
     &                 (t(i,j,k,nnew,itemp).gt.t_fr)) THEN
                gamma_k=pfac*(t_fr-t(i,j,k,nnew,itemp))/                &
     &                       (Lhat+t(i,j,k,nnew,itemp)*(1.0_r8-r)+      &
     &                        0.0543_r8*t(i,j,k,nnew,isalt))
                delta_wfr=gamma_k*Hz(i,j,k)*(rho0+rho(i,j,k))*orhoi
                IF ((Fi(i,j,icW_fr)+delta_wfr).gt.0.0_r8) THEN
                  Fi(i,j,icW_fr)=Fi(i,j,icW_fr)+delta_wfr
                ELSE
                  gamma_k=-Fi(i,j,icW_fr)*                              &
     &                    IceRho(ng)/(Hz(i,j,k)*(rho0+rho(i,j,k)))
                  Fi(i,j,icW_fr)=0.0_r8
                ENDIF
                t(i,j,k,nnew,itemp)=t(i,j,k,nnew,itemp)+                &
     &                              gamma_k*                            &
     &                              (Lhat+                              &
     &                               t(i,j,k,nnew,itemp)*(1.0_r8-r))
              END IF
            END IF
          END DO
!
!  Compute rate of frazil ice growth.
!
          Fi(i,j,icW_fr)=Fi(i,j,icW_fr)/dt(ng)
!
!  If negative frazil growth, set it to zero.  It implies that the melt
!  effect exceeded the freezing effect (SMD).
!
          IF (Fi(i,j,icW_fr).lt.0.0_r8) THEN
            Fi(i,j,icW_fr)=0.0_r8
          END IF
        END DO
      END DO
!
!  Set lateral boundary conditions.
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fi(:,:,icW_fr))

# ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Fi(:,:,icW_fr))
# endif
!
      RETURN
      END SUBROUTINE ice_frazil_tile
!
!***********************************************************************
      FUNCTION freezing_point (S, Z) RESULT (FP)
!***********************************************************************
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: S, Z
!
!  Local variable declarations.
!
      real(r8) :: FP
!
!-----------------------------------------------------------------------
!  Freezing point temperature of sea water.
!-----------------------------------------------------------------------
!
!! Gill (1982).
!!
!!    FP = S*(-0.0575_r8 + 1.710523E-3_r8*SQRT(S) - 2.154996E-4_r8*S) + &
!!   &     0.000753_r8*Z
!!
!! Steele et al. (1989).
!!
!!    FP = -0.0543_r8*S + 0.000759_r8*Z
!!
      FP = -0.0543_r8*S
!
      RETURN
      END FUNCTION freezing_point
#endif
      END MODULE ice_frazil_mod

