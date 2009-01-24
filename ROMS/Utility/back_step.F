#include "cppdefs.h"
      MODULE back_step_mod
#ifdef BACKGROUND
!
!svn $Id: back_step.F 733 2008-09-07 01:56:45Z jcwarner $
!=================================================== Andrew M. Moore ===
!  Copyright (c) ROMS/TOMS Adjoint Group                               !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine computes descent algorithm step size due to the     !
!  background cost function in v-space.  It computes the numerator     !
!  and denominator norms separately:                                   !
!                                                                      !
!    StepSizeR = - (StepTopBck + StepTopObs) /                         !
!                  (StepBotBck + StepBotObs)                           !
!                                                                      !
!  and                                                                 !
!                                                                      !
!    StepTopBck = deltaV(m-1) * [deltaV(m) - deltaV(m-1)]              !
!                                                                      !
!    StepBotBck = [deltaV(m) - deltaV(m-1)] ** 2                       !
!                                                                      !
!    StepTopObs = [deltaX(m-1) - d] * 1/O * [deltaX(m) - deltaX(m-1)]  !
!                                                                      !
!    StepBotBcK = [deltaX(m) - deltaX(m-1)] ** 2                       !
!                                                                      !
!  where                                                               !
!                                                                      !
!    m         Minimization (descent algorithm) iteration              !
!    deltaV    Increments due to the background (v-space)              !
!    deltaX    Increments due to the observations (x-space)            !
!    O         Observation error covariance                            !
!                                                                      !
!  The StepTopObs and StepBotObs norms are computed in routine         !
!  "obs_cost".  Recall that the  increments are evaluated with         !
!  the tangent linear model.                                           !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE back_step (ng, tile, Lold, Lnew)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Lold, Lnew
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL back_step_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Lnew, Lold,                                  &
# ifdef SOLVE3D
     &                     OCEAN(ng) % tl_t,                            &
     &                     OCEAN(ng) % tl_u,                            &
     &                     OCEAN(ng) % tl_v,                            &
# else
     &                     OCEAN(ng) % tl_ubar,                         &
     &                     OCEAN(ng) % tl_vbar,                         &
# endif
     &                     OCEAN(ng) % tl_zeta)
      RETURN
      END SUBROUTINE back_step
!
!***********************************************************************
      SUBROUTINE back_step_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           Lold, Lnew,                            &
# ifdef SOLVE3D
     &                           tl_t,                                  &
     &                           tl_u,                                  &
     &                           tl_v,                                  &
# else
     &                           tl_ubar,                               &
     &                           tl_vbar,                               &
# endif
     &                           tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lold, Lnew
!
# ifdef ASSUMED_SHAPE
#  ifdef SOLVE3D
      real(r8), intent(in) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(in) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: tl_zeta(LBi:,LBj:,:)
# else
#  ifdef SOLVE3D
      real(r8), intent(in) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(in) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(in) :: tl_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: NSUB, i, j
# ifdef SOLVE3D
      integer :: itrc, k
# endif

      real(r8) ::  cff, my_StepBot, my_StepTop

# ifdef DISTRIBUTE
      real(r8), dimension(2) :: buffer
      character (len=3), dimension(2) :: op_handle
# endif

# include "set_bounds.h"
!
!----------------------------------------------------------------------
!  Compute background contributions to descent alogorithm step size.
!----------------------------------------------------------------------
!
      my_StepTop=0.0_r8
      my_StepBot=0.0_r8
!
!  Free-surface contribution.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          cff=tl_zeta(i,j,Lnew)-tl_zeta(i,j,Lold)
          my_StepTop=my_StepTop+0.5_r8*cff*tl_zeta(i,j,Lold)
          my_StepBot=my_StepBot+0.5_r8*cff*cff
        END DO
      END DO

# if !defined SOLVE3D
!
!  2D momentum contribution.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          cff=tl_ubar(i,j,Lnew)-tl_ubar(i,j,Lold)
          my_StepTop=my_StepTop+0.5_r8*cff*tl_ubar(i,j,Lold)
          my_StepBot=my_StepBot+0.5_r8*cff*cff
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          cff=tl_vbar(i,j,Lnew)-tl_vbar(i,j,Lold)
          my_StepTop=my_StepTop+0.5_r8*cff*tl_vbar(i,j,Lold)
          my_StepBot=my_StepBot+0.5_r8*cff*cff
        END DO
      END DO
# else
!
!  3D momentum contribution.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            cff=tl_u(i,j,k,Lnew)-tl_u(i,j,k,Lold)
            my_StepTop=my_StepTop+0.5_r8*cff*tl_u(i,j,k,Lold)
            my_StepBot=my_StepBot+0.5_r8*cff*cff
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            cff=tl_v(i,j,k,Lnew)-tl_v(i,j,k,Lold)
            my_StepTop=my_StepTop+0.5_r8*cff*tl_v(i,j,k,Lold)
            my_StepBot=my_StepBot+0.5_r8*cff*cff
          END DO
        END DO
      END DO
!
!  Tracers contribution.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              cff=tl_t(i,j,k,Lnew,itrc)-tl_t(i,j,k,Lold,itrc)
              my_StepTop=my_StepTop+0.5_r8*cff*tl_t(i,j,k,Lold,itrc)
              my_StepBot=my_StepBot+0.5_r8*cff*cff
            END DO
          END DO
        END DO
      END DO
# endif
!
!-----------------------------------------------------------------------
! Compute global background cost function.
!-----------------------------------------------------------------------
!
      IF (SOUTH_WEST_CORNER.and.                                        &
     &    NORTH_EAST_CORNER) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
!$OMP CRITICAL (BACKCOST)
      IF (tile_count.eq.0) THEN
        StepBotBck=my_StepBot
        StepTopBck=my_StepTop
      ELSE
        StepBotBck=StepBotBck+my_StepBot
        StepTopBck=StepTopBck+my_StepTop
      END IF
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
# ifdef DISTRIBUTE
        buffer(1)=StepBotBck
        buffer(2)=StepTopBck
        op_handle(1)='SUM'
        op_handle(2)='SUM'
        CALL mp_reduce (ng, iTLM, 2, buffer, op_handle)
        StepBotBck=buffer(1)
        StepTopBck=buffer(2)
# endif
      END IF
!$OMP END CRITICAL (BACKCOST)

      RETURN
      END SUBROUTINE back_step_tile
#endif
      END MODULE back_step_mod
