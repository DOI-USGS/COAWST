#include "cppdefs.h"
#if defined FOUR_DVAR && defined OBSERVATIONS && \
    !(defined WEAK_CONSTRAINT || defined IOM)
      SUBROUTINE obs_cost (ng, model)
!
!svn $Id: obs_cost.F 503 2008-01-10 00:11:51Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  If minimization first pass, this routine computes the observation   !
!  cost function (Jo) as the misfit (squared difference) between the   !
!  model and observations.                                             !
!                                                                      !
!  If conventional strong contraint 4DVAR:                             !
!                                                                      !
!         Jo = 1/2 transpose(H X - Xo) * O^(-1) * (H X - Xo)           !
!                                                                      !
!  or if incremental strong contraint 4DVAR:                           !
!                                                                      !
!         Jo = 1/2 transpose(H deltaX - d) * O^(-1) * (H deltaX - d)   !
!                                                                      !
!  where                                                               !
!                                                                      !
!          d = Xo - H Xb                                               !
!                                                                      !
!         d  : innovation vector                                       !
!         H  : observation operator (linearized if incremental)        !
!       H Xb : background at observation points previous forecast)     !
!         Xo : observations vector                                     !
!       H X  : nonlinear model at observation points                   !
!  H deltaX  : increment at observation point                          !
!         O  : observations error covariance                           !
!                                                                      !
!  If minimization second pass (m=2), it computes conjugate gradient   !
!  step size numerator and denominator factors due to observations:    !
!                                                                      !
!  If conventional strong contraint 4DVAR:                             !
!                                                                      !
!    StepTopObs = 1/2 * (H X(m-1) - Xo) * O^(-1) * (H X(m) - H X(m-1)) !
!                                                                      !
!    StepBotObs = 1/2 * O^(-1) * (H X(m) - H X(m-1)) ** 2              ! 
!                                                                      !
!  or if incremental strong contraint 4DVAR:                           !
!                                                                      !
!    StepTopObs = 1/2 * (H deltaX(m-1) - d) * O^(-1) *                 !
!                       (H deltaX(m) - H deltaX(m-1))                  !
!                                                                      !
!    StepBotObs = 1/2 * O^(-1) * (H deltaX(m) - H deltaX(m-1)) ** 2    ! 
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer :: NSUB, iobs, ivar

      real(r8) ::  cff
# if defined IS4DVAR_OLD || defined S4DVAR
      real(r8) ::  my_StepBot, my_StepTop
# endif
      real(r8), dimension(0:NstateVar(ng)) :: my_ObsCost
!
!-----------------------------------------------------------------------
!  Compute observation misfit cost function (ObsCost).
!-----------------------------------------------------------------------
!
      DO ivar=0,NstateVar(ng)
        my_ObsCost(ivar)=0.0_r8
      END DO
# if defined IS4DVAR
      DO iobs=1,Nobs(ng)
        ivar=ObsType(iobs)
        cff=0.5_r8*ObsScale(iobs)*ObsErr(iobs)*                         &
     &      (NLmodVal(iobs)+TLmodVal(iobs)-ObsVal(iobs))**2
        my_ObsCost(0)=my_ObsCost(0)+cff
        my_ObsCost(ivar)=my_ObsCost(ivar)+cff
      END DO
# else
      IF (Ipass.eq.1) THEN
        DO iobs=1,Nobs(ng)
          ivar=ObsType(iobs)
#  if defined S4DVAR
          cff=0.5_r8*ObsScale(iobs)*ObsErr(iobs)*                       &
     &        (NLmodVal(iobs,Ipass)-ObsVal(iobs))**2
#  elif defined IS4DVAR_OLD
          cff=0.5_r8*ObsScale(iobs)*ObsErr(iobs)*                       &
     &        (NLmodVal(iobs)+TLmodVal(iobs,Ipass)-ObsVal(iobs))**2
#  endif
          my_ObsCost(0)=my_ObsCost(0)+cff
          my_ObsCost(ivar)=my_ObsCost(ivar)+cff
        END DO
      END IF
# endif
# if defined S4DVAR || defined IS4DVAR_OLD
!
!-----------------------------------------------------------------------
!  Compute Conjugate gradient algorithm step size dot products,
!  StepTopObs, StepBotObs.
!-----------------------------------------------------------------------
!
      IF (Ipass.eq.2) THEN
        my_StepTop=0.0_r8
        my_StepBot=0.0_r8
        DO iobs=1,Nobs(ng)
#  if defined S4DVAR
          cff=0.5_r8*ObsScale(iobs)*ObsErr(iobs)
          my_StepTop=my_StepTop+                                        &
     &               cff*                                               &
     &               (NLmodVal(iobs,Ipass-1)-ObsVal(iobs))*             &
     &               (NLmodVal(iobs,Ipass)-NLmodVal(iobs,Ipass-1))
          my_StepBot=my_StepBot+                                        &
     &               cff*                                               &
     &               (NLmodVal(iobs,Ipass)-NLmodVal(iobs,Ipass-1))**2
#  elif defined IS4DVAR_OLD
          cff=0.5_r8*ObsScale(iobs)*ObsErr(iobs)
          my_StepTop=my_StepTop+                                        &
     &               cff*                                               &
     &               (NLmodVal(iobs)+TLmodVal(iobs,Ipass-1)-            &
     &                ObsVal(iobs))*                                    &
     &               (TLmodVal(iobs,Ipass)-TLmodVal(iobs,Ipass-1))
          my_StepBot=my_StepBot+                                        &
     &               cff*                                               &
     &               (TLmodVal(iobs,Ipass)-TLmodVal(iobs,Ipass-1))**2
#  endif 
        END DO
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Load global values.  Notice that there is not need for a global
!  reduction here since all the threads have the same copy of all
!  the vectors used.
!-----------------------------------------------------------------------
!
      NSUB=NtileX(ng)*NtileE(ng)
!$OMP CRITICAL (COST_FUN)
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
        DO ivar=0,NstateVar(ng)
          FOURDVAR(ng)%ObsCost(ivar)=FOURDVAR(ng)%ObsCost(ivar)+        &
     &                               my_ObsCost(ivar)
        END DO
# ifdef IS4DVAR_OLD
        StepBotObs=StepBotObs+my_StepBot
        StepTopObs=StepTopObs+my_StepTop
# endif
      END IF
!$OMP END CRITICAL (COST_FUN)
!
!  If start of minimization, set cost function scales used to report
!  normalized values.
!
# ifdef IS4DVAR_OLD
      IF ((Nrun.eq.1).and.(Ipass.eq.1)) THEN
# else
      IF (Nrun.eq.1) THEN
# endif
        DO ivar=0,NstateVar(ng)
          FOURDVAR(ng)%CostNorm(ivar)=FOURDVAR(ng)%ObsCost(ivar)
        END DO
      END IF

      RETURN
      END SUBROUTINE obs_cost
#else
      SUBROUTINE obs_cost
      RETURN
      END SUBROUTINE obs_cost
#endif
