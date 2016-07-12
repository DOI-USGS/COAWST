#include "cppdefs.h"
#if defined FOUR_DVAR && defined OBSERVATIONS
      SUBROUTINE obs_cost (ng, model)
!
!svn $Id: obs_cost.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================

# ifdef WEAK_CONSTRAINT
#  if defined TL_W4DVAR || defined W4DVAR || defined W4DVAR_SENSITIVITY
!                                                                      !
!  This routine computes the data penalty function directly in during  !
!  runs of the representer model:                                      !
!                                                                      !
#  else
!                                                                      !
!  This routine computes the data penalty function directly in during  !
!  runs of the nonlinear model:                                        !
!                                                                      !
#  endif
!         Jdata = transpose(H X - Xo) * O^(-1) * (H X - Xo)            !
!                                                                      !
!         H  : observation operator (linearized if incremental)        !
!         Xo : observations vector                                     !
!       H X  : representer model at observation points                 !
!         O  : observations error covariance                           !
# else
!                                                                      !
!  This routine computes the observation cost function (Jo) as the     !
!  misfit (squared difference) between the model and observations.     !
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
# endif
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

      real(r8) ::  cff, cff1

      real(r8), dimension(0:NstateVar(ng)) :: my_ObsCost
# if defined DATALESS_LOOPS && \
    (defined TL_W4DVAR || defined W4DVAR || defined W4DVAR_SENSITIVITY)
      real(r8), dimension(0:NstateVar(ng)) :: my_ObsCost1
# endif
!
!-----------------------------------------------------------------------
!  Compute observation misfit cost function (ObsCost).
!-----------------------------------------------------------------------

# if defined TL_W4DVAR || defined W4DVAR || defined W4DVAR_SENSITIVITY
!
!  Compute data penalty function.
!
      IF (model.eq.iRPM) THEN
        DO ivar=0,NstateVar(ng)
          my_ObsCost(ivar)=0.0_r8
#  ifdef DATALESS_LOOPS
          my_ObsCost1(ivar)=0.0_r8
#  endif
        END DO
        DO iobs=NstrObs(ng),NendObs(ng)
          ivar=ObsType(iobs)
          IF (ObsErr(iobs).ne.0.0_r8) THEN
            cff=ObsScale(iobs)*(TLmodVal(iobs)-ObsVal(iobs))**2/        &
     &          ObsErr(iobs)
#  ifdef DATALESS_LOOPS
            cff1=ObsScale(iobs)*(NLmodVal(iobs)-ObsVal(iobs))**2/       &
     &           ObsErr(iobs)
#  endif
          END IF
          my_ObsCost(0)=my_ObsCost(0)+cff
          my_ObsCost(ivar)=my_ObsCost(ivar)+cff
#  ifdef DATALESS_LOOPS
          my_ObsCost1(0)=my_ObsCost1(0)+cff1
          my_ObsCost1(ivar)=my_ObsCost1(ivar)+cff1
#  endif
        END DO
      END IF

# elif defined TL_W4DPSAS          || defined W4DPSAS || \
       defined W4DPSAS_SENSITIVITY
!
!  Compute nonlinear model data penalty function.
!
      IF (model.eq.iNLM) THEN
        DO ivar=0,NstateVar(ng)
          my_ObsCost(ivar)=0.0_r8
        END DO
        DO iobs=NstrObs(ng),NendObs(ng)
          ivar=ObsType(iobs)
          IF (ObsErr(iobs).NE.0.0_r8) THEN
            cff=ObsScale(iobs)*(NLmodVal(iobs)-ObsVal(iobs))**2/        &
     &          ObsErr(iobs)
          END IF
          my_ObsCost(0)=my_ObsCost(0)+cff
          my_ObsCost(ivar)=my_ObsCost(ivar)+cff
        END DO
      END IF
# else
!
!  Compute tangent linear model cost function.
!
      IF (model.eq.iTLM) THEN
        DO ivar=0,NstateVar(ng)
          my_ObsCost(ivar)=0.0_r8
        END DO
        DO iobs=1,Nobs(ng)
          ivar=ObsType(iobs)
          cff=0.5_r8*ObsScale(iobs)*ObsErr(iobs)*                       &
     &        (NLmodVal(iobs)+TLmodVal(iobs)-ObsVal(iobs))**2
          my_ObsCost(0)=my_ObsCost(0)+cff
          my_ObsCost(ivar)=my_ObsCost(ivar)+cff
        END DO
!
!  Compute nonlinear model cost function.
!
      ELSE IF (model.eq.iNLM) THEN
        DO ivar=0,NstateVar(ng)
          my_ObsCost(ivar)=0.0_r8
        END DO
        DO iobs=1,Nobs(ng)
          ivar=ObsType(iobs)
          cff=0.5_r8*ObsScale(iobs)*ObsErr(iobs)*                       &
     &        (NLmodVal(iobs)-ObsVal(iobs))**2
          my_ObsCost(0)=my_ObsCost(0)+cff
          my_ObsCost(ivar)=my_ObsCost(ivar)+cff
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
# ifdef DISTRIBUTE
      NSUB=1                             ! distributed-memory
# else
      NSUB=NtileX(ng)*NtileE(ng)
# endif
!$OMP CRITICAL (COST_FUN)
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
# if defined TL_W4DVAR || defined W4DVAR || defined W4DVAR_SENSITIVITY
        IF (model.eq.iRPM) THEN
          DO ivar=0,NstateVar(ng)
            FOURDVAR(ng)%DataPenalty(ivar)=my_ObsCost(ivar)+            &
     &                                    FOURDVAR(ng)%DataPenalty(ivar)
#  ifdef DATALESS_LOOPS
            FOURDVAR(ng)%NLPenalty(ivar)=my_ObsCost1(ivar)+             &
     &                                   FOURDVAR(ng)%NLPenalty(ivar)
#  endif
          END DO
        END IF
# elif defined TL_W4DPSAS          || defined W4DPSAS || \
       defined W4DPSAS_SENSITIVITY
        IF (model.eq.iNLM) THEN
          DO ivar=0,NstateVar(ng)
            FOURDVAR(ng)%NLPenalty(ivar)=my_ObsCost(ivar)+              &
     &                                   FOURDVAR(ng)%NLPenalty(ivar)
          END DO
        END IF
# else
        IF (model.eq.iTLM) THEN
          DO ivar=0,NstateVar(ng)
            FOURDVAR(ng)%ObsCost(ivar)=FOURDVAR(ng)%ObsCost(ivar)+      &
     &                                 my_ObsCost(ivar)
          END DO
        ELSE IF (model.eq.iNLM) THEN
          DO ivar=0,NstateVar(ng)
            FOURDVAR(ng)%NLobsCost(ivar)=FOURDVAR(ng)%NLobsCost(ivar)+  &
     &                                   my_ObsCost(ivar)
          END DO
        END IF
# endif
      END IF
!$OMP END CRITICAL (COST_FUN)

# ifndef WEAK_CONSTRAINT
!
!  If start of minimization, set cost function scales used to report
!  normalized values.
!
      IF ((Nrun.eq.1).and.(model.eq.iTLM)) THEN
        DO ivar=0,NstateVar(ng)
          FOURDVAR(ng)%CostNorm(ivar)=FOURDVAR(ng)%ObsCost(ivar)
        END DO
      END IF
!
!  Save initial inner loop cost function.
!
      IF ((inner.eq.0).and.(model.eq.iTLM)) THEN
        FOURDVAR(ng)%Cost0(outer)=FOURDVAR(ng)%ObsCost(0)
      END IF
# endif

      RETURN
      END SUBROUTINE obs_cost
#else
      SUBROUTINE obs_cost
      RETURN
      END SUBROUTINE obs_cost
#endif
