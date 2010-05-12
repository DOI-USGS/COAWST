#include "cppdefs.h"

#ifdef SENSITIVITY_4DVAR
      SUBROUTINE ad_congrad (ng, model, outLoop, innLoop, NinnLoop,     &
     &                       Lcgini)
!
!svn $Id: ad_congrad.F 429 2009-12-20 17:30:26Z arango $
!=================================================== Andrew M. Moore ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group      Hernan G. Arango   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Weak Constraint 4-Dimensional Variational (4DVar) Pre-conditioned   !
!                     Conjugate Gradient Algorithm                     !

# if defined TL_W4DVAR || defined W4DVAR_SENSITIVITY
!                                                                      !
!  The indirect representer method solves the system:                  !
!                                                                      !
!              (R_n + Cobs) * Beta_n = h_n                             !
!                                                                      !
!              h_n = Xo - H * X_n                                      !
!                                                                      !
!  where R_n is the representer matrix, Cobs is the observation-error  !
!  covariance,  Beta_n  are the representer coefficients,  h_n is the  !
!  misfit between observations (Xo) and model (H*X_n),  and  H is the  !
!  linearized observation operator. Here, _n denotes iteration.        !
!                                                                      !
!  This system does not need to be solved explicitly by inverting the  !
!  symmetric stabilized representer matrix, P_n:                       !
!                                                                      !
!              P_n = R_n + Cobs                                        !
!                                                                      !
!  but by computing the action of P_n on any vector PSI, such that     !
!                                                                      !
!              P_n * PSI = R_n * PSI + Cobs * PSI                      !
!                                                                      !
!  The representer matrix is not explicitly computed but evaluated by  !
!  one integration backward of the adjoint model  and one integration  !
!  forward of the tangent linear model for any forcing vector PSI.     !
!                                                                      !
!  Notice that "ObsScale" vector is used for screenning observations.  !
!  This scale is one (zero) for good (bad) observations.               !
!                                                                      !
!  Currently, parallelization of this algorithm is not needed because  !
!  each parallel node has a full copy of the assimilation vectors.     !
!                                                                      !
!  This code solves Ax=b by minimizing the cost function               !
!  0.5*x*A*x-x*b, assuming an initial guess of x=x0. In this case the  !
!  gradient is Ax-b and the Hessian is A.                              !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Chua, B. S. and A. F. Bennett,  2001:  An inverse ocean modeling  !
!      sytem, Ocean Modelling, 3, 137-165.                             !

# elif defined TL_W4DPSAS || defined W4DPSAS_SENSITIVITY
!                                                                      !
!  Solve the system (following Courtier, 1997):                        !
!                                                                      !
!              (H M_n B (M_n)' H' + Cobs) * w_n = d_n                  !
!                                                                      !
!              d_n = yo - H * Xb_n                                     !
!                                                                      !
!  where  M_n is the tangent linear model matrix and  _n denotes a     !
!  sequence of outer-loop estimates, Cobs is the observation-error     !
!  covariance,  B is the background error  covariance,  d_n is the     !
!  misfit between observations (yo) and model (H * Xb_n), and H is     !
!  the linearized observation operator.                                !
!                                                                      !
!  The analysis increment is:                                          !
!                                                                      !
!             dx_n = B M' H' w_n                                       !
!                                                                      !
!  so that Xa = Xb + dx_n.                                             !
!                                                                      !
!  The system does not need to be  solved explicitly  by inverting     !
!  the symmetric matrix, P_n:                                          !
!                                                                      !
!              P_n = H M_n B (M_n)' H' + Cobs                          !
!                                                                      !
!  but by computing the action of P_n on any vector PSI, such that     !
!                                                                      !
!              P_n * PSI =  H M_n B (M_n)' H' * PSI + Cobs * PSI       !
!                                                                      !
!  The (H M_n B (M_n)' H') matrix is not  explicitly computed  but     !
!  evaluated by  one integration backward of the adjoint model and     !
!  one  integration  forward of the  tangent linear model  for any     !
!  forcing vector PSI.                                                 !
!                                                                      !
!  A preconditioned conjugate gradient algorithm is used to compute    !
!  an approximation PSI for w_n.                                       !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Courtier, P., 1997: Dual formulation of four-dimensional          !
!      variational assimilation, Quart. J. Roy. Meteor. Soc.,          !
!      123, 2449-2461.                                                 !
# endif
!                                                                      !
!  Lanczos Algorithm Reference:                                        !
!                                                                      !
!    Fisher, M., 1998: Minimization Algorithms for Variational Data    !
!      Assimilation. In Seminar on Recent Developments in Numerical    !
!      Methods for Atmospheric Modelling, 1998.                        !
!                                                                      !
!    Tchimanga, J., S. Gratton, A.T. Weaver, and A. Sartenaer, 2008:   !
!      Limited-memory preconditioners, with application to incremental !
!      four-dimensional variational ocean data assimilation, Q.J.R.    !
!      Meteorol. Soc., 134, 753-771.                                   !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti, mp_bcastl
# endif
      implicit none
!
!  Imported variable declarations
!
      integer, intent(in) :: ng, model, outLoop, innLoop, NinnLoop
      logical, intent(in) :: Lcgini
!
!  Local variable declarations.
!
      logical :: Ltrans

      integer :: i, j, iobs, ivec, Lscale, info

      real(r8) :: dla, zbet
      real(r8) :: adfac, ad_dla

      real(r8), dimension(NinnLoop) :: zu, zgam
      real(r8), dimension(NinnLoop) :: ad_zu, ad_zrhs
      real(r8), dimension(Ndatum(ng)) :: pgrad, zt, pgrad_S
      real(r8), dimension(Ndatum(ng)) :: ad_px, ad_pgrad, ad_zt
!
!-----------------------------------------------------------------------
!  Weak constraint 4DVar conjugate gradient, Lanczos-based algorithm.
!-----------------------------------------------------------------------
!
!  This conjugate gradient algorithm is not run in parallel since the
!  weak constraint is done in observation space. Mostly all the
!  variables are 1D arrays. Therefore, in parallel applications (only
!  distributed-memory is possible) the master node does the computations
!  and then broadcasts results to remaining nodes in the communicator.
!
      Ltrans=.FALSE.

      MASTER_THREAD : IF (Master) THEN
!
!  Initialize internal parameters.  This is needed here for output
!  reasons.
!
!
! Clear arrays.
!
        ADmodVal=0.0_r8
        ad_pgrad=0.0_r8
        ad_px=0.0_r8
        ad_zu=0.0_r8
        ad_zrhs=0.0_r8
        ad_zt=0.0_r8
!
!  Initialize cg_Gnorm. The adjoint of reprecond is not available.
!
        DO i=1,outLoop
          cg_Gnorm(i)=cg_Gnorm_v(i)
        END DO
!
        IF (innLoop.eq.0) THEN

# ifdef W4DPSAS
!
!  If this is the first inner-loop, save NLmodVal in BCKmodVal.
!
          DO iobs=1,Ndatum(ng)
            BCKmodVal(iobs)=NLmodVal(iobs)
          END DO
# endif
!
          IF ((outLoop.eq.1).or.(.not.LhotStart)) THEN
!
            DO iobs=1,Ndatum(ng)
!>            tl_cg_QG(1)=tl_cg_QG(1)+                                  &
!>   &                    tl_zcglwk(iobs,1)*zgrad0(iobs,outLoop)+       &
!>   &                    zcglwk(iobs,1,outLoop)*tl_zgrad0(iobs)
!>
              ad_zcglwk(iobs,1)=ad_zcglwk(iobs,1)+                      &
     &                          zgrad0(iobs,outLoop)*ad_cg_QG(1)
              ad_zgrad0(iobs)=ad_zgrad0(iobs)+                          &
     &                        zcglwk(iobs,1,outLoop)*ad_cg_QG(1)
            END DO
!>          tl_cg_QG(1)=0.0_r8
!>
            ad_cg_QG(1)=0.0_r8
!
! Convert ADmodVal from v-space to x-space.
!
            DO iobs=1,Ndatum(ng)
              IF (ObsErr(iobs).NE.0.0_r8) THEN
!<>             tl_ADmodVal(iobs)=tl_ADmodVal(iobs)/SQRT(ObsErr(iobs))
!<>             ad_ADmodVal(iobs)=ad_ADmodVal(iobs)/SQRT(ObsErr(iobs))
!>
                TLmodVal(iobs)=TLmodVal(iobs)/SQRT(ObsErr(iobs))
              END IF
            END DO
!
!  If preconditioning, convert ADmodVal from y-space to v-space.
!
!>          IF (Lprecond.and.(outLoop.gt.1)) THEN
!>            Lscale=2                 ! SQRT spectral LMP
!>            Ltrans=.FALSE.
!>            CALL rprecond (ng, Lscale, Ltrans, outLoop, NinnLoop,     &
!>   &                       ADmodVal)
!>          END IF
!>
            DO iobs=1,Ndatum(ng)
!<>           tl_ADmodVal(iobs)=tl_zcglwk(iobs,1)
!<>           ad_zcglwk(iobs,1)=ad_zcglwk(iobs,1)+ad_ADmodVal(iobs)
!>
              ad_zcglwk(iobs,1)=ad_zcglwk(iobs,1)+TLmodVal(iobs)
!<>           ad_ADmodVal(iobs)=0.0_r8
              TLmodVal(iobs)=0.0_r8

!>            tl_zcglwk(iobs,1)=(tl_pgrad(iobs)-                        &
!>   &                           tl_Gnorm(outLoop)*                     &
!>   &                           zcglwk(iobs,1,outLoop))/               &
!>   &                          cg_Gnorm(outLoop)
!>
              adfac=ad_zcglwk(iobs,1)/cg_Gnorm(outLoop)
              ad_cg_Gnorm(outLoop)=ad_cg_Gnorm(outLoop)-                &
     &                             zcglwk(iobs,1,outLoop)*adfac
              ad_pgrad(iobs)=adfac
              ad_zcglwk(iobs,1)=0.0_r8
            END DO
!>          tl_cg_Gnorm(outLoop)=0.5_r8*tl_cg_Gnorm(outLoop)/           &
!>   &                           cg_Gnorm(outLoop)
!>
            ad_cg_Gnorm(outLoop)=0.5_r8*ad_cg_Gnorm(outLoop)/           &
     &                           cg_Gnorm(outLoop)

            DO iobs=1,Ndatum(ng)
!>            tl_cg_Gnorm(outLoop)=tl_cg_Gnorm(outLoop)+                &
!>   &                             2.0_r8*tl_zgrad0(iobs)*              &
!>   &                             zgrad0(iobs,outLoop)
!>
              ad_zgrad0(iobs)=ad_zgrad0(iobs)+                          &
     &                        2.0_r8*ad_cg_Gnorm(outLoop)*              &
     &                        zgrad0(iobs,outLoop)
!>            tl_zgrad0(iobs)=tl_pgrad(iobs)
!>
              ad_pgrad(iobs)=ad_pgrad(iobs)+ad_zgrad0(iobs)
              ad_zgrad0(iobs)=0.0_r8
            END DO
!>          tl_cg_Gnorm(outLoop)=0.0_r8
!>
            ad_cg_Gnorm(outLoop)=0.0_r8
!
!  If preconditioning, convert pgrad from v-space to y-space.
!
!>          IF (Lprecond.and.(outLoop.gt.1)) THEN
!>            Lscale=2                 ! SQRT spectral LMP
!>            Ltrans=.TRUE.
!>            CALL rprecond (ng, Lscale, Ltrans, outLoop, NinnLoop,     &
!>   &                       pgrad)
!>          END IF

            DO iobs=1,Ndatum(ng)
!
! Convert pgrad from x-space to v-space.
!
              IF (ObsErr(iobs).NE.0.0_r8) THEN
!>              tl_pgrad(iobs)=tl_pgrad(iobs)/SQRT(ObsErr(iobs))
!>
                ad_pgrad(iobs)=ad_pgrad(iobs)/SQRT(ObsErr(iobs))
              END IF
# ifdef W4DPSAS
!>            tl_pgrad(iobs)=-ObsScale(iobs)*tl_ObsVal(iobs)
!>
              ad_ObsVal(iobs)=ad_ObsVal(iobs)-ObsScale(iobs)*           &
     &                        ad_pgrad(iobs)
              ad_pgrad(iobs)=0.0_r8
# else
!<>           tl_pgrad(iobs)=-ObsScale(iobs)*                           &
!<>  &                       (tl_ObsVal(iobs)-tl_TLmodVal(iobs))
!<>           ad_TLmodVal(iobs)=ad_TLmodVal(iobs)+                      &
!<>  &                          ObsScale(iobs)*ad_pgrad(iobs)
!>
              ADmodVal(iobs)=ADmodVal(iobs)+                            &
     &                       ObsScale(iobs)*ad_pgrad(iobs)
              ad_ObsVal(iobs)=ad_ObsVal(iobs)-ObsScale(iobs)*           &
     &                        ad_pgrad(iobs)
              ad_pgrad(iobs)=0.0_r8
# endif
            END DO

# ifdef W4DPSAS_SENSITIVITY
!
!  For observation sensitivity, copy ADmodVal into TLmodVal.
!
!           DO iobs=1,Ndatum(ng)
!             TLmodVal(iobs)=ADmodVal(iobs)
!           END DO
# endif
          ELSE

            IF (Lcgini) THEN
!
!  Convert ADmodVal from v-space to x-space and cg_innov (the
!  contribution to the starting gradient) from x-space to v-space.
!
              DO iobs=1,Ndatum(ng)
                IF (ObsErr(iobs).NE.0.0_r8) THEN
!>                tl_cg_innov(iobs)=tl_cg_innov(iobs)/SQRT(ObsErr(iobs))
!>
                  ad_cg_innov(iobs)=ad_cg_innov(iobs)/SQRT(ObsErr(iobs))
!<>               tl_ADmodVal(iobs)=tl_ADmodVal(iobs)/SQRT(ObsErr(iobs))
!<>               ad_ADmodVal(iobs)=ad_ADmodVal(iobs)/SQRT(ObsErr(iobs))
!>
                  TLmodVal(iobs)=TLmodVal(iobs)/SQRT(ObsErr(iobs))
                END IF
              END DO
!
!  When outer>1 we need to evaluate Ax0 so for inner=0 we use
!  cg_pxsave in v-space as the starting vector.
!
              DO iobs=1,Ndatum(ng)
# ifdef W4DPSAS
!>              tl_cg_innov(iobs)=-ObsScale(iobs)*tl_ObsVal(iobs)
!>
                ad_ObsVal(iobs)=ad_ObsVal(iobs)-ObsScale(iobs)*         &
     &                          ad_cg_innov(iobs)
                ad_cg_innov(iobs)=0.0_r8
# else
!<>             tl_cg_innov(iobs)=-ObsScale(iobs)*                      &
!<>  &                            (tl_ObsVal(iobs)-tl_TLmodVal(iobs))
!<>             ad_TLmodVal(iobs)=ad_TLmodVal(iobs)+                    &
!<>  &                            ObsScale(iobs)*ad_cg_innov(iobs)
!>
                ADmodVal(iobs)=ADmodVal(iobs)+                          &
     &                         ObsScale(iobs)*ad_cg_innov(iobs)
                ad_ObsVal(iobs)=ad_ObsVal(iobs)-ObsScale(iobs)*         &
     &                          ad_cg_innov(iobs)
                ad_cg_innov(iobs)=0.0_r8
# endif
!<>             tl_ADmodVal(iobs)=tl_cg_pxsave(iobs)
!<>             ad_cg_pxsave(iobs)=ad_cg_pxsave(iobs)+ad_ADmodVal(iobs)
!>
                ad_cg_pxsave(iobs)=ad_cg_pxsave(iobs)+TLmodVal(iobs)
!<>             ad_ADmodVal(iobs)=0.0_r8
                TLmodVal(iobs)=0.0_r8
              END DO

            ELSE

              DO iobs=1,Ndatum(ng)
!>              tl_cg_QG(1)=tl_cg_QG(1)+                                &
!>   &                      tl_zcglwk(iobs,1)*zgrad0(iobs,outLoop)+     &
!>   &                      zcglwk(iobs,1,outLoop)*tl_zgrad0(iobs)
!>
                ad_zgrad0(iobs)=ad_zgrad0(iobs)+                        &
     &                          zcglwk(iobs,1,outLoop)*ad_cg_QG(1)
                ad_zcglwk(iobs,1)=ad_zcglwk(iobs,1)+                    &
     &                            zgrad0(iobs,outLoop)*ad_cg_QG(1)
              END DO
!>            tl_cg_QG(1)=0.0_r8
!>
              ad_cg_QG(1)=0.0_r8

              DO iobs=1,Ndatum(ng)
                IF (ObsErr(iobs).NE.0.0_r8) THEN
!<>               tl_ADmodVal(iobs)=tl_ADmodVal(iobs)/SQRT(ObsErr(iobs))
!<>               ad_ADmodVal(iobs)=ad_ADmodVal(iobs)/SQRT(ObsErr(iobs))
!>
                  TLmodVal(iobs)=TLmodVal(iobs)/SQRT(ObsErr(iobs))
                END IF
              END DO
!
!  If preconditioning, convert ADmodVal from y-space to v-space.
!
!>            IF (Lprecond.and.(outLoop.gt.1)) THEN
!>              Lscale=2                 ! SQRT spectral LMP
!>              Ltrans=.FALSE.
!>              CALL rprecond (ng, Lscale, Ltrans, outLoop, NinnLoop,   &
!>   &                         ADmodVal)
!>            END IF
!>
              DO iobs=1,Ndatum(ng)
!<>             tl_ADmodVal(iobs)=tl_zcglwk(iobs,1)
!<>             ad_zcglwk(iobs,1)=ad_zcglwk(iobs,1)+ad_ADmodVal(iobs)
!>
                ad_zcglwk(iobs,1)=ad_zcglwk(iobs,1)+TLmodVal(iobs)
!<>             ad_ADmodVal(iobs)=0.0_r8
                TLmodVal(iobs)=0.0_r8
!>              tl_zcglwk(iobs,1)=(tl_pgrad(iobs)-                      &
!>   &                             tl_cg_Gnorm(outLoop)*                &
!>   &                             zcglwk(iobs,1,outLoop))/             &
!>   &                            cg_Gnorm(outLoop)
!>
                adfac=ad_zcglwk(iobs,1)/cg_Gnorm(outLoop)
                ad_cg_Gnorm(outLoop)=ad_cg_Gnorm(outLoop)-              &
     &                               zcglwk(iobs,1,outLoop)*adfac
                ad_pgrad(iobs)=ad_pgrad(iobs)+adfac
                ad_zcglwk(iobs,1)=0.0_r8
              END DO
!>            tl_cg_Gnorm(outLoop)=0.5_r8*tl_cg_Gnorm(outLoop)/         &
!>   &                             cg_Gnorm(outLoop)
!>
              ad_cg_Gnorm(outLoop)=0.5_r8*ad_cg_Gnorm(outLoop)/         &
     &                             cg_Gnorm(outLoop)

              DO iobs=1,Ndatum(ng)
!>              tl_cg_Gnorm(outLoop)=tl_cg_Gnorm(outLoop)+              &
!>   &                               2.0_r8*tl_zgrad0(iobs)*            &
!>   &                               zgrad0(iobs,outLoop)
!>
                ad_zgrad0(iobs)=ad_zgrad0(iobs)+                        &
     &                          2.0_r8*zgrad0(iobs,outLoop)*            &
     &                          ad_cg_Gnorm(outLoop)
!>              tl_zgrad0(iobs)=tl_pgrad(iobs)
!>
                ad_pgrad(iobs)=ad_pgrad(iobs)+ad_zgrad0(iobs)
                ad_zgrad0(iobs)=0.0_r8
              END DO
!>            tl_cg_Gnorm(outLoop)=0.0_r8
!>
              ad_cg_Gnorm(outLoop)=0.0_r8
!
!  If preconditioning, convert pgrad from v-space to y-space.
!
!>            IF (Lprecond.and.(outLoop.gt.1)) THEN
!>              Lscale=2                 ! SQRT spectral LMP
!>              Ltrans=.TRUE.
!>              CALL rprecond (ng, Lscale, Ltrans, outLoop, NinnLoop,   &
!>   &                         pgrad)
!>            END IF

              DO iobs=1,Ndatum(ng)
!
!  Add I*x0=cg_pxsave contribution to the gradient and the term
!  -b=cg_innov (everything is in v-space at this point).
!
!>              tl_pgrad(iobs)=tl_pgrad(iobs)+ObsScale(iobs)*           &
!>   &                         (tl_cg_pxsave(iobs)+tl_cg_innov(iobs))
!>
                adfac=ObsScale(iobs)*ad_pgrad(iobs)
                ad_cg_innov(iobs)=ad_cg_innov(iobs)+adfac
                ad_cg_pxsave(iobs)=ad_cg_pxsave(iobs)+adfac
!
!  Convert gradient contribution from x-space to v-space.
!
                IF (ObsErr(iobs).NE.0.0_r8) THEN
!>                tl_pgrad(iobs)=tl_pgrad(iobs)/SQRT(ObsErr(iobs))
!>
                  ad_pgrad(iobs)=ad_pgrad(iobs)/SQRT(ObsErr(iobs))
                END IF
!<>             tl_pgrad(iobs)=ObsScale(iobs)*tl_TLmodVal(iobs)
!<>             ad_TLmodVal(iobs)=ad_TLmodVal(iobs)+                    &
!<>  &                            ObsScale(iobs)*ad_pgrad(iobs)
!>
                ADmodVal(iobs)=ADmodVal(iobs)+                          &
     &                         ObsScale(iobs)*ad_pgrad(iobs)
                ad_pgrad(iobs)=0.0_r8
              END DO
            END IF
          END IF

        ELSE
!
!  Convert ADmodVal from v-space to x-space.
!
          DO iobs=1,Ndatum(ng)
            IF (ObsErr(iobs).NE.0.0_r8) THEN
!<>           tl_ADmodVal(iobs)=tl_ADmodVal(iobs)/SQRT(ObsErr(iobs))
!<>           ad_ADmodVal(iobs)=ad_ADmodVal(iobs)/SQRT(ObsErr(iobs))
!>
              TLmodVal(iobs)=TLmodVal(iobs)/SQRT(ObsErr(iobs))
            END IF
          END DO

          IF (innLoop.eq.NinnLoop) THEN

            IF (LhotStart) THEN
              DO iobs=1,Ndatum(ng)
!<>             tl_cg_pxsave(iobs)=tl_ADmodVal(iobs)
!<>             ad_ADmodVal(iobs)=ad_ADmodVal(iobs)+ad_cg_pxsave(iobs)
!>
                TLmodVal(iobs)=TLmodVal(iobs)+ad_cg_pxsave(iobs)
                ad_cg_pxsave(iobs)=0.0_r8
              END DO
              DO iobs=1,Ndatum(ng)
!<>             tl_ADmodVal(iobs)=tl_ADmodVal(iobs)+tl_cg_pxsave(iobs)
!<>             ad_cg_pxsave(iobs)=ad_ADmodVal(iobs)
!>
                ad_cg_pxsave(iobs)=TLmodVal(iobs)
              END DO
            END IF
!
!  Note: px is already in v-space so there is no need to convert
!  if preconditioning.
!
            DO iobs=1,Ndatum(ng)
!<>           tl_ADmodVal(iobs)=tl_px(iobs)
!<>           ad_px(iobs)=ad_px(iobs)+ad_ADmodVal(iobs)
!>
              ad_px(iobs)=ad_px(iobs)+TLmodVal(iobs)
!<>           ad_ADmodVal(iobs)=0.0_r8
              TLmodVal(iobs)=0.0_r8
            END DO
          ELSE
!
!  If preconditioning, convert ADmodVal from y-space to v-space.
!
!>          IF (Lprecond.and.(outLoop.gt.1)) THEN
!>            Lscale=2                 ! SQRT spectral LMP
!>            Ltrans=.FALSE.
!>            CALL rprecond (ng, Lscale, Ltrans, outLoop, NinnLoop,     &
!>   &                       ADmodVal)
!>          END IF
!>
            DO iobs=1,Ndatum(ng)
!<>           tl_ADmodVal(iobs)=tl_zcglwk(iobs,innLoop+1)
!<>           ad_zcglwk(iobs,innLoop+1)=ad_zcglwk(iobs,innLoop+1)+      &
!<>  &                                  ad_ADmodVal(iobs)
!>
              ad_zcglwk(iobs,innLoop+1)=ad_zcglwk(iobs,innLoop+1)+      &
     &                                  TLmodVal(iobs)
!<>           ad_ADmodVal(iobs)=0.0_r8
              TLmodVal(iobs)=0.0_r8
            END DO
          END IF

!
!  Calculate the new solution based upon the new, orthonormalized
!  Lanczos vector. First, the tridiagonal system is solved by
!  decomposition and forward/back substitution.
!
          IF (innLoop.eq.NinnLoop) THEN

!
!  Compute zu and zgam first.
!
            zbet=cg_delta(1,outLoop)
            zu(1)=-cg_QG(1,outLoop)/zbet
            DO ivec=2,innLoop
              zgam(ivec)=cg_beta(ivec,outLoop)/zbet
              zbet=cg_delta(ivec,outLoop)-                              &
     &             cg_beta(ivec,outLoop)*zgam(ivec)
              zu(ivec)=(-cg_QG(ivec,outLoop)-cg_beta(ivec,outLoop)*     &
     &                  zu(ivec-1))/zbet
            END DO
            DO ivec=innLoop-1,1,-1
              zu(ivec)=zu(ivec)-zgam(ivec+1)*zu(ivec+1)
            END DO

            DO iobs=1,Ndatum(ng)
              DO ivec=1,innLoop
!>              tl_px(iobs)=tl_px(iobs)+                                &
!>   &                      tl_zcglwk(iobs,ivec)*zu(ivec)+              &
!>   &                      zcglwk(iobs,ivec,outLoop)*tl_zu(ivec)
!>
                ad_zu(ivec)=ad_zu(ivec)+                                &
     &                      zcglwk(iobs,ivec,outLoop)*ad_px(iobs)
                ad_zcglwk(iobs,ivec)=ad_zcglwk(iobs,ivec)+              &
     &                               zu(ivec)*ad_px(iobs)
              END DO
!>            tl_px(iobs)=0.0_r8
!>
              ad_px(iobs)=0.0_r8
            END DO

            IF (NinnLoop.eq.1) THEN
              zu(1)=-cg_QG(1,outLoop)/cg_delta(1,outLoop)
!>            tl_zu(1)=tl_zrhs(1)/cg_delta(1,outLoop)
!>
              ad_zrhs(1)=ad_zrhs(1)+ad_zu(1)/cg_delta(1,outLoop)
              ad_zu(1)=0.0_r8
!>            tl_zrhs(1)=-tl_cg_QG(1)-tl_cg_delta(1)*zu(1)
!>
              ad_cg_QG(1)=ad_cg_QG(1)-ad_zrhs(1)
              ad_cg_delta(1)=ad_cg_delta(1)-zu(1)*ad_zrhs(1)
              ad_zrhs(1)=0.0_r8
            ELSE
!>            DO ivec=innLoop-1,1,-1
!>
              DO ivec=1,innLoop-1
!>              tl_zu(ivec)=tl_zu(ivec)-zgam(ivec+1)*tl_zu(ivec+1)
!>
                ad_zu(ivec+1)=ad_zu(ivec+1)-zgam(ivec+1)*ad_zu(ivec)
              END DO
!>            DO ivec=2,innLoop
!>
              DO ivec=innLoop,2,-1
                zbet=cg_delta(ivec,outLoop)-                            &
     &               cg_beta(ivec,outLoop)*zgam(ivec)
!>              tl_zu(ivec)=(tl_zrhs(ivec)-cg_beta(ivec,outLoop)*       &
!>   &                       tl_zu(ivec-1))/zbet
!>
                adfac=ad_zu(ivec)/zbet
                ad_zu(ivec-1)=ad_zu(ivec-1)-                            &
     &                        cg_beta(ivec,outLoop)*adfac
                ad_zrhs(ivec)=ad_zrhs(ivec)+adfac
                ad_zu(ivec)=0.0_r8
              END DO
              zbet=cg_delta(1,outLoop)
!>            tl_zu(1)=tl_zrhs(1)/zbet
!>
              ad_zrhs(1)=ad_zrhs(1)+ad_zu(1)/zbet
              ad_zu(1)=0.0_r8
!>            tl_zrhs(innLoop)=-tl_cg_QG(innLoop)-                      &
!>   &                         tl_cg_beta(innLoop)*zu(innLoop-1)-       &
!>   &                         tl_cg_delta(innLoop)*zu(innLoop)
!>
              ad_cg_QG(innLoop)=ad_cg_QG(innLoop)-ad_zrhs(innLoop)
              ad_cg_beta(innLoop)=ad_cg_beta(innLoop)-                  &
     &                            zu(innLoop-1)*ad_zrhs(innLoop)
              ad_cg_delta(innLoop)=ad_cg_delta(innLoop)-                &
     &                             zu(innLoop)*ad_zrhs(innLoop)
              ad_zrhs(innLoop)=0.0_r8

              DO ivec=2,innLoop-1
!>              tl_zrhs(ivec)=-tl_cg_QG(ivec)-                          &
!>   &                        tl_cg_beta(ivec)*zu(ivec-1)-              &
!>   &                        tl_cg_delta(ivec)*zu(ivec)-               &
!>   &                        tl_cg_beta(ivec+1)*zu(ivec+1)
!>
                ad_cg_QG(ivec)=ad_cg_QG(ivec)-ad_zrhs(ivec)
                ad_cg_beta(ivec)=ad_cg_beta(ivec)-                      &
     &                           zu(ivec-1)*ad_zrhs(ivec)
                ad_cg_delta(ivec)=ad_cg_delta(ivec)-                    &
     &                            zu(ivec)*ad_zrhs(ivec)
                ad_cg_beta(ivec+1)=ad_cg_beta(ivec+1)-                  &
     &                             zu(ivec+1)*ad_zrhs(ivec)
                ad_zrhs(ivec)=0.0_r8
              END DO
!>            tl_zrhs(1)=-tl_cg_QG(1)-                                  &
!>   &                   tl_cg_delta(1)*zu(1)-                          &
!>   &                   tl_cg_beta(2)*zu(2)
!>
              ad_cg_QG(1)=ad_cg_QG(1)-ad_zrhs(1)
              ad_cg_delta(1)=ad_cg_delta(1)-zu(1)*ad_zrhs(1)
              ad_cg_beta(2)=ad_cg_beta(2)-zu(2)*ad_zrhs(1)
              ad_zrhs(1)=0.0_r8
            END IF
          END IF

          DO iobs=1,Ndatum(ng)
!>          tl_cg_QG(innLoop+1)=tl_cg_QG(innLoop+1)+                    &
!>   &                          tl_zcglwk(iobs,innLoop+1)*              &
!>   &                          zgrad0(iobs,outLoop)+                   &
!>   &                          zcglwk(iobs,innLoop+1,outLoop)*         &
!>   &                          tl_zgrad0(iobs)
!>
            ad_zgrad0(iobs)=ad_zgrad0(iobs)+                            &
     &                      zcglwk(iobs,innLoop+1,outLoop)*             &
     &                      ad_cg_QG(innLoop+1)
            ad_zcglwk(iobs,innLoop+1)=ad_zcglwk(iobs,innLoop+1)+        &
     &                                zgrad0(iobs,outLoop)*             &
     &                                ad_cg_QG(innLoop+1)
          END DO
!>        tl_cg_QG(innLoop+1)=0.0_r8
!>
          ad_cg_QG(innLoop+1)=0.0_r8
          DO iobs=1,Ndatum(ng)
!>          tl_zcglwk(iobs,innLoop+1)=(tl_pgrad(iobs)-                  &
!>   &                                 tl_cg_beta(innLoop+1)*           &
!>   &                                 zcglwk(iobs,innLoop+1,outLoop))/ &
!>   &                                cg_beta(innLoop+1,outLoop)
!>
            adfac=ad_zcglwk(iobs,innLoop+1)/cg_beta(innLoop+1,outLoop)
            ad_cg_beta(innLoop+1)=ad_cg_beta(innLoop+1)-                &
     &                            zcglwk(iobs,innLoop+1,outLoop)*adfac
            ad_pgrad(iobs)=ad_pgrad(iobs)+adfac
            ad_zcglwk(iobs,innLoop+1)=0.0_r8
          END DO
!>        tl_cg_beta(innLoop+1)=0.5_r8*tl_cg_beta(innLoop+1)/           &
!>   &                          cg_beta(innLoop+1,outLoop)
!>
          ad_cg_beta(innLoop+1)=0.5_r8*ad_cg_beta(innLoop+1)/           &
     &                          cg_beta(innLoop+1,outLoop)
!
!  Compute appropriate value of pgrad.
!
          DO iobs=1,Ndatum(ng)
            pgrad(iobs)=ObsScale(iobs)*TLmodVal_S(iobs,innLoop,outLoop)
!
!  Convert gradient contribution from x-space to v-space.
!
            IF (ObsErr(iobs).NE.0.0_r8) THEN
              pgrad(iobs)=pgrad(iobs)/SQRT(ObsErr(iobs))
            END IF
          END DO
          DO iobs=1,Ndatum(ng)
            zt(iobs)=zcglwk(iobs,innLoop,outLoop)
          END DO
!
!  If preconditioning, convert zt from y-space to v-space.
!
!>        IF (Lprecond.and.(outLoop.gt.1)) THEN
!>          Lscale=2                 ! SQRT spectral LMP
!>          Ltrans=.FALSE.
!>          CALL rprecond (ng, Lscale, Ltrans, outLoop, NinnLoop, zt)
!>        END IF
!>
          DO iobs=1,Ndatum(ng)
            pgrad(iobs)=pgrad(iobs)+ObsScale(iobs)*zt(iobs)
          END DO
!
!  If preconditioning, convert pgrad from v-space to y-space.
!
!>        IF (Lprecond.and.(outLoop.gt.1)) THEN
!>          Lscale=2                 ! SQRT spectral LMP
!>          Ltrans=.TRUE.
!>          CALL rprecond (ng, Lscale, Ltrans, outLoop, NinnLoop, pgrad)
!>        END IF
!>
!
!  Compute the new Lanczos vector.
!
          DO iobs=1,Ndatum(ng)
            pgrad(iobs)=pgrad(iobs)-                                    &
     &                  cg_delta(innLoop,outLoop)*                      &
     &                  zcglwk(iobs,innLoop,outLoop)
          END DO
          IF (innLoop.gt.1) THEN
            DO iobs=1,Ndatum(ng)
              pgrad(iobs)=pgrad(iobs)-                                  &
     &                    cg_beta(innLoop,outLoop)*                     &
     &                    zcglwk(iobs,innLoop-1,outLoop)
            END DO
          END IF
!
!  Orthonormalize against previous Lanczos vectors.
!
          DO ivec=innLoop,1,-1
            DO iobs=1,Ndatum(ng)
              pgrad(iobs)=pgrad(iobs)-                                  &
     &                    cg_dla(ivec,outLoop)*                         &
     &                    zcglwk(iobs,ivec,outLoop)
            END DO
          END DO

          DO iobs=1,Ndatum(ng)
!>          tl_cg_beta(innLoop+1)=tl_cg_beta(innLoop+1)+                &
!>   &                            2.0_r8*tl_pgrad(iobs)*pgrad(iobs)
!>
            ad_pgrad(iobs)=ad_pgrad(iobs)+                              &
     &                     2.0_r8*pgrad(iobs)*ad_cg_beta(innLoop+1)
          END DO
!>        tl_cg_beta(innLoop+1)=0.0_r8
!>
          ad_cg_beta(innLoop+1)=0.0_r8
!
!  Compute appropriate value of pgrad.
!
          DO iobs=1,Ndatum(ng)
            pgrad(iobs)=ObsScale(iobs)*TLmodVal_S(iobs,innLoop,outLoop)
!
!  Convert gradient contribution from x-space to v-space.
!
            IF (ObsErr(iobs).NE.0.0_r8) THEN
              pgrad(iobs)=pgrad(iobs)/SQRT(ObsErr(iobs))
            END IF
          END DO
          DO iobs=1,Ndatum(ng)
            zt(iobs)=zcglwk(iobs,innLoop,outLoop)
          END DO
!
!  If preconditioning, convert zt from y-space to v-space.
!
!>        IF (Lprecond.and.(outLoop.gt.1)) THEN
!>          Lscale=2                 ! SQRT spectral LMP
!>          Ltrans=.FALSE.
!>          CALL rprecond (ng, Lscale, Ltrans, outLoop, NinnLoop, zt)
!>        END IF
!>
          DO iobs=1,Ndatum(ng)
            pgrad(iobs)=pgrad(iobs)+ObsScale(iobs)*zt(iobs)
          END DO
!
!  If preconditioning, convert pgrad from v-space to y-space.
!
!>        IF (Lprecond.and.(outLoop.gt.1)) THEN
!>          Lscale=2                 ! SQRT spectral LMP
!>          Ltrans=.TRUE.
!>          CALL rprecond (ng, Lscale, Ltrans, outLoop, NinnLoop, pgrad)
!>        END IF
!>
!
!  Compute the new Lanczos vector.
!
          DO iobs=1,Ndatum(ng)
            pgrad(iobs)=pgrad(iobs)-                                    &
     &                  cg_delta(innLoop,outLoop)*                      &
     &                  zcglwk(iobs,innLoop,outLoop)
          END DO
          IF (innLoop.gt.1) THEN
            DO iobs=1,Ndatum(ng)
              pgrad(iobs)=pgrad(iobs)-                                  &
     &                    cg_beta(innLoop,outLoop)*                     &
     &                    zcglwk(iobs,innLoop-1,outLoop)
            END DO
          END IF
!
!  Save the intermediate value of pgrad.
!
          DO iobs=1,Ndatum(ng)
            pgrad_S(iobs)=pgrad(iobs)
          END DO
!
!  Orthonormalize against previous Lanczos vectors.
!
!>        DO ivec=innLoop,1,-1
!>
          DO ivec=1,innLoop

            DO iobs=1,Ndatum(ng)
              pgrad(iobs)=pgrad_S(iobs)
            END DO
            DO i=ivec+1,innLoop
              DO iobs=1,Ndatum(ng)
                pgrad(iobs)=pgrad(iobs)-                                &
     &                      cg_dla(i,outLoop)*                          &
     &                      zcglwk(iobs,i,outLoop)
              END DO
            END DO

            ad_dla=0.0_r8
            DO iobs=1,Ndatum(ng)
!>            tl_pgrad(iobs)=tl_pgrad(iobs)-                            &
!>   &                       cg_dla(ivec,outLoop)*tl_zcglwk(iobs,ivec)- &
!>   &                       tl_dla*zcglwk(iobs,ivec,outLoop)
!>
              ad_zcglwk(iobs,ivec)=ad_zcglwk(iobs,ivec)-                &
     &                             cg_dla(ivec,outLoop)*ad_pgrad(iobs)
              ad_dla=ad_dla-zcglwk(iobs,ivec,outLoop)*ad_pgrad(iobs)
            END DO
            DO iobs=1,Ndatum(ng)
!>            tl_dla=tl_dla+                                            &
!>   &               tl_pgrad(iobs)*zcglwk(iobs,ivec,outLoop)+          &
!>   &               pgrad(iobs)*tl_zcglwk(iobs,ivec)
!>
              ad_zcglwk(iobs,ivec)=ad_zcglwk(iobs,ivec)+                &
     &                             pgrad(iobs)*ad_dla
              ad_pgrad(iobs)=ad_pgrad(iobs)+                            &
     &                       zcglwk(iobs,ivec,outLoop)*ad_dla
            END DO
!>          tl_dla=0.0_r8
!>
            ad_dla=0.0_r8
          END DO

          IF (innLoop.gt.1) THEN
            DO iobs=1,Ndatum(ng)
!>            tl_pgrad(iobs)=tl_pgrad(iobs)-                            &
!>   &                       tl_cg_beta(innLoop)*                       &
!>   &                       zcglwk(iobs,innLoop-1,outLoop)-            &
!>   &                       cg_beta(innLoop,outLoop)*                  &
!>   &                       tl_zcglwk(iobs,innLoop-1)
!>
              ad_zcglwk(iobs,innLoop-1)=ad_zcglwk(iobs,innLoop-1)-      &
     &                                  cg_beta(innLoop,outLoop)*       &
     &                                  ad_pgrad(iobs)
              ad_cg_beta(innLoop)=ad_cg_beta(innLoop)-                  &
     &                            zcglwk(iobs,innLoop-1,outLoop)*       &
     &                            ad_pgrad(iobs)
            END DO
          END IF
          DO iobs=1,Ndatum(ng)
!>          tl_pgrad(iobs)=tl_pgrad(iobs)-                              &
!>   &                     tl_cg_delta(innLoop)*                        &
!>   &                     zcglwk(iobs,innLoop,outLoop)-                &
!>   &                     cg_delta(innLoop,outLoop)*                   &
!>   &                     tl_zcglwk(iobs,innLoop)
!>
            ad_zcglwk(iobs,innLoop)=ad_zcglwk(iobs,innLoop)-            &
     &                              cg_delta(innLoop,outLoop)*          &
     &                              ad_pgrad(iobs)
            ad_cg_delta(innLoop)=ad_cg_delta(innLoop)-                  &
     &                           zcglwk(iobs,innLoop,outLoop)*          &
     &                           ad_pgrad(iobs)
          END DO
!
!  Compute appropriate value of pgrad.
!
          DO iobs=1,Ndatum(ng)
            pgrad(iobs)=ObsScale(iobs)*TLmodVal_S(iobs,innLoop,outLoop)
!
!  Convert gradient contribution from x-space to v-space.
!
            IF (ObsErr(iobs).NE.0.0_r8) THEN
              pgrad(iobs)=pgrad(iobs)/SQRT(ObsErr(iobs))
            END IF
          END DO
          DO iobs=1,Ndatum(ng)
            zt(iobs)=zcglwk(iobs,innLoop,outLoop)
          END DO
!
!  If preconditioning, convert zt from y-space to v-space.
!
!>        IF (Lprecond.and.(outLoop.gt.1)) THEN
!>          Lscale=2                 ! SQRT spectral LMP
!>          Ltrans=.FALSE.
!>          CALL rprecond (ng, Lscale, Ltrans, outLoop, NinnLoop, zt)
!>        END IF
!>
          DO iobs=1,Ndatum(ng)
            pgrad(iobs)=pgrad(iobs)+ObsScale(iobs)*zt(iobs)
          END DO
!
!  If preconditioning, convert pgrad from v-space to y-space.
!
!>        IF (Lprecond.and.(outLoop.gt.1)) THEN
!>          Lscale=2                 ! SQRT spectral LMP
!>          Ltrans=.TRUE.
!>          CALL rprecond (ng, Lscale, Ltrans, outLoop, NinnLoop, pgrad)
!>        END IF
!>
          DO iobs=1,Ndatum(ng)
!>          tl_cg_delta(innLoop)=tl_cg_delta(innLoop)+                  &
!>   &                           tl_zcglwk(iobs,innLoop)*pgrad(iobs)+   &
!>   &                           zcglwk(iobs,innLoop,outLoop)*          &
!>   &                           tl_pgrad(iobs)
!>
            ad_pgrad(iobs)=ad_pgrad(iobs)+                              &
     &                     zcglwk(iobs,innLoop,outLoop)*                &
     &                     ad_cg_delta(innLoop)
            ad_zcglwk(iobs,innLoop)=ad_zcglwk(iobs,innLoop)+            &
     &                              pgrad(iobs)*ad_cg_delta(innLoop)
          END DO
!>        tl_cg_delta(innLoop)=0.0_r8
!>
          ad_cg_delta(innLoop)=0.0_r8
!
!  If preconditioning, convert pgrad from v-space to y-space.
!
!>        IF (Lprecond.and.(outLoop.gt.1)) THEN
!>          Lscale=2                 ! SQRT spectral LMP
!>          Ltrans=.TRUE.
!>          CALL rprecond (ng, Lscale, Ltrans, outLoop, NinnLoop, pgrad)
!>        END IF
!>
          DO iobs=1,Ndatum(ng)
!>          tl_pgrad(iobs)=tl_pgrad(iobs)+ObsScale(iobs)*tl_zt(iobs)
!>
            ad_zt(iobs)=ad_zt(iobs)+ObsScale(iobs)*ad_pgrad(iobs)
          END DO
!
!  If preconditioning, convert zt from y-space to v-space.
!
!>        IF (Lprecond.and.(outLoop.gt.1)) THEN
!>          Lscale=2                 ! SQRT spectral LMP
!>          Ltrans=.FALSE.
!>          CALL rprecond (ng, Lscale, Ltrans, outLoop, NinnLoop, zt)
!>        END IF
!>
          DO iobs=1,Ndatum(ng)
!>          tl_zt(iobs)=tl_zcglwk(iobs,innLoop)
!>
            ad_zcglwk(iobs,innLoop)=ad_zcglwk(iobs,innLoop)+ad_zt(iobs)
            ad_zt(iobs)=0.0_r8
          END DO
          DO iobs=1,Ndatum(ng)
            IF (ObsErr(iobs).NE.0.0_r8) THEN
!>            tl_pgrad(iobs)=tl_pgrad(iobs)/SQRT(ObsErr(iobs))
!>
              ad_pgrad(iobs)=ad_pgrad(iobs)/SQRT(ObsErr(iobs))
            END IF
!>          tl_pgrad(iobs)=ObsScale(iobs)*tl_TLmodVal(iobs)
!<>         ad_TLmodVal(iobs)=ad_TLmodVal(iobs)+ObsScale(iobs)*         &
!<>  &                        ad_pgrad(iobs)
!>
            ADmodVal(iobs)=ADmodVal(iobs)+ObsScale(iobs)*ad_pgrad(iobs)
            ad_pgrad(iobs)=0.0_r8
          END DO

        END IF
      END IF MASTER_THREAD

# ifdef DISTRIBUTE
!
!  Broadcast new solution to other nodes.
!
      CALL mp_bcasti (ng, model, exit_flag)
      CALL mp_bcastf (ng, model, ADmodVal)
      CALL mp_bcastf (ng, model, ad_ObsVal)
      CALL mp_bcastf (ng, model, TLmodVal)
      CALL mp_bcasti (ng, model, info)
      CALL mp_bcastf (ng, model, ad_cg_beta(:))
      CALL mp_bcastf (ng, model, ad_cg_delta(:))
      CALL mp_bcastf (ng, model, ad_zcglwk(:,:))
      CALL mp_bcastf (ng, model, ad_cg_QG(:))
      CALL mp_bcastf (ng, model, ad_cg_Gnorm(:))
      CALL mp_bcastf (ng, model, ad_cg_pxsave(:))
      CALL mp_bcastf (ng, model, ad_cg_innov(:))
      CALL mp_bcastf (ng, model, ad_zgrad0(:))
# endif

      RETURN
      END SUBROUTINE ad_congrad
#else
      SUBROUTINE ad_congrad
      RETURN
      END SUBROUTINE ad_congrad
#endif

