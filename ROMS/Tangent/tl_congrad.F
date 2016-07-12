#include "cppdefs.h"

#if defined SENSITIVITY_4DVAR || defined TL_W4DPSAS || defined TL_W4DVAR
      SUBROUTINE tl_congrad (ng, model, outLoop, innLoop, NinnLoop,     &
     &                       Lcgini)
!
!svn $Id: tl_congrad.F 795 2016-05-11 01:42:43Z arango $
!=================================================== Andrew M. Moore ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group      Hernan G. Arango   !
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
      logical, intent(in) :: Lcgini
      integer, intent(in) :: ng, model, outLoop, innLoop, NinnLoop
!
!  Local variable declarations.
!
      logical :: Ltrans

      integer :: i, j, iobs, ivec, Lscale, info

      real(r8) :: dla, zbet
      real(r8) :: tl_dla
# ifdef MINRES
      real(r8) :: zsum, zck, zgk
      real(r8) :: tl_zsum, tl_zck, tl_zgk
# endif

      real(r8), dimension(NinnLoop) :: zu, zgam
      real(r8), dimension(NinnLoop) :: tl_zu, tl_zrhs
      real(r8), dimension(Ndatum(ng)) :: pgrad, zt
      real(r8), dimension(Ndatum(ng)) :: tl_px, tl_pgrad, tl_zt
# ifdef MINRES
      real(r8), dimension(innLoop,innLoop) :: ztriT, zLT, zLTt
      real(r8), dimension(innLoop,innLoop) :: tl_ztriT, tl_zLT
      real(r8), dimension(innLoop,innLoop) :: tl_zLTt
      real(r8), dimension(innLoop) :: tau, zwork1, ze, zeref
      real(r8), dimension(innLoop) :: tl_tau, tl_zwork1, tl_ze, tl_zeref
# endif
!
!=======================================================================
!  Weak constraint 4DVar conjugate gradient, Lanczos-based algorithm.
!=======================================================================
!
!  This conjugate gradient algorithm is not run in parallel since the
!  weak constraint is done in observation space. Mostly all the
!  variables are 1D arrays. Therefore, in parallel applications (only
!  distributed-memory is possible) the master node does the computations
!  and then broadcasts results to remaining nodes in the communicator.
!
!  This version of congrad solves A(x+x0)=b for x, by minimizing
!  J=0.5*x'Ax-x'(b+Ax0), where x0 is a first-guess estimate of the
!  representer coefficients from the previous outer-loop.
!  For the first outer-loop, x0=0. In the code x0=cg_pxsave and
!  x=px.
!
      Ltrans=.FALSE.

      MASTER_THREAD : IF (Master) THEN
!
!  Initialize cg_Gnorm. The TL of precond is not available.
!
        DO i=1,outLoop
          cg_Gnorm(i)=cg_Gnorm_v(i)
        END DO
!
!  Initialize internal parameters.  This is needed here for output
!  reasons.
!
        IF (innLoop.eq.0) THEN

# if defined W4DPSAS || defined TL_W4DPSAS
!
!  If this is the first inner-loop, save NLmodVal in BCKmodVal.
!
          DO iobs=1,Ndatum(ng)
            BCKmodVal(iobs)=NLmodVal(iobs)
          END DO
# endif
!
!  If this is the first outer-loop, clear the solution vector px.
!
          IF ((outLoop.eq.1).or.(.not.LhotStart)) THEN
!
!  For the first outer-loop, x0=0.
!
            DO iobs=1,Ndatum(ng)
              tl_px(iobs)=0.0_r8
              tl_cg_pxsave(iobs)=0.0_r8
            END DO
!
!  If this is the first pass of the inner loop, set up the vectors and
!  store the first guess. The initial starting guess is assumed to be
!  zero in which case the gradient is just: -(obs-model).
!  A first-level preconditioning is applied using the inverse of the
!  observation error standard deviations (i.e. sqrt(ObsErr)).
!
            DO iobs=1,Ndatum(ng)
# if defined W4DPSAS || defined TL_W4DPSAS
!>            pgrad(iobs)=-ObsScale(iobs)*                              &
!>   &                    (ObsVal(iobs)-BCKmodVal(iobs))
!>
              tl_pgrad(iobs)=-ObsScale(iobs)*tl_ObsVal(iobs)
# else
!>            pgrad(iobs)=-ObsScale(iobs)*                              &
!>   &                    (ObsVal(iobs)-TLmodVal(iobs))
!<>           tl_pgrad(iobs)=-ObsScale(iobs)*                           &
!<>  &                       (tl_ObsVal(iobs)-tl_TLmodVal(iobs))
!>
              tl_pgrad(iobs)=-ObsScale(iobs)*                           &
     &                       (tl_ObsVal(iobs)-TLmodVal(iobs))
# endif
!
! Convert pgrad from x-space to v-space.
!
              IF (ObsErr(iobs).NE.0.0_r8) THEN
!>              pgrad(iobs)=pgrad(iobs)/SQRT(ObsErr(iobs))
!>
                tl_pgrad(iobs)=tl_pgrad(iobs)/SQRT(ObsErr(iobs))
              END IF
!>            vgrad0(iobs)=pgrad(iobs)
!>
            END DO
!
!  If preconditioning, convert pgrad from v-space to y-space.
!
!>          IF (Lprecond.and.(outLoop.gt.1)) THEN
!>            Lscale=2                 ! SQRT spectral LMP
!>            Ltrans=.TRUE.
!>            CALL rprecond (ng, Lscale, Ltrans, outLoop, NinnLoop,     &
!>   &                       pgrad)
!>          END IF
!>
!>          cg_Gnorm(outLoop)=0.0_r8
!>
            tl_cg_Gnorm(outLoop)=0.0_r8
!>          vgnorm(outLoop)=0.0_r8
!>
            DO iobs=1,Ndatum(ng)
!>            zgrad0(iobs,outLoop)=pgrad(iobs)
!>
              tl_zgrad0(iobs)=tl_pgrad(iobs)
!>            cg_Gnorm(outLoop)=cg_Gnorm(outLoop)+                      &
!>   &                          zgrad0(iobs)*zgrad0(iobs)
!>
              tl_cg_Gnorm(outLoop)=tl_cg_Gnorm(outLoop)+                &
     &                             2.0_r8*tl_zgrad0(iobs)*              &
     &                             zgrad0(iobs,outLoop)
!>            vgnorm(outLoop)=vgnorm(outLoop)+vgrad0(iobs)*vgrad0(iobs)
!>
            END DO
!>          cg_Gnorm(outLoop)=SQRT(cg_Gnorm(outLoop))
!>
            tl_cg_Gnorm(outLoop)=0.5_r8*tl_cg_Gnorm(outLoop)/           &
     &                           cg_Gnorm(outLoop)
!>          vgnorm(outLoop)=SQRT(vgnorm(outLoop))
!>
            DO iobs=1,Ndatum(ng)
!>            zcglwk(iobs,1,outLoop)=pgrad(iobs)/cg_Gnorm(outLoop)
!>
              tl_zcglwk(iobs,1)=(tl_pgrad(iobs)-                        &
     &                           tl_cg_Gnorm(outLoop)*                  &
     &                           zcglwk(iobs,1,outLoop))/               &
     &                          cg_Gnorm(outLoop)
!>            ADmodVal(iobs)=zcglwk(iobs,1,outLoop)
!<>           tl_ADmodVal(iobs)=tl_zcglwk(iobs,1)
!>
              ADmodVal(iobs)=tl_zcglwk(iobs,1)
            END DO
!
!  If preconditioning, convert ADmodVal from y-space to v-space.
!
!>          IF (Lprecond.and.(outLoop.gt.1)) THEN
!>            Lscale=2                 ! SQRT spectral LMP
!>            Ltrans=.FALSE.
!>            CALL rprecond (ng, Lscale, Ltrans, outLoop, NinnLoop,     &
!>     &                     ADmodVal)
!>          END IF
!
! Convert ADmodVal from v-space to x-space.
!
            DO iobs=1,Ndatum(ng)
              IF (ObsErr(iobs).NE.0.0_r8) THEN
!>              ADmodVal(iobs)=ADmodVal(iobs)/SQRT(ObsErr(iobs))
!<>             tl_ADmodVal(iobs)=tl_ADmodVal(iobs)/SQRT(ObsErr(iobs))
!>
                ADmodVal(iobs)=ADmodVal(iobs)/SQRT(ObsErr(iobs))
              END IF
            END DO

!>          cg_QG(1,outLoop)=0.0_r8
!>
            tl_cg_QG(1)=0.0_r8
            DO iobs=1,Ndatum(ng)
!>            cg_QG(1,outLoop)=cg_QG(1,outLoop)+                        &
!>   &                         zcglwk(iobs,1,outLoop)*zgrad0(iobs)
!>
              tl_cg_QG(1)=tl_cg_QG(1)+                                  &
     &                    tl_zcglwk(iobs,1)*zgrad0(iobs,outLoop)+       &
     &                    zcglwk(iobs,1,outLoop)*tl_zgrad0(iobs)
            END DO

          ELSE
            IF (Lcgini) THEN
!
!  When outer>1 we need to evaluate Ax0 so for inner=0 we use
!  cg_pxsave in v-space as the starting vector.
!
              DO iobs=1,Ndatum(ng)
!>              ADmodVal(iobs)=cg_pxsave(iobs)
!<>             tl_ADmodVal(iobs)=tl_cg_pxsave(iobs)
!>
                ADmodVal(iobs)=tl_cg_pxsave(iobs)
# if defined W4DPSAS || defined TL_W4DPSAS
!>              cg_innov(iobs)=-ObsScale(iobs)*                         &
!>   &                         (ObsVal(iobs)-BCKmodVal(iobs))
!<>             tl_cg_innov(iobs)=0.0_r8
!>
                tl_cg_innov(iobs)=-ObsScale(iobs)*tl_ObsVal(iobs)
# else
!>              cg_innov(iobs)=-ObsScale(iobs)*                         &
!>   &                  (ObsVal(iobs)-TLmodVal(iobs))
!<>             tl_cg_innov(iobs)=-ObsScale(iobs)*                      &
!<>  &                            (tl_ObsVal(iobs)-tl_TLmodVal(iobs))
!>
                tl_cg_innov(iobs)=-ObsScale(iobs)*                      &
     &                            (tl_ObsVal(iobs)-TLmodVal(iobs))
# endif
              END DO
!
!  Convert ADmodVal from v-space to x-space and cg_innov (the
!  contribution to the starting gradient) from x-space to v-space.
!
              DO iobs=1,Ndatum(ng)
                IF (ObsErr(iobs).NE.0.0_r8) THEN
!>                ADmodVal(iobs)=ADmodVal(iobs)/SQRT(ObsErr(iobs))
!<>               tl_ADmodVal(iobs)=tl_ADmodVal(iobs)/SQRT(ObsErr(iobs))
!>
                  ADmodVal(iobs)=ADmodVal(iobs)/SQRT(ObsErr(iobs))
!>                cg_innov(iobs)=cg_innov(iobs)/SQRT(ObsErr(iobs))
!>
                  tl_cg_innov(iobs)=tl_cg_innov(iobs)/SQRT(ObsErr(iobs))
                END IF
              END DO

            ELSE

              DO iobs=1,Ndatum(ng)
!
!  Convert gradient contribution from x-space to v-space.
!
!>              pgrad(iobs)=ObsScale(iobs)*TLmodVal(iobs)
!>              tl_pgrad(iobs)=ObsScale(iobs)*tl_TLmodVal(iobs)
!>
                tl_pgrad(iobs)=ObsScale(iobs)*TLmodVal(iobs)
                IF (ObsErr(iobs).NE.0.0_r8) THEN
!>                pgrad(iobs)=pgrad(iobs)/SQRT(ObsErr(iobs))
!>
                  tl_pgrad(iobs)=tl_pgrad(iobs)/SQRT(ObsErr(iobs))
                END IF
!
!  Add I*x0=cg_pxsave contribution to the gradient and the term
!  -b=cg_innov (everything is in v-space at this point).
!
!>              pgrad(iobs)=pgrad(iobs)+ObsScale(iobs)*                 &
!>   &                      (cg_pxsave(iobs)+cg_innov(iobs))
!>
                tl_pgrad(iobs)=tl_pgrad(iobs)+ObsScale(iobs)*           &
     &                         (tl_cg_pxsave(iobs)+tl_cg_innov(iobs))
!>              vgrad0(iobs)=pgrad(iobs)
!>
              END DO
!
!  If preconditioning, convert pgrad from v-space to y-space.
!
!>            IF (Lprecond.and.(outLoop.gt.1)) THEN
!>              Lscale=2                 ! SQRT spectral LMP
!>              Ltrans=.TRUE.
!>              CALL rprecond(ng, Lscale, Ltrans, outLoop, NinnLoop,    &
!>   &                        pgrad)
!>            END IF
!>
!>            cg_Gnorm(outLoop)=0.0_r8
!>
              tl_cg_Gnorm(outLoop)=0.0_r8
!>            vgnorm(outLoop)=0.0_r8
!>
              DO iobs=1,Ndatum(ng)
!>              zgrad0(iobs,outLoop)=pgrad(iobs)
!>
                tl_zgrad0(iobs)=tl_pgrad(iobs)
!>              cg_Gnorm(outLoop)=cg_Gnorm(outLoop)+                    &
!>   &                            zgrad0(iobs,outLoop)*                 &
!>   &                            zgrad0(iobs,outLoop)
!>
                tl_cg_Gnorm(outLoop)=tl_cg_Gnorm(outLoop)+              &
     &                               2.0_r8*tl_zgrad0(iobs)*            &
     &                               zgrad0(iobs,outLoop)
!>              vgnorm(outLoop)=vgnorm(outLoop)+                        &
!>   &                          vgrad0(iobs)*vgrad0(iobs)
!>
              END DO
!>            cg_Gnorm(outLoop)=SQRT(cg_Gnorm(outLoop))
!>
              tl_cg_Gnorm(outLoop)=0.5_r8*tl_cg_Gnorm(outLoop)/         &
     &                             cg_Gnorm(outLoop)
!>            vgnorm(outLoop)=SQRT(vgnorm(outLoop))
!>
              DO iobs=1,Ndatum(ng)
!>              zcglwk(iobs,1,outLoop)=pgrad(iobs)/cg_Gnorm(outLoop)
!>
                tl_zcglwk(iobs,1)=(tl_pgrad(iobs)-                      &
     &                             tl_cg_Gnorm(outLoop)*                &
     &                             zcglwk(iobs,1,outLoop))/             &
     &                            cg_Gnorm(outLoop)
!>              ADmodVal(iobs)=zcglwk(iobs,1,outLoop)
!<>             tl_ADmodVal(iobs)=tl_zcglwk(iobs,1)
!>
                ADmodVal(iobs)=tl_zcglwk(iobs,1)
              END DO
!
!  If preconditioning, convert ADmodVal from y-space to v-space.
!
!>            IF (Lprecond.and.(outLoop.gt.1)) THEN
!>              Lscale=2                 ! SQRT spectral LMP
!>              Ltrans=.FALSE.
!>              CALL rprecond(ng, Lscale, Ltrans, outLoop, NinnLoop,     &
!>   &                        ADmodVal)
!>            END IF
!>
              DO iobs=1,Ndatum(ng)
                IF (ObsErr(iobs).NE.0.0_r8) THEN
!>                ADmodVal(iobs)=ADmodVal(iobs)/SQRT(ObsErr(iobs))
!<>               tl_ADmodVal(iobs)=tl_ADmodVal(iobs)/SQRT(ObsErr(iobs))
!>
                  ADmodVal(iobs)=ADmodVal(iobs)/SQRT(ObsErr(iobs))
                END IF
              END DO
!>            cg_QG(1,outLoop)=0.0_r8
!>
              tl_cg_QG(1)=0.0_r8
              DO iobs=1,Ndatum(ng)
!>              cg_QG(1,outLoop)=cg_QG(1,outLoop)+                      &
!>   &                           zcglwk(iobs,1,outLoop)*                &
!>   &                           zgrad0(iobs,outLoop)
!>
                tl_cg_QG(1)=tl_cg_QG(1)+                                &
     &                      tl_zcglwk(iobs,1)*zgrad0(iobs,outLoop)+     &
     &                      zcglwk(iobs,1,outLoop)*tl_zgrad0(iobs)
              END DO

            END IF

          END IF

        ELSE
!
!  After the initialization, for every other inner loop, calculate a
!  new Lanczos vector, store it in the matrix, and update search.
!
          DO iobs=1,Ndatum(ng)
            pgrad(iobs)=ObsScale(iobs)*TLmodVal_S(iobs,innLoop,outLoop)
!<>         tl_pgrad(iobs)=ObsScale(iobs)*tl_TLmodVal(iobs)
            tl_pgrad(iobs)=ObsScale(iobs)*TLmodVal(iobs)
!
!  Convert gradient contribution from x-space to v-space.
!
            IF (ObsErr(iobs).NE.0.0_r8) THEN
              pgrad(iobs)=pgrad(iobs)/SQRT(ObsErr(iobs))
              tl_pgrad(iobs)=tl_pgrad(iobs)/SQRT(ObsErr(iobs))
            END IF
          END DO

          DO iobs=1,Ndatum(ng)
            zt(iobs)=zcglwk(iobs,innLoop,outLoop)
            tl_zt(iobs)=tl_zcglwk(iobs,innLoop)
          END DO
!
!  If preconditioning, convert zt from y-space to v-space.
!
!>        IF (Lprecond.and.(outLoop.gt.1)) THEN
!>          Lscale=2                 ! SQRT spectral LMP
!>          Ltrans=.FALSE.
!>          CALL rprecond(ng, Lscale, Ltrans, outLoop, NinnLoop, zt)
!>        END IF
!>
          DO iobs=1,Ndatum(ng)
            pgrad(iobs)=pgrad(iobs)+ObsScale(iobs)*zt(iobs)
            tl_pgrad(iobs)=tl_pgrad(iobs)+ObsScale(iobs)*tl_zt(iobs)
          END DO
!
!  If preconditioning, convert pgrad from v-space to y-space.
!
!>        IF (Lprecond.and.(outLoop.gt.1)) THEN
!>          Lscale=2                 ! SQRT spectral LMP
!>          Ltrans=.TRUE.
!>          CALL rprecond(ng, Lscale, Ltrans, outLoop, NinnLoop, pgrad)
!>        END IF
!>
!>        cg_delta(innLoop,outLoop)=0.0_r8
!>
          tl_cg_delta(innLoop)=0.0_r8
          DO iobs=1,Ndatum(ng)
!>          cg_delta(innLoop,outLoop)=cg_delta(innLoop,outLoop)+        &
!>   &                                zcglwk(iobs,innLoop,outLoop)*     &
!>   &                                pgrad(iobs)
!>
            tl_cg_delta(innLoop)=tl_cg_delta(innLoop)+                  &
     &                           tl_zcglwk(iobs,innLoop)*pgrad(iobs)+   &
     &                           zcglwk(iobs,innLoop,outLoop)*          &
     &                           tl_pgrad(iobs)
          END DO
!
!  Exit, if not a positive definite algorithm.
!
!>        IF (cg_delta(innLoop,outLoop).le.0.0_r8) THEN
!>          WRITE (stdout,20) cg_delta(innLoop,outLoop), outLoop,       &
!>   &                        innLoop
!>          exit_flag=8
!>          GO TO 10
!>        END IF
!>
!
!  Compute the new Lanczos vector.
!
          DO iobs=1,Ndatum(ng)
            pgrad(iobs)=pgrad(iobs)-                                    &
                        cg_delta(innLoop,outLoop)*                      &
     &                  zcglwk(iobs,innLoop,outLoop)
            tl_pgrad(iobs)=tl_pgrad(iobs)-                              &
     &                     tl_cg_delta(innLoop)*                        &
     &                     zcglwk(iobs,innLoop,outLoop)-                &
     &                     cg_delta(innLoop,outLoop)*                   &
     &                     tl_zcglwk(iobs,innLoop)
          END DO
          IF (innLoop.gt.1) THEN
            DO iobs=1,Ndatum(ng)
              pgrad(iobs)=pgrad(iobs)-                                  &
     &                    cg_beta(innLoop,outLoop)*                     &
     &                    zcglwk(iobs,innLoop-1,outLoop)
              tl_pgrad(iobs)=tl_pgrad(iobs)-                            &
     &                       tl_cg_beta(innLoop)*                       &
     &                       zcglwk(iobs,innLoop-1,outLoop)-            &
     &                       cg_beta(innLoop,outLoop)*                  &
     &                       tl_zcglwk(iobs,innLoop-1)
            END DO
          END IF
!
!  Orthonormalize against previous Lanczos vectors.
!
          DO ivec=innLoop,1,-1
!>          cg_dla(ivec,outLoop)=0.0_r8
!>
            tl_dla=0.0_r8
            DO iobs=1,Ndatum(ng)
!>            cg_dla(ivec,outLoop)=cg_dla(ivec,outLoop)+                &
!>   &                             pgrad(iobs)*                         &
!>   &                             zcglwk(iobs,ivec,outLoop)
!>
              tl_dla=tl_dla+                                            &
     &               tl_pgrad(iobs)*zcglwk(iobs,ivec,outLoop)+          &
     &               pgrad(iobs)*tl_zcglwk(iobs,ivec)
            END DO
            DO iobs=1,Ndatum(ng)
              pgrad(iobs)=pgrad(iobs)-                                  &
     &                    cg_dla(ivec,outLoop)*                         &
     &                    zcglwk(iobs,ivec,outLoop)
              tl_pgrad(iobs)=tl_pgrad(iobs)-                            &
     &                       cg_dla(ivec,outLoop)*tl_zcglwk(iobs,ivec)- &
     &                       tl_dla*zcglwk(iobs,ivec,outLoop)
            END DO
          END DO
!
!>        cg_beta(innLoop+1,outLoop)=0.0_r8
!>
          tl_cg_beta(innLoop+1)=0.0_r8
          DO iobs=1,Ndatum(ng)
!>          cg_beta(innLoop+1,outLoop)=cg_beta(innLoop+1,outLoop)+      &
!>   &                                 pgrad(iobs)*pgrad(iobs)
!>
            tl_cg_beta(innLoop+1)=tl_cg_beta(innLoop+1)+                &
     &                            2.0_r8*tl_pgrad(iobs)*pgrad(iobs)
          END DO
!>        cg_beta(innLoop+1,outLoop)=SQRT(cg_beta(innLoop+1,outLoop))
!>
          tl_cg_beta(innLoop+1)=0.5_r8*tl_cg_beta(innLoop+1)/           &
     &                          cg_beta(innLoop+1,outLoop)
!
          DO iobs=1,Ndatum(ng)
!>          zcglwk(iobs,innLoop+1,outLoop)=pgrad(iobs)/                 &
!>   &                                     cg_beta(innLoop+1,outLoop)
!>
            tl_zcglwk(iobs,innLoop+1)=(tl_pgrad(iobs)-                  &
     &                                 tl_cg_beta(innLoop+1)*           &
     &                                 zcglwk(iobs,innLoop+1,outLoop))/ &
     &                                cg_beta(innLoop+1,outLoop)
          END DO
!
!>        cg_QG(innLoop+1,outLoop)=0.0_r8
!>
          tl_cg_QG(innLoop+1)=0.0_r8
          DO iobs=1,Ndatum(ng)
!>          cg_QG(innLoop+1,outLoop)=cg_QG(innLoop+1,outLoop)+          &
!>   &                               zcglwk(iobs,innLoop+1,outLoop)*
!>   &                               zgrad0(iobs)
!>
            tl_cg_QG(innLoop+1)=tl_cg_QG(innLoop+1)+                    &
     &                          tl_zcglwk(iobs,innLoop+1)*              &
     &                          zgrad0(iobs,outLoop)+                   &
     &                          zcglwk(iobs,innLoop+1,outLoop)*         &
     &                          tl_zgrad0(iobs)
          END DO
          IF (innLoop.eq.NinnLoop) THEN
# ifdef MINRES
!
!  Use the minimum residual method as described by Paige and Saunders
!  ("Sparse Indefinite Systems of Linear Equations", 1975, SIAM Journal
!  on Numerical Analysis, 617-619). Specifically we refer to equations
!  6.10 and 6.11 of this paper.
!
!  Perform a LQ factorization of the tridiagonal matrix.
!
          ztriT=0.0_r8
          tl_ztriT=0.0_r8
          DO i=1,innLoop
            ztriT(i,i)=cg_delta(i,outLoop)
            tl_ztriT(i,i)=tl_cg_delta(i)
          END DO
          DO i=1,innLoop-1
            ztriT(i,i+1)=cg_beta(i+1,outLoop)
            tl_ztriT(i,i+1)=tl_cg_beta(i+1)
          END DO
          DO i=2,innLoop
            ztriT(i,i-1)=cg_beta(i,outLoop)
            tl_ztriT(i,i-1)=tl_cg_beta(i)
          END DO
!
!  Note: tl_sqlq also computes the LQ factorization of ztriT.
!
          CALL tl_sqlq(innLoop, ztriT, tl_ztriT, tau, tl_tau, zwork1,   &
     &                 tl_zwork1)
!
!   Isolate L=zLT and its transpose.
!
          zLT=0.0_r8
          tl_zLT=0.0_r8
          zLTt=0.0_r8
          tl_zLTt=0.0_r8
          DO i=1,innLoop
            DO j=1,i
              zLT(i,j)=ztriT(i,j)
              tl_zLT(i,j)=tl_ztriT(i,j)
            END DO
          END DO
          DO j=1,innLoop
            DO i=1,innLoop
              zLTt(i,j)=zLT(j,i)
              tl_zLTt(i,j)=tl_zLT(j,i)
            END DO
          END DO
!
!   Now form ze=beta1*Q*e1.
!
          ze=0.0_r8
          tl_ze=0.0_r8
          DO i=1,innLoop
            ze(i)=-cg_QG(i,outLoop)
            tl_ze(i)=-tl_cg_QG(i)
          END DO
          DO i=1,innLoop
            DO j=1,innLoop
              zeref(j)=0.0_r8
              tl_zeref(j)=0.0_r8
            END DO
            zeref(i)=1.0_r8
            tl_zeref(i)=0.0_r8
            DO j=i+1,innLoop
              zeref(j)=ztriT(i,j)
              tl_zeref(j)=tl_ztriT(i,j)
            END DO
            zsum=0.0_r8
            tl_zsum=0.0_r8
            DO j=1,innLoop
              zsum=zsum+ze(j)*zeref(j)
              tl_zsum=tl_zsum+tl_ze(j)*zeref(j)+ze(j)*tl_zeref(j)
            END DO
            DO j=1,innLoop
              ze(j)=ze(j)-tau(i)*zsum*zeref(j)
              tl_ze(j)=tl_ze(j)-tl_tau(i)*zsum*zeref(j)-                &
     &                          tau(i)*tl_zsum*zeref(j)-                &
     &                          tau(i)*zsum*tl_zeref(j)
            END DO
          END DO
!
!   Now form ze=D*ze (using equation 5.6 and 6.5 also).
!
          zgk=SQRT(zLT(innLoop,innLoop)*zLT(innLoop,innLoop)+           &
     &        cg_beta(innLoop+1,outLoop)*cg_beta(innLoop+1,outLoop))
          IF (zgk.GT.0.0_r8) THEN
            tl_zgk=(tl_zLT(innLoop,innLoop)*zLT(innLoop,innLoop)+       &
     &              tl_cg_beta(innLoop+1)*cg_beta(innLoop+1,outLoop))/  &
     &             zgk
          ELSE
            tl_zgk=0.0_r8
          ENDIF
          zck=zLT(innLoop,innLoop)/zgk
          tl_zck=tl_zLT(innLoop,innLoop)/zgk-tl_zgk*zck/zgk
          ze(innLoop)=zck*ze(innLoop)
          tl_ze(innLoop)=tl_zck*ze(innLoop)+zck*tl_ze(innLoop)
!
!   Now compute tl_ze=inv(L')*(tl_ze-tl_L'*ze).
!
!   First solve for ze=inv(L')*ze.
!
          DO j=innLoop,1,-1
            ze(j)=ze(j)/zLTt(j,j)
            DO i=1,j-1
              ze(i)=ze(i)-ze(j)*zLTt(i,j)
            END DO
          END DO
!
!   Next compute tl_rhs=tl_L'*ze then subtract from tl_ze.
!
          DO i=1,innLoop
            tl_zrhs(i)=0.0_r8
            DO j=1,innLoop
              tl_zrhs(i)=tl_zrhs(i)+tl_zLTt(i,j)*ze(j)
            END DO
            tl_ze(i)=tl_ze(i)-tl_zrhs(i)
          END DO
!
!   Now solve the linear triangular system.
!
          DO j=innLoop,1,-1
            tl_ze(j)=tl_ze(j)/zLTt(j,j)
            DO i=1,j-1
              tl_ze(i)=tl_ze(i)-tl_ze(j)*zLTt(i,j)
            END DO
          END DO
!
!   Copy the solution ze into zu.
!
          DO i=1,innLoop
            zu(i)=ze(i)
            tl_zu(i)=tl_ze(i)
          END DO
# else
!
!  Calculate the new solution based upon the new, orthonormalized
!  Lanczos vector. First, the tridiagonal system is solved by
!  decomposition and forward/backward substitution.
!
            IF (NinnLoop.eq.1) THEN
              zu(1)=-cg_QG(1,outLoop)/cg_delta(1,outLoop)
              tl_zrhs(1)=-tl_cg_QG(1)-tl_cg_delta(1)*zu(1)
              tl_zu(1)=tl_zrhs(1)/cg_delta(1,outLoop)
            ELSE
!
!  Compute zu first.
!
              zbet=cg_delta(1,outLoop)
              zu(1)=-cg_QG(1,outLoop)/zbet
              DO ivec=2,innLoop
                zgam(ivec)=cg_beta(ivec,outLoop)/zbet
                zbet=cg_delta(ivec,outLoop)-                            &
     &               cg_beta(ivec,outLoop)*zgam(ivec)
                zu(ivec)=(-cg_QG(ivec,outLoop)-cg_beta(ivec,outLoop)*   &
     &                    zu(ivec-1))/zbet
              END DO
              DO ivec=innLoop-1,1,-1
                zu(ivec)=zu(ivec)-zgam(ivec+1)*zu(ivec+1)
              END DO
!
!  Now compute tl_zrhs.
!
              tl_zrhs(1)=-tl_cg_QG(1)-                                  &
     &                   tl_cg_delta(1)*zu(1)-                          &
     &                   tl_cg_beta(2)*zu(2)
              DO ivec=2,innLoop-1
                tl_zrhs(ivec)=-tl_cg_QG(ivec)-                          &
     &                        tl_cg_beta(ivec)*zu(ivec-1)-              &
     &                        tl_cg_delta(ivec)*zu(ivec)-               &
     &                        tl_cg_beta(ivec+1)*zu(ivec+1)
              END DO
              tl_zrhs(innLoop)=-tl_cg_QG(innLoop)-                      &
     &                         tl_cg_beta(innLoop)*zu(innLoop-1)-       &
     &                         tl_cg_delta(innLoop)*zu(innLoop)
!
!  Now solve the TL tridiagonal system A*dx=b-dA*x
!
              zbet=cg_delta(1,outLoop)
              tl_zu(1)=tl_zrhs(1)/zbet
              DO ivec=2,innLoop
                zgam(ivec)=cg_beta(ivec,outLoop)/zbet
                zbet=cg_delta(ivec,outLoop)-                            &
     &               cg_beta(ivec,outLoop)*zgam(ivec)
                tl_zu(ivec)=(tl_zrhs(ivec)-cg_beta(ivec,outLoop)*       &
     &                       tl_zu(ivec-1))/zbet
              END DO

              DO ivec=innLoop-1,1,-1
!>              zu(ivec)=zu(ivec)-zgam(ivec+1)*zu(ivec+1)
!>
                tl_zu(ivec)=tl_zu(ivec)-zgam(ivec+1)*tl_zu(ivec+1)
              END DO
            END IF

!>          DO iobs=1,Ndatum(ng)
!>            zw(iobs)=zgrad0(iobs)+                                    &
!>   &                 cg_beta(innLoop+1,outLoop)*                      &
!>   &                 zcglwk(iobs,innLoop+1,outLoop)*zwork(innLoop,3)
!>          END DO
# endif

            DO iobs=1,Ndatum(ng)
!>            px(iobs)=0.0_r8
!>
              tl_px(iobs)=0.0_r8
              DO ivec=1,innLoop
!>              px(iobs)=px(iobs)+                                      &
!>   &                   zcglwk(iobs,ivec,outLoop)*zu(ivec)
!>
                tl_px(iobs)=tl_px(iobs)+                                &
     &                      tl_zcglwk(iobs,ivec)*zu(ivec)+              &
     &                      zcglwk(iobs,ivec,outLoop)*tl_zu(ivec)
!>              zw(iobs)=zw(iobs)-                                      &
!>   &                   zcglwk(iobs,ivec,outLoop)*cg_QG(ivec,outLoop)
!>
              END DO
            END DO
!
!  If preconditioning, convert px from y-space to v-space.
!  We will always keep px in v-space.
!
!>          IF (Lprecond.and.(outLoop.gt.1)) THEN
!>            Lscale=2                 ! SQRT spectral LMP
!>            Ltrans=.FALSE.
!>            CALL rprecond(ng, Lscale, Ltrans, outLoop, NinnLoop, px)
!>          END IF
!>
          END IF
!
!  Put the new trial solution into the adjoint vector for the next loop
!  Put the final solution into the adjoint vector when converged
!  of on the final inner-loop.
!
          IF ((innLoop.eq.NinnLoop)) THEN
!
!  Note: px is already in v-space so there is no need to convert
!  if preconditioning. cg_pxsave is also in v-space.
!
            DO iobs=1,Ndatum(ng)
!>            ADmodVal(iobs)=px(iobs)
!<>           tl_ADmodVal(iobs)=tl_px(iobs)
!>
              ADmodVal(iobs)=tl_px(iobs)
            END DO
            IF (LhotStart) THEN
              DO iobs=1,Ndatum(ng)
!>              ADmodVal(iobs)=ADmodVal(iobs)+cg_pxsave(iobs)
!<>             tl_ADmodVal(iobs)=tl_ADmodVal(iobs)+tl_cg_pxsave(iobs)
!>
                ADmodVal(iobs)=ADmodVal(iobs)+tl_cg_pxsave(iobs)
              END DO
              DO iobs=1,Ndatum(ng)
!>              cg_pxsave(iobs)=ADmodVal(iobs)
!<>             tl_cg_pxsave(iobs)=tl_ADmodVal(iobs)
!>
                tl_cg_pxsave(iobs)=ADmodVal(iobs)
              END DO
            END IF
          ELSE
            DO iobs=1,Ndatum(ng)
!>            ADmodVal(iobs)=zcglwk(iobs,innLoop+1,outLoop)
!<>           tl_ADmodVal(iobs)=tl_zcglwk(iobs,innLoop+1)
!>
              ADmodVal(iobs)=tl_zcglwk(iobs,innLoop+1)
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
          END IF
!
!  Convert ADmodVal from v-space to x-space.
!
          DO iobs=1,Ndatum(ng)
            IF (ObsErr(iobs).NE.0.0_r8) THEN
!>            ADmodVal(iobs)=ADmodVal(iobs)/SQRT(ObsErr(iobs))
!<>           tl_ADmodVal(iobs)=tl_ADmodVal(iobs)/SQRT(ObsErr(iobs))
!>
              ADmodVal(iobs)=ADmodVal(iobs)/SQRT(ObsErr(iobs))
            END IF
          END DO
        END IF

      END IF MASTER_THREAD

# ifdef DISTRIBUTE
!
!  Broadcast new solution to other nodes.
!
      CALL mp_bcasti (ng, model, exit_flag)
      CALL mp_bcastf (ng, model, ADmodVal)
      CALL mp_bcastf (ng, model, tl_cg_QG(:))
      CALL mp_bcastf (ng, model, tl_cg_Gnorm(:))
      CALL mp_bcastf (ng, model, tl_cg_pxsave(:))
      CALL mp_bcastf (ng, model, tl_cg_innov(:))
      CALL mp_bcastf (ng, model, tl_cg_beta(:))
      CALL mp_bcastf (ng, model, tl_cg_delta(:))
      CALL mp_bcastf (ng, model, tl_zcglwk(:,:))
# endif
!
      RETURN
      END SUBROUTINE tl_congrad
#else
      SUBROUTINE tl_congrad
      RETURN
      END SUBROUTINE tl_congrad
#endif
