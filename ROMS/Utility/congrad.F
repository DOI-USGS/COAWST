#include "cppdefs.h"

#if defined WEAK_CONSTRAINT && !defined RPCG
      SUBROUTINE congrad (ng, model, outLoop, innLoop, NinnLoop,        &
     &                    Lcgini)
!
!svn $Id: congrad.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Weak Constraint 4-Dimensional Variational (4DVar) Pre-conditioned   !
!                     Conjugate Gradient Algorithm                     !

# if defined TL_W4DVAR || defined W4DVAR || defined W4DVAR_SENSITIVITY
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

# elif defined TL_W4DPSAS || \
       defined W4DPSAS    || defined W4DPSAS_SENSITIVITY
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
      integer :: i, ic, j, iobs, ivec, Lscale, info

      real(r8) :: zbet, eps, preducv, preducy
# ifdef MINRES
      real(r8) :: zsum, zck, zgk
# endif
      real(r8) :: Jopt, Jf, Jmod, Jdata, Jb, Jobs, Jact, cff

      real(r8), dimension(NinnLoop) :: zu, zgam
      real(r8), dimension(NinnLoop) :: ztemp1, ztemp2, zu1, zu2
      real(r8), dimension(Ndatum(ng)) :: px, pgrad, zw, zt
      real(r8), dimension(Ndatum(ng)) :: zhv, zht, zd
      real(r8), dimension(Ninner,3) :: zwork
      real(r8), dimension(2*(NinnLoop-1)) :: work
      real(r8), dimension(Ninner,Ninner) :: zgv
# ifdef MINRES
      real(r8), dimension(innLoop,innLoop) :: ztriT, zLT, zLTt
      real(r8), dimension(innLoop) :: tau, zwork1, ze, zeref
# endif

      character (len=13) :: string
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
!
      MASTER_THREAD : IF (Master) THEN
!
!  Initialize internal parameters.  This is needed here for output
!  reasons.
!
        Jopt=0.0_r8
        Jf=0.0_r8
        Jdata=0.0_r8
        Jmod=0.0_r8
        Jb=0.0_r8
        Jobs=0.0_r8
        Jact=0.0_r8
        eps=0.0_r8
        preducv=0.0_r8
        preducy=0.0_r8
!
!-----------------------------------------------------------------------
!  Initialization step, innloop = 0.
!-----------------------------------------------------------------------
!
        MINIMIZE : IF (innLoop.eq.0) THEN

# if defined W4DPSAS             || defined TL_W4DPSAS || \
     defined W4DPSAS_SENSITIVITY
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
              px(iobs)=0.0_r8
              cg_pxsave(iobs)=0.0_r8
            END DO
!
!  If this is the first pass of the inner loop, set up the vectors and
!  store the first guess. The initial starting guess is assumed to be
!  zero in which case the gradient is just: -(obs-model).
!  A first-level preconditioning is applied using the inverse of the
!  observation error standard deviations (i.e. sqrt(ObsErr)).
!
            DO iobs=1,Ndatum(ng)
# if defined W4DPSAS             || defined TL_W4DPSAS || \
     defined W4DPSAS_SENSITIVITY
              pgrad(iobs)=-ObsScale(iobs)*                              &
     &                    (ObsVal(iobs)-BCKmodVal(iobs))
# else
              pgrad(iobs)=-ObsScale(iobs)*                              &
     &                    (ObsVal(iobs)-TLmodVal(iobs))
# endif
!
! Convert pgrad from x-space to v-space.
!
              IF (ObsErr(iobs).NE.0.0_r8) THEN
                pgrad(iobs)=pgrad(iobs)/SQRT(ObsErr(iobs))
              END IF
              vgrad0(iobs)=pgrad(iobs)
              vgrad0s(iobs)=-vgrad0(iobs)
            END DO
!
!  If preconditioning, convert pgrad from v-space to y-space.
!
            IF (Lprecond.and.(outLoop.gt.1)) THEN
              Lscale=2                 ! SQRT spectral LMP
              Ltrans=.TRUE.
              CALL rprecond (ng, Lscale, Ltrans, outLoop, NinnLoop,     &
     &                       pgrad)
            END IF
!
            cg_Gnorm_y(outLoop)=0.0_r8
            cg_Gnorm_v(outLoop)=0.0_r8
            DO iobs=1,Ndatum(ng)
              zgrad0(iobs,outLoop)=pgrad(iobs)
              cg_Gnorm_y(outLoop)=cg_Gnorm_y(outLoop)+                  &
     &                            zgrad0(iobs,outLoop)*                 &
     &                            zgrad0(iobs,outLoop)
              cg_Gnorm_v(outLoop)=cg_Gnorm_v(outLoop)+                  &
     &                            vgrad0(iobs)*vgrad0(iobs)
            END DO
            cg_Gnorm_y(outLoop)=SQRT(cg_Gnorm_y(outLoop))
            cg_Gnorm_v(outLoop)=SQRT(cg_Gnorm_v(outLoop))
            DO iobs=1,Ndatum(ng)
              zcglwk(iobs,1,outLoop)=pgrad(iobs)/cg_Gnorm_y(outLoop)
              ADmodVal(iobs)=zcglwk(iobs,1,outLoop)
            END DO
!
!  If preconditioning, convert ADmodVal from y-space to v-space.
!
            IF (Lprecond.and.(outLoop.gt.1)) THEN
              Lscale=2                 ! SQRT spectral LMP
              Ltrans=.FALSE.
              CALL rprecond (ng, Lscale, Ltrans, outLoop, NinnLoop,     &
     &                       ADmodVal)
            END IF
!
! Convert ADmodVal from v-space to x-space.
!
            DO iobs=1,Ndatum(ng)
              IF (ObsErr(iobs).NE.0.0_r8) THEN
                ADmodVal(iobs)=ADmodVal(iobs)/SQRT(ObsErr(iobs))
              END IF
            END DO

            cg_QG(1,outLoop)=0.0_r8
            DO iobs=1,Ndatum(ng)
               cg_QG(1,outLoop)=cg_QG(1,outLoop)+                       &
     &                          zcglwk(iobs,1,outLoop)*                 &
     &                          zgrad0(iobs,outLoop)
            END DO

            preducv=1.0_r8
            preducy=1.0_r8
!
! Compute the actual penalty function terms.
!
            Jb=0.0_r8
            Jobs=0.0_r8
            DO iobs=1,Ndatum(ng)
              Jobs=Jobs+vgrad0s(iobs)*vgrad0s(iobs)
            END DO
            Jact=Jb+Jobs

          ELSE

            IF (Lcgini) THEN
!
!  When outer>1 we need to evaluate Ax0 so for inner=0 we use
!  cg_pxsave in v-space as the starting vector.
!
              DO iobs=1,Ndatum(ng)
                ADmodVal(iobs)=cg_pxsave(iobs)
# if defined W4DPSAS             || defined TL_W4DPSAS || \
     defined W4DPSAS_SENSITIVITY
                cg_innov(iobs)=-ObsScale(iobs)*                         &
     &                         (ObsVal(iobs)-BCKmodVal(iobs))
# else
                cg_innov(iobs)=-ObsScale(iobs)*                         &
     &                         (ObsVal(iobs)-TLmodVal(iobs))
# endif
              END DO
!
!  Convert ADmodVal from v-space to x-space and cg_innov (the
!  contribution to the starting gradient) from x-space to v-space.
!
              DO iobs=1,Ndatum(ng)
                IF (ObsErr(iobs).NE.0.0_r8) THEN
                  ADmodVal(iobs)=ADmodVal(iobs)/SQRT(ObsErr(iobs))
                  cg_innov(iobs)=cg_innov(iobs)/SQRT(ObsErr(iobs))
                END IF
              END DO

            ELSE

              DO iobs=1,Ndatum(ng)
!
!  Convert gradient contribution from x-space to v-space.
!
                pgrad(iobs)=ObsScale(iobs)*TLmodVal(iobs)
                gdgw(iobs)=pgrad(iobs)
                IF (ObsErr(iobs).NE.0.0_r8) THEN
                  pgrad(iobs)=pgrad(iobs)/SQRT(ObsErr(iobs))
                  vgrad0s(iobs)=pgrad(iobs)+cg_innov(iobs)+             &
     &                                      cg_pxsave(iobs)
                END IF
!
!  Add I*x0=cg_pxsave contribution to the gradient and the term
!  -b=cg_innov (everything is in v-space at this point).
!
                pgrad(iobs)=pgrad(iobs)+ObsScale(iobs)*                 &
     &                      (cg_pxsave(iobs)+cg_innov(iobs))
                vgrad0(iobs)=pgrad(iobs)
              END DO
!
!  If preconditioning, convert pgrad from v-space to y-space.
!
              IF (Lprecond.and.(outLoop.gt.1)) THEN
                Lscale=2                 ! SQRT spectral LMP
                Ltrans=.TRUE.
                CALL rprecond (ng, Lscale, Ltrans, outLoop, NinnLoop,   &
     &                         pgrad)
              END IF

              cg_Gnorm_y(outLoop)=0.0_r8
              cg_Gnorm_v(outLoop)=0.0_r8
              DO iobs=1,Ndatum(ng)
                zgrad0(iobs,outLoop)=pgrad(iobs)
                cg_Gnorm_y(outLoop)=cg_Gnorm_y(outLoop)+                &
     &                              zgrad0(iobs,outLoop)*               &
     &                              zgrad0(iobs,outLoop)
                cg_Gnorm_v(outLoop)=cg_Gnorm_v(outLoop)+                &
     &                              vgrad0(iobs)*vgrad0(iobs)
              END DO
              cg_Gnorm_y(outLoop)=SQRT(cg_Gnorm_y(outLoop))
              cg_Gnorm_v(outLoop)=SQRT(cg_Gnorm_v(outLoop))
!
              DO iobs=1,Ndatum(ng)
                zcglwk(iobs,1,outLoop)=pgrad(iobs)/cg_Gnorm_y(outLoop)
                ADmodVal(iobs)=zcglwk(iobs,1,outLoop)
              END DO
!
!  If preconditioning, convert ADmodVal from y-space to v-space.
!
              IF (Lprecond.and.(outLoop.gt.1)) THEN
                Lscale=2                 ! SQRT spectral LMP
                Ltrans=.FALSE.
                CALL rprecond (ng, Lscale, Ltrans, outLoop, NinnLoop,   &
     &                         ADmodVal)
              END IF
!
              DO iobs=1,Ndatum(ng)
                IF (ObsErr(iobs).NE.0.0_r8) THEN
                  ADmodVal(iobs)=ADmodVal(iobs)/SQRT(ObsErr(iobs))
                END IF
              END DO
              cg_QG(1,outLoop)=0.0_r8
              DO iobs=1,Ndatum(ng)
                cg_QG(1,outLoop)=cg_QG(1,outLoop)+                      &
     &                           zcglwk(iobs,1,outLoop)*                &
     &                           zgrad0(iobs,outLoop)
              END DO

            END IF

          END IF

          preducv=1.0_r8
          preducy=1.0_r8
!
!-----------------------------------------------------------------------
!  Minimization step, innLoop > 0.
!-----------------------------------------------------------------------
!
        ELSE
!
!  After the initialization, for every other inner loop, calculate a
!  new Lanczos vector, store it in the matrix, and update search.
!
          DO iobs=1,Ndatum(ng)
            pgrad(iobs)=ObsScale(iobs)*TLmodVal(iobs)
# if defined SENSITIVITY_4DVAR || \
     defined TL_W4DPSAS        || defined TL_W4DVAR || \
     defined W4DPSAS           || defined W4DVAR
            TLmodVal_S(iobs,innLoop,outLoop)=TLmodVal(iobs)
# endif
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
          IF (Lprecond.and.(outLoop.gt.1)) THEN
            Lscale=2                 ! SQRT spectral LMP
            Ltrans=.FALSE.
            CALL rprecond (ng, Lscale, Ltrans, outLoop, NinnLoop, zt)
          END IF
!
          DO iobs=1,Ndatum(ng)
            pgrad(iobs)=pgrad(iobs)+ObsScale(iobs)*zt(iobs)
          END DO
!
!  If preconditioning, convert pgrad from v-space to y-space.
!
          IF (Lprecond.and.(outLoop.gt.1)) THEN
            Lscale=2                 ! SQRT spectral LMP
            Ltrans=.TRUE.
            CALL rprecond (ng, Lscale, Ltrans, outLoop, NinnLoop, pgrad)
          END IF
!
          cg_delta(innLoop,outLoop)=0.0_r8
          DO iobs=1,Ndatum(ng)
            cg_delta(innLoop,outLoop)=cg_delta(innLoop,outLoop)+        &
     &                                zcglwk(iobs,innLoop,outLoop)*     &
     &                                pgrad(iobs)
          END DO
!
!  Exit, if not a positive definite algorithm.
!
          IF (cg_delta(innLoop,outLoop).le.0.0_r8) THEN
            WRITE (stdout,20) cg_delta(innLoop,outLoop), outLoop,       &
     &                        innLoop
            exit_flag=8
            GO TO 10
          END IF
!
!  Compute the new Lanczos vector.
!
          DO iobs=1,Ndatum(ng)
            pgrad(iobs)=pgrad(iobs)-cg_delta(innLoop,outLoop)*          &
     &                  zcglwk(iobs,innLoop,outLoop)
          END DO
          IF (innLoop.gt.1) THEN
            DO iobs=1,Ndatum(ng)
              pgrad(iobs)=pgrad(iobs)-cg_beta(innLoop,outLoop)*         &
     &                                 zcglwk(iobs,innLoop-1,outLoop)
            END DO
          END IF
!
!  Orthonormalize against previous Lanczos vectors.
!
          DO ivec=innLoop,1,-1
            cg_dla(ivec,outLoop)=0.0_r8
            DO iobs=1,Ndatum(ng)
              cg_dla(ivec,outLoop)=cg_dla(ivec,outLoop)+pgrad(iobs)*    &
     &                             zcglwk(iobs,ivec,outLoop)
            END DO
            DO iobs=1,Ndatum(ng)
              pgrad(iobs)=pgrad(iobs)-                                  &
     &                    cg_dla(ivec,outLoop)*                         &
     &                    zcglwk(iobs,ivec,outLoop)
            END DO
          END DO
!
          cg_beta(innLoop+1,outLoop)=0.0_r8
          DO iobs=1,Ndatum(ng)
            cg_beta(innLoop+1,outLoop)=cg_beta(innLoop+1,outLoop)+      &
     &                                 pgrad(iobs)*pgrad(iobs)
          END DO
          cg_beta(innLoop+1,outLoop)=SQRT(cg_beta(innLoop+1,outLoop))
!
          DO iobs=1,Ndatum(ng)
            zcglwk(iobs,innLoop+1,outLoop)=pgrad(iobs)/                 &
     &                                     cg_beta(innLoop+1,outLoop)
          END DO
!
          cg_QG(innLoop+1,outLoop)=0.0_r8
          DO iobs=1,Ndatum(ng)
            cg_QG(innLoop+1,outLoop)=cg_QG(innLoop+1,outLoop)+          &
     &                               zcglwk(iobs,innLoop+1,outLoop)*    &
     &                               zgrad0(iobs,outLoop)
          END DO
# ifdef MINRES
!
!  Use the minimum residual method as recommended by El Akkraoui and
!  Gauthier (2010, QJRMS, 136, 107-115) and described by Paige and
!  Saunders ("Sparse Indefinite Systems of Linear Equations", 1975,
!  SIAM Journal on Numerical Analysis, 617-619). Specifically we refer
!  to equations 6.10 and 6.11 of this paper.
!
!  Perform a LQ factorization of the tridiagonal matrix.
!
          ztriT=0.0_r8
          DO i=1,innLoop
            ztriT(i,i)=cg_delta(i,outLoop)
          END DO
          DO i=1,innLoop-1
            ztriT(i,i+1)=cg_beta(i+1,outLoop)
          END DO
          DO i=2,innLoop
            ztriT(i,i-1)=cg_beta(i,outLoop)
          END DO
          CALL sqlq(innLoop, ztriT, tau, zwork1)
!
!   Isolate L=zLT and its transpose.
!
          zLT=0.0_r8
          zLTt=0.0_r8
          DO i=1,innLoop
            DO j=1,i
              zLT(i,j)=ztriT(i,j)
            END DO
          END DO
          DO j=1,innLoop
            DO i=1,innLoop
              zLTt(i,j)=zLT(j,i)
            END DO
          END DO
!
!   Now form ze=beta1*Q*e1.
!
          ze=0.0_r8
          DO i=1,innLoop
            ze(i)=-cg_QG(i,outLoop)
          END DO
          DO i=1,innLoop
            DO j=1,innLoop
              zeref(j)=0.0_r8
            END DO
            zeref(i)=1.0_r8
            DO j=i+1,innLoop
              zeref(j)=ztriT(i,j)
            END DO
            zsum=0.0_r8
            DO j=1,innLoop
              zsum=zsum+ze(j)*zeref(j)
            END DO
            DO j=1,innLoop
              ze(j)=ze(j)-tau(i)*zsum*zeref(j)
            END DO
          END DO
!
!   Now form ze=D*ze (using equation 5.6 and 6.5 also).
!
          zgk=SQRT(zLT(innLoop,innLoop)*zLT(innLoop,innLoop)+           &
     &        cg_beta(innLoop+1,outLoop)*cg_beta(innLoop+1,outLoop))
          zck=zLT(innLoop,innLoop)/zgk
          ze(innLoop)=zck*ze(innLoop)
!
!   Now compute ze=inv(L')*ze.
!
          DO j=innLoop,1,-1
            ze(j)=ze(j)/zLTt(j,j)
            DO i=1,j-1
              ze(i)=ze(i)-ze(j)*zLTt(i,j)
            END DO
          END DO
!
!   Copy the solution ze into zwork.
!
          DO i=1,innLoop
            zwork(i,3)=ze(i)
          END DO
# else
!
!  Calculate the new solution based upon the new, orthonormalized
!  Lanczos vector. First, the tridiagonal system is solved by
!  decomposition and forward/back substitution.
!
          zbet=cg_delta(1,outLoop)
          zu(1)=-cg_QG(1,outLoop)/zbet
          DO ivec=2,innLoop
            zgam(ivec)=cg_beta(ivec,outLoop)/zbet
            zbet=cg_delta(ivec,outLoop)-cg_beta(ivec,outLoop)*zgam(ivec)
            zu(ivec)=(-cg_QG(ivec,outLoop)-cg_beta(ivec,outLoop)*       &
     &                zu(ivec-1))/zbet
          END DO
          zwork(innLoop,3)=zu(innLoop)

          DO ivec=innLoop-1,1,-1
            zu(ivec)=zu(ivec)-zgam(ivec+1)*zu(ivec+1)
            zwork(ivec,3)=zu(ivec)
          END DO
# endif
          DO iobs=1,Ndatum(ng)
            zw(iobs)=zgrad0(iobs,outLoop)+                              &
     &               cg_beta(innLoop+1,outLoop)*                        &
     &               zcglwk(iobs,innLoop+1,outLoop)*                    &
     &               zwork(innLoop,3)
          END DO

          DO iobs=1,Ndatum(ng)
            px(iobs)=0.0_r8
            DO ivec=1,innLoop
              px(iobs)=px(iobs)+zcglwk(iobs,ivec,outLoop)*zwork(ivec,3)
              zw(iobs)=zw(iobs)-zcglwk(iobs,ivec,outLoop)*              &
     &                                       cg_QG(ivec,outLoop)
            END DO
          END DO
!
!  If preconditioning, convert px from y-space to v-space.
!  We will always keep px in v-space.
!
          IF (Lprecond.and.(outLoop.gt.1)) THEN
            Lscale=2                 ! SQRT spectral LMP
            Ltrans=.FALSE.
            CALL rprecond (ng, Lscale, Ltrans, outLoop, NinnLoop, px)
          END IF
!
!  Compute the reduction in the gradient Ax-b (in y-space).
!
          preducy=0.0_r8
          DO iobs=1,Ndatum(ng)
            preducy=preducy+zw(iobs)*zw(iobs)
          END DO
          preducy=SQRT(preducy)/cg_Gnorm_y(outLoop)
!
!  Compute the reduction in the gradient Ax-b (in v-space).
!
          IF (Lprecond.and.(outLoop.gt.1)) THEN
            Lscale=-2                 ! SQRT spectral LMP
            Ltrans=.TRUE.
            CALL rprecond (ng, Lscale, Ltrans, outLoop, NinnLoop, zw)
          END IF
!
          preducv=0.0_r8
          DO iobs=1,Ndatum(ng)
            preducv=preducv+zw(iobs)*zw(iobs)
          END DO
          preducv=SQRT(preducv)/cg_Gnorm_v(outLoop)
!
!  Estimate the residual (in y-space) by:
!
!    ||Ax-b|| = cg_beta(k+1,outLoop)*zwork(innLoop,3)
!
          eps=ABS(cg_beta(innLoop+1,outLoop)*zwork(innLoop,3))
!
!  Estimate the various penalty function components based on
!  Bennett (2002), page 22 and pages 42-44:
!
!     Jopt  = value of the penalty function when true px
!             has been identified.
!     Jf    = initial data misfit for the first-guess
!     Jdata = data misfit for the state estimate
!     Jmod  = the model penalty function (the sum of the
!             dynamical, initial and boundary penalties)
!     Jb    = the ACTUAL model penalty function.
!     Jobs  = the ACTUAL data penalty function.
!     Jact  = the ACTUAL total penalty function.
!
!  NOTE:  Unlike IS4DVAR where we compute the actual cost function
!         for each iteration that decreases with increasing
!         number of iterations, the various penalty terms
!         represent the optimal values (i.e. the values when
!         optimal state estimate has been identified). Therefore
!         the various penalty terms will only be correct when the
!         optimal state has been found. Therefore, as the
!         W4DVAR proceeds, Jopt, Jdata and Jmod will asymptote
!         to their final values and WILL NOT necessary decrease
!         with increasing iteration number.
!
          IF (outLoop.eq.1.or.(.not.LhotStart)) THEN
            DO iobs=1,Ndatum(ng)
              Jopt=Jopt-px(iobs)*vgrad0(iobs)
              IF (ObsErr(iobs).ne.0.0_r8) THEN
                Jf=Jf+vgrad0(iobs)*vgrad0(iobs)
              END IF
              Jdata=Jdata+px(iobs)*px(iobs)
            END DO
            Jmod=Jopt-Jdata
          ELSE
            DO iobs=1,Ndatum(ng)
              Jopt=Jopt-(px(iobs)+cg_pxsave(iobs))*cg_innov(iobs)
              IF (ObsErr(iobs).ne.0.0_r8) THEN
                Jf=Jf+cg_innov(iobs)*cg_innov(iobs)
              END IF
              Jdata=Jdata+                                              &
     &              (px(iobs)+cg_pxsave(iobs))*                         &
     &              (px(iobs)+cg_pxsave(iobs))
            END DO
            Jmod=Jopt-Jdata
          END IF
!
!  Compute the actual penalty function terms.
!
          DO iobs=1,Ndatum(ng)
            zd(iobs)=vgrad0s(iobs)*SQRT(ObsErr(iobs))
          END DO
          DO iobs=1,Ndatum(ng)
            IF (ObsErr(iobs).ne.0.0_r8) THEN
              zhv(iobs)=zd(iobs)/SQRT(ObsErr(iobs))
            END IF
          END DO
          DO ivec=1,innLoop
            ztemp1(ivec)=0.0_r8
            DO iobs=1,Ndatum(ng)
              ztemp1(ivec)=ztemp1(ivec)+                                &
     &                     zcglwk(iobs,ivec,outLoop)*zhv(iobs)
            END DO
          END DO
# ifdef MINRES
!
!   Now form ze=Q*ztemp1.
!
          ze=0.0_r8
          DO i=1,innLoop
            ze(i)=ztemp1(i)
          END DO
          DO i=1,innLoop
            DO j=1,innLoop
              zeref(j)=0.0_r8
            END DO
            zeref(i)=1.0_r8
            DO j=i+1,innLoop
              zeref(j)=ztriT(i,j)
            END DO
            zsum=0.0_r8
            DO j=1,innLoop
              zsum=zsum+ze(j)*zeref(j)
            END DO
            DO j=1,innLoop
              ze(j)=ze(j)-tau(i)*zsum*zeref(j)
            END DO
          END DO
!
!   Now form ze=D*ze (using equation 5.6 and 6.5 also).
!
          zgk=SQRT(zLT(innLoop,innLoop)*zLT(innLoop,innLoop)+           &
     &        cg_beta(innLoop+1,outLoop)*cg_beta(innLoop+1,outLoop))
          zck=zLT(innLoop,innLoop)/zgk
          ze(innLoop)=zck*ze(innLoop)
!
!   Now compute ze=inv(L')*ze.
!
          DO j=innLoop,1,-1
            ze(j)=ze(j)/zLTt(j,j)
            DO i=1,j-1
              ze(i)=ze(i)-ze(j)*zLTt(i,j)
            END DO
          END DO
!
!   Copy the solution ze into ztemp2.
!
          DO i=1,innLoop
            ztemp2(i)=ze(i)
          END DO
!
# else
          zbet=cg_delta(1,outLoop)
          zu1(1)=ztemp1(1)/zbet
          DO ivec=2,innLoop
            zgam(ivec)=cg_beta(ivec,outLoop)/zbet
            zbet=cg_delta(ivec,outLoop)-                                &
     &           cg_beta(ivec,outLoop)*zgam(ivec)
            zu1(ivec)=(ztemp1(ivec)-                                    &
     &                 cg_beta(ivec,outLoop)*zu1(ivec-1))/zbet
          END DO
          ztemp2(innLoop)=zu1(innLoop)

          DO ivec=innLoop-1,1,-1
            zu1(ivec)=zu1(ivec)-zgam(ivec+1)*zu1(ivec+1)
            ztemp2(ivec)=zu1(ivec)
          END DO
# endif
          DO iobs=1,Ndatum(ng)
            zhv(iobs)=0.0_r8
          END DO
          DO iobs=1,Ndatum(ng)
            DO ivec=1,innLoop
              zhv(iobs)=zhv(iobs)+                                      &
     &                  ztemp2(ivec)*TLmodVal_S(iobs,ivec,outLoop)
            END DO
          END DO
!
!  Compute Jb.
!
          DO ivec=1,innLoop
            ztemp1(ivec)=0.0_r8
            DO iobs=1,Ndatum(ng)
              IF (ObsErr(iobs).NE.0.0_r8) THEN
                ztemp1(ivec)=ztemp1(ivec)+                              &
     &                       zcglwk(iobs,ivec,outLoop)*zhv(iobs)/       &
     &                       SQRT(ObsErr(iobs))
              END IF
            END DO
          END DO
# ifdef MINRES
!
!   Now form ze=Q*ztemp1.
!
          ze=0.0_r8
          DO i=1,innLoop
            ze(i)=ztemp1(i)
          END DO
          DO i=1,innLoop
            DO j=1,innLoop
              zeref(j)=0.0_r8
            END DO
            zeref(i)=1.0_r8
            DO j=i+1,innLoop
              zeref(j)=ztriT(i,j)
            END DO
            zsum=0.0_r8
            DO j=1,innLoop
              zsum=zsum+ze(j)*zeref(j)
            END DO
            DO j=1,innLoop
              ze(j)=ze(j)-tau(i)*zsum*zeref(j)
            END DO
          END DO
!
!   Now form ze=D*ze (using equation 5.6 and 6.5 also).
!
          zgk=SQRT(zLT(innLoop,innLoop)*zLT(innLoop,innLoop)+           &
     &        cg_beta(innLoop+1,outLoop)*cg_beta(innLoop+1,outLoop))
          zck=zLT(innLoop,innLoop)/zgk
          ze(innLoop)=zck*ze(innLoop)
!
!   Now compute ze=inv(L')*ze.
!
          DO j=innLoop,1,-1
            ze(j)=ze(j)/zLTt(j,j)
            DO i=1,j-1
              ze(i)=ze(i)-ze(j)*zLTt(i,j)
            END DO
          END DO
!
!   Copy the solution ze into ztemp2.
!
          DO i=1,innLoop
            ztemp2(i)=ze(i)
          END DO
!
# else
          zbet=cg_delta(1,outLoop)
          zu2(1)=ztemp1(1)/zbet
          DO ivec=2,innLoop
            zgam(ivec)=cg_beta(ivec,outLoop)/zbet
            zbet=cg_delta(ivec,outLoop)-                                &
     &           cg_beta(ivec,outLoop)*zgam(ivec)
            zu2(ivec)=(ztemp1(ivec)-                                    &
     &                 cg_beta(ivec,outLoop)*zu2(ivec-1))/zbet
          END DO
          ztemp2(innLoop)=zu2(innLoop)
          DO ivec=innLoop-1,1,-1
            zu2(ivec)=zu2(ivec)-zgam(ivec+1)*zu2(ivec+1)
            ztemp2(ivec)=zu2(ivec)
          END DO
# endif
          DO iobs=1,Ndatum(ng)
            zht(iobs)=0.0_r8
          END DO
          DO iobs=1,Ndatum(ng)
            DO ivec=1,innLoop
              IF (ObsErr(iobs).NE.0.0_r8) THEN
                zht(iobs)=zht(iobs)+                                    &
     &                    ztemp2(ivec)*zcglwk(iobs,ivec,outLoop)/       &
     &                    SQRT(ObsErr(iobs))
              END IF
            END DO
          END DO
          Jb=0.0_r8
          DO iobs=1,Ndatum(ng)
            IF (ObsErr(iobs).ne.0.0_r8) THEN
              cff=cg_pxsave(iobs)/SQRT(ObsErr(iobs))
              Jb=Jb+zd(iobs)*zht(iobs)-2.0_r8*cff*zhv(iobs)+            &
     &           cff*gdgw(iobs)
            END IF
          END DO
!
!  Compute Jobs.
!
          Jobs=0.0_r8
          IF (outLoop.eq.1.or.(.not.LhotStart)) THEN
            DO iobs=1,Ndatum(ng)
              IF (ObsErr(iobs).NE.0.0_r8) THEN
                cff=zhv(iobs)-zd(iobs)
                Jobs=Jobs+cff*cff/ObsErr(iobs)
              END IF
            END DO
          ELSE
            DO iobs=1,Ndatum(ng)
              IF (ObsErr(iobs).NE.0.0_r8) THEN
                cff=gdgw(iobs)-zhv(iobs)+cg_innov(iobs)*SQRT(ObsErr(iobs))
                Jobs=Jobs+cff*cff/ObsErr(iobs)
              END IF
            END DO
          END IF
          Jact=Jb+Jobs
!
!  Put the new trial solution into the adjoint vector for the next
!  loop or put the final solution into the adjoint vector on the
!  final inner-loop.
!
          IF (innLoop.eq.NinnLoop) THEN
!
!  Note: px is already in v-space so there is no need to convert
!        if preconditioning; cg_pxsave is also in v-space.
!
            DO iobs=1,Ndatum(ng)
              ADmodVal(iobs)=px(iobs)
            END DO
            IF (LhotStart) THEN
              DO iobs=1,Ndatum(ng)
                ADmodVal(iobs)=ADmodVal(iobs)+cg_pxsave(iobs)
              END DO
              DO iobs=1,Ndatum(ng)
                cg_pxsave(iobs)=ADmodVal(iobs)
              END DO
            END IF
          ELSE
            DO iobs=1,Ndatum(ng)
              ADmodVal(iobs)=zcglwk(iobs,innLoop+1,outLoop)
            END DO
!
!  If preconditioning, convert ADmodVal from y-space to v-space.
!
            IF (Lprecond.and.(outLoop.gt.1)) THEN
              Lscale=2                 ! SQRT spectral LMP
              Ltrans=.FALSE.
              CALL rprecond (ng, Lscale, Ltrans, outLoop, NinnLoop,     &
     &                       ADmodVal)
            END IF
          END IF
!
!  Convert ADmodVal from v-space to x-space.
!
          DO iobs=1,Ndatum(ng)
            IF (ObsErr(iobs).NE.0.0_r8) THEN
              ADmodVal(iobs)=ADmodVal(iobs)/SQRT(ObsErr(iobs))
            END IF
          END DO

        END IF MINIMIZE
!
!-----------------------------------------------------------------------
!  Report minimization progress.
!-----------------------------------------------------------------------
!
        WRITE (stdout,30) Ndatum(ng),                                   &
     &                    outLoop, innLoop, eps,                        &
     &                    outLoop, innLoop, preducy,                    &
     &                    outLoop, innLoop, preducv,                    &
     &                    outLoop, innLoop, Jf,                         &
     &                    outLoop, innLoop, Jdata,                      &
     &                    outLoop, innLoop, Jmod,                       &
     &                    outLoop, innLoop, Jopt,                       &
     &                    outLoop, innLoop, Jb,                         &
     &                    outLoop, innLoop, Jobs,                       &
     &                    outLoop, innLoop, Jact,                       &
     &                    outLoop, innLoop
        DO ivec=1,innLoop
          WRITE (stdout,40) ivec, cg_delta(ivec,outLoop),               &
     &                      cg_beta(ivec,outLoop), zwork(ivec,3)
        END DO
        WRITE (stdout,50) innLoop+1, cg_beta(innLoop+1,outLoop)
!
!  If preconditioning, report the Ritz eigenvalues used.
!
        IF ((LhessianEV.or.Lprecond).and.(outLoop.gt.1)) THEN
          IF (NritzEV.gt.0) THEN
            WRITE (stdout,60) outLoop, innLoop, nConvRitz
            ic=0
            DO ivec=1,NinnLoop
              IF (ivec.gt.(NinnLoop-nConvRitz)) THEN
                string='         '
                ic=ic+1
                WRITE (stdout,70) ivec, cg_Ritz(ivec,outLoop-1),        &
     &                            cg_RitzErr(ivec,outLoop-1),           &
     &                            TRIM(ADJUSTL(string)),                &
     &                            'Used=', ic
              ELSE
                string='             '
                WRITE (stdout,80) ivec, cg_Ritz(ivec,outLoop-1),        &
     &                            cg_RitzErr(ivec,outLoop-1),           &
     &                            TRIM(ADJUSTL(string))
              END IF
            END DO
          ELSE
            WRITE (stdout,90) outLoop, innLoop, RitzMaxErr
            ic=0
            DO ivec=1,NinnLoop
              IF (cg_RitzErr(ivec,outLoop-1).le.RitzMaxErr) THEN
                string='converged'
                ic=ic+1
                WRITE (stdout,70) ivec, cg_Ritz(ivec,outLoop-1),        &
     &                            cg_RitzErr(ivec,outLoop-1),           &
     &                            TRIM(ADJUSTL(string)),                &
     &                            'Good=', ic
              ELSE
                string='not converged'
                WRITE (stdout,80) ivec, cg_Ritz(ivec,outLoop-1),        &
     &                            cg_RitzErr(ivec,outLoop-1),           &
     &                            TRIM(ADJUSTL(string))
              END IF
            END DO
          END IF
        END IF
!
!-----------------------------------------------------------------------
!  If last inner loop, innLoop = NinnLoop.
!-----------------------------------------------------------------------
!
!  Compute the eigenvalues and eigenvectors of the tridiagonal matrix.
!
        IF (innLoop.eq.NinnLoop) THEN
          IF (LhessianEV.or.Lprecond) THEN
            DO ivec=1,innLoop
              cg_Ritz(ivec,outLoop)=cg_delta(ivec,outLoop)
            END DO
            DO ivec=1,innLoop-1
              zwork(ivec,1)=cg_beta(ivec+1,outLoop)
            END DO
!
!  Use the LAPACK routine DSTEQR to compute the eigenvectors and
!  eigenvalues of the tridiagonal matrix. If applicable, the
!  eigenpairs are computed by master thread only.
!
            CALL DSTEQR ('I', innLoop, cg_Ritz(1,outLoop), zwork(1,1),  &
     &                   zgv, Ninner, work, info)
            IF (info.ne.0) THEN
              PRINT *,' CONGRAD - Error in DSTEQR: info = ',info
              exit_flag=8
              GO TO 10
            END IF
!
            DO j=1,Ninner
              DO i=1,Ninner
                cg_zv(i,j,outLoop)=zgv(i,j)
              END DO
            END DO
!
!  Estimate the Ritz value error bounds.
!
            DO ivec=1,innLoop
              cg_RitzErr(ivec,outLoop)=ABS(cg_beta(innLoop+1,outLoop)*  &
     &                                 cg_zv(innLoop,ivec,outLoop))
            END DO
!
!  Check for exploding or negative Ritz values.
!
            DO ivec=1,innLoop
              IF (cg_RitzErr(ivec,outLoop).lt.0.0_r8) THEN
                PRINT *,'CONGRAD - negative Ritz value found.'
                exit_flag=8
                GO TO 10
              END IF
            END DO
!
!  Change the error bounds for the acceptable eigenvectors.
!
            RitzMaxErr=HevecErr
            DO ivec=1,NinnLoop
              cg_RitzErr(ivec,outLoop)=cg_RitzErr(ivec,outLoop)/        &
     &                                 cg_Ritz(NinnLoop,outLoop)
            END DO
!
!  Report Eigenvalues and their relative accuracy.
!
            WRITE (stdout,100) outLoop, InnLoop, RitzMaxErr
            ic=0
            DO ivec=1,NinnLoop
              IF (cg_RitzErr(ivec,outLoop).le.RitzMaxErr) THEN
                string='converged'
                ic=ic+1
                WRITE (stdout,70) ivec, cg_Ritz(ivec,outLoop),          &
     &                            cg_RitzErr(ivec,outLoop),             &
     &                            TRIM(ADJUSTL(string)),                &
     &                            'Good=', ic
              ELSE
                string='not converged'
                WRITE (stdout,80) ivec, cg_Ritz(ivec,outLoop),          &
     &                            cg_RitzErr(ivec,outLoop),             &
     &                            TRIM(ADJUSTL(string))
              END IF
            END DO
!
!  Calculate the converged eigenvectors (vcglev) of the [R_n + Cobs].
!
            CALL RPevecs (ng, outLoop, NinnLoop)
          END IF
        END IF

      END IF MASTER_THREAD
!
!-----------------------------------------------------------------------
!  Come here if not a possite definite algorithm.
!-----------------------------------------------------------------------
!
 10   CONTINUE
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, model, exit_flag)
# endif
      IF (exit_flag.ne.NoError) RETURN

# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Broadcast new solution to other nodes.
!-----------------------------------------------------------------------
!
      CALL mp_bcastf (ng, model, ADmodVal)
#  if defined SENSITIVITY_4DVAR || \
      defined TL_W4DPSAS        || defined TL_W4DVAR || \
      defined W4DPSAS           || defined W4DVAR
      CALL mp_bcastf (ng, model, TLmodVal_S(:,:,outLoop))
#  endif
      CALL mp_bcasti (ng, model, info)
      CALL mp_bcastf (ng, model, cg_beta(:,outLoop))
      CALL mp_bcastf (ng, model, cg_QG(:,outLoop))
      CALL mp_bcastf (ng, model, cg_delta(:,outLoop))
      CALL mp_bcastf (ng, model, cg_dla(:,outLoop))
      CALL mp_bcastf (ng, model, cg_Gnorm_y(:))
      CALL mp_bcastf (ng, model, cg_Gnorm_v(:))
      CALL mp_bcastf (ng, model, cg_pxsave(:))
      CALL mp_bcastf (ng, model, cg_innov(:))
      CALL mp_bcastf (ng, model, zgrad0(:,outLoop))
      CALL mp_bcastf (ng, model, zcglwk(:,:,outLoop))
      IF ((LhessianEV.or.Lprecond).and.(innLoop.eq.NinnLoop)) THEN
        CALL mp_bcastf (ng, model, cg_Ritz(:,outLoop))
        CALL mp_bcastf (ng, model, cg_RitzErr(:,outLoop))
        CALL mp_bcastf (ng, model, cg_zv(:,:,outLoop))
        CALL mp_bcastf (ng, model, vcglev(:,:,outLoop))
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Write out conjugate gradient vectors into 4DVAR NetCDF file.
!-----------------------------------------------------------------------
!
      CALL cg_write (ng, model, innLoop, outLoop,                       &
     &               Jf, Jdata, Jmod, Jopt, Jb, Jobs, Jact,             &
     &               preducv, preducy)
!
 20   FORMAT (/,' CONGRAD - Fatal error, not possitive definite',       &
     &        ' algorithm:',/,                                          &
     &        /,11x,'cg_delta = ',1p,e15.8,0p,3x,'(',i3.3,', ',i3.3,')')
 30   FORMAT (/,' CONGRAD - Conjugate Gradient Information:',/,         &
     &        /,11x,'Ndatum     = ',i7.7,/,/,                           &
     &        1x,'(',i3.3,',',i3.3,'): ',                               &
     &        'Residual estimate,                      eps = ',         &
     &        1p,e14.7,/,1x,'(',i3.3,',',i3.3,'): ',                    &
     &        'Reduction in gradient norm,  Greduc y-space = ',         &
     &        1p,e14.7,/,1x,'(',i3.3,',',i3.3,'): ',                    &
     &        'Reduction in gradient norm,  Greduc v-space = ',         &
     &        1p,e14.7,/,1x,'(',i3.3,',',i3.3,'): ',                    &
     &        'First guess initial data misfit,         Jf = ',         &
     &        1p,e14.7,/,1x,'(',i3.3,',',i3.3,'): ',                    &
     &        'State estimate data misfit,           Jdata = ',         &
     &        1p,e14.7,/,1x,'(',i3.3,',',i3.3,'): ',                    &
     &        'Model penalty function,                Jmod = ',         &
     &        1p,e14.7,/,1x,'(',i3.3,',',i3.3,'): ',                    &
     &        'Optimal penalty function,              Jopt = ',         &
     &        1p,e14.7,/,/,1x,'(',i3.3,',',i3.3,'): ',                  &
     &        'Actual Model penalty function,           Jb = ',         &
     &        1p,e14.7,/,1x,'(',i3.3,',',i3.3,'): ',                    &
     &        'Actual data penalty function,          Jobs = ',         &
     &        1p,e14.7,/,/,1x,'(',i3.3,',',i3.3,'): ',                  &
     &        'Actual total penalty function,         Jact = ',         &
     &        1p,e14.7,/,/,1x,'(',i3.3,',',i3.3,'): ',                  &
     &        'Lanczos vectors - cg_delta, cg_beta, zwork:',/)
 40   FORMAT (6x,i3.3,4x,3(1p,e15.8,5x))
 50   FORMAT (6x,i3.3,24x,1p,e15.8)
 60   FORMAT (/,1x,'(',i3.3,',',i3.3,'): ',                             &
     &        'Ritz eigenvalues used in preconditioning, ',             &
     &        'nConvRitz = ',i3.3,/)
 70   FORMAT (6x,i3.3,2x,1p,e14.7,2x,1p,e14.7,2x,a,5x,'('a,i3.3,')')
 80   FORMAT (6x,i3.3,2x,1p,e14.7,2x,1p,e14.7,2x,a)
 90   FORMAT (/,1x,'(',i3.3,',',i3.3,'): ',                             &
     &        'Ritz eigenvalues used in preconditioning, ',             &
     &        'RitzMaxErr = ',1p,e12.5,/)
100   FORMAT (/,1x,'(',i3.3,',',i3.3,'): ',                             &
     &        'New Ritz eigenvalues and their accuracy,  ',             &
     &        'RitzMaxErr = ',1p,e12.5,/)

      RETURN
      END SUBROUTINE congrad

      SUBROUTINE RPevecs (ng, outLoop, NinnLoop)
!
!=======================================================================
!                                                                      !
!  This routine computes the converged eigenvectors of (R_n+Cobs).     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
!
      implicit none
!
!  Imported variable declarations
!
      integer, intent(in) :: ng, outLoop, NinnLoop
!
!  Local variable declarations.
!
      integer :: iobs, ivec, jn, jk, ij, ingood

      real(r8) :: dla
!
!-----------------------------------------------------------------------
!  Compute and store converged eigenvectors of (R_n+Cobs).
!-----------------------------------------------------------------------
!
!  Count and collect the converged eigenvalues.
!
      ingood=0
      DO ivec=NinnLoop,1,-1
        ingood=ingood+1
        Ritz(ingood)=cg_Ritz(ivec,outLoop)
      END DO
!
!  Calculate the converged eigenvectors of (R_n+Cobs).
!
      ij=0
      DO jn=NinnLoop,1,-1
        ij=ij+1
        DO iobs=1,Ndatum(ng)
          vcglev(iobs,ij,outLoop)=0.0_r8
        END DO
        DO jk=1,NinnLoop
          DO iobs=1,Ndatum(ng)
            vcglev(iobs,ij,outLoop)=vcglev(iobs,ij,outLoop)+            &
     &                              zcglwk(iobs,jk,outLoop)*            &
     &                              cg_zv(jk,jn,outLoop)
          END DO
        END DO
        DO jk=1,ij-1
          dla=0.0_r8
          DO iobs=1,Ndatum(ng)
            dla=dla+vcglev(iobs,jk,outLoop)*vcglev(iobs,ij,outLoop)
          END DO
          DO iobs=1,Ndatum(ng)
            vcglev(iobs,ij,outLoop)=vcglev(iobs,ij,outLoop)-            &
     &                              dla*vcglev(iobs,jk,outLoop)
          END DO
        END DO
        dla=0.0_r8
        DO iobs=1,Ndatum(ng)
          dla=dla+vcglev(iobs,ij,outLoop)*vcglev(iobs,ij,outLoop)
        END DO
        DO iobs=1,Ndatum(ng)
          vcglev(iobs,ij,outLoop)=vcglev(iobs,ij,outLoop)/SQRT(dla)
        END DO
      END DO

      RETURN
      END SUBROUTINE RPevecs

      SUBROUTINE rprecond (ng, Lscale, Ltrans, outLoop, NinnLoop, py)
!
!=======================================================================
!                                                                      !
!  This routine is the preconditioner in observation space.            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Ltrans
      integer, intent(in) :: ng, Lscale, outLoop, NinnLoop
      real(r8), intent(inout) :: py(Ndatum(ng))
!
!  Local variable declarations.
!
      integer :: is, ie, inc, iss, i
      integer :: nol, nols, nole, ninc
      integer :: ingood
      integer :: iobs, nvec

      real(r8) :: dla, fac, facritz
!
!  Set the do-loop indices for the sequential preconditioner
!  loop.
!
      IF (Ltrans) THEN
        nols=outLoop-1
        nole=1
        ninc=-1
      ELSE
        nols=1
        nole=outLoop-1
        ninc=1
      END IF
!
!  Apply the preconditioners derived from all previous outer-loops
!  sequentially.
!
      DO nol=nols,nole,ninc
!
!  Determine the number of Ritz vectors to use.
!  For Lritz=.TRUE., choose HvecErr to be larger.
!
        ingood=0
        DO i=1,NinnLoop
          IF (cg_RitzErr(i,nol).le.RitzMaxErr) THEN
            ingood=ingood+1
          END IF
        END DO
        IF (NritzEV.gt.0) THEN
          nConvRitz=NritzEV
          ingood=nConvRitz
        ELSE
          nConvRitz=ingood
        END IF
!
        IF (Lscale.gt.0) THEN
          is=1
          ie=ingood
          inc=1
        ELSE
          is=ingood
          ie=1
          inc=-1
        END IF
!
        IF (Ltrans) THEN
          iss=is
          is=ie
          ie=iss
          inc=-inc
        END IF
!
!  Lscale determines the form of the preconditioner:
!
!     1= LMP (NOT CODED YET)
!    -1= Inverse LMP (NOT CODED YET)
!     2= Square root LMP
!    -2= Inverse square root LMP
!
        DO nvec=is,ie,inc
!
          fac=0.0_r8
!
          IF (Lritz) THEN
!
!  Note that the primitive Ritz vectors cg_zv are stored in order of
!  ascending Ritz values.
!
            facritz=cg_beta(NinnLoop+1,nol)*                            &
     &              cg_zv(NinnLoop,NinnLoop+1-nvec,nol)
!
!  Compute dot-product of py with q_k+1 from previous outer-loop.
!
            IF (.not.Ltrans) THEN
              dla=0.0_r8
              DO iobs=1,Ndatum(ng)
                dla=dla+zcglwk(iobs,NinnLoop+1,nol)*py(iobs)
              END DO
              facritz=facritz*dla
            END IF

          END IF
!
!  Compute the dot-product of py with the Ritz vectors from previous
!  outer-loop.
!
!  Note that vcglev are arranged in order of descending Ritz values,
!  while the Ritz values themselves that are stored in cg_Ritz are in
!  ascending order.
!
          dla=0.0_r8
          DO iobs=1,Ndatum(ng)
            dla=dla+vcglev(iobs,nvec,nol)*py(iobs)
          END DO
!
!  NOTE: Lscale=1 and Lscale=-1 cases are not complete yet.
!
          IF (Lscale.eq.-1) THEN
            fac=(cg_Ritz(NinnLoop+1-nvec,nol)-1.0_r8)*dla
          ELSE IF (Lscale.eq.1) THEN
            fac=(1.0_r8/cg_Ritz(NinnLoop+1-nvec,nol)-1.0_r8)*dla
          ELSE IF (Lscale.eq.-2) THEN
            fac=(SQRT(cg_Ritz(NinnLoop+1-nvec,nol))-1.0_r8)*dla
          ELSE IF (Lscale.eq.2) THEN
            fac=(1.0_r8/SQRT(cg_Ritz(NinnLoop+1-nvec,nol))-1.0_r8)*dla
          END IF
!
          IF (.not.Ltrans) THEN
            IF (Lritz.and.(Lscale.eq.-2)) THEN
              fac=fac+facritz/SQRT(cg_Ritz(NinnLoop+1-nvec,nol))
            END IF
            IF (Lritz.and.(Lscale.eq.2)) THEN
              fac=fac-facritz/cg_Ritz(NinnLoop+1-nvec,nol)
            END IF
          END IF
!
          DO iobs=1,Ndatum(ng)
             py(iobs)=py(iobs)+fac*vcglev(iobs,nvec,nol)
          END DO
!
          IF (Lritz.and.Ltrans) THEN

            IF (Lscale.eq.2) THEN
              fac=-facritz*dla/cg_Ritz(NinnLoop+1-nvec,nol)
            END IF
            IF (Lscale.eq.-2) THEN
              fac=facritz*dla/SQRT(cg_Ritz(NinnLoop+1-nvec,nol))
            END IF

            DO iobs=1,Ndatum(ng)
               py(iobs)=py(iobs)+fac*zcglwk(iobs,NinnLoop+1,nol)
            END DO

          END IF

        END DO
      END DO

      END SUBROUTINE rprecond

      SUBROUTINE cg_write (ng, model, innLoop, outLoop,                 &
     &                     Jf, Jdata, Jmod, Jopt, Jb, Jobs, Jact,       &
     &                     preducv, preducy)
!
!=======================================================================
!                                                                      !
!  This routine writes conjugate gradient vectors into 4DVAR NetCDF    !
!  for restart purposes.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      Use mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations
!
      integer, intent(in) :: ng, model, innLoop, outLoop

      real(r8), intent(in) :: Jf, Jdata, Jmod, Jopt, preducv, preducy
      real(r8), intent(in) :: Jb, Jobs, Jact
!
!  Local variable declarations.
!
      integer :: status
!
      SourceFile='congrad.F, cg_write'
!
!-----------------------------------------------------------------------
!  Write out conjugate gradient variables.
!-----------------------------------------------------------------------
!
!  Write out outer and inner iteration.
!
      CALL netcdf_put_ivar (ng, model, DAV(ng)%name, 'outer',           &
     &                      outer, (/0/), (/0/),                        &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, DAV(ng)%name, 'inner',           &
     &                      inner, (/0/), (/0/),                        &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out number of converged Ritz eigenvalues.
!
      IF (innLoop.eq.Ninner) THEN
        CALL netcdf_put_ivar (ng, model, DAV(ng)%name, 'nConvRitz',     &
     &                        nConvRitz, (/outLoop/), (/1/),            &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out converged Ritz eigenvalues.
!
      IF (innLoop.eq.Ninner) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'Ritz',          &
     &                        Ritz, (/1,outLoop/), (/nConvRitz,1/),     &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out conjugate gradient norm.
!
      IF (innLoop.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_beta',       &
     &                        cg_beta, (/1,1/), (/Ninner+1,Nouter/),    &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out Lanczos algorithm coefficients.
!
      IF (innLoop.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_delta',      &
     &                        cg_delta, (/1,1/), (/Ninner,Nouter/),     &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write normalization coefficients for Lanczos vectos.
!
      IF (innLoop.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_dla',        &
     &                        cg_dla, (/1,1/), (/Ninner,Nouter/),       &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Initial gradient normalization factor.
!
      IF (innLoop.eq.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_Gnorm_y',    &
     &                        cg_Gnorm_y, (/1/), (/Nouter/),            &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN

        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_Gnorm_v',    &
     &                        cg_Gnorm_v, (/1/), (/Nouter/),            &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Lanczos vector normalization factor.
!
      IF (innLoop.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_QG',         &
     &                        cg_QG, (/1,1/), (/Ninner+1,Nouter/),      &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Reduction in the gradient norm.
!
      IF (innLoop.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_Greduc_y',   &
     &                        preducy, (/innLoop,outLoop/), (/1,1/),    &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN

        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_Greduc_v',   &
     &                        preducv, (/innLoop,outLoop/), (/1,1/),    &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Eigenvalues of Lanczos recurrence relationship.
!
      IF (innLoop.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_Ritz',       &
     &                        cg_Ritz, (/1,1/), (/Ninner,Nouter/),      &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Eigenvalues relative error.
!
      IF (innLoop.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_RitzErr',    &
     &                        cg_RitzErr, (/1,1/), (/Ninner,Nouter/),   &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Eigenvectors of Lanczos recurrence relationship.
!
      IF (innLoop.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_zv',         &
     &                        cg_zv(:,:,outLoop),                       &
     &                        (/1,1,outLoop/), (/Ninner,Ninner,1/),     &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  First guess initial data misfit.
!
      CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'Jf',              &
     &                      Jf, (/innLoop+1,outLoop/), (/1,1/),         &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  State estimate data misfit.
!
      CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'Jdata',           &
     &                      Jdata, (/innLoop+1,outLoop/), (/1,1/),      &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Model penalty function.
!
      CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'Jmod',            &
     &                      Jmod, (/innLoop+1,outLoop/), (/1,1/),       &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Optimal penalty function.
!
      CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'Jopt',            &
     &                      Jopt, (/innLoop+1,outLoop/), (/1,1/),       &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Actual model penalty function.
!
      CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'Jb',              &
     &                      Jb, (/innLoop+1,outLoop/), (/1,1/),         &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Actual data penalty function.
!
      CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'Jobs',            &
     &                      Jobs, (/innLoop+1,outLoop/), (/1,1/),       &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Actual total penalty function.
!
      CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'Jact',            &
     &                      Jact, (/innLoop+1,outLoop/), (/1,1/),       &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Initial gradient for minimization.
!
      IF (innLoop.eq.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'zgrad0',        &
     &                        zgrad0(:,outLoop),                        &
     &                        (/1,outLoop/), (/Ndatum(ng),1/),          &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Lanczos vectors.
!
      CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'zcglwk',          &
     &                      zcglwk(:,innLoop+1,outLoop),                &
     &                      (/1,innLoop+1,outLoop/), (/Ndatum(ng),1,1/),&
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Saved TLmodVal_S.
!
      IF (innLoop.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'TLmodVal_S',    &
     &                        TLmodVal_S(:,innLoop,outLoop),            &
     &                        (/1,innLoop,outLoop/), (/Ndatum(ng),1,1/),&
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize model/observation NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, model, DAV(ng)%name, DAV(ng)%ncid)

      RETURN
      END SUBROUTINE cg_write
#else
      SUBROUTINE congrad
      RETURN
      END SUBROUTINE congrad
#endif
