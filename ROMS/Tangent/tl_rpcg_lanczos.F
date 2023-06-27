#include "cppdefs.h"

#if (defined SENSITIVITY_4DVAR || \
     defined TL_RBL4DVAR       || \
     defined TL_R4DVAR)        && defined RPCG

      SUBROUTINE tl_rpcg_lanczos (ng, model, outLoop, innLoop,          &
     &                            NinnLoop, Lcgini)
!
!git $Id$
!svn $Id: tl_rpcg_lanczos.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group       Selime Gurol      !
!    Licensed under a MIT/X style license            Andrew M. Moore   !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Weak Constraint 4-Dimensional Variational (4DVar) Restricted        !
!      Pre-conditioned Conjugate Gradient (RPCG) Algorithm             !

# if defined R4DVAR_ANA_SENSITIVITY || defined TL_R4DVAR
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
!  gradient is Ax-b and the stabilized representer matrix is A.        !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Chua, B. S. and A. F. Bennett,  2001:  An inverse ocean modeling  !
!      sytem, Ocean Modelling, 3, 137-165.                             !

# elif defined RBL4DVAR_ANA_SENSITIVITY || \
       defined RBL4DVAR_FCT_SENSITIVITY || \
       defined TL_RBL4DVAR
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
!  Reference:                                                          !
!                                                                      !
!    Courtier, P., 1997: Dual formulation of four-dimensional          !
!      variational assimilation, Quart. J. Roy. Meteor. Soc.,          !
!      123, 2449-2461.                                                 !
# endif
!                                                                      !
!  The restricted preconditioned conjugate gradient (RPCG) algorithm   !
!  is used to compute the solution. Specfically, the following system  !
!  is solved:                                                          !
!                                                                      !
!              (Cobs^-1 H M_n B (M_n)' H'+I)*w=Cobs^-1 d_n             !
!                                                                      !
!  All inner-products are defined with respect to a norm involving     !
!  H M_n B (M_n)' H', and the convergence of J in model space is       !
!  guaranteed to converge at the same rate as I4D-Var in primal space. !
!                                                                      !
!  RPCG Reference:                                                     !
!                                                                      !
!    Gratton, S. and J. Tshimanga, 2009: An observation-space          !
!      formulation of variational assimilation using a restricted      !
!      preconditioned conjugate gradient algorithm. QJRMS, 135,        !
!      1573-1585.
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
      logical :: Ltrans, Laug
      integer :: i, ic, j, iobs, ivec, Lscale, info

      real(r8) :: zbet, eps, preducv, preducy
      real(r8) :: Jopt, Jf, Jmod, Jdata, Jb, Jobs, Jact, cff
      real(r8) :: tl_dla, fact

      real(r8), dimension(NinnLoop) :: zu, zgam, zfact, tl_zfact
      real(r8), dimension(NinnLoop) :: tl_zu, tl_zrhs
      real(r8), dimension(Ndatum(ng)+1) :: px, pgrad
      real(r8), dimension(Ndatum(ng)+1) :: tl_px, tl_pgrad
      real(r8), dimension(Ninner,3) :: zwork, tl_zwork

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
!  Set vgrad0 based on the current outer-loop.
!
        DO iobs=1,Ndatum(ng)+1
          vgrad0(iobs)=zgrad0(iobs,outLoop)
        END DO
!
!  Multiply TLmodVal by ObsScale to effectively remove any rejected
!  observations from the analysis.
!
        IF (innLoop.gt.0) THEN
          DO iobs=1,Ndatum(ng)
!^          TLmodVal(iobs)=ObsScale(iobs)*TLmodVal(iobs)
!^
            TLmodVal(iobs)=ObsScale(iobs)*                              &
     &                     TLmodVal(iobs)+                              &
     &                     tl_ObsScale(iobs)*                           &
     &                     TLmodVal_S(iobs,innLoop,outLoop)
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  Initialization step, innloop = 0.
!-----------------------------------------------------------------------
!
        MINIMIZE : IF (innLoop.eq.0) THEN

# if defined RBL4DVAR_ANA_SENSITIVITY || \
     defined RBL4DVAR_FCT_SENSITIVITY || \
     defined TL_RBL4DVAR              || \
     defined TL_R4DVAR
!
          DO iobs=1,Ndatum(ng)
            BCKmodVal(iobs)=NLmodVal(iobs)
            tl_BCKmodVal(iobs)=0.0_r8
!!
!! THIS IS A TEMPORARY FIX.
!!
!!AMM       tl_ObsVal(iobs)=ObsVal(iobs)
!!AMM       tl_ObsErr(iobs)=ObsErr(iobs)
!!
          END DO
# endif
!
!  For the first outer-loop, x0=0.
!
          IF (outLoop.eq.1) THEN
            DO iobs=1,Ndatum(ng)+1
!^            px(iobs)=0.0_r8
!^
              tl_px(iobs)=0.0_r8
            END DO
            DO iobs=1,Ndatum(ng)
!^            Hbk(iobs,outLoop)=0.0_r8
!^
              tl_Hbk(iobs)=0.0_r8
            END DO
          ELSE
!!
!! NOTE: CODE WON'T WORK FOR R4D-Var AT THE MOMENT SINCE THIS IS
!! THE WRONG TLmodVal HERE!!!
!!
            DO iobs=1,Ndatum(ng)
!^            Hbk(iobs,outLoop)=-TLmodVal(iobs)
!^            tl_Hbk(iobs)=-tl_TLmodVal(iobs)
!^
              tl_Hbk(iobs)=-TLmodVal(iobs)
            END DO
!^          IF (Lprecond) THEN
!^            DO iobs=1,Ndatum(ng)+1
!^              FOURDVAR(ng)%cg_Dpxsum(iobs)=FOURDVAR(ng)%cg_pxsum(iobs)
!^            END DO
!^            Lscale=1                 ! Spectral LMP
!^            Laug=.FALSE.
!^            CALL rprecond (ng, Lscale, Laug, outLoop, NinnLoop-1,     &
!^   &                       FOURDVAR(ng)%cg_Dpxsum)
!^            DO iobs=1,Ndatum(ng)+1
!^              FOURDVAR(ng)%cg_Dpxsum(iobs)=                           &
!^   &                                 -FOURDVAR(ng)%cg_Dpxsum(iobs)+
!^   &                                  FOURDVAR(ng)%cg_pxsum(iobs)
!^            END DO
!^          END IF
          END IF
!
!  If this is the first pass of the inner loop, set up the vectors and
!  store the first guess. The initial starting guess is assumed to be
!  zero in which case the residual is just: (obs-model).
!
          DO iobs=1,Ndatum(ng)
!
!  Selime code:  r=b.
!
# if defined RBL4DVAR                 || \
     defined RBL4DVAR_ANA_SENSITIVITY || \
     defined RBL4DVAR_FCT_SENSITIVITY || \
     defined TL_RBL4DVAR
!^          pgrad(iobs)=ObsScale(iobs)*                                 &
!^   &                  (ObsVal(iobs)-BCKmodVal(iobs))
!^
            tl_pgrad(iobs)=ObsScale(iobs)*                              &
     &                     (tl_ObsVal(iobs)-                            &
     &                      tl_BCKmodVal(iobs))+                        &
     &                      tl_ObsScale(iobs)*                          &
     &                     (ObsVal(iobs)-BCKmodVal(iobs))
# else
!!
!! NOTE: THIS CODE IS NOT CORRECT YET FOR R4D-Var
!!
!^          pgrad(iobs)=ObsScale(iobs)*                                 &
!^   &                  (ObsVal(iobs)-TLmodVal(iobs))
!^          pgrad(iobs)=ObsScale(iobs)*                                 &
!^   &                  (ObsVal(iobs)-TLmodVal_S(iobs,innLoop,outLoop))
!^          tl_pgrad(iobs)=ObsScale(iobs)*                              &
!^   &                     (tl_ObsVal(iobs)-tl_TLmodVal(iobs))
!^
            tl_pgrad(iobs)=ObsScale(iobs)*                              &
     &                     (tl_ObsVal(iobs)-                            &
     &                      TLmodVal(iobs))+                            &
     &                     tl_ObsScale(iobs)*                           &
     &                     (ObsVal(iobs)-                               &
     &                      TLmodVal_S(iobs,innLoop,outLoop))
# endif
            IF (ObsErr(iobs).NE.0.0_r8) THEN
!^             pgrad(iobs)=pgrad(iobs)/ObsErr(iobs)
!^
               pgrad(iobs)=vgrad0(iobs)
               tl_pgrad(iobs)=tl_pgrad(iobs)/ObsErr(iobs)-              &
     &                        tl_ObsErr(iobs)*pgrad(iobs)/ObsErr(iobs)
            END IF
!^          vgrad0(iobs)=pgrad(iobs)
!^
            tl_vgrad0(iobs)=tl_pgrad(iobs)
!^          vgrad0s(iobs)=vgrad0(iobs)
!^
            tl_vgrad0s(iobs)=tl_vgrad0(iobs)
          END DO
!
!  AUGMENTED.
!
!  Augment the residual vector.
!
          pgrad(Ndatum(ng)+1)=1.0_r8
          tl_pgrad(Ndatum(ng)+1)=0.0_r8
!^        vgrad0(Ndatum(ng)+1)=pgrad(Ndatum(ng)+1)
!^
          tl_vgrad0(Ndatum(ng)+1)=tl_pgrad(Ndatum(ng)+1)
!^        vgrad0s(Ndatum(ng)+1)=vgrad0(Ndatum(ng)+1)
!^
          tl_vgrad0s(Ndatum(ng)+1)=tl_vgrad0(Ndatum(ng)+1)
!
!  If preconditioning, convert pgrad from v-space to y-space.
!
!  Selime code: z1=G*r.
!
!^        IF (Lprecond.and.(outLoop.gt.1)) THEN
!^          Lscale=1                 ! Spectral LMP
!^          Laug=.TRUE.
!^          CALL rprecond (ng, Lscale, Laug, outLoop, NinnLoop-1,       &
!^   &                     pgrad)
!^        END IF
!^
!
!  Selime code: zgrad0=z1.
!
!  Selime code: vgrad0=r.
!
          DO iobs=1,Ndatum(ng)+1
!^          zgrad0(iobs,outLoop)=pgrad(iobs)
!^
            tl_zgrad0(iobs)=tl_pgrad(iobs)
          END DO
!
!  Selime code: ADmodVal=z1.
!
          DO iobs=1,Ndatum(ng)
!^          ADmodVal(iobs)=zgrad0(iobs,outLoop)
!^          tl_ADmodVal(iobs)=tl_zgrad0(iobs,outLoop)
!^
            ADmodVal(iobs)=tl_zgrad0(iobs)
          END DO
!
          preducv=1.0_r8
          preducy=1.0_r8
!
!-----------------------------------------------------------------------
!  Minimization step, innLoop > 0.
!-----------------------------------------------------------------------
!
        ELSE
!
!  COMPLETE THE CALCULATION OF PREVIOUS LANCZOS VECTOR.
!
!  Selime code: t1=GDG'*z1=TLmodVal.
!
!
!  Complete the computation of the Lanczos vector from previous
!  inner-loop.
!
          IF (innLoop.EQ.1) THEN
!^          cg_Gnorm_v(outLoop)=0.0_r8
!^
            tl_cg_Gnorm_v(outLoop)=0.0_r8
            DO iobs=1,Ndatum(ng)
!^            cg_Gnorm_v(outLoop)=cg_Gnorm_v(outLoop)+                  &
!^   &                            TLmodVal(iobs)*vgrad0(iobs)
!^
              tl_cg_Gnorm_v(outLoop)=tl_cg_Gnorm_v(outLoop)+            &
     &                               TLmodVal(iobs)*                    &
     &                               vgrad0(iobs)+                      &
     &                               TLmodVal_S(iobs,innLoop,outLoop)*  &
     &                               tl_vgrad0(iobs)
            END DO
!
!  AUGMENTED. Here Hbk=H(x0-xb). Hbk=0 and Jb0=0 when outer=1.
!  Jb0 is computed using sum_i (G'*w)_i, i.e. the sum of the outputs
!  of the adjoint model at the end of AD_LOOP2.
!
            DO iobs=1,Ndatum(ng)
!^            cg_Gnorm_v(outLoop)=cg_Gnorm_v(outLoop)+                  &
!^   &                            Hbk(iobs,outLoop)*                    &
!^   &                            zgrad0(Ndatum(ng)+1,outLoop)*         &
!^   &                            vgrad0(iobs)+                         &
!^   &                            Hbk(iobs,outLoop)*                    &
!^   &                            vgrad0(Ndatum(ng)+1)*                 &
!^   &                            zgrad0(iobs,outLoop)
!^
              tl_cg_Gnorm_v(outLoop)=tl_cg_Gnorm_v(outLoop)+            &
     &                               tl_Hbk(iobs)*                      &
     &                               zgrad0(Ndatum(ng)+1,outLoop)*      &
     &                               vgrad0(iobs)+                      &
     &                               Hbk(iobs,outLoop)*                 &
     &                               tl_zgrad0(Ndatum(ng)+1)*           &
     &                               vgrad0(iobs)+                      &
     &                               Hbk(iobs,outLoop)*                 &
     &                               zgrad0(Ndatum(ng)+1,outLoop)*      &
     &                               tl_vgrad0(iobs)+                   &
     &                               tl_Hbk(iobs)*                      &
     &                               vgrad0(Ndatum(ng)+1)*              &
     &                               zgrad0(iobs,outLoop)+              &
     &                               Hbk(iobs,outLoop)*                 &
     &                               tl_vgrad0(Ndatum(ng)+1)*           &
     &                               zgrad0(iobs,outLoop)+              &
     &                               Hbk(iobs,outLoop)*                 &
     &                               vgrad0(Ndatum(ng)+1)*              &
     &                               tl_zgrad0(iobs)
            END DO
!^          cg_Gnorm_v(outLoop)=cg_Gnorm_v(outLoop)+                    &
!^   &                          Jb0(outLoop-1)*                         &
!^   &                          vgrad0(Ndatum(ng)+1)*                   &
!^   &                          zgrad0(Ndatum(ng)+1,outLoop)
!^
            tl_cg_Gnorm_v(outLoop)=tl_cg_Gnorm_v(outLoop)+              &
     &                             tl_Jb0(outLoop-1)*                   &
     &                             vgrad0(Ndatum(ng)+1)*                &
     &                             zgrad0(Ndatum(ng)+1,outLoop)+        &
     &                             Jb0(outLoop-1)*                      &
     &                             tl_vgrad0(Ndatum(ng)+1)*             &
     &                             zgrad0(Ndatum(ng)+1,outLoop)+        &
     &                             Jb0(outLoop-1)*                      &
     &                             vgrad0(Ndatum(ng)+1)*                &
     &                             tl_zgrad0(Ndatum(ng)+1)
!
!^          cg_Gnorm_v(outLoop)=SQRT(cg_Gnorm_v(outLoop))
!^
            tl_cg_Gnorm_v(outLoop)=0.5_r8*tl_cg_Gnorm_v(outLoop)/       &
     &                             cg_Gnorm_v(outLoop)
!
!  Selime code: beta_first=sqrt(r't1)=cg_Gnorm_v.
!                 cg_QG(1)=beta_first
!
!^          cg_QG(1,outLoop)=cg_Gnorm_v(outLoop)
!^
            tl_cg_QG(1)=tl_cg_Gnorm_v(outLoop)
!
!  Normalize TLmodVal here and save in TLmodVal_S.
!
!  AUGMENTED.
!
            DO iobs=1,Ndatum(ng)+1
!
!  Selime code: v1=r/beta_first=zcglwk(1).
!
!^            zcglwk(iobs,1,outLoop)=vgrad0(iobs)/                      &
!^   &                               cg_Gnorm_v(outLoop)
!^
              tl_zcglwk(iobs,1)=tl_vgrad0(iobs)/                        &
     &                          cg_Gnorm_v(outLoop)-                    &
     &                          tl_cg_Gnorm_v(outLoop)*                 &
     &                          zcglwk(iobs,1,outLoop)/                 &
     &                          cg_Gnorm_v(outLoop)
!
!  Selime code: z1=z1/beta_first=vcglwk(1).
!
!^            vcglwk(iobs,1,outLoop)=zgrad0(iobs,outLoop)/              &
!^   &                               cg_Gnorm_v(outLoop)
!^
              tl_vcglwk(iobs,1)=tl_zgrad0(iobs)/                        &
     &                          cg_Gnorm_v(outLoop)-                    &
     &                          tl_cg_Gnorm_v(outLoop)*                 &
     &                          vcglwk(iobs,1,outLoop)/                 &
     &                          cg_Gnorm_v(outLoop)
            END DO
!
!  AUGMENTED.
!
            DO iobs=1,Ndatum(ng)
!
!  Selime code: t1=t1/beta_first=TLmodVal
!
!^            TLmodVal_S(iobs,innLoop,outLoop)=TLmodVal(iobs)
!^
              tl_TLmodVal_S(iobs,innLoop,outLoop)=TLmodVal(iobs)
!^            TLmodVal(iobs)=TLmodVal(iobs)/                            &
!^   &                       cg_Gnorm_v(outLoop)
!^
              TLmodVal(iobs)=TLmodVal(iobs)/                            &
     &                       cg_Gnorm_v(outLoop)-                       &
     &                       tl_cg_Gnorm_v(outLoop)*                    &
     &                       TLmodVal_S(iobs,innLoop,outLoop)/          &
     &                       (cg_Gnorm_v(outLoop)*                      &
     &                        cg_Gnorm_v(outLoop))
            END DO
!
          ELSE
!
!  First grab the non-normalized Lanczos vector.
!
!  AUGMENTED.
!
            DO iobs=1,Ndatum(ng)+1
!
!  Selime code: w=pgrad, t1=TLmodVal.
!
!
!  Need to multiply zcglwk by cg_beta when recomputing pgrad
!  since zcglwk is the normalized Lanczos vector.
!
              pgrad(iobs)=zcglwk(iobs,innLoop,outLoop)*                 &
     &                    cg_beta(innLoop,outLoop)
              tl_pgrad(iobs)=tl_zcglwk(iobs,innLoop)
            END DO
!
!  Selime's code:  beta=pgrad*HBH'*pgrad=cg_beta
!
!^          cg_beta(innLoop,outLoop)=0.0_r8
!^
            tl_cg_beta(innLoop)=0.0_r8
            DO iobs=1,Ndatum(ng)
!^            cg_beta(innLoop,outLoop)=cg_beta(innLoop,outLoop)+        &
!^   &                                 pgrad(iobs)*                     &
!^   &                                 TLmodVal(iobs)
!^
              tl_cg_beta(innLoop)=tl_cg_beta(innLoop)+                  &
     &                            tl_pgrad(iobs)*                       &
     &                            TLmodVal_S(iobs,innLoop,outLoop)+     &
     &                            pgrad(iobs)*                          &
     &                            TLmodVal(iobs)
            END DO
!
!  AUGMENTED.
!
            DO iobs=1,Ndatum(ng)
!^            cg_beta(innLoop,outLoop)=cg_beta(innLoop,outLoop)+        &
!^   &                                 Hbk(iobs,outLoop)*               &
!^   &                                 vcglwk(Ndatum(ng)+1,innLoop,     &
!^   &                                        outLoop)*                 &
!^   &                                 pgrad(iobs)+                     &
!^   &                                 Hbk(iobs,outLoop)*               &
!^   &                                 pgrad(Ndatum(ng)+1)*             &
!^   &                                 vcglwk(iobs,innLoop,outLoop)
!^
              tl_cg_beta(innLoop)=tl_cg_beta(innLoop)+                  &
     &                            tl_Hbk(iobs)*                         &
     &                            vcglwk(Ndatum(ng)+1,innLoop,outLoop)* &
     &                            pgrad(iobs)+                          &
     &                            Hbk(iobs,outLoop)*                    &
     &                            tl_vcglwk(Ndatum(ng)+1,innLoop)*      &
     &                            pgrad(iobs)+                          &
     &                            Hbk(iobs,outLoop)*                    &
     &                            vcglwk(Ndatum(ng)+1,innLoop,outLoop)* &
     &                            tl_pgrad(iobs)+                       &
     &                            tl_Hbk(iobs)*                         &
     &                            pgrad(Ndatum(ng)+1)*                  &
     &                            vcglwk(iobs,innLoop,outLoop)+         &
     &                            Hbk(iobs,outLoop)*                    &
     &                            tl_pgrad(Ndatum(ng)+1)*               &
     &                            vcglwk(iobs,innLoop,outLoop)+         &
     &                            Hbk(iobs,outLoop)*                    &
     &                            pgrad(Ndatum(ng)+1)*                  &
     &                            tl_vcglwk(iobs,innLoop)
            END DO
!^          cg_beta(innLoop,outLoop)=cg_beta(innLoop,outLoop)+          &
!^   &                               Jb0(outLoop-1)*                    &
!^   &                               pgrad(Ndatum(ng)+1)*               &
!^   &                               vcglwk(Ndatum(ng)+1,innLoop,       &
!^   &                                      outLoop)
!^
            tl_cg_beta(innLoop)=tl_cg_beta(innLoop)+                    &
     &                          tl_Jb0(outLoop-1)*                      &
     &                          pgrad(Ndatum(ng)+1)*                    &
     &                          vcglwk(Ndatum(ng)+1,innLoop,outLoop)+   &
     &                          Jb0(outLoop-1)*                         &
     &                          tl_pgrad(Ndatum(ng)+1)*                 &
     &                          vcglwk(Ndatum(ng)+1,innLoop,outLoop)+   &
     &                          Jb0(outLoop-1)*                         &
     &                          pgrad(Ndatum(ng)+1)*                    &
     &                          tl_vcglwk(Ndatum(ng)+1,innLoop)
!
!^          cg_beta(innLoop,outLoop)=SQRT(cg_beta(innLoop,outLoop))
!^
            tl_cg_beta(innLoop)=0.5_r8*tl_cg_beta(innLoop)/             &
     &                           cg_beta(innLoop,outLoop)
!
!  Normalize TLmodVal here and save in TLmodVal_S.
!
!  AUGMENTED.
!
            DO iobs=1,Ndatum(ng)+1
!
!  Selime code: v1=w/beta=zcglwk
!
!^            zcglwk(iobs,innLoop,outLoop)=pgrad(iobs)/                 &
!^   &                                     cg_beta(innLoop,outLoop)
!^
              tl_zcglwk(iobs,innLoop)=tl_pgrad(iobs)/                   &
     &                                       cg_beta(innLoop,outLoop)-  &
     &                                       tl_cg_beta(innLoop)*       &
     &                                zcglwk(iobs,innLoop,outLoop)/     &
     &                                cg_beta(innLoop,outLoop)
!
!   Selime code: z1=w/beta=vcglwk
!
!^            vcglwk(iobs,innLoop,outLoop)=vcglwk(iobs,innLoop,         &
!^   &                                            outLoop)/             &
!^   &                                     cg_beta(innLoop,outLoop)
!^
              tl_vcglwk(iobs,innLoop)=tl_vcglwk(iobs,innLoop)/          &
     &                                cg_beta(innLoop,outLoop)-         &
     &                                tl_cg_beta(innLoop)*              &
     &                                vcglwk(iobs,innLoop,outLoop)/     &
     &                                cg_beta(innLoop,outLoop)
            END DO
!
!  Selime code: t1=t1/beta=TLmodVal
!
!  AUGMENTED.
!
            DO iobs=1,Ndatum(ng)
!^            TLmodVal_S(iobs,innLoop,outLoop)=TLmodVal(iobs)
!^
              tl_TLmodVal_S(iobs,innLoop,outLoop)=TLmodVal(iobs)
!^            TLmodVal(iobs)=TLmodVal(iobs)/                            &
!^   &                       cg_beta(innLoop,outLoop)
              TLmodVal(iobs)=TLmodVal(iobs)/cg_beta(innLoop,outLoop)-   &
     &                       tl_cg_beta(innLoop)*                       &
     &                       TLmodVal_S(iobs,innLoop,outLoop)/          &
     &                       (cg_beta(innLoop,outLoop)*                 &
     &                        cg_beta(innLoop,outLoop))
            END DO
!
!  Selime code:  cg_QG=zgrad0*HBH'*zcglwk
!
!^          cg_QG(innLoop,outLoop)=0.0_r8
!^
            tl_cg_QG(innLoop)=0.0_r8
            DO iobs=1,Ndatum(ng)
!^            cg_QG(innLoop,outLoop)=cg_QG(innLoop,outLoop)+            &
!^   &                               TLmodVal(iobs)*                    &
!^   &                               vgrad0(iobs)
!^
              tl_cg_QG(innLoop)=tl_cg_QG(innLoop)+                      &
     &                          TLmodVal(iobs)*                         &
     &                          vgrad0(iobs)+                           &
     &                          TLmodVal_S(iobs,innLoop,outLoop)*       &
     &                          tl_vgrad0(iobs)/                        &
     &                          cg_beta(innLoop,outLoop)
            END DO
!
!  AUGMENTED.
!
            DO iobs=1,Ndatum(ng)
!^            cg_QG(innLoop,outLoop)=cg_QG(innLoop,outLoop)+            &
!^   &                               Hbk(iobs,outLoop)*                 &
!^   &                               vcglwk(Ndatum(ng)+1,innLoop,       &
!^   &                                      outLoop)*                   &
!^   &                               vgrad0(iobs)+                      &
!^   &                               Hbk(iobs,outLoop)*                 &
!^   &                               vcglwk(iobs,innLoop,outLoop)*      &
!^   &                               vgrad0(Ndatum(ng)+1)
!^
              tl_cg_QG(innLoop)=tl_cg_QG(innLoop)+                      &
     &                          tl_Hbk(iobs)*                           &
     &                          vcglwk(Ndatum(ng)+1,innLoop,outLoop)*   &
     &                          vgrad0(iobs)+                           &
     &                          Hbk(iobs,outLoop)*                      &
     &                          tl_vcglwk(Ndatum(ng)+1,innLoop)*        &
     &                          vgrad0(iobs)+                           &
     &                          Hbk(iobs,outLoop)*                      &
     &                          vcglwk(Ndatum(ng)+1,innLoop,outLoop)*   &
     &                          tl_vgrad0(iobs)+                        &
     &                          tl_Hbk(iobs)*                           &
     &                          vcglwk(iobs,innLoop,outLoop)*           &
     &                          vgrad0(Ndatum(ng)+1)+                   &
     &                          Hbk(iobs,outLoop)*                      &
     &                          tl_vcglwk(iobs,innLoop)*                &
     &                          vgrad0(Ndatum(ng)+1)+                   &
     &                          Hbk(iobs,outLoop)*                      &
     &                          vcglwk(iobs,innLoop,outLoop)*           &
     &                          tl_vgrad0(Ndatum(ng)+1)
            END DO
!^          cg_QG(innLoop,outLoop)=cg_QG(innLoop,outLoop)+              &
!^   &                             Jb0(outLoop-1)*                      &
!^   &                             vcglwk(Ndatum(ng)+1,innLoop,         &
!^   &                                    outLoop)*                     &
!^   &                             vgrad0(Ndatum(ng)+1)
!^
            tl_cg_QG(innLoop)=tl_cg_QG(innLoop)+                        &
     &                        tl_Jb0(outLoop-1)*                        &
     &                        vcglwk(Ndatum(ng)+1,innLoop,outLoop)*     &
     &                        vgrad0(Ndatum(ng)+1)+                     &
     &                        Jb0(outLoop-1)*                           &
     &                        tl_vcglwk(Ndatum(ng)+1,innLoop)*          &
     &                        vgrad0(Ndatum(ng)+1)+                     &
     &                        Jb0(outLoop-1)*                           &
     &                        vcglwk(Ndatum(ng)+1,innLoop,outLoop)*     &
     &                        tl_vgrad0(Ndatum(ng)+1)
!
!!
!! TEST     IF (innLoop.gt.1) cg_QG(innLoop,outLoop)=0.0_r8
!!
            IF (innLoop.eq.NinnLoop) THEN
!
!  Calculate the new solution based upon the new, orthonormalized
!  Lanczos vector. First, the tridiagonal system is solved by
!  decomposition and forward/back substitution.
!
              IF (NinnLoop.eq.1) THEN
                print *,'Illegal configuration!'
                print *,'Ninner must be ge 2'
                STOP
              END IF
              IF (NinnLoop.eq.2) THEN
                zu(1)=cg_QG(1,outLoop)/cg_delta(1,outLoop)
                tl_zrhs(1)=tl_cg_QG(1)-tl_cg_delta(1)*zu(1)
                tl_zu(1)=tl_zrhs(1)/cg_delta(1,outLoop)
                zwork(1,3)=zu(1)
                tl_zwork(1,3)=tl_zu(1)
              ELSE
!
!  Compute zu first.
!
                zbet=cg_delta(1,outLoop)
                zu(1)=cg_QG(1,outLoop)/zbet
                DO ivec=2,innLoop-1
                  zgam(ivec)=cg_beta(ivec,outLoop)/zbet
                  zbet=cg_delta(ivec,outLoop)-                          &
     &                 cg_beta(ivec,outLoop)*zgam(ivec)
                  zu(ivec)=(cg_QG(ivec,outLoop)-                        &
     &                      cg_beta(ivec,outLoop)*                      &
     &                      zu(ivec-1))/zbet
                END DO
                zwork(innLoop-1,3)=zu(innLoop-1)
!
                DO ivec=innLoop-2,1,-1
                  zu(ivec)=zu(ivec)-zgam(ivec+1)*zu(ivec+1)
                  zwork(ivec,3)=zu(ivec)
                END DO
!
!  Now compute tl_zrhs.
!
                tl_zrhs(1)=tl_cg_QG(1)-                                 &
     &                     tl_cg_delta(1)*zu(1)-                        &
     &                     tl_cg_beta(2)*zu(2)
                DO ivec=2,innLoop-2
                  tl_zrhs(ivec)=tl_cg_QG(ivec)-                         &
     &                          tl_cg_beta(ivec)*zu(ivec-1)-            &
     &                          tl_cg_delta(ivec)*zu(ivec)-             &
     &                          tl_cg_beta(ivec+1)*zu(ivec+1)
                END DO
                tl_zrhs(innLoop-1)=tl_cg_QG(innLoop-1)-                 &
     &                             tl_cg_beta(innLoop-1)*zu(innLoop-2)- &
     &                             tl_cg_delta(innLoop-1)*zu(innLoop-1)
!
!  Now solve the TL tridiagonal system A*dx=b-dA*x
!
                zbet=cg_delta(1,outLoop)
                tl_zu(1)=tl_zrhs(1)/zbet
                DO ivec=2,innLoop-1
                  zgam(ivec)=cg_beta(ivec,outLoop)/zbet
                  zbet=cg_delta(ivec,outLoop)-                          &
     &                 cg_beta(ivec,outLoop)*zgam(ivec)
                  tl_zu(ivec)=(tl_zrhs(ivec)-                           &
     &                         cg_beta(ivec,outLoop)*                   &
     &                         tl_zu(ivec-1))/zbet
                END DO
                tl_zwork(innLoop-1,3)=tl_zu(innLoop-1)
!
                DO ivec=innLoop-2,1,-1
!^                zu(ivec)=zu(ivec)-zgam(ivec+1)*zu(ivec+1)
!^
                  tl_zu(ivec)=tl_zu(ivec)-zgam(ivec+1)*tl_zu(ivec+1)
                  tl_zwork(ivec,3)=tl_zu(ivec)
                END DO
              END IF
!
!  AUGMENTED.
!
              DO iobs=1,Ndatum(ng)+1
                px(iobs)=0.0_r8
                tl_px(iobs)=0.0_r8
                DO ivec=1,innLoop-1
                  px(iobs)=px(iobs)+                                    &
     &                     zcglwk(iobs,ivec,outLoop)*zwork(ivec,3)
                  tl_px(iobs)=tl_px(iobs)+                              &
     &                        tl_zcglwk(iobs,ivec)*                     &
     &                        zwork(ivec,3)+                            &
     &                        zcglwk(iobs,ivec,outLoop)*                &
     &                        tl_zwork(ivec,3)
                END DO
              END DO
!
!  If preconditioning, convert px from y-space to v-space.
!  We will always keep px in v-space.
!
!^            IF (Lprecond.and.(outLoop.gt.1)) THEN
!^              Lscale=1                 ! Spectral LMP
!^              Laug=.TRUE.
!^              CALL rprecond (ng, Lscale, Laug, outLoop, NinnLoop-1,   &
!^   &                         px)
!^            END IF

           END IF

         END IF
!
!  NEXT CALCULATION
!
!  After the initialization, for every other inner loop, calculate a
!  new Lanczos vector, store it in the matrix, and update search.
!
!  AUGMENTED.
!
          IF (innLoop.eq.1) THEN
            fact=1.0_r8/cg_Gnorm_v(outLoop)
          ELSE
            fact=1.0_r8/cg_beta(innLoop,outLoop)
          END IF
!
          DO iobs=1,Ndatum(ng)
!
!  Selime code: w=t1.
!
!^          pgrad(iobs)=TLmodVal(iobs)+                                 &
!^   &                  Hbk(iobs,outLoop)*                              &
!^   &                  vcglwk(Ndatum(ng)+1,innLoop,outLoop)
!^
            pgrad(iobs)=TLmodVal_S(iobs,innLoop,outLoop)*               &
     &                  fact+                                           &
     &                  Hbk(iobs,outLoop)*                              &
     &                  vcglwk(Ndatum(ng)+1,innLoop,outLoop)
            tl_pgrad(iobs)=TLmodVal(iobs)+                              &
     &                     tl_Hbk(iobs)*                                &
     &                     vcglwk(Ndatum(ng)+1,innLoop,outLoop)+        &
     &                     Hbk(iobs,outLoop)*                           &
     &                     tl_vcglwk(Ndatum(ng)+1,innLoop)
          END DO
!
          DO iobs=1,Ndatum(ng)
!
!  Convert gradient contribution from x-space to v-space.
!
            IF (ObsErr(iobs).ne.0.0_r8) THEN
!
!  Selime code: w=t1/R. No augmented term here since Raug^-1=[R^-1 0].
!
              pgrad(iobs)=pgrad(iobs)/ObsErr(iobs)
              tl_pgrad(iobs)=tl_pgrad(iobs)/ObsErr(iobs)-               &
     &                       tl_ObsErr(iobs)*pgrad(iobs)/ObsErr(iobs)
            END IF
          END DO
          pgrad(Ndatum(ng)+1)=0.0_r8
          tl_pgrad(Ndatum(ng)+1)=0.0_r8
!
!  Selime code: w=t1/R+z1
!
          DO iobs=1,Ndatum(ng)
            pgrad(iobs)=pgrad(iobs)+                                    &
     &                  vcglwk(iobs,innLoop,outLoop)
            tl_pgrad(iobs)=tl_pgrad(iobs)+                              &
     &                     tl_vcglwk(iobs,innLoop)
          END DO
          pgrad(Ndatum(ng)+1)=pgrad(Ndatum(ng)+1)+                      &
     &                        vcglwk(Ndatum(ng)+1,innLoop,outLoop)
          tl_pgrad(Ndatum(ng)+1)=tl_pgrad(Ndatum(ng)+1)+                &
     &                           tl_vcglwk(Ndatum(ng)+1,innLoop)
!
          IF (innLoop.gt.1) THEN
!
!  Selime code: w=w-v0*beta
!
!  AUGMENTED.
!
            DO iobs=1,Ndatum(ng)+1
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
!^        cg_delta(innLoop,outLoop)=0.0_r8
!^
          tl_cg_delta(innLoop)=0.0_r8
!
!  Selime code: alpha=w*t1=cg_delta
!
!  AUGMENTED.
!
          IF (innLoop.eq.1) THEN
            fact=1.0_r8/cg_Gnorm_v(outLoop)
          ELSE
            fact=1.0_r8/cg_beta(innLoop,outLoop)
          END IF
!
          DO iobs=1,Ndatum(ng)
!^          cg_delta(innLoop,outLoop)=cg_delta(innLoop,outLoop)+        &
!^   &                                TLmodVal(iobs)*                   &
!^   &                                pgrad(iobs)+                      &
!^   &                                Hbk(iobs,outLoop)*                &
!^   &                                vcglwk(Ndatum(ng)+1,innLoop,      &
!^   &                                       outLoop)*                  &
!^   &                                pgrad(iobs)+                      &
!^   &                                Hbk(iobs,outLoop)*                &
!^   &                                vcglwk(iobs,innLoop,outLoop)*     &
!^   &                                pgrad(Ndatum(ng)+1)
!^
            tl_cg_delta(innLoop)=tl_cg_delta(innLoop)+                  &
     &                           TLmodVal(iobs)*                        &
     &                           pgrad(iobs)+                           &
     &                           TLmodVal_S(iobs,innLoop,outLoop)*      &
     &                           tl_pgrad(iobs)*                        &
     &                           fact+                                  &
     &                           tl_Hbk(iobs)*                          &
     &                           vcglwk(Ndatum(ng)+1,innLoop,outLoop)*  &
     &                           pgrad(iobs)+                           &
     &                           Hbk(iobs,outLoop)*                     &
     &                           tl_vcglwk(Ndatum(ng)+1,innLoop)*       &
     &                           pgrad(iobs)+                           &
     &                           Hbk(iobs,outLoop)*                     &
     &                           vcglwk(Ndatum(ng)+1,innLoop,outLoop)*  &
     &                           tl_pgrad(iobs)+                        &
     &                           tl_Hbk(iobs)*                          &
     &                           vcglwk(iobs,innLoop,outLoop)*          &
     &                           pgrad(Ndatum(ng)+1)+                   &
     &                           Hbk(iobs,outLoop)*                     &
     &                           tl_vcglwk(iobs,innLoop)*               &
     &                           pgrad(Ndatum(ng)+1)+                   &
     &                           Hbk(iobs,outLoop)*                     &
     &                           vcglwk(iobs,innLoop,outLoop)*          &
     &                           tl_pgrad(Ndatum(ng)+1)
          END DO
!^        cg_delta(innLoop,outLoop)=cg_delta(innLoop,outLoop)+          &
!^   &                              Jb0(outLoop-1)*                     &
!^   &                              vcglwk(Ndatum(ng)+1,innLoop,        &
!^   &                                     outLoop)*                    &
!^   &                              pgrad(Ndatum(ng)+1)
!^
          tl_cg_delta(innLoop)=tl_cg_delta(innLoop)+                    &
     &                         tl_Jb0(outLoop-1)*                       &
     &                         vcglwk(Ndatum(ng)+1,innLoop,outLoop)*    &
     &                         pgrad(Ndatum(ng)+1)+                     &
     &                         Jb0(outLoop-1)*                          &
     &                         tl_vcglwk(Ndatum(ng)+1,innLoop)*         &
     &                         pgrad(Ndatum(ng)+1)+                     &
     &                         Jb0(outLoop-1)*                          &
     &                         vcglwk(Ndatum(ng)+1,innLoop,outLoop)*    &
     &                         tl_pgrad(Ndatum(ng)+1)
!
!  Selime code: w=w-alpha*v1
!
!  AUGMENTED.
!
          DO iobs=1,Ndatum(ng)+1
            pgrad(iobs)=pgrad(iobs)-                                    &
     &                  cg_delta(innLoop,outLoop)*                      &
     &                  zcglwk(iobs,innLoop,outLoop)
            tl_pgrad(iobs)=tl_pgrad(iobs)-                              &
     &                     tl_cg_delta(innLoop)*                        &
     &                     zcglwk(iobs,innLoop,outLoop)-                &
     &                     cg_delta(innLoop,outLoop)*                   &
     &                     tl_zcglwk(iobs,innLoop)
          END DO
!
!  Orthonormalize against previous directions.
!
!  Identify the appropriate Lanczos vector normalization coefficient.
!
          DO ivec=1,innLoop
            IF (ivec.eq.1) THEN
              zfact(ivec)=1.0_r8/cg_Gnorm_v(outLoop)
              tl_zfact(ivec)=-tl_cg_Gnorm_v(outLoop)*                   &
     &                       zfact(ivec)/                               &
     &                       cg_Gnorm_v(outLoop)
            ELSE
              zfact(ivec)=1.0_r8/cg_beta(ivec,outLoop)
              tl_zfact(ivec)=-tl_cg_beta(ivec)*                         &
     &                       zfact(ivec)/                               &
     &                       cg_beta(ivec,outLoop)
            ENDIF
          END DO
!
!  Selime code: cg_dla=pgrad*HBH'*v1=pgrad*TLmodVal_S/zfact
!
!  AUGMENTED.
!
          DO ivec=innLoop,1,-1
!^          cg_dla(ivec,outLoop)=0.0_r8
!^
            tl_dla=0.0_r8
            DO iobs=1,Ndatum(ng)
!^            cg_dla(ivec,outLoop)=cg_dla(ivec,outLoop)+                &
!^   &                             pgrad(iobs)*                         &
!^   &                             TLmodVal_S(iobs,ivec,outLoop)*       &
!^   &                             zfact(ivec)+                         &
!^   &                             pgrad(iobs)*                         &
!^   &                             Hbk(iobs,outLoop)*                   &
!^   &                             vcglwk(Ndatum(ng)+1,ivec,outLoop)+   &
!^   &                             Hbk(iobs,outLoop)*                   &
!^   &                             vcglwk(iobs,ivec,outLoop)*           &
!^   &                             pgrad(Ndatum(ng)+1)
!^
              tl_dla=tl_dla+                                            &
     &               tl_pgrad(iobs)*                                    &
     &               TLmodVal_S(iobs,ivec,outLoop)*                     &
     &               zfact(ivec)+                                       &
     &               pgrad(iobs)*                                       &
     &               tl_TLmodVal_S(iobs,ivec,outLoop)*                  &
     &               zfact(ivec)+                                       &
     &               pgrad(iobs)*                                       &
     &               TLmodVal_S(iobs,ivec,outLoop)*                     &
     &               tl_zfact(ivec)+                                    &
     &               tl_pgrad(iobs)*                                    &
     &               Hbk(iobs,outLoop)*                                 &
     &               vcglwk(Ndatum(ng)+1,ivec,outLoop)+                 &
     &               pgrad(iobs)*                                       &
     &               tl_Hbk(iobs)*                                      &
     &               vcglwk(Ndatum(ng)+1,ivec,outLoop)+                 &
     &               pgrad(iobs)*                                       &
     &               Hbk(iobs,outLoop)*                                 &
     &               tl_vcglwk(Ndatum(ng)+1,ivec)+                      &
     &               tl_Hbk(iobs)*                                      &
     &               vcglwk(iobs,ivec,outLoop)*                         &
     &               pgrad(Ndatum(ng)+1)+                               &
     &               Hbk(iobs,outLoop)*                                 &
     &               tl_vcglwk(iobs,ivec)*                              &
     &               pgrad(Ndatum(ng)+1)+                               &
     &               Hbk(iobs,outLoop)*                                 &
     &               vcglwk(iobs,ivec,outLoop)*                         &
     &               tl_pgrad(Ndatum(ng)+1)
            END DO
!^          cg_dla(ivec,outLoop)=cg_dla(ivec,outLoop)+                  &
!^   &                           Jb0(outLoop-1)*                        &
!^   &                           vcglwk(Ndatum(ng)+1,ivec,outLoop)*     &
!^   &                           pgrad(Ndatum(ng)+1)
!^
            tl_dla=tl_dla+                                              &
     &             tl_Jb0(outLoop-1)*                                   &
     &             vcglwk(Ndatum(ng)+1,ivec,outLoop)*                   &
     &             pgrad(Ndatum(ng)+1)+                                 &
     &             Jb0(outLoop-1)*                                      &
     &             tl_vcglwk(Ndatum(ng)+1,ivec)*                        &
     &             pgrad(Ndatum(ng)+1)+                                 &
     &             Jb0(outLoop-1)*                                      &
     &             vcglwk(Ndatum(ng)+1,ivec,outLoop)*                   &
     &             tl_pgrad(Ndatum(ng)+1)
            DO iobs=1,Ndatum(ng)+1
              pgrad(iobs)=pgrad(iobs)-                                  &
     &                    cg_dla(ivec,outLoop)*                         &
     &                    vcglwk(iobs,ivec,outLoop)
              tl_pgrad(iobs)=tl_pgrad(iobs)-tl_dla*                     &
     &                       vcglwk(iobs,ivec,outLoop)-                 &
     &                       cg_dla(ivec,outLoop)*                      &
     &                       tl_vcglwk(iobs,ivec)
            END DO
          END DO
!
!  Save the non-normalized Lanczos vector. zcglwk is used as temporary
!  storage. The calculation will be completed at the start of the next
!  inner-loop.
!
!  Selime code: w=zcglwk.
!
!  AUGMENTED.
!
          DO iobs=1,Ndatum(ng)+1
!^          zcglwk(iobs,innLoop+1,outLoop)=pgrad(iobs)
!^
            tl_zcglwk(iobs,innLoop+1)=tl_pgrad(iobs)
          END DO
!
!  If preconditioning, convert pgrad from v-space to y-space.
!
!
!  Selime code: z1=Gw.
!
!^        IF (Lprecond.and.(outLoop.gt.1)) THEN
!^          Lscale=1                 ! Spectral LMP
!^          Laug=.TRUE.
!^          CALL rprecond (ng, Lscale, Laug, outLoop, NinnLoop-1, pgrad)
!^        END IF
!
!  Selime code: z1=vcglwk.
!
!  AUGMENTED.
!
          DO iobs=1,Ndatum(ng)+1
!^          vcglwk(iobs,innLoop+1,outLoop)=pgrad(iobs)
!^
            tl_vcglwk(iobs,innLoop+1)=tl_pgrad(iobs)
          END DO
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
!^            ADmodVal(iobs)=px(iobs)
!^            tl_ADmodVal(iobs)=tl_px(iobs)
!^
              ADmodVal(iobs)=tl_px(iobs)
            END DO
            DO iobs=1,Ndatum(ng)+1
!^            FOURDVAR(ng)%cg_pxsum(iobs)=FOURDVAR(ng)%cg_pxsum(iobs)+  &
!^   &                                    FOURDVAR(ng)%cg_pxsave(iobs)
!^            FOURDVAR(ng)%cg_pxsave(iobs)=px(iobs)
!^
              FOURDVAR(ng)%tl_cg_pxsave(iobs)=tl_px(iobs)
!AMM         FOURDVAR(ng)%cg_pxsum(iobs)=FOURDVAR(ng)%cg_pxsum(iobs)+   &
!AMM &                                   px(iobs)
            END DO
          ELSE
!
!  Selime code: z1=ADmodVal.
!
            DO iobs=1,Ndatum(ng)
!^            ADmodVal(iobs)=vcglwk(iobs,innLoop+1,outLoop)
!^            tl_ADmodVal(iobs)=tl_vcglwk(iobs,innLoop+1)
!^
              ADmodVal(iobs)=tl_vcglwk(iobs,innLoop+1)
            END DO
!
          END IF
!
        END IF MINIMIZE

      END IF MASTER_THREAD

# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Broadcast new solution to other nodes.
!-----------------------------------------------------------------------
!
      CALL mp_bcastf (ng, model, ADmodVal)
      CALL mp_bcastf (ng, model, vgrad0(:))
      CALL mp_bcastf (ng, model, tl_ObsVal)
      CALL mp_bcastf (ng, model, tl_ObsErr)
      CALL mp_bcastf (ng, model, tl_ObsScale)
#  if defined RBL4DVAR          || defined R4DVAR      || \
      defined SENSITIVITY_4DVAR || defined TL_RBL4DVAR || \
      defined TL_R4DVAR
      CALL mp_bcastf (ng, model, tl_Hbk(:))
      CALL mp_bcastf (ng, model, tl_Jb0(outLoop-1))
#  endif
      CALL mp_bcasti (ng, model, info)
      CALL mp_bcastf (ng, model, tl_cg_beta(:))
      CALL mp_bcastf (ng, model, tl_cg_QG(:))
      CALL mp_bcastf (ng, model, tl_cg_delta(:))
      CALL mp_bcastf (ng, model, tl_zgrad0(:))
      CALL mp_bcastf (ng, model, tl_zcglwk(:,:))
      CALL mp_bcastf (ng, model, tl_vcglwk(:,:))
      CALL mp_bcastf (ng, model, FOURDVAR(ng)%tl_cg_pxsave(:))
# endif
!
      RETURN
      END SUBROUTINE tl_rpcg_lanczos

#else
      SUBROUTINE tl_rpcg_lanczos
      RETURN
      END SUBROUTINE tl_rpcg_lanczos
#endif
