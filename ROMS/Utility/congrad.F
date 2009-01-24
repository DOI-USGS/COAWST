#include "cppdefs.h"

#if defined W4DVAR || defined W4DPSAS
      SUBROUTINE congrad (ng, outer, inner, Ninner, converged)
!
!svn $Id: congrad.F 503 2008-01-10 00:11:51Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Weak Constraint 4-Dimensional Variational (4DVar) Pre-conditioned   !
!                     Conjugate Gradient Algorithm                     !
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
!  0.5*x*A*x-x*b, assuming an initial guess of x=0. In this case the   !
!  gradient is Ax-b and the Hessian is A.                              !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Chua, B. S. and A. F. Bennett,  2001:  An inverse ocean modeling  !
!      sytem, Ocean Modelling, 3, 137-165.                             !
!                                                                      !
!  Lanczos Algorithm Reference:                                        !
!                                                                      !
!    Fisher, M., 1998: Minimization Algorithms for Variational Data    !
!      Assimilation. In Seminar on Recent Developments in Numerical    !
!      Methods for Atmospheric Modelling, 1998.                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar

# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti, mp_bcastl
# endif
      implicit none
!
!  Imported variable declarations
!
      integer, intent(in) :: ng, outer, inner, Ninner
      logical, intent(out) :: converged
!
!  Local variable declarations.
!
      integer :: iobs, jinner
      real(r8), dimension(Ninner) :: zu, zgam
      real(r8), dimension(Ndatum(ng)) :: px, pgrad, zw
      real(r8), dimension(Ninner,3) :: zwork
      real(r8) :: dla, zbet, eps, preduc, Jo


!
!-----------------------------------------------------------------------
!  Since this is not run in parallel, we should be careful about who
!  is calculating these vectors. Only Master will run this so that
!  all nodes have computed their respective solutions.
!-----------------------------------------------------------------------
!

      IF (Master) THEN

!
!-----------------------------------------------------------------------
!  If this is the first call to the gradient algorithm, we will
!  set up the vectors and store the first guess.
!  The initial starting guess is assumed to be zero in which case
!  the gradient is just -(obs-modval).
!-----------------------------------------------------------------------
!
        IF (inner.eq.0) THEN
          zgnorm=0.0
          DO iobs=1,Ndatum(ng)
# ifdef W4DPSAS
            pgrad(iobs)=-1.0_r8*(ObsScale(iobs)*                        &
     &                 (ObsVal(iobs)-NLmodVal(iobs)))
# else
            pgrad(iobs)=-1.0_r8*(ObsScale(iobs)*                        &
     &                 (ObsVal(iobs)-TLmodVal(iobs)))
# endif
            zgrad0(iobs)=pgrad(iobs)
            zgnorm=zgnorm+zgrad0(iobs)*zgrad0(iobs)
          END DO
          zgnorm=SQRT(zgnorm)

          DO iobs=1,Ndatum(ng)
            zcglwk(iobs,1)=pgrad(iobs)/zgnorm
            ADmodVal(iobs)=zcglwk(iobs,1)
          END DO

          zqg0(1)=0.0_r8
          DO iobs=1,Ndatum(ng)
            zqg0(1)=zqg0(1)+zcglwk(iobs,1)*zgrad0(iobs)
          END DO

        ELSE

!
!-----------------------------------------------------------------------
!  After the initialization, for every other inner value, we will
!  calculate a new Lanczos vector, store it in the matrix, and update
!  our search
!-----------------------------------------------------------------------
!
          zdelta(inner)=0.0
          DO iobs=1,Ndatum(ng)
            pgrad(iobs)=ObsScale(iobs)*                                 &
     &                 (TLmodVal(iobs)+ObsErr(iobs)*zcglwk(iobs,inner))
            zdelta(inner)=zdelta(inner)+zcglwk(iobs,inner)*pgrad(iobs)
          END DO

! HERNAN -- Do something here. If zdelta is < 0 then we are in deep doo-doo
          IF (zdelta(inner).le.0.0) THEN
            print *,'zdelta not positive definite!'
            print *,'inner=',inner,' zdelta=',zdelta(inner)
            STOP
          END IF
!
!    Compute the new Lanczos vector.
!
          DO iobs=1,Ndatum(ng)
            pgrad(iobs)=pgrad(iobs)-zdelta(inner)*zcglwk(iobs,inner)
          END DO

          IF (inner.gt.1) THEN
            DO iobs=1,Ndatum(ng)
              pgrad(iobs)=pgrad(iobs)-zbeta(inner)*zcglwk(iobs,inner-1)
            END DO
          END IF

!
!    Orthonormalize against previous Lanczos vectors.
!
          DO jinner=inner,1,-1
            dla=0.0_r8
            DO iobs=1,Ndatum(ng)
              dla=dla+pgrad(iobs)*zcglwk(iobs,jinner)
            END DO
            DO iobs=1,Ndatum(ng)
              pgrad(iobs)=pgrad(iobs)-dla*zcglwk(iobs,jinner)
            END DO
          END DO

          zbeta(inner+1)=0.0_r8
          DO iobs=1,Ndatum(ng)
            zbeta(inner+1)=zbeta(inner+1)+pgrad(iobs)*pgrad(iobs)
          END DO
          zbeta(inner+1)=SQRT(zbeta(inner+1))

          DO iobs=1,Ndatum(ng)
            zcglwk(iobs,inner+1)=pgrad(iobs)/zbeta(inner+1)
          END DO

          zqg0(inner+1)=0.0_r8
          DO iobs=1,Ndatum(ng)
            zqg0(inner+1)=zqg0(inner+1)+zcglwk(iobs,inner+1)*zgrad0(iobs)
          END DO


!
!-----------------------------------------------------------------------
!  Calculate the new solution based upon the new, orthonormalized
!  Lanczos vector. First, the tridiagonal system is solved by
!  decomposition and forward/back substitution.
!-----------------------------------------------------------------------
!
          zbet=zdelta(1)
          zu(1)=-zqg0(1)/zbet
          DO jinner=2,inner
            zgam(jinner)=zbeta(jinner)/zbet
            zbet=zdelta(jinner)-zbeta(jinner)*zgam(jinner)
            zu(jinner)=(-zqg0(jinner)-zbeta(jinner)*zu(jinner-1))/zbet
          END DO
          zwork(inner,3)=zu(inner)

          DO jinner=inner-1,1,-1
            zu(jinner)=zu(jinner)-zgam(jinner+1)*zu(jinner+1)
            zwork(jinner,3)=zu(jinner)
          END DO

          DO iobs=1,Ndatum(ng)
            zw(iobs)=zgrad0(iobs)+zbeta(inner+1)*zcglwk(iobs,inner+1)   &
     &                       *zwork(inner,3)
          END DO

          DO iobs=1,Ndatum(ng)
            px(iobs)=0.0
            DO jinner=1,inner
              px(iobs)=px(iobs)+zcglwk(iobs,jinner)*zwork(jinner,3)
              zw(iobs)=zw(iobs)-zcglwk(iobs,jinner)*zqg0(jinner)
            END DO
          END DO
!
!   Compute the reduction in the gradient Ax-b.
!
          preduc=0.0_r8
          DO iobs=1,Ndatum(ng)
            preduc=preduc+zw(iobs)*zw(iobs)
          END DO
          preduc=SQRT(preduc)/zgnorm

!
!   Estimate the residual by: ||Ax-b||=zbeta(k+1)
!
          eps=ABS(zbeta(inner+1))

!
!   Estimate the cost function
!
          Jo=0.0_r8
          DO iobs=1,Ndatum(ng)
            Jo=Jo-px(iobs)*zgrad0(iobs)
          END DO
! Temporary printing to track the progress
  print *,'zdelta=',(zdelta(jinner),jinner=1,inner)
  print *,'zbeta=',(zbeta(jinner),jinner=1,inner+1)
  print *,'zwork=',(zwork(jinner,3),jinner=1,inner)
  print *,'eps=',eps
  print *,'preduc=',preduc
  print *,'Jo=',Jo
  print *,'px=',(px(iobs),iobs=1,Ndatum(ng))

!
!   Check the convergence criteria
!-----------------------------------------------------------------------
!  The Lanczos algorithm should not be run for more loops than the
!  size of the data. This will cause problems. If we have performed
!  more loops than we have datum, then we set converged to TRUE so 
!  that the inner loop stops. 
!-----------------------------------------------------------------------
!
          converged=.FALSE.
          IF (preduc.le.CGeps.or.inner.eq.Ndatum(ng)) THEN
            converged=.TRUE.
          END IF

!
!   Put the new trial solution into the adjoint vector for the next loop
!   Put the final solution into the adjoint vector when converged
!   of on the final inner-loop.
!
          IF ((inner.eq.Ninner).or.converged) THEN
            DO iobs=1,Ndatum(ng)
              ADmodVal(iobs)=px(iobs)
            END DO
          ELSE
            DO iobs=1,Ndatum(ng)
              ADmodVal(iobs)=zcglwk(iobs,inner+1)
            END DO
          END IF
        END IF
! DEBUG
!       ADmodVal=(/1.0, 0.0, 0.0/)
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcastf (ng, iTLM, ADmodVal, Mobs)
# endif

!
!  Write out conjugate gradient vectors into 4DVAR NetCDF file.
!
!      CALL cg_write (ng, converged)
!
!  Report minimization parameters.
!
!      IF (Master) THEN
!        PRINT 10, outer, inner, Ndatum(ng), converged,                  &
!     &            eps,      outer, inner,                               &
!     &            CGeps,    outer, inner,                               &
!     &            zbeta,    outer, inner,                               &
!     &            cg_sigma, outer, inner,                               &
!     &            cg_rnorm, outer, inner,                               &
!     &            MINVAL(cg_x), MAXVAL(cg_x),                           &
!     &            MINVAL(cg_z), MAXVAL(cg_z),                           &
!     &            MINVAL(cg_p), MAXVAL(cg_p),                           &
!     &            MINVAL(cg_v), MAXVAL(cg_v),                           &
!     &            MINVAL(cg_s), MAXVAL(cg_s),                           &
!     &            MINVAL(cg_r), MAXVAL(cg_r)
! 10     FORMAT (/,' CONGRAD - Conjugate Gradient Information: ',/,      &
!     &          /,11x,'Outer Loop = ',i3.3,                             &
!     &          /,11x,'Inner Loop = ',i3.3,                             &
!     &          /,11x,'Ndatum     = ',i7.7,                             &
!     &          /,11x,'converged  = ',l1,/,                             &
!     &          /,11x,'cg_eps   = ',1p,e15.8,0p,                        &
!     &            3x,'(',i3.3,', ',i3.3,')',                            &
!     &          /,11x,'cg_tol   = ',1p,e15.8,0p,                        &
!     &            3x,'(',i3.3,', ',i3.3,')',                            &
!     &          /,11x,'cg_gamma = ',1p,e15.8,0p,                        &
!     &            3x,'(',i3.3,', ',i3.3,')',                            &
!     &          /,11x,'cg_sigma = ',1p,e15.8,0p,                        &
!     &            3x,'(',i3.3,', ',i3.3,')',                            &
!     &          /,11x,'cg_rnorm = ',1p,e15.8,0p,                        &
!     &            3x,'(',i3.3,', ',i3.3,')',/,                          &
!     &          /,11x,'Min cg_x = ',1p,e15.8,0p,2x,                     &
!     &                'Max cg_x = ',1p,e15.8,0p,                        &
!     &          /,11x,'Min cg_z = ',1p,e15.8,0p,2x,                     &
!     &                'Max cg_z = ',1p,e15.8,0p,                        &
!     &          /,11x,'Min cg_p = ',1p,e15.8,0p,2x,                     &
!     &                'Max cg_p = ',1p,e15.8,0p,                        &
!     &          /,11x,'Min cg_v = ',1p,e15.8,0p,2x,                     &
!     &                'Max cg_v = ',1p,e15.8,0p,                        &
!     &          /,11x,'Min cg_s = ',1p,e15.8,0p,2x,                     &
!     &                'Max cg_s = ',1p,e15.8,0p,                        &
!     &          /,11x,'Min cg_r = ',1p,e15.8,0p,2x,                     &
!     &                'Max cg_r = ',1p,e15.8,0p)
!      END IF

      RETURN
      END SUBROUTINE congrad

      SUBROUTINE cg_write (ng, converged)
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
      logical, intent(in) :: converged

      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: status, varid
      integer :: start(2), total(2)
!
!-----------------------------------------------------------------------
!  Write out conjugate gradient vectors.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
!
!  Write out outer and inner iteration.
!
        status=nf90_inq_varid(ncMODid(ng), 'outer', varid)
        status=nf90_put_var(ncMODid(ng), varid, outer)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'outer', TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF

        status=nf90_inq_varid(ncMODid(ng), 'inner', varid)
        status=nf90_put_var(ncMODid(ng), varid, inner)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'inner', TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Write out representer coefficients.
!
        start(1)=1
        total(1)=Ndatum(ng)
        status=nf90_put_var(ncMODid(ng), modVid(idRepC,ng),             &
     &                      ADmodVal, start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRepC)), TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Write out conjugate gradient vectors.
!
!  HERNAN -- Here we need to write out zdelta, zbeta, and zcglwk
!
       END IF
!
!-----------------------------------------------------------------------
!  Synchronize observations NetCDF file to disk.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
        status=nf90_sync(ncMODid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

  10  FORMAT (/,' CG_WRITE - error while writing variable: ',a,/,       &
     &        12x,'into NetCDF file: ',a)
  20  FORMAT (/,' CG_WRITE - unable to synchronize 4DVAR',              &
     &        1x,'NetCDF file to disk.')

      END SUBROUTINE cg_write
#else
      SUBROUTINE congrad
      RETURN
      END SUBROUTINE congrad
#endif
