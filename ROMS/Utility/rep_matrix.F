#include "cppdefs.h"
#if defined WEAK_CONSTRAINT && defined OBS_IMPACT
      SUBROUTINE rep_matrix (ng, model, outLoop, NinnLoop)
!
!svn $Id: rep_matrix.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine estimates the TRANSPOSE of the inverse of the          !
!  stabilized representer matrix using the Lanczos vectors from        !
!  the inner-loops.                                                    !
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
      USE distribute_mod, ONLY : mp_bcastf
# endif
!
      implicit none
!
!  Imported variable declarations
!
      integer, intent(in) :: ng, model, outLoop, NinnLoop
!
!  Local variable declarations.
!
      integer :: i, j, iobs, ivec, innLoop
# ifdef MINRES
      integer :: info
      integer, dimension(NinnLoop) :: ipiv
# endif

      real(r8) :: zbet

      real(r8), dimension(NinnLoop) :: zu, zgam, zt
# ifdef MINRES
      real(r8) :: zsum, zck, zgk
      real(r8), dimension(NinnLoop,NinnLoop) :: ztriT, zLT, zLTt
      real(r8), dimension(NinnLoop) :: tau, zwork1, zeref
# endif
# ifdef RPCG
      real(r8), dimension(NinnLoop) :: zfact
# endif
!
!-----------------------------------------------------------------------
!  Compute observation impact on the weak constraint 4DVAR data
!  assimilation systems.
!-----------------------------------------------------------------------
!
      MASTER_THREAD : IF (Master) THEN

        DO iobs=1,Ndatum(ng)
          ad_ObsVal(iobs)=0.0_r8
        END DO
        DO innLoop=1,NinnLoop
          zt(innLoop)=0.0_r8
        END DO
# ifdef RPCG
        DO innLoop=1,NinnLoop
          IF (innLoop.eq.1) THEN
            zfact(innLoop)=1.0_r8/cg_QG(1,outLoop)
          ELSE
            zfact(innLoop)=1.0_r8/cg_beta(innLoop,outLoop)
          ENDIF
        END DO
# endif
!
!  Multiply TLmodVal by the tranpose matrix of Lanczos vectors.
!  Note that the factor of 1/SQRT(ObsErr) is added to convert to
!  x-space.
!
        DO innLoop=1,NinnLoop
          DO iobs=1,Ndatum(ng)
            IF (ObsErr(iobs).NE.0.0_r8) THEN
# ifdef RPCG
              zt(innLoop)=zt(innLoop)+ObsScale(iobs)*                   &
     &                    zcglwk(iobs,innLoop,outLoop)*TLmodVal(iobs)
# else
              zt(innLoop)=zt(innLoop)+ObsScale(iobs)*                   &
     &                    zcglwk(iobs,innLoop,outLoop)*TLmodVal(iobs)/  &
     &                    SQRT(ObsErr(iobs))
# endif
            END IF
          END DO
        END DO

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
        DO i=1,NinnLoop
          ztriT(i,i)=cg_delta(i,outLoop)
        END DO
        DO i=1,NinnLoop-1
          ztriT(i,i+1)=cg_beta(i+1,outLoop)
        END DO
        DO i=2,NinnLoop
          ztriT(i,i-1)=cg_beta(i,outLoop)
        END DO
        CALL sqlq (NinnLoop, ztriT, tau, zwork1)
!
!   Isolate L=zLT and its tranpose.
!
        zLT=0.0_r8
        zLTt=0.0_r8
        DO i=1,NinnLoop
          DO j=1,i
            zLT(i,j)=ztriT(i,j)
          END DO
        END DO
        DO j=1,NinnLoop
          DO i=1,NinnLoop
            zLTt(i,j)=zLT(j,i)
          END DO
        END DO
!
!  Compute zck.
!
        zgk=SQRT(zLT(NinnLoop,NinnLoop)*zLT(NinnLoop,NinnLoop)+         &
     &      cg_beta(NinnLoop+1,outLoop)*cg_beta(NinnLoop+1,outLoop))
        zck=zLT(NinnLoop,NinnLoop)/zgk
!
!  Now form inv(L)*zt - NOTE: we use L not L', so we use the
!  adjoint of the original solver.
!
        DO j=1,NinnLoop
          DO i=j-1,1,-1
            zt(j)=zt(j)-zt(i)*zLTt(i,j)
          END DO
          zt(j)=zt(j)/zLTt(j,j)
        END DO
!
!   Now form zt=D*zt.
!
        zt(NinnLoop)=zck*zt(NinnLoop)
!
!   Now form zt=Q'*zt.
!
        DO i=NinnLoop,1,-1
          DO j=1,NinnLoop
            zeref(j)=0.0_r8
          END DO
          zeref(i)=1.0_r8
          DO j=i+1,NinnLoop
            zeref(j)=ztriT(i,j)
          END DO
          zsum=0.0_r8
          DO j=1,NinnLoop
            zsum=zsum+zt(j)*zeref(j)
          END DO
          DO j=1,NinnLoop
            zt(j)=zt(j)-tau(i)*zsum*zeref(j)
          END DO
        END DO
!
!   Copy the solution zt into zu.
!
        DO i=1,NinnLoop
          zu(i)=zt(i)
        END DO
# else
!
!  Now multiply the result by the inverse tridiagonal matrix.
!
        zbet=cg_delta(1,outLoop)
        zu(1)=zt(1)/zbet
        DO ivec=2,NinnLoop
          zgam(ivec)=cg_beta(ivec,outLoop)/zbet
          zbet=cg_delta(ivec,outLoop)-cg_beta(ivec,outLoop)*zgam(ivec)
          zu(ivec)=(zt(ivec)-cg_beta(ivec,outLoop)*                     &
     &              zu(ivec-1))/zbet
        END DO

        DO ivec=NinnLoop-1,1,-1
          zu(ivec)=zu(ivec)-zgam(ivec+1)*zu(ivec+1)
        END DO
# endif
!
!  Finally multiply by the matrix of Lanczos vactors.
# ifndef RPCG
!  Note that the factor of 1/SQRT(ObsErr) is added to covert to
!  x-space.
# endif
!
        DO iobs=1,Ndatum(ng)
          DO innLoop=1,NinnLoop
            IF (ObsErr(iobs).NE.0.0_r8) THEN
# ifdef RPCG
              ad_ObsVal(iobs)=ad_ObsVal(iobs)+                          &
     &                        ObsScale(iobs)*                           &
     &                        TLmodVal_S(iobs,innLoop,outLoop)*         &
     &                        zu(innLoop)*zfact(innLoop)/ObsErr(iobs)
#  else
              ad_ObsVal(iobs)=ad_ObsVal(iobs)+                          &
     &                        ObsScale(iobs)*                           &
     &                        zcglwk(iobs,innLoop,outLoop)*zu(innLoop)/ &
     &                        SQRT(ObsErr(iobs))
#  endif
            END IF
          END DO
        END DO

      END IF MASTER_THREAD

# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Broadcast new solution to other nodes.
!-----------------------------------------------------------------------
!
      CALL mp_bcastf (ng, model, ad_ObsVal)
# endif
!
      RETURN
      END SUBROUTINE rep_matrix
#else
      SUBROUTINE rep_matrix
      RETURN
      END SUBROUTINE rep_matrix
#endif
