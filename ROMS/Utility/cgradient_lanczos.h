!
!svn $Id: cgradient_lanczos.h 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module minimizes a quadratic cost function using the conjugate !
!  gradient algorithm proposed by Mike Fisher (ECMWF).                 !
!                                                                      !
!  These routines exploit the  close connection  between the conjugate !
!  gradient minimization and the Lanczos algorithm:                    !
!                                                                      !
!    q(k) = g(k) / ||g(k)||                                            !
!                                                                      !
!  If we eliminate the  descent directions and multiply by the Hessian !
!  matrix, we get the Lanczos recurrence relationship:                 !
!                                                                      !
!    H q(k+1) = Gamma(k+1) q(k+2) + Delta(k+1) q(k+1) + Gamma(k) q(k)  !
!                                                                      !
!  with                                                                !
!                                                                      !
!    Delta(k+1) = (1 / Alpha(k+1)) + (Beta(k+1) / Alpha(k))            !
!                                                                      !
!    Gamma(k) = - SQRT(Beta(k+1)) / Alpha(k)                           !
!                                                                      !
!  since the gradient and Lanczos vectors  are mutually orthogonal the !
!  recurrence maybe written in matrix form as:                         !
!                                                                      !
!    H Q(k) = Q(k) T(k) + Gamma(k) q(k+1) e'(k)                        !
!                                                                      !
!  with                                                                !
!                                                                      !
!            { q(1), q(2), q(3), ..., q(k) }                           !
!    Q(k) =  {  .     .     .          .   }                           !
!            {  .     .     .          .   }                           !
!            {  .     .     .          .   }                           !
!                                                                      !
!            { Delta(1)  Gamma(1)                                }     !
!            { Gamma(1)  Delta(2)  Gamma(2)                      }     !
!            {         .         .         .                     }     !
!    T(k) =  {          .         .         .                    }     !
!            {           .         .         .                   }     !
!            {              Gamma(k-2)   Delta(k-1)   Gamma(k-1) }     !
!            {                           Gamma(k-1)   Delta(k)   }     !
!                                                                      !
!    e'(k) = { 0, ...,0, 1 }                                           !
!                                                                      !
!  The eigenvalues of  T(k)  and the vectors formed by  Q(k)*T(k)  are !
!  approximations to the eigenvalues and eigenvectors of the  Hessian. !
!  They can be used for pre-conditioning.                              !
!                                                                      !
!  The tangent linear model conditions and associated adjoint in terms !
!  of the Lanzos algorithm are:                                        !
!                                                                      !
!    X(k) = X(0) + Q(k) Z(k)                                           !
!                                                                      !
!    T(k) Z(k) = - transpose[Q(k0)] g(0)                               !
!                                                                      !
!  where                                                               !
!                                                                      !
!    k           Inner loop iteration                                  !
!    Alpha(k)    Conjugate gradient coefficient                        !
!    Beta(k)     Conjugate gradient coefficient                        !
!    Delta(k)    Lanczos algorithm coefficient                         !
!    Gamma(k)    Lanczos algorithm coefficient                         !
!    H           Hessian matrix                                        !
!    Q(k)        Matrix of orthonormal Lanczos vectors                 !
!    T(k)        Symmetric, tri-diagonal matrix                        !
!    Z(k)        Eigenvectors of Q(k)*T(k)                             !
!    e'(k)       Tansposed unit vector                                 !
!    g(k)        Gradient vectors (adjoint solution: GRAD(J))          !
!    q(k)        Lanczos vectors                                       !
!    <...>       Dot product                                           !
!    ||...||     Euclidean norm, ||g(k)|| = SQRT( <g(k),g(k)> )        !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Fisher, M., 1997: Efficient Minimization of Quadratic Penalty     !
!      funtions, unpublish manuscript, 1-14.                           !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: cgradient

      CONTAINS
!
!***********************************************************************
      SUBROUTINE cgradient (ng, tile, model, innLoop, outLoop)
!***********************************************************************
!
      USE mod_param
#ifdef SOLVE3D
      USE mod_coupling
#endif
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
#endif
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, innLoop, outLoop
!
!  Local variable declarations.
!
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, model, 36)
#endif
      CALL cgradient_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Lold(ng), Lnew(ng),                          &
     &                     innLoop, outLoop,                            &
#ifdef MASKING
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
#endif
#ifdef ADJUST_WSTRESS
     &                     FORCES(ng) % tl_ustr,                        &
     &                     FORCES(ng) % tl_vstr,                        &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                     FORCES(ng) % tl_tflux,                       &
# endif
     &                     OCEAN(ng) % tl_t,                            &
     &                     OCEAN(ng) % tl_u,                            &
     &                     OCEAN(ng) % tl_v,                            &
#else
     &                     OCEAN(ng) % tl_ubar,                         &
     &                     OCEAN(ng) % tl_vbar,                         &
#endif
     &                     OCEAN(ng) % tl_zeta,                         &
#ifdef ADJUST_WSTRESS
     &                     FORCES(ng) % d_sustr,                        &
     &                     FORCES(ng) % d_svstr,                        &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                     FORCES(ng) % d_stflx,                        &
# endif
     &                     OCEAN(ng) % d_t,                             &
     &                     OCEAN(ng) % d_u,                             &
     &                     OCEAN(ng) % d_v,                             &
#else
     &                     OCEAN(ng) % d_ubar,                          &
     &                     OCEAN(ng) % d_vbar,                          &
#endif
     &                     OCEAN(ng) % d_zeta,                          &
#ifdef ADJUST_WSTRESS
     &                     FORCES(ng) % ad_ustr,                        &
     &                     FORCES(ng) % ad_vstr,                        &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                     FORCES(ng) % ad_tflux,                       &
# endif
     &                     OCEAN(ng) % ad_t,                            &
     &                     OCEAN(ng) % ad_u,                            &
     &                     OCEAN(ng) % ad_v,                            &
#else
     &                     OCEAN(ng) % ad_ubar,                         &
     &                     OCEAN(ng) % ad_vbar,                         &
#endif
     &                     OCEAN(ng) % ad_zeta)
#ifdef PROFILE
      CALL wclock_on (ng, model, 36)
#endif
      RETURN
      END SUBROUTINE cgradient
!
!***********************************************************************
      SUBROUTINE cgradient_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           Lold, Lnew,                            &
     &                           innLoop, outLoop,                      &
#ifdef MASKING
     &                           rmask, umask, vmask,                   &
#endif
#ifdef ADJUST_WSTRESS
     &                           tl_ustr, tl_vstr,                      &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                           tl_tflux,                              &
# endif
     &                           tl_t, tl_u, tl_v,                      &
#else
     &                           tl_ubar, tl_vbar,                      &
#endif
     &                           tl_zeta,                               &
#ifdef ADJUST_WSTRESS
     &                           d_sustr, d_svstr,                      &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                           d_stflx,                               &
# endif
     &                           d_t, d_u, d_v,                         &
#else
     &                           d_ubar, d_vbar,                        &
#endif
     &                           d_zeta,                                &
#ifdef ADJUST_WSTRESS
     &                           ad_ustr, ad_vstr,                      &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                           ad_tflux,                              &
# endif
     &                           ad_t, ad_u, ad_v,                      &
#else
     &                           ad_ubar, ad_vbar,                      &
#endif
     &                           ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf, mp_bcastf_m, mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lold, Lnew
      integer, intent(in) :: innLoop, outLoop
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: d_sustr(LBi:,LBj:,:)
      real(r8), intent(inout) :: d_svstr(LBi:,LBj:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: d_stflx(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(inout) :: d_t(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: d_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: d_v(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: d_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: d_vbar(LBi:,LBj:)
# endif
      real(r8), intent(inout) :: d_zeta(LBi:,LBj:)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: d_sustr(LBi:UBi,LBj:UBj,Nfrec(ng))
      real(r8), intent(inout) :: d_svstr(LBi:UBi,LBj:UBj,Nfrec(ng))
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: d_stflx(LBi:UBi,LBj:UBj,               &
     &                                   Nfrec(ng),NT(ng))
#  endif
      real(r8), intent(inout) :: d_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: d_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: d_v(LBi:UBi,LBj:UBj,N(ng))
# else
      real(r8), intent(inout) :: d_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: d_vbar(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(inout) :: d_zeta(LBi:UBi,LBj:UBj)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: Linp, Lout, Lwrk, Lwrk1, i, ic
      integer :: info, ingood, itheta1

      real(r8) :: norm, zbeta, ztheta1

      real(r8), dimension(0:NstateVar(ng)) :: Adjust
      real(r8), dimension(0:NstateVar(ng)) :: dot_old, dot_new
      real(r8), dimension(0:NstateVar(ng)) :: old_dot, new_dot
      real(r8), dimension(2*Ninner-2) :: work

      character (len=13) :: string

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize trial step size.
!-----------------------------------------------------------------------
!
      cg_tau(innLoop,outLoop)=CGstepI            
      IF (innLoop.eq.0) THEN
        ingood=0                  ! number of good eigenpairs
        DO i=0,NstateVar(ng)
          dot_old(i)=0.0_r8
          dot_new(i)=0.0_r8
          old_dot(i)=0.0_r8
          new_dot(i)=0.0_r8
          FOURDVAR(ng)%CostGradDot(i)=0.0_r8
        END DO
      END IF
      IF (Master) WRITE (stdout,10)
 10   FORMAT (/,' <<<< Descent Algorithm >>>>',/)
!
!  Estimate the Hessian and save the starting vector in ad_*(Lold).
!
      IF (innLoop.gt.0) THEN
        Lwrk=2
        Linp=1
        Lout=2
        CALL hessian (ng, tile, model,                                  &
     &                LBi, UBi, LBj, UBj,                               &
     &                Linp, Lout, Lwrk,                                 &
     &                innLoop, outLoop,                                 &
#ifdef MASKING
     &                rmask, umask, vmask,                              &
#endif
#ifdef ADJUST_WSTRESS
     &                ad_ustr, ad_vstr,                                 &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                ad_tflux,                                         &
# endif
     &                ad_t, ad_u, ad_v,                                 &
#else
     &                ad_ubar, ad_vbar,                                 &
#endif
     &                ad_zeta,                                          &
#ifdef ADJUST_WSTRESS
     &                tl_ustr, tl_vstr,                                 &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                tl_tflux,                                         &
# endif
     &                tl_t, tl_u, tl_v,                                 &
#else
     &                tl_ubar, tl_vbar,                                 &
#endif
     &                tl_zeta)
!
!  Check for positive Hessian, J''.
!
        IF (cg_delta(innLoop,outLoop).le.0.0_r8) THEN
          PRINT *,'CG_DELTA not positive'
          PRINT *, 'CG_DELTA = ', innLoop, cg_delta(innLoop,outLoop)
          STOP
        END IF
      END IF
!
!  Apply the Lanczos recurrence and orthonormalize.
!
      Linp=1
      Lout=2
      Lwrk=2
      CALL lanczos (ng, tile, model,                                    &
     &              LBi, UBi, LBj, UBj,                                 &
     &              Linp, Lout, Lwrk,                                   &
     &              innLoop, outLoop,                                   &
#ifdef MASKING
     &              rmask, umask, vmask,                                &
#endif
#ifdef ADJUST_WSTRESS
     &              tl_ustr, tl_vstr,                                   &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &              tl_tflux,                                           &
# endif
     &              tl_t, tl_u, tl_v,                                   &
#else
     &              tl_ubar, tl_vbar,                                   &
#endif
     &              tl_zeta,                                            &
#ifdef ADJUST_WSTRESS
     &              ad_ustr, ad_vstr,                                   &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &              ad_tflux,                                           &
# endif
     &              ad_t, ad_u, ad_v,                                   &
#else
     &              ad_ubar, ad_vbar,                                   &
#endif
     &              ad_zeta)
!
!  Compute new direction, d(k+1).
!
      CALL new_direction (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Linp, Lout,                                   &
#ifdef MASKING
     &                    rmask, umask, vmask,                          &
#endif
#ifdef ADJUST_WSTRESS
     &                    ad_ustr, ad_vstr,                             &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                    ad_tflux,                                     &
# endif
     &                    ad_t, ad_u, ad_v,                             &
#else
     &                    ad_ubar, ad_vbar,                             &
#endif
     &                    ad_zeta,                                      &
#ifdef ADJUST_WSTRESS
     &                    d_sustr, d_svstr,                             &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                    d_stflx,                                      &
# endif
     &                    d_t, d_u, d_v,                                &
#else
     &                    d_ubar, d_vbar,                               &
#endif
     &                    d_zeta)
!
!-----------------------------------------------------------------------
!  Calculate the reduction in the gradient norm by solving a 
!  tridiagonal system.
!-----------------------------------------------------------------------
!
      IF (innLoop.gt.1) THEN
!
!  Decomposition and forward substitution.
!
        zbeta=cg_delta(1,outLoop)
        cg_zu(1,outLoop)=-cg_QG(1,outLoop)/zbeta
        DO i=2,innLoop
          cg_gamma(i,outLoop)=cg_beta(i,outLoop)/zbeta
          zbeta=cg_delta(i,outLoop)-                                    &
     &          cg_beta(i,outLoop)*cg_gamma(i,outLoop)
          cg_zu(i,outLoop)=(-cg_QG(i,outLoop)-                          &
     &                      cg_beta(i,outLoop)*cg_zu(i-1,outLoop))/zbeta
        END DO
!
!  Back substitution.
!
        cg_Tmatrix(innLoop,3)=cg_zu(innLoop,outLoop)
        DO i=innLoop-1,1,-1
          cg_zu(i,outLoop)=cg_zu(i,outLoop)-                            &
     &                     cg_gamma(i+1,outLoop)*cg_zu(i+1,outLoop)
          cg_Tmatrix(i,3)=cg_zu(i,outLoop)
        END DO
!
!  Compute gradient norm using ad*(:,:,1) and tl_*(:,:,2) as temporary
!  storage.
!
        Linp=1
        Lout=2
        Lwrk=2
        CALL new_gradient (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Linp, Lout, Lwrk,                            &
     &                     innLoop, outLoop,                            &
#ifdef MASKING
     &                     rmask, umask, vmask,                         &
#endif
#ifdef ADJUST_WSTRESS
     &                     tl_ustr, tl_vstr,                            &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                     tl_tflux,                                    &
# endif
     &                     tl_t, tl_u, tl_v,                            &
#else
     &                     tl_ubar, tl_vbar,                            &
#endif
     &                     tl_zeta,                                     &
#ifdef ADJUST_WSTRESS
     &                     ad_ustr, ad_vstr,                            &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                     ad_tflux,                                    &
# endif
     &                     ad_t, ad_u, ad_v,                            &
#else
     &                     ad_ubar, ad_vbar,                            &
#endif
     &                     ad_zeta)
      END IF
!
!-----------------------------------------------------------------------
!  Determine the eigenvalues and eigenvectors of the tridiagonal matrix.
!  These will be used on the last inner-loop to compute the eigenvectors
!  of the Hessian.
!-----------------------------------------------------------------------
!
      IF (LhessianEV.and.(innLoop.gt.0)) THEN
        DO i=1,innLoop
          cg_Ritz(i,outLoop)=cg_delta(i,outLoop)
        END DO
        DO i=1,innLoop-1
          cg_Tmatrix(i,1)=cg_beta(i+1,outLoop)
        END DO
!
!  Use the LAPACK routine DSTEQR to compute the eigenvectors and
!  eigenvalues of the tridiagonal matrix. If applicable, the 
!  eigenpairs is computed by master thread only. Notice that on
!  exit, the matrix cg_Tmatrix is destroyed.
!
        IF (Master) THEN
          CALL DSTEQR ('I', innLoop, cg_Ritz(1,outLoop), cg_Tmatrix,    &
     &                 cg_zv, Ninner, work, info)
        END IF
#ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iTLM, info, 1)
#endif
        IF (info.ne.0) THEN
          PRINT *,'Error in DSTEQR: info=',info
          STOP
        END IF
#ifdef DISTRIBUTE
        CALL mp_bcastf (ng, iTLM, cg_Ritz(:,outLoop), Ninner)
        CALL mp_bcastf_m (ng, ITLM, cg_zv, Ninner, Ninner)
#endif
!
!  Estimate the Ritz value error bounds.
!
        DO i=1,innLoop
          cg_RitzErr(i,outLoop)=ABS(cg_beta(innLoop+1,outLoop)*         &
     &                              cg_zv(innLoop,i))
        END DO
!
!  Check for exploding or negative Ritz values.
!
        DO i=1,innLoop
          IF (cg_Ritz(i,outLoop).lt.0.0_r8) THEN
            PRINT *,'negative Ritz value found'
            STOP
          END IF
        END DO
!
!  Count the converged eigenvectors.
!
        ingood=0
        RitzMaxErr=GradErr*cg_Ritz(innLoop,outLoop)
        DO i=1,innLoop
          IF (cg_RitzErr(i,outLoop).le.RitzMaxErr) THEN
            ingood=ingood+1
          END IF
        END DO
!
!  Deal with newly converged eigenvector and save leading converged
!  eigenvector for explosion test.
!
        IF (ingood.gt.0) THEN
          DO i=innLoop,1,-1
            IF (cg_RitzErr(i,outLoop).le.RitzMaxErr) THEN
              ztheta1=cg_Ritz(i,outLoop)
              itheta1=i
              EXIT
            END IF
          END DO
        END IF
!
!  Calculate the converged eigenvectors of the Hessian.
!
        IF (innLoop.eq.Ninner-1) THEN
          RitzMaxErr=HevecErr
          DO i=1,innLoop
            cg_RitzErr(i,outLoop)=cg_RitzErr(i,outLoop)/                &
     &                            cg_Ritz(i,outLoop)
          END DO
          PRINT *, 'Entering hessian_evecs', MyRank
          Lwrk=2
          Linp=1
          Lout=2
          CALL hessian_evecs (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Linp, Lout, Lwrk,                         &
     &                        innLoop, outLoop,                         &
# ifdef MASKING
     &                        rmask, umask, vmask,                      &
# endif
# ifdef ADJUST_WSTRESS
     &                        tl_ustr, tl_vstr,                         &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                        tl_tflux,                                 &
#  endif
     &                        tl_t, tl_u, tl_v,                         &
# else
     &                        tl_ubar, tl_vbar,                         &
# endif
     &                        tl_zeta,                                  &
# ifdef ADJUST_WSTRESS
     &                        ad_ustr, ad_vstr,                         &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                        ad_tflux,                                 &
#  endif
     &                        ad_t, ad_u, ad_v,                         &
# else
     &                        ad_ubar, ad_vbar,                         &
# endif
     &                        ad_zeta)

          PRINT *, 'Exiting hessian_evecs', MyRank

          IF (Master.and.(nConvRitz.eq.0)) THEN
            PRINT *,' No converged Hesssian eigenvectors found.'
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Set TLM initial conditions for next inner loop, X(k+1).
!-----------------------------------------------------------------------
!
!   X(k+1) = tau(k+1) * d(k+1)
!
!   For the Lanczos algorithm, X(Linp) is ALWAYS the starting TL initial 
!   condition which for IS4DVAR is zero.
!
      Linp=1
      Lout=2
      CALL tl_new_state (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Linp, Lout,                                    &
     &                   innLoop, outLoop,                              &
     &                   cg_tau(innLoop,outLoop),                       &
#ifdef MASKING
     &                   rmask, umask, vmask,                           &
#endif
#ifdef ADJUST_WSTRESS
     &                   d_sustr, d_svstr,                              &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                   d_stflx,                                       &
# endif
     &                   d_t, d_u, d_v,                                 &
#else
     &                   d_ubar, d_vbar,                                &
#endif
     &                   d_zeta,                                        &
#ifdef ADJUST_WSTRESS
     &                   tl_ustr, tl_vstr,                              &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                   tl_tflux,                                      &
# endif
     &                   tl_t, tl_u, tl_v,                              &
#else
     &                   tl_ubar, tl_vbar,                              &
#endif
     &                   tl_zeta,                                       &
#ifdef ADJUST_WSTRESS
     &                   ad_ustr, ad_vstr,                              &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                   ad_tflux,                                      &
# endif
     &                   ad_t, ad_u, ad_v,                              &
#else
     &                   ad_ubar, ad_vbar,                              &
#endif
     &                   ad_zeta)
!
!-----------------------------------------------------------------------
!  Write out conjugate gradient information into NetCDF file.
!-----------------------------------------------------------------------
!
      CALL cg_write (ng, innLoop, outLoop)
!
!  Report algorithm parameters.
!
      IF (Master) THEN
        IF (inner.eq.0) THEN
          WRITE (stdout,20) outLoop, innLoop,                           &
     &                      cg_Gnorm(outLoop)
 20       FORMAT (1x,'(',i3.3,',',i3.3,'): ',                           &
     &            'Gnorm  = ',1p,e14.7)
        END IF
        IF (innLoop.gt.0) THEN
          WRITE (stdout,30) outLoop, innLoop,                           &
     &                      cg_Greduc(innLoop,outLoop),                 &
     &                      outLoop, innLoop,                           &
     &                      cg_delta(innLoop,outLoop)
 30       FORMAT (1x,'(',i3.3,',',i3.3,'): ',                           &
     &            'Greduc = ',1p,e14.7,/,                               &
     &            1x,'(',i3.3,',',i3.3,'): ',                           &
     &            'delta  = ',1p,e14.7)
          WRITE (stdout,40) RitzMaxErr
 40       FORMAT (/,'Ritz Eigenvalues and relative accuracy: ',         &
     &             'RitzMaxErr = ',1p,e14.7,/)
          ic=0
          DO i=1,innLoop
            IF (cg_RitzErr(i,outLoop).le.RitzMaxErr) THEN
              string='converged'
              ic=ic+1
              WRITE (stdout,50) i, cg_Ritz(i,outLoop),                  &
     &                          cg_RitzErr(i,outLoop),                  &
     &                          TRIM(ADJUSTL(string)), ic
 50           FORMAT(5x,i3.3,2x,1p,e14.7,2x,1p,e14.7,2x,a,2x,           &
     &               '(Good='i3.3,')')
            ELSE
              string='not converged'
              WRITE (stdout,60) i, cg_Ritz(i,outLoop),                  &
     &                          cg_RitzErr(i,outLoop),                  &
     &                          TRIM(ADJUSTL(string))
 60           FORMAT(5x,i3.3,2x,1p,e14.7,2x,1p,e14.7,2x,a)
            END IF
          END DO
        END IF
      END IF

      RETURN 
      END SUBROUTINE cgradient_tile

!
!***********************************************************************
      SUBROUTINE tl_new_state (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Linp, Lout,                              &
     &                         innLoop, outLoop,                        &
     &                         alphaK,                                  &
#ifdef MASKING
     &                         rmask, umask, vmask,                     &
#endif
#ifdef ADJUST_WSTRESS
     &                         d_sustr, d_svstr,                        &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                         d_stflx,                                 &
# endif
     &                         d_t, d_u, d_v,                           &
#else
     &                         d_ubar, d_vbar,                          &
#endif
     &                         d_zeta,                                  &
#ifdef ADJUST_WSTRESS
     &                         tl_ustr, tl_vstr,                        &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                         tl_tflux,                                &
# endif
     &                         tl_t, tl_u, tl_v,                        &
#else
     &                         tl_ubar, tl_vbar,                        &
#endif
     &                         tl_zeta,                                 &
#ifdef ADJUST_WSTRESS
     &                         ad_ustr, ad_vstr,                        &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                         ad_tflux,                                &
# endif
     &                         ad_t, ad_u, ad_v,                        &
#else
     &                         ad_ubar, ad_vbar,                        &
#endif
     &                         ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_fourdvar
      USE mod_ncparam
      USE mod_scalars
      USE mod_iounits
!
      USE state_addition_mod, ONLY : state_addition
      USE state_copy_mod, ONLY : state_copy
      USE state_initialize_mod, ONLY : state_initialize
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Linp, Lout
      integer, intent(in) :: innLoop, outLoop

      real(r8), intent(in) :: alphaK
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: d_sustr(LBi:,LBj:,:)
      real(r8), intent(in) :: d_svstr(LBi:,LBj:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(in) :: d_stflx(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(in) :: d_t(LBi:,LBj:,:,:)
      real(r8), intent(in) :: d_u(LBi:,LBj:,:)
      real(r8), intent(in) :: d_v(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: d_ubar(LBi:,LBj:)
      real(r8), intent(in) :: d_vbar(LBi:,LBj:)
# endif
      real(r8), intent(in) :: d_zeta(LBi:,LBj:)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: d_sustr(LBi:UBi,LBj:UBj,Nfrec(ng))
      real(r8), intent(inout) :: d_svstr(LBi:UBi,LBj:UBj,Nfrec(ng))
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: d_stflx(LBi:UBi,LBj:UBj,               &
     &                                   Nfrec(ng),NT(ng))
#  endif
      real(r8), intent(inout) :: d_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: d_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: d_v(LBi:UBi,LBj:UBj,N(ng))
# else
      real(r8), intent(inout) :: d_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: d_vbar(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(inout) :: d_zeta(LBi:UBi,LBj:UBj)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: i, j, k, lstr, rec
#ifdef SOLVE3D
      integer :: itrc
#endif

      real(r8) :: fac, fac1, fac2

      character (len=80) :: ncname

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute new starting tangent linear state vector, X(k+1).
!-----------------------------------------------------------------------
!
      IF (innLoop.ne.Ninner-1) THEN
!
!  Free-surface.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            tl_zeta(i,j,Lout)=alphaK*d_zeta(i,j)
#ifdef MASKING
            tl_zeta(i,j,Lout)=tl_zeta(i,j,Lout)*rmask(i,j)
#endif
          END DO
        END DO
#ifndef SOLVE3D
!
!  2D momentum.
!
        DO j=JstrR,JendR
          DO i=Istr,IendR
            tl_ubar(i,j,Lout)=alphaK*d_ubar(i,j)
# ifdef MASKING
            tl_ubar(i,j,Lout)=tl_ubar(i,j,Lout)*umask(i,j)
# endif
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            tl_vbar(i,j,Lout)=alphaK*d_vbar(i,j)
# ifdef MASKING
            tl_vbar(i,j,Lout)=tl_vbar(i,j,Lout)*vmask(i,j)
# endif
          END DO
        END DO
#endif
#ifdef ADJUST_WSTRESS
!
!  Surface momentum stress.
!
        DO k=1,Nfrec(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              tl_ustr(i,j,k,Lout)=alphaK*d_sustr(i,j,k)
# ifdef MASKING
              tl_ustr(i,j,k,Lout)=tl_ustr(i,j,k,Lout)*umask(i,j)
# endif
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              tl_vstr(i,j,k,Lout)=alphaK*d_svstr(i,j,k)
# ifdef MASKING
              tl_vstr(i,j,k,Lout)=tl_vstr(i,j,k,Lout)*vmask(i,j)
# endif
            END DO
          END DO
        END DO
#endif
#ifdef SOLVE3D
!
!  3D momentum.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              tl_u(i,j,k,Lout)=alphaK*d_u(i,j,k)
# ifdef MASKING
              tl_u(i,j,k,Lout)=tl_u(i,j,k,Lout)*umask(i,j)
# endif
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              tl_v(i,j,k,Lout)=alphaK*d_v(i,j,k)
# ifdef MASKING
              tl_v(i,j,k,Lout)=tl_v(i,j,k,Lout)*vmask(i,j)
# endif
            END DO
          END DO
        END DO
!
!  Tracers.
!
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                tl_t(i,j,k,Lout,itrc)=alphaK*d_t(i,j,k,itrc)
# ifdef MASKING
                tl_t(i,j,k,Lout,itrc)=tl_t(i,j,k,Lout,itrc)*rmask(i,j)
# endif
              END DO
            END DO
          END DO
        END DO
# ifdef ADJUST_STFLUX
!
!  Surface tracers flux.
!
        DO itrc=1,NT(ng)
          DO k=1,Nfrec(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                tl_tflux(i,j,k,Lout,itrc)=alphaK*d_stflx(i,j,k,itrc)
#  ifdef MASKING
                tl_tflux(i,j,k,Lout,itrc)=tl_tflux(i,j,k,Lout,itrc)*    &
     &                                    rmask(i,j)
#  endif
              END DO
            END DO
          END DO
        END DO
# endif
#endif
!
!-----------------------------------------------------------------------
!  If last inner-loop, compute the tangent linear model initial
!  conditions from the Lanczos algorithm.  Compute the actual final
!  value of the cost function. Use adjoint state arrays, index Linp,
!  as temporary storage.
!-----------------------------------------------------------------------
!
      ELSE
!
!  Clear the adjoint working arrays (index Linp) since the tangent
!  linear model initial condition on the first inner-loop is zero:
!
!    ad_var(Linp) = fac
!
        fac=0.0_r8

        CALL state_initialize (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Linp, fac,                               &
#ifdef MASKING
     &                         rmask, umask, vmask,                     &
#endif
#ifdef ADJUST_WSTRESS
     &                         ad_ustr, ad_vstr,                        &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                         ad_tflux,                                &
# endif
     &                         ad_t, ad_u, ad_v,                        &
#else
     &                         ad_ubar, ad_vbar,                        &
#endif
     &                         ad_zeta)
!
!  Read in each previous gradient state solutions, g(0) to g(k), and
!  compute its associated dot against current g(k+1).
!
        DO rec=1,innLoop
!
!  Determine adjoint file to process.
!
          IF (ndefADJ(ng).gt.0) THEN
            lstr=LEN_TRIM(ADJbase(ng))
            WRITE (ncname,10) ADJbase(ng)(1:lstr-3), rec
 10         FORMAT (a,'_',i3.3,'.nc')
          ELSE
            ncname=ADJname(ng)
          END IF
!
!  Read gradient solution and load it into TANGENT LINEAR STATE ARRAYS
!  at index Lout.
!
          CALL read_state (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Lout, rec,                                   &
     &                     ndefADJ(ng), ncADJid(ng), ncname,            &
#ifdef MASKING
     &                     rmask, umask, vmask,                         &
#endif
#ifdef ADJUST_WSTRESS
     &                     tl_ustr, tl_vstr,                            &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                     tl_tflux,                                    &
# endif
     &                     tl_t, tl_u, tl_v,                            &
#else
     &                     tl_ubar, tl_vbar,                            &
#endif
     &                     tl_zeta)
!
!  Sum all previous normalized gradients:
!
!    ad_var(Linp) = fac1 * ad_var(Linp) + fac2 * tl_var(Lout)
!     
          fac1=1.0_r8
          fac2=cg_zu(rec,outLoop)

          CALL state_addition (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Linp, Lout, Linp, fac1, fac2,            &
#ifdef MASKING
     &                         rmask, umask, vmask,                     &
#endif
#ifdef ADJUST_WSTRESS
     &                         ad_ustr, tl_ustr,                        &
     &                         ad_vstr, tl_vstr,                        &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                         ad_tflux, tl_tflux,                      &
# endif
     &                         ad_t, tl_t,                              &
     &                         ad_u, tl_u,                              &
     &                         ad_v, tl_v,                              &
#else
     &                         ad_ubar, tl_ubar,                        &
     &                         ad_vbar, tl_vbar,                        &
#endif
     &                         ad_zeta, tl_zeta)
        END DO
!
!  Load new tangent linear model initial conditions to respective state
!  arrays, index Lout:
!
!    tl_var(Lout) = ad_var(Linp)
!
        CALL state_copy (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Linp, Lout,                                    &
#ifdef ADJUST_WSTRESS
     &                   tl_ustr, ad_ustr,                              &
     &                   tl_vstr, ad_vstr,                              &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                   tl_tflux, ad_tflux,                            &
# endif
     &                   tl_t, ad_t,                                    &
     &                   tl_u, ad_u,                                    &
     &                   tl_v, ad_v,                                    &
#else
     &                   tl_ubar, ad_ubar,                              &
     &                   tl_vbar, ad_vbar,                              &
#endif
     &                   tl_zeta, ad_zeta)
      END IF

      RETURN
      END SUBROUTINE tl_new_state
!
!***********************************************************************
      SUBROUTINE read_state (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Lwrk, rec,                                 &
     &                       ndef, ncfileid, ncname,                    &
#ifdef MASKING
     &                       rmask, umask, vmask,                       &
#endif
#ifdef ADJUST_WSTRESS
     &                       s_ustr, s_vstr,                            &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                       s_tflux,                                   &
# endif
     &                       s_t, s_u, s_v,                             &
#else
     &                       s_ubar, s_vbar,                            &
#endif
     &                       s_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lwrk, rec, ndef

      integer, intent(inout) :: ncfileid

      character (len=*), intent(in) :: ncname
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: s_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: s_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: s_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: s_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: s_zeta(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: s_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s_tflux(LBi:UBi,LBj:UBj,               &
     &                                   Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: s_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: s_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: s_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: s_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: s_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: s_zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: i, j, k
#ifdef SOLVE3D
      integer :: itrc
#endif
      integer :: gtype, ncid, status
      integer, dimension(NV) :: vid
      integer, dimension(4) :: Vsize

      integer :: nf_fread2d
#ifdef SOLVE3D
      integer :: nf_fread3d
#endif

      real(r8) :: Fmin, Fmax, scale

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Read in requested model state record. Load data into state array
!  index Lwrk.
!-----------------------------------------------------------------------
!
!  Determine file and variables ids.
!
      IF (ndef.gt.0) THEN
        IF (InpThread) THEN
          status=nf90_open(TRIM(ncname), nf90_nowrite, ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(ncname)
            exit_flag=2
            ioerror=status
          END IF            
          ncfileid=ncid
        END IF
#ifdef DISTRIBUTE
        CALL mp_bcasti (ng, model, ncfileid, 1)
        CALL mp_bcasti (ng, model, exit_flag, 1)
        IF (exit_flag.ne.NoError) RETURN
#endif
      ELSE
        ncid=ncfileid
      END IF
      IF (InpThread) THEN
#ifndef SOLVE3D
        status=nf90_inq_varid(ncid, TRIM(Vname(1,idUbar)), vid(idUbar))
        status=nf90_inq_varid(ncid, TRIM(Vname(1,idVbar)), vid(idVbar))
#endif
        status=nf90_inq_varid(ncid, TRIM(Vname(1,idFsur)), vid(idFsur))
#ifdef ADJUST_WSTRESS
        status=nf90_inq_varid(ncid, TRIM(Vname(1,idUsms)), vid(idUsms))
        status=nf90_inq_varid(ncid, TRIM(Vname(1,idVsms)), vid(idVsms))
#endif
#ifdef SOLVE3D
        status=nf90_inq_varid(ncid, TRIM(Vname(1,idUvel)), vid(idUvel))
        status=nf90_inq_varid(ncid, TRIM(Vname(1,idVvel)), vid(idVvel))
        DO itrc=1,NT(ng)
          status=nf90_inq_varid(ncid, TRIM(Vname(1,idTvar(itrc))),      &
     &                          vid(idTvar(itrc)))
# ifdef ADJUST_STFLUX
          status=nf90_inq_varid(ncid, TRIM(Vname(1,idTsur(itrc))),      &
     &                          vid(idTsur(itrc)))
# endif
        END DO
#endif
      END IF
      DO i=1,4
        Vsize(i)=0
      END DO
      scale=1.0_r8
!
!  Read in free-surface.
!
      gtype=r2dvar
      status=nf_fread2d(ng, iTLM, ncid, vid(idFsur), rec, gtype,        &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  scale, Fmin, Fmax,                              &
#ifdef MASKING
     &                  rmask(LBi,LBj),                                 &
#endif
     &                  s_zeta(LBi,LBj,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idFsur)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#ifndef SOLVE3D
!
!  Read in 2D momentum.
!
      gtype=u2dvar
      status=nf_fread2d(ng, iTLM, ncid, vid(idUbar), rec, gtype,        &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  umask(LBi,LBj),                                 &
# endif
     &                  s_ubar(LBi,LBj,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUbar)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

      gtype=v2dvar
      status=nf_fread2d(ng, iTLM, ncid, vid(idVbar), rec, gtype,        &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  vmask(LBi,LBj),                                 &
# endif
     &                  s_vbar(LBi,LBj,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVbar)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#endif

#ifdef ADJUST_WSTRESS
!
!  Read surface momentum stress.
!
      gtype=u3dvar
      status=nf_fread3d(ng, iTLM, ncid, vid(idUsms), rec, gtype,        &
     &                  Vsize, LBi, UBi, LBj, UBj, 1, Nfrec(ng),        &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  umask(LBi,LBj),                                 &
# endif
     &                  s_ustr(LBi,LBj,1,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUsms)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

      gtype=v3dvar
      status=nf_fread3d(ng, iTLM, ncid, vid(idVsms), rec, gtype,        &
     &                  Vsize, LBi, UBi, LBj, UBj, 1, Nfrec(ng),        &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  vmask(LBi,LBj),                                 &
# endif
     &                  s_vstr(LBi,LBj,1,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVsms)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#endif

#ifdef SOLVE3D
!
!  Read in 3D momentum.
!
      gtype=u3dvar
      status=nf_fread3d(ng, iTLM, ncid, vid(idUvel), rec, gtype,        &
     &                  Vsize, LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  umask(LBi,LBj),                                 &
# endif
     &                  s_u(LBi,LBj,1,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUvel)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

      gtype=v3dvar
      status=nf_fread3d(ng, iTLM, ncid, vid(idVvel), rec, gtype,        &
     &                  Vsize, LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  vmask(LBi,LBj),                                 &
# endif
     &                  s_v(LBi,LBj,1,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVvel)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Read in tracers.
!
      gtype=r3dvar
      DO itrc=1,NT(ng)
        status=nf_fread3d(ng, iTLM, ncid, vid(idTvar(itrc)), rec,       &
     &                    gtype, Vsize, LBi, UBi, LBj, UBj, 1, N(ng),   &
     &                    scale, Fmin, Fmax,                            &
# ifdef MASKING
     &                    rmask(LBi,LBj),                               &
# endif
     &                    s_t(LBi,LBj,1,Lwrk,itrc))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTvar(itrc))), rec,         &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO

# ifdef ADJUST_STFLUX
!
!  Read in surface tracers flux.
!
      gtype=r3dvar
      DO itrc=1,NT(ng)
        status=nf_fread3d(ng, iTLM, ncid, vid(idTsur(itrc)), rec,       &
     &                    gtype, Vsize, LBi, UBi, LBj, UBj, 1,Nfrec(ng),&
     &                    scale, Fmin, Fmax,                            &
#  ifdef MASKING
     &                    rmask(LBi,LBj),                               &
#  endif
     &                    s_tflux(LBi,LBj,1,Lwrk,itrc))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTsur(itrc))), rec,         &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
# endif
#endif
!
!  If multiple files, close current file.
!
      IF (InpThread.and.(ndef.gt.0)) THEN
        status=nf90_close(ncid)
      END IF
!
 10   FORMAT (' READ_STATE - unable to open NetCDF file: ',a)
 20   FORMAT (' READ_STATE - error while reading variable: ',a,2x,      &
     &        'at time record = ',i3,/,14x,'in NetCDF file: ',a)

      RETURN
      END SUBROUTINE read_state

!
!***********************************************************************
      SUBROUTINE new_direction (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Lold, Lnew,                             &
#ifdef MASKING
     &                          rmask, umask, vmask,                    &
#endif
#ifdef ADJUST_WSTRESS
     &                          ad_ustr, ad_vstr,                       &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                          ad_tflux,                               &
# endif
     &                          ad_t, ad_u, ad_v,                       &
#else
     &                          ad_ubar, ad_vbar,                       &
#endif
     &                          ad_zeta,                                &
#ifdef ADJUST_WSTRESS
     &                          d_sustr, d_svstr,                       &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                          d_stflx,                                &
# endif
     &                          d_t, d_u, d_v,                          &
#else
     &                          d_ubar, d_vbar,                         &
#endif
     &                          d_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_scalars
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lold, Lnew
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: d_sustr(LBi:,LBj:,:)
      real(r8), intent(inout) :: d_svstr(LBi:,LBj:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: d_stflx(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(inout) :: d_t(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: d_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: d_v(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: d_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: d_vbar(LBi:,LBj:)
# endif
      real(r8), intent(inout) :: d_zeta(LBi:,LBj:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBI,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: d_sustr(LBi:UBi,LBj:UBj,Nfrec(ng))
      real(r8), intent(inout) :: d_svstr(LBi:UBI,LBj:UBj,Nfrec(ng))
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: d_stflx(LBi:UBi,LBj:UBj,               &
     &                                   Nfrec(ng),NT(ng))
#  endif
      real(r8), intent(inout) :: d_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: d_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: d_v(LBi:UBi,LBj:UBj,N(ng))
# else
      real(r8), intent(inout) :: d_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: d_vbar(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(inout) :: d_zeta(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: i, j, k
#ifdef SOLVE3D
      integer :: itrc
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute new conjugate descent direction, d(k+1). Notice that the old
!  descent direction is overwritten.
!-----------------------------------------------------------------------
!
!  Free-sruface.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          d_zeta(i,j)=ad_zeta(i,j,Lnew)
#ifdef MASKING
          d_zeta(i,j)=d_zeta(i,j)*rmask(i,j)
#endif
        END DO
      END DO
#ifndef SOLVE3D
!
!  2D momentum.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          d_ubar(i,j)=ad_ubar(i,j,Lnew)
# ifdef MASKING
          d_ubar(i,j)=d_ubar(i,j)*umask(i,j)
# endif
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          d_vbar(i,j)=ad_vbar(i,j,Lnew)
# ifdef MASKING
          d_vbar(i,j)=d_vbar(i,j)*vmask(i,j)
# endif
        END DO
      END DO
#endif
#ifdef ADJUST_WSTRESS
!
!  Surface momentum stress.
!
      DO k=1,Nfrec(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            d_sustr(i,j,k)=ad_ustr(i,j,k,Lnew)
# ifdef MASKING
            d_sustr(i,j,k)=d_sustr(i,j,k)*umask(i,j)
# endif
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            d_svstr(i,j,k)=ad_vstr(i,j,k,Lnew)
# ifdef MASKING
            d_svstr(i,j,k)=d_svstr(i,j,k)*vmask(i,j)
# endif
          END DO
        END DO
      END DO
#endif
#ifdef SOLVE3D
!
!  3D momentum.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            d_u(i,j,k)=ad_u(i,j,k,Lnew)
# ifdef MASKING
            d_u(i,j,k)=d_u(i,j,k)*umask(i,j)
# endif
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            d_v(i,j,k)=ad_v(i,j,k,Lnew)
# ifdef MASKING
            d_v(i,j,k)=d_v(i,j,k)*vmask(i,j)
# endif
          END DO
        END DO
      END DO
!
!  Tracers.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              d_t(i,j,k,itrc)=ad_t(i,j,k,Lnew,itrc)
# ifdef MASKING
              d_t(i,j,k,itrc)=d_t(i,j,k,itrc)*rmask(i,j)
# endif
            END DO
          END DO
        END DO
      END DO
# ifdef ADJUST_STFLUX
!
!  Surface tracers flux.
!
      DO itrc=1,NT(ng)
        DO k=1,Nfrec(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              d_stflx(i,j,k,itrc)=ad_tflux(i,j,k,Lnew,itrc)
#  ifdef MASKING
              d_stflx(i,j,k,itrc)=d_stflx(i,j,k,itrc)*rmask(i,j)
#  endif
            END DO
          END DO
        END DO
      END DO
# endif
#endif

      RETURN
      END SUBROUTINE new_direction
!
!***********************************************************************
      SUBROUTINE hessian (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Lold, Lnew, Lwrk,                             &
     &                    innLoop, outLoop,                             &
#ifdef MASKING
     &                    rmask, umask, vmask,                          &
#endif
#ifdef ADJUST_WSTRESS
     &                    ad_ustr, ad_vstr,                             &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                    ad_tflux,                                     &
# endif
     &                    ad_t, ad_u, ad_v,                             &
#else
     &                    ad_ubar, ad_vbar,                             &
#endif
     &                    ad_zeta,                                      &
#ifdef ADJUST_WSTRESS
     &                    tl_ustr, tl_vstr,                             &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                    tl_tflux,                                     &
# endif
     &                    tl_t, tl_u, tl_v,                             &
#else
     &                    tl_ubar, tl_vbar,                             &
#endif
     &                    tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE state_dotprod_mod, ONLY : state_dotprod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lold, Lnew, Lwrk
      integer, intent(in) :: innLoop, outLoop
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: i, j, k, lstr
#ifdef SOLVE3D
      integer :: itrc
#endif
      real(r8) :: fac

      real(r8), dimension(0:NstateVar(ng)) :: dot

      character (len=80) :: ncname

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Estimate the Hessian.
!-----------------------------------------------------------------------
!
!  Need to multiply the adjoint state arrays (index Lold) by zgnorm to
!  convert back to the non-normalized gradient. Here, the tangent linear
!  state arrays (index Lold) contain the background cost function
!  gradient.
!
      fac=1.0_r8/cg_tau(innLoop,outLoop)
!
!  Free-surface.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ad_zeta(i,j,Lnew)=fac*(ad_zeta(i,j,Lnew)+                     &
     &                           tl_zeta(i,j,Lold)-                     &
     &                           ad_zeta(i,j,Lold)*                     &
     &                           cg_Gnorm(outLoop))
#ifdef MASKING
          ad_zeta(i,j,Lnew)=ad_zeta(i,j,Lnew)*rmask(i,j)
#endif
        END DO
      END DO
#ifndef SOLVE3D
!
!  2D momentum.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ad_ubar(i,j,Lnew)=fac*(ad_ubar(i,j,Lnew)+                     &
     &                           tl_ubar(i,j,Lold)-                     &
     &                           ad_ubar(i,j,Lold)*                     &
     &                           cg_Gnorm(outLoop))
# ifdef MASKING
          ad_ubar(i,j,Lnew)=ad_ubar(i,j,Lnew)*umask(i,j)
# endif
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          ad_vbar(i,j,Lnew)=fac*(ad_vbar(i,j,Lnew)+                     &
     &                           tl_vbar(i,j,Lold)-                     &
     &                           ad_vbar(i,j,Lold)*                     &
     &                           cg_Gnorm(outLoop))
# ifdef MASKING
          ad_vbar(i,j,Lnew)=ad_vbar(i,j,Lnew)*vmask(i,j)
# endif
        END DO
      END DO
#endif
#ifdef ADJUST_WSTRESS
!
!  Surface momentum stress.
!
      DO k=1,Nfrec(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_ustr(i,j,k,Lnew)=fac*(ad_ustr(i,j,k,Lnew)+               &
     &                               tl_ustr(i,j,k,Lold)-               &
     &                               ad_ustr(i,j,k,Lold)*               &
     &                               cg_Gnorm(outLoop))
# ifdef MASKING
            ad_ustr(i,j,k,Lnew)=ad_ustr(i,j,k,Lnew)*umask(i,j)
# endif
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_vstr(i,j,k,Lnew)=fac*(ad_vstr(i,j,k,Lnew)+               &
     &                               tl_vstr(i,j,k,Lold)-               &
     &                               ad_vstr(i,j,k,Lold)*               &
     &                               cg_Gnorm(outLoop))
# ifdef MASKING
            ad_vstr(i,j,k,Lnew)=ad_vstr(i,j,k,Lnew)*vmask(i,j)
# endif
          END DO
        END DO
      END DO
#endif
#ifdef SOLVE3D
!
!  3D momentum.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_u(i,j,k,Lnew)=fac*(ad_u(i,j,k,Lnew)+                     &
     &                            tl_u(i,j,k,Lold)-                     &
     &                            ad_u(i,j,k,Lold)*                     &
     &                            cg_Gnorm(outLoop))
# ifdef MASKING
            ad_u(i,j,k,Lnew)=ad_u(i,j,k,Lnew)*umask(i,j)
# endif
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_v(i,j,k,Lnew)=fac*(ad_v(i,j,k,Lnew)+                     &
     &                            tl_v(i,j,k,Lold)-                     &
     &                            ad_v(i,j,k,Lold)*                     &
     &                            cg_Gnorm(outLoop))
# ifdef MASKING
            ad_v(i,j,k,Lnew)=ad_v(i,j,k,Lnew)*vmask(i,j)
# endif
          END DO
        END DO
      END DO
!
!  Tracers.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              ad_t(i,j,k,Lnew,itrc)=fac*(ad_t(i,j,k,Lnew,itrc)+         &
     &                                   tl_t(i,j,k,Lold,itrc)-         &
     &                                   ad_t(i,j,k,Lold,itrc)*         &
     &                                   cg_Gnorm(outLoop))
# ifdef MASKING
              ad_t(i,j,k,Lnew,itrc)=ad_t(i,j,k,Lnew,itrc)*rmask(i,j)
# endif
            END DO
          END DO
        END DO
      END DO

# ifdef ADJUST_STFLUX
!
!  Surface tracers flux.
!
      DO itrc=1,NT(ng)
        DO k=1,Nfrec(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              ad_tflux(i,j,k,Lnew,itrc)=fac*(ad_tflux(i,j,k,Lnew,itrc)+ &
     &                                       tl_tflux(i,j,k,Lold,itrc)- &
     &                                       ad_tflux(i,j,k,Lold,itrc)* &
     &                                       cg_Gnorm(outLoop))
#  ifdef MASKING
              ad_tflux(i,j,k,Lnew,itrc)=ad_tflux(i,j,k,Lnew,itrc)*      &
     &                                  rmask(i,j)
#  endif
            END DO
          END DO
        END DO
      END DO
# endif
#endif
!
!-----------------------------------------------------------------------
!  Compute norm Delta(k) as the dot product between the new gradient
!  and current iteration gradient solution.
!-----------------------------------------------------------------------
!
!  Determine gradient file to process.
!
      IF (ndefADJ(ng).gt.0) THEN
        lstr=LEN_TRIM(ADJbase(ng))
        WRITE (ncname,10) ADJbase(ng)(1:lstr-3), innLoop
 10     FORMAT (a,'_',i3.3,'.nc')
      ELSE
        ncname=ADJname(ng)
      END IF
!
!  Read current gradient solution into tangent linear state array,
!  index Lwrk.
!
      CALL read_state (ng, tile, model,                                 &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 Lwrk, innLoop,                                   &
     &                 ndefADJ(ng), ncADJid(ng), ncname,                &
#ifdef MASKING
     &                 rmask, umask, vmask,                             &
#endif
#ifdef ADJUST_WSTRESS
     &                 tl_ustr, tl_vstr,                                &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                 tl_tflux,                                        &
# endif
     &                 tl_t, tl_u, tl_v,                                &
#else
     &                 tl_ubar, tl_vbar,                                &
#endif
     &                 tl_zeta)
!
!  Compute current iteration norm Delta(k) used to compute tri-diagonal
!  matrix T(k) in the Lanczos recurrence.
!
      CALL state_dotprod (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NstateVar(ng), dot(0:),                       &
#ifdef MASKING
     &                    rmask, umask, vmask,                          &
#endif
#ifdef ADJUST_WSTRESS
     &                    ad_ustr(:,:,:,Lnew), tl_ustr(:,:,:,Lwrk),     &
     &                    ad_vstr(:,:,:,Lnew), tl_vstr(:,:,:,Lwrk),     &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                    ad_tflux(:,:,:,Lnew,:),                       &
     &                    tl_tflux(:,:,:,Lwrk,:),                       &
# endif
     &                    ad_t(:,:,:,Lnew,:), tl_t(:,:,:,Lwrk,:),       &
     &                    ad_u(:,:,:,Lnew), tl_u(:,:,:,Lwrk),           &
     &                    ad_v(:,:,:,Lnew), tl_v(:,:,:,Lwrk),           &
#else
     &                    ad_ubar(:,:,Lnew), tl_ubar(:,:,Lwrk),         &
     &                    ad_vbar(:,:,Lnew), tl_vbar(:,:,Lwrk),         &
#endif
     &                    ad_zeta(:,:,Lnew), tl_zeta(:,:,Lwrk))

      cg_delta(innLoop,outLoop)=dot(0)

      RETURN
      END SUBROUTINE hessian
!
!***********************************************************************
      SUBROUTINE lanczos (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Lold, Lnew, Lwrk,                             &
     &                    innLoop, outLoop,                             &
#ifdef MASKING
     &                    rmask, umask, vmask,                          &
#endif
#ifdef ADJUST_WSTRESS
     &                    tl_ustr, tl_vstr,                             &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                    tl_tflux,                                     &
# endif
     &                    tl_t, tl_u, tl_v,                             &
#else
     &                    tl_ubar, tl_vbar,                             &
#endif
     &                    tl_zeta,                                      &
#ifdef ADJUST_WSTRESS
     &                    ad_ustr, ad_vstr,                             &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                    ad_tflux,                                     &
# endif
     &                    ad_t, ad_u, ad_v,                             &
#else
     &                    ad_ubar, ad_vbar,                             &
#endif
     &                    ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE state_addition_mod, ONLY : state_addition
      USE state_dotprod_mod, ONLY : state_dotprod
      USE state_scale_mod, ONLY : state_scale
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lold, Lnew, Lwrk
      integer, intent(in) :: innLoop, outLoop
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: i, j, lstr, rec
#ifdef SOLVE3D
      integer :: itrc, k
#endif
      real(r8) :: fac, fac1, fac2

      real(r8), dimension(0:NstateVar(ng)) :: dot
      real(r8), dimension(0:Ninner) :: DotProd, dot_new, dot_old

      character (len=80) :: ncname

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Calculate the new Lanczos vector, q(k+1) using reccurence equation
!  for the gradient vectors:
!
!     H q(k+1) = Gamma(k+1) q(k+2) + Delta(k+1) q(k+1) + Gamma(k) q(k)
!
!  where  Gamma(k) = - SQRT ( Beta(k+1) ) / Alpha(k)
!-----------------------------------------------------------------------
!
!  At this point, the previous orthonormal Lanczos vector is still in
!  tangent linear state arrays (index Lwrk).
!
      IF (innLoop.gt.0) THEN
!
!  Compute new Lanczos vector:
!
!    ad_var(Lnew) = fac1 * ad_var(Lnew) + fac2 * tl_var(Lwrk)
!     
        fac1=1.0_r8
        fac2=-cg_delta(innLoop,outLoop)

        CALL state_addition (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Lnew, Lwrk, Lnew, fac1, fac2,              &
#ifdef MASKING
     &                       rmask, umask, vmask,                       &
#endif
#ifdef ADJUST_WSTRESS
     &                       ad_ustr, tl_ustr,                          &
     &                       ad_vstr, tl_vstr,                          &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                       ad_tflux, tl_tflux,                        &
# endif
     &                       ad_t, tl_t,                                &
     &                       ad_u, tl_u,                                &
     &                       ad_v, tl_v,                                &
#else
     &                       ad_ubar, tl_ubar,                          &
     &                       ad_vbar, tl_vbar,                          &
#endif
     &                       ad_zeta, tl_zeta)
      END IF
!
!  Substract previous orthonormal Lanczos vector.
!
      IF (innLoop.gt.1) THEN
!
!  Determine adjoint file to process.
!
        IF (ndefADJ(ng).gt.0) THEN
          lstr=LEN_TRIM(ADJbase(ng))
          WRITE (ncname,10) ADJbase(ng)(1:lstr-3), innLoop-1
 10       FORMAT (a,'_',i3.3,'.nc')
        ELSE
          ncname=ADJname(ng)
        END IF
!
!  Read in the previous (innLoop-1) orthonormal Lanczos vector.
!
        CALL read_state (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Lwrk, innLoop-1,                               &
     &                   ndefADJ(ng), ncADJid(ng), ncname,              &
#ifdef MASKING
     &                   rmask, umask, vmask,                           &
#endif
#ifdef ADJUST_WSTRESS
     &                   tl_ustr, tl_vstr,                              &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                   tl_tflux,                                      &
# endif
     &                   tl_t, tl_u, tl_v,                              &
#else
     &                   tl_ubar, tl_vbar,                              &
#endif
     &                   tl_zeta)
!
!  Substract previous orthonormal Lanczos vector:
!
!    ad_var(Lnew) = fac1 * ad_var(Lnew) + fac2 * tl_var(Lwrk)
!     
        fac1=1.0_r8
        fac2=-cg_beta(innLoop,outLoop)

        CALL state_addition (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Lnew, Lwrk, Lnew, fac1, fac2,              &
#ifdef MASKING
     &                       rmask, umask, vmask,                       &
#endif
#ifdef ADJUST_WSTRESS
     &                       ad_ustr, tl_ustr,                          &
     &                       ad_vstr, tl_vstr,                          &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                       ad_tflux, tl_tflux,                        &
# endif
     &                       ad_t, tl_t,                                &
     &                       ad_u, tl_u,                                &
     &                       ad_v, tl_v,                                &
#else
     &                       ad_ubar, tl_ubar,                          &
     &                       ad_vbar, tl_vbar,                          &
#endif
     &                       ad_zeta, tl_zeta)
      END IF
!
!-----------------------------------------------------------------------
!  Orthogonalize current gradient, q(k+1), against all previous
!  gradients (reverse order) using Gramm-Schmidt procedure.
!-----------------------------------------------------------------------
!
!  We can overwrite adjoint arrays at index Lnew each time around the
!  the following loop because the preceding gradient vectors that we
!  read are orthogonal to each other. The reversed order of the loop
!  is important for the Lanczos vector calculations.
!
      DO rec=innLoop,1,-1
!
!  Determine adjoint file to process.
!
        IF (ndefADJ(ng).gt.0) THEN
          lstr=LEN_TRIM(ADJbase(ng))
          WRITE (ncname,10) ADJbase(ng)(1:lstr-3), rec
        ELSE
          ncname=ADJname(ng)
        END IF
!
!  Read in each previous gradient state solutions, G(0) to G(k), and
!  compute its associated dot angaint curret G(k+1). Each gradient
!  solution is loaded into TANGENT LINEAR STATE ARRAYS at index Lwrk.
!
        CALL read_state (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Lwrk, rec,                                     &
     &                   ndefADJ(ng), ncADJid(ng), ncname,              &
#ifdef MASKING
     &                   rmask, umask, vmask,                           &
#endif
#ifdef ADJUST_WSTRESS
     &                   tl_ustr, tl_vstr,                              &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                   tl_tflux,                                      &
# endif
     &                   tl_t, tl_u, tl_v,                              &
#else
     &                   tl_ubar, tl_vbar,                              &
#endif
     &                   tl_zeta)
!
!  Compute dot product <q(k+1), q(rec)>.
!
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), dot(0:),                     &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef ADJUST_WSTRESS
     &                      ad_ustr(:,:,:,Lnew), tl_ustr(:,:,:,Lwrk),   &
     &                      ad_vstr(:,:,:,Lnew), tl_vstr(:,:,:,Lwrk),   &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      ad_tflux(:,:,:,Lnew,:),                     &
     &                      tl_tflux(:,:,:,Lwrk,:),                     &
# endif
     &                      ad_t(:,:,:,Lnew,:), tl_t(:,:,:,Lwrk,:),     &
     &                      ad_u(:,:,:,Lnew), tl_u(:,:,:,Lwrk),         &
     &                      ad_v(:,:,:,Lnew), tl_v(:,:,:,Lwrk),         &
#else
     &                      ad_ubar(:,:,Lnew), tl_ubar(:,:,Lwrk),       &
     &                      ad_vbar(:,:,Lnew), tl_vbar(:,:,Lwrk),       &
#endif
     &                      ad_zeta(:,:,Lnew), tl_zeta(:,:,Lwrk))
!
!  Compute Gramm-Schmidt scaling coefficient.
!
        DotProd(rec)=dot(0)
!
!  Gramm-Schmidt orthonormalization, free-surface.
!
!    ad_var(Lnew) = fac1 * ad_var(Lnew) + fac2 * tl_var(Lwrk)
!     
        fac1=1.0_r8
        fac2=-DotProd(rec)

        CALL state_addition (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Lnew, Lwrk, Lnew, fac1, fac2,              &
#ifdef MASKING
     &                       rmask, umask, vmask,                       &
#endif
#ifdef ADJUST_WSTRESS
     &                       ad_ustr, tl_ustr,                          &
     &                       ad_vstr, tl_vstr,                          &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                       ad_tflux, tl_tflux,                        &
# endif
     &                       ad_t, tl_t,                                &
     &                       ad_u, tl_u,                                &
     &                       ad_v, tl_v,                                &
#else
     &                       ad_ubar, tl_ubar,                          &
     &                       ad_vbar, tl_vbar,                          &
#endif
     &                       ad_zeta, tl_zeta)
      END DO
!
!-----------------------------------------------------------------------
!  Normalize current orthogonal gradient vector.
!-----------------------------------------------------------------------
!
      CALL state_dotprod (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NstateVar(ng), dot(0:),                       &
#ifdef MASKING
     &                    rmask, umask, vmask,                          &
#endif
#ifdef ADJUST_WSTRESS
     &                    ad_ustr(:,:,:,Lnew), ad_ustr(:,:,:,Lnew),     &
     &                    ad_vstr(:,:,:,Lnew), ad_vstr(:,:,:,Lnew),     &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                    ad_tflux(:,:,:,Lnew,:),                       &
     &                    ad_tflux(:,:,:,Lnew,:),                       &
# endif
     &                    ad_t(:,:,:,Lnew,:), ad_t(:,:,:,Lnew,:),       &
     &                    ad_u(:,:,:,Lnew), ad_u(:,:,:,Lnew),           &
     &                    ad_v(:,:,:,Lnew), ad_v(:,:,:,Lnew),           &
#else
     &                    ad_ubar(:,:,Lnew), ad_ubar(:,:,Lnew),         &
     &                    ad_vbar(:,:,Lnew), ad_vbar(:,:,Lnew),         &
#endif
     &                    ad_zeta(:,:,Lnew), ad_zeta(:,:,Lnew))
!
!  Compute normalization factor.
!
      IF (innLoop.eq.0) THEN
        cg_Gnorm(outLoop)=SQRT(dot(0))
      ELSE
        cg_beta(innLoop+1,outLoop)=SQRT(dot(0))
      END IF
!
!  Normalize gradient: ad_var(Lnew) = fac * ad_var(Lnew)
!
      fac=1.0_r8/SQRT(dot(0))

      CALL state_scale (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Lnew, Lnew, fac,                                &
#ifdef MASKING
     &                  rmask, umask, vmask,                            &
#endif
#ifdef ADJUST_WSTRESS
     &                  ad_ustr, ad_vstr,                               &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                  ad_tflux,                                       &
# endif
     &                  ad_t, ad_u, ad_v,                               &
#else
     &                  ad_ubar, ad_vbar,                               &
#endif
     &                  ad_zeta)
!
!-----------------------------------------------------------------------
!  Compute dot product of new Lanczos vector with gradient.
!-----------------------------------------------------------------------
!
      IF (innLoop.eq.0) THEN
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), dot(0:),                     &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef ADJUST_WSTRESS
     &                      ad_ustr(:,:,:,Lnew), ad_ustr(:,:,:,Lnew),   &
     &                      ad_vstr(:,:,:,Lnew), ad_vstr(:,:,:,Lnew),   &
#endif 
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      ad_tflux(:,:,:,Lnew,:),                     &
     &                      ad_tflux(:,:,:,Lnew,:),                     &
# endif 
     &                      ad_t(:,:,:,Lnew,:), ad_t(:,:,:,Lnew,:),     &
     &                      ad_u(:,:,:,Lnew), ad_u(:,:,:,Lnew),         &
     &                      ad_v(:,:,:,Lnew), ad_v(:,:,:,Lnew),         &
#else
     &                      ad_ubar(:,:,Lnew), ad_ubar(:,:,Lnew),       &
     &                      ad_vbar(:,:,Lnew), ad_vbar(:,:,Lnew),       &
#endif
     &                      ad_zeta(:,:,Lnew), ad_zeta(:,:,Lnew))
      ELSE
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), dot(0:),                     &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef ADJUST_WSTRESS
     &                      ad_ustr(:,:,:,Lold), ad_ustr(:,:,:,Lnew),   &
     &                      ad_vstr(:,:,:,Lold), ad_vstr(:,:,:,Lnew),   &
#endif 
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      ad_tflux(:,:,:,Lold,:),                     &
     &                      ad_tflux(:,:,:,Lnew,:),                     &
# endif 
     &                      ad_t(:,:,:,Lold,:), ad_t(:,:,:,Lnew,:),     &
     &                      ad_u(:,:,:,Lold), ad_u(:,:,:,Lnew),         &
     &                      ad_v(:,:,:,Lold), ad_v(:,:,:,Lnew),         &
#else
     &                      ad_ubar(:,:,Lold), ad_ubar(:,:,Lnew),       &
     &                      ad_vbar(:,:,Lold), ad_vbar(:,:,Lnew),       &
#endif
     &                      ad_zeta(:,:,Lold), ad_zeta(:,:,Lnew))
      ENDIF
!
!  Need to multiply dot(0) by zgnorm because the gradient (index Lold)
!  has been normalized.
!
      cg_QG(innLoop+1,outLoop)=cg_Gnorm(outLoop)*dot(0)

#ifdef TEST_ORTHOGONALIZATION
!
!-----------------------------------------------------------------------
!  Test orthogonal properties of the new gradient.
!-----------------------------------------------------------------------
!
      DO rec=innLoop,1,-1
!
!  Determine adjoint file to process.
!
        IF (ndefADJ(ng).gt.0) THEN
          lstr=LEN_TRIM(ADJbase(ng))
          WRITE (ncname,10) ADJbase(ng)(1:lstr-3), rec
        ELSE
          ncname=ADJname(ng)
        END IF
!
!  Read in each previous gradient state solutions, q(0) to q(k), and
!  compute its associated dot angaint orthogonalized q(k+1). Again, 
!  each gradient solution is loaded into TANGENT LINEAR STATE ARRAYS
!  at index Lwrk.
!
        CALL read_state (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Lwrk, rec,                                     &
     &                   ndefADJ(ng), ncADJid(ng), ncname,              &
# ifdef MASKING
     &                   rmask, umask, vmask,                           &
# endif
# ifdef ADJUST_WSTRESS
     &                   tl_ustr, tl_vstr,                              &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                   tl_tflux,                                      &
#  endif
     &                   tl_t, tl_u, tl_v,                              &
# else
     &                   tl_ubar, tl_vbar,                              &
# endif
     &                   tl_zeta)
!
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), dot(0:),                     &
# ifdef MASKING
     &                      rmask, umask, vmask,                        &
# endif
# ifdef ADJUST_WSTRESS
     &                      ad_ustr(:,:,:,Lnew), tl_ustr(:,:,:,Lwrk),   &
     &                      ad_vstr(:,:,:,Lnew), tl_vstr(:,:,:,Lwrk),   &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                      ad_tflux(:,:,:,Lnew,:),                     &
     &                      tl_tflux(:,:,:,Lwrk,:),                     &
#  endif
     &                      ad_t(:,:,:,Lnew,:), tl_t(:,:,:,Lwrk,:),     &
     &                      ad_u(:,:,:,Lnew), tl_u(:,:,:,Lwrk),         &
     &                      ad_v(:,:,:,Lnew), tl_v(:,:,:,Lwrk),         &
# else
     &                      ad_ubar(:,:,Lnew), tl_ubar(:,:,Lwrk),       &
     &                      ad_vbar(:,:,Lnew), tl_vbar(:,:,Lwrk),       &
# endif
     &                      ad_zeta(:,:,Lnew), tl_zeta(:,:,Lwrk))
        dot_new(rec)=dot(0)
      END DO
!
!  Report dot products. If everything is working correctly, at the
!  end of the orthogonalization dot_new(rec) << dot_old(rec).
!
      IF (Master) THEN
        WRITE (stdout,20) outLoop, innLoop
        DO rec=innLoop,1,-1
          WRITE (stdout,30) DotProd(rec), rec-1
        END DO
        WRITE (stdout,*) ' '
        DO rec=innLoop,1,-1
          WRITE (stdout,40) innLoop, rec-1, dot_new(rec),               &
     &                      rec-1, rec-1, dot_old(rec)
        END DO
 20     FORMAT (/,1x,'(',i3.3,',',i3.3,'): ',                           &
     &          'Gramm-Schmidt Orthogonalization:',/)
 30     FORMAT (12x,'Orthogonalization Factor = ',1p,e19.12,3x,         &
     &          '(Iter=',i3.3,')')
 40     FORMAT (2x,'Ortho Test: ',                                      &
     &          '<G(',i3.3,'),G(',i3.3,')> = ',1p,e15.8,1x,             &
     &          '<G(',i3.3,'),G(',i3.3,')> = ',1p,e15.8)
      END IF
#endif

      RETURN
      END SUBROUTINE lanczos
!
!***********************************************************************
      SUBROUTINE new_gradient (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Lold, Lnew, Lwrk,                        &
     &                         innLoop, outLoop,                        &
#ifdef MASKING
     &                         rmask, umask, vmask,                     &
#endif
#ifdef ADJUST_WSTRESS
     &                         tl_ustr, tl_vstr,                        &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                         tl_tflux,                                &
# endif
     &                         tl_t, tl_u, tl_v,                        &
#else
     &                         tl_ubar, tl_vbar,                        &
#endif
     &                         tl_zeta,                                 &
#ifdef ADJUST_WSTRESS
     &                         ad_ustr, ad_vstr,                        &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                         ad_tflux,                                &
# endif
     &                         ad_t, ad_u, ad_v,                        &
#else
     &                         ad_ubar, ad_vbar,                        &
#endif
     &                         ad_zeta)
!
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE state_addition_mod, ONLY : state_addition
      USE state_dotprod_mod, ONLY : state_dotprod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lold, Lnew, Lwrk
      integer, intent(in) :: innLoop, outLoop
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: i, j, lstr, rec
#ifdef SOLVE3D
      integer :: itrc, k
#endif
      real(r8) :: fac1, fac2

      real(r8), dimension(0:NstateVar(ng)) :: dot
      real(r8), dimension(0:Ninner) :: DotProd, dot_new, dot_old

      character (len=80) :: ncname

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Computes the gradient of the cost function at the new point.
!-----------------------------------------------------------------------
!
!  Need to multiply the gradient (index Lold) by cg_Gnorm because it has
!  been normalized:
!
!    ad_var(Lold) = fac1 * ad_var(Lold) + fac2 * ad_var(Lnew)
!     
      fac1=cg_Gnorm(outLoop)
      fac2=cg_beta(innLoop+1,outLoop)*cg_Tmatrix(innLoop,3)

      CALL state_addition (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Lold, Lnew, Lold, fac1, fac2,                &
#ifdef MASKING
     &                     rmask, umask, vmask,                         &
#endif
#ifdef ADJUST_WSTRESS
     &                     ad_ustr, ad_ustr,                            &
     &                     ad_vstr, ad_vstr,                            &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                     ad_tflux, ad_tflux,                          &
# endif
     &                     ad_t, ad_t,                                  &
     &                     ad_u, ad_u,                                  &
     &                     ad_v, ad_v,                                  &
#else
     &                     ad_ubar, ad_ubar,                            &
     &                     ad_vbar, ad_vbar,                            &
#endif
     &                     ad_zeta, ad_zeta)
!
!  Adjust gradient against all previous gradients
!
      DO rec=1,innLoop
!
!  Determine adjoint file to process.
!
        IF (ndefADJ(ng).gt.0) THEN
          lstr=LEN_TRIM(ADJbase(ng))
          WRITE (ncname,10) ADJbase(ng)(1:lstr-3), rec
 10       FORMAT (a,'_',i3.3,'.nc')
        ELSE
          ncname=ADJname(ng)
        END IF
!
!  Read in each previous gradient state solutions, G(0) to G(k), and
!  compute its associated dot angaint curret G(k+1). Each gradient
!  solution is loaded into TANGENT LINEAR STATE ARRAYS at index Lwrk.
!
        CALL read_state (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Lwrk, rec,                                     &
     &                   ndefADJ(ng), ncADJid(ng), ncname,              &
#ifdef MASKING
     &                   rmask, umask, vmask,                           &
#endif
#ifdef ADJUST_WSTRESS
     &                   tl_ustr, tl_vstr,                              &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                   tl_tflux,                                      &
# endif
     &                   tl_t, tl_u, tl_v,                              &
#else
     &                   tl_ubar, tl_vbar,                              &
#endif
     &                   tl_zeta)
!
!  In this expression for FAC2, the term cg_QG gives the contribution
!  to the gradient of Jo, and the term cg_Tmatrix gives the contribution
!  of Jb:
!
!    ad_var(Lold) = fac1 * ad_var(Lold) + fac2 * tl_var(Lwrk)
!     
        fac1=1.0_r8
        fac2=-(cg_Tmatrix(rec,3)+cg_QG(rec,outLoop))

        CALL state_addition (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Lold, Lwrk, Lold, fac1, fac2,              &
#ifdef MASKING
     &                       rmask, umask, vmask,                       &
#endif
#ifdef ADJUST_WSTRESS
     &                       ad_ustr, tl_ustr,                          &
     &                       ad_vstr, tl_vstr,                          &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                       ad_tflux, tl_tflux,                        &
# endif
     &                       ad_t, tl_t,                                &
     &                       ad_u, tl_u,                                &
     &                       ad_v, tl_v,                                &
#else
     &                       ad_ubar, tl_ubar,                          &
     &                       ad_vbar, tl_vbar,                          &
#endif
     &                       ad_zeta, tl_zeta)
      END DO
!
!  Compute excess cost function.
!
      CALL state_dotprod (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NstateVar(ng), dot(0:),                       &
#ifdef MASKING
     &                    rmask, umask, vmask,                          &
#endif
#ifdef ADJUST_WSTRESS
     &                    ad_ustr(:,:,:,Lold), ad_ustr(:,:,:,Lold),     &
     &                    ad_vstr(:,:,:,Lold), ad_vstr(:,:,:,Lold),     &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                    ad_tflux(:,:,:,Lold,:),                       &
     &                    ad_tflux(:,:,:,Lold,:),                       &
# endif
     &                    ad_t(:,:,:,Lold,:), ad_t(:,:,:,Lold,:),       &
     &                    ad_u(:,:,:,Lold), ad_u(:,:,:,Lold),           &
     &                    ad_v(:,:,:,Lold), ad_v(:,:,:,Lold),           &
#else
     &                    ad_ubar(:,:,Lold), ad_ubar(:,:,Lold),         &
     &                    ad_vbar(:,:,Lold), ad_vbar(:,:,Lold),         &
#endif
     &                    ad_zeta(:,:,Lold), ad_zeta(:,:,Lold))

      cg_Greduc(innLoop,outLoop)=SQRT(dot(0))/cg_Gnorm(outLoop)

      RETURN
      END SUBROUTINE new_gradient
!
!***********************************************************************
      SUBROUTINE hessian_evecs (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Lold, Lnew, Lwrk,                       &
     &                          innLoop, outLoop,                       &
#ifdef MASKING
     &                          rmask, umask, vmask,                    &
#endif
#ifdef ADJUST_WSTRESS
     &                          tl_ustr, tl_vstr,                       &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                          tl_tflux,                               &
# endif
     &                          tl_t, tl_u, tl_v,                       &
#else
     &                          tl_ubar, tl_vbar,                       &
#endif
     &                          tl_zeta,                                &
#ifdef ADJUST_WSTRESS
     &                          ad_ustr, ad_vstr,                       &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                          ad_tflux,                               &
# endif
     &                          ad_t, ad_u, ad_v,                       &
#else
     &                          ad_ubar, ad_vbar,                       &
#endif
     &                          ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE state_addition_mod, ONLY : state_addition
      USE state_dotprod_mod, ONLY : state_dotprod
      USE state_initialize_mod, ONLY : state_initialize
      USE state_scale_mod, ONLY : state_scale
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lold, Lnew, Lwrk
      integer, intent(in) :: innLoop, outLoop
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: i, ingood, j, lstr, rec, nvec, status, varid
#ifdef SOLVE3D
      integer :: itrc, k
#endif
      integer :: start(4), total(4)

      real(r8) :: fac, fac1, fac2

      real(r8), dimension(Ninner) :: RitzErr

      real(r8), dimension(0:NstateVar(ng)) :: dot
      real(r8), dimension(0:Ninner) :: DotProd, dot_new, dot_old

      character (len=80) :: ncname

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Calculate converged eigenvectors of the Hessian.
!-----------------------------------------------------------------------
!
!  Count and collect the converged eigenvalues.
!
      ingood=0
      DO i=innLoop,1,-1
        IF (cg_RitzErr(i,outLoop).le.RitzMaxErr) THEN
          ingood=ingood+1
          Ritz(ingood)=cg_Ritz(i,outLoop)
          RitzErr(ingood)=cg_RitzErr(i,outLoop)
        END IF
      END DO
      nConvRitz=ingood
!
!  Write out number of converged eigenvalues.
!      
      IF (OutThread) THEN
        status=nf90_inq_varid(ncHSSid(ng), 'nConvRitz', varid)
        status=nf90_put_var(ncHSSid(ng), varid, nConvRitz)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'nConvRitz', TRIM(HSSname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  First, premultiply the converged eigenvectors of the tridiagonal
!  matrix T(k) by the matrix of Lanczos vectors Q(k).  Use tangent
!  linear (index Lwrk) and adjoint (index Lold) state arrays as
!  temporary storage.
!-----------------------------------------------------------------------
!
      COLUMNS : DO nvec=innLoop,1,-1
        BOUNDED : IF (cg_RitzErr(nvec,outLoop).le.RitzMaxErr) THEN
!
!  Initialize adjoint state arrays: ad_var(Lold) = fac
!
          fac=0.0_r8

          CALL state_initialize (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           Lold, fac,                             &
#ifdef MASKING
     &                           rmask, umask, vmask,                   &
#endif
#ifdef ADJUST_WSTRESS
     &                           ad_ustr, ad_vstr,                      &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                           ad_tflux,                              &
# endif
     &                           ad_t, ad_u, ad_v,                      &
#else
     &                           ad_ubar, ad_vbar,                      &
#endif
     &                           ad_zeta)
!
!  Compute Hessian eigenvectors.
!
          ROWS : DO rec=1,innLoop
!
!  Determine adjoint file to process.
!
            IF (ndefADJ(ng).gt.0) THEN
              lstr=LEN_TRIM(ADJbase(ng))
              WRITE (ncname,20) ADJbase(ng)(1:lstr-3), rec
            ELSE
              ncname=ADJname(ng)
            END IF
!
!  Read gradient solution and load it into TANGENT LINEAR STATE ARRAYS
!  at index Lwrk.
!
            CALL read_state (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Lwrk, rec,                                 &
     &                       ndefADJ(ng), ncADJid(ng), ncname,          &
#ifdef MASKING
     &                       rmask, umask, vmask,                       &
#endif
#ifdef ADJUST_WSTRESS
     &                       tl_ustr, tl_vstr,                          &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                       tl_tflux,                                  &
# endif
     &                       tl_t, tl_u, tl_v,                          &
#else
     &                       tl_ubar, tl_vbar,                          &
#endif
     &                       tl_zeta)
!
!  Compute Hessian eigenvectors:
!
!    ad_var(Lold) = fac1 * ad_var(Lold) + fac2 * tl_var(Lwrk)
!     
            fac1=1.0_r8
            fac2=cg_zv(rec,nvec)

            CALL state_addition (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           Lold, Lwrk, Lold, fac1, fac2,          &
#ifdef MASKING
     &                           rmask, umask, vmask,                   &
#endif
#ifdef ADJUST_WSTRESS
     &                           ad_ustr, tl_ustr,                      &
     &                           ad_vstr, tl_vstr,                      &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                           ad_tflux, tl_tflux,                    &
# endif
     &                           ad_t, tl_t,                            &
     &                           ad_u, tl_u,                            &
     &                           ad_v, tl_v,                            &
#else
     &                           ad_ubar, tl_ubar,                      &
     &                           ad_vbar, tl_vbar,                      &
#endif
     &                           ad_zeta, tl_zeta)
          END DO ROWS
!
!  Write eigenvectors into Hessian NetCDF.
!
          LwrtState2d(ng)=.TRUE.
          CALL wrt_hessian (ng, Lold, Lold) 
          LwrtState2d(ng)=.FALSE.
          IF (exit_flag.ne.NoERRor) RETURN

        END IF BOUNDED

      END DO COLUMNS
!
!-----------------------------------------------------------------------
!  Second, orthonormalize the converged Hessian vectors against each
!  other. Use tangent linear state arrays (index Lwrk) as temporary
!  storage.
!-----------------------------------------------------------------------
!
!  In the following we use index Lnew adjoint state arrays as temporary
!  storage because at this point we are done with the inner loops and do
!  not need the Lanczos vector stored in it.
!
      DO nvec=1,ingood
!
!  Read in just computed Hessian eigenvectors into adjoint state array
!  index Lold.
!
        CALL read_state (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Lold, nvec,                                    &
     &                   0, ncHSSid(ng), HSSname(ng),                   &
#ifdef MASKING
     &                   rmask, umask, vmask,                           &
#endif
#ifdef ADJUST_WSTRESS
     &                   ad_ustr, ad_vstr,                              &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                   ad_tflux,                                      &
# endif
     &                   ad_t, ad_u, ad_v,                              &
#else
     &                   ad_ubar, ad_vbar,                              &
#endif
     &                   ad_zeta)
!
!  Initialize adjoint state arrays index Lnew with just read Hessian
!  vector in index Lold (initialize the summation):
!
!    ad_var(Lnew) = fac * ad_var(Lold)
!
        fac=1.0_r8

        CALL state_scale (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Lold, Lnew, fac,                              &
#ifdef MASKING
     &                    rmask, umask, vmask,                          &
#endif
#ifdef ADJUST_WSTRESS
     &                    ad_ustr, ad_vstr,                             &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                    ad_tflux,                                     &
# endif
     &                    ad_t, ad_u, ad_v,                             &
#else
     &                    ad_ubar, ad_vbar,                             &
#endif
     &                    ad_zeta)
!
!  Orthogonalize Hessian eigenvectors against each other.
!
        DO rec=1,nvec-1
!
!  Read in gradient just computed Hessian eigenvectors into tangent
!  linear state array index Lwrk.
!
          CALL read_state (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Lwrk, rec,                                   &
     &                     0, ncHSSid(ng), HSSname(ng),                 &
#ifdef MASKING
     &                     rmask, umask, vmask,                         &
#endif
#ifdef ADJUST_WSTRESS
     &                     tl_ustr, tl_vstr,                            &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                     tl_tflux,                                    &
# endif
     &                     tl_t, tl_u, tl_v,                            &
#else
     &                     tl_ubar, tl_vbar,                            &
#endif
     &                     tl_zeta)
!
!  Compute dot product.
!
          CALL state_dotprod (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NstateVar(ng), dot(0:),                   &
#ifdef MASKING
     &                        rmask, umask, vmask,                      &
#endif
#ifdef ADJUST_WSTRESS
     &                        ad_ustr(:,:,:,Lold), tl_ustr(:,:,:,Lwrk), &
     &                        ad_vstr(:,:,:,Lold), ad_vstr(:,:,:,Lwrk), &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                        ad_tflux(:,:,:,Lold,:),                   &
     &                        tl_tflux(:,:,:,Lwrk,:),                   &
# endif
     &                        ad_t(:,:,:,Lold,:), tl_t(:,:,:,Lwrk,:),   &
     &                        ad_u(:,:,:,Lold), tl_u(:,:,:,Lwrk),       &
     &                        ad_v(:,:,:,Lold), tl_v(:,:,:,Lwrk),       &
#else
     &                        ad_ubar(:,:,Lold), tl_ubar(:,:,Lwrk),     &
     &                        ad_vbar(:,:,Lold), tl_vbar(:,:,Lwrk),     &
#endif
     &                        ad_zeta(:,:,Lold), tl_zeta(:,:,Lwrk))
!
!  Orthogonalize Hessian eigenvectors:
!
!    ad_var(Lnew) = fac1 * ad_var(Lnew) + fac2 * tl_var(Lwrk)
!
          fac1=1.0_r8
          fac2=-dot(0)

          CALL state_addition (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Lnew, Lwrk, Lnew, fac1, fac2,            &
#ifdef MASKING
     &                         rmask, umask, vmask,                     &
#endif
#ifdef ADJUST_WSTRESS
     &                         ad_ustr, tl_ustr,                        &
     &                         ad_vstr, tl_vstr,                        &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                         ad_tflux, tl_tflux,                      &
# endif
     &                         ad_t, tl_t,                              &
     &                         ad_u, tl_u,                              &
     &                         ad_v, tl_v,                              &
#else
     &                         ad_ubar, tl_ubar,                        &
     &                         ad_vbar, tl_vbar,                        &
#endif
     &                         ad_zeta, tl_zeta)
        END DO
!
!  Compute normalization factor.
!
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), dot(0:),                     &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef ADJUST_WSTRESS
     &                      ad_ustr(:,:,:,Lnew), ad_ustr(:,:,:,Lnew),   &
     &                      ad_vstr(:,:,:,Lnew), ad_vstr(:,:,:,Lnew),   &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      ad_tflux(:,:,:,Lnew,:),                     &
     &                      ad_tflux(:,:,:,Lnew,:),                     &
# endif
     &                      ad_t(:,:,:,Lnew,:), ad_t(:,:,:,Lnew,:),     &
     &                      ad_u(:,:,:,Lnew), ad_u(:,:,:,Lnew),         &
     &                      ad_v(:,:,:,Lnew), ad_v(:,:,:,Lnew),         &
#else
     &                      ad_ubar(:,:,Lnew), ad_ubar(:,:,Lnew),       &
     &                      ad_vbar(:,:,Lnew), ad_vbar(:,:,Lnew),       &
#endif
     &                      ad_zeta(:,:,Lnew), ad_zeta(:,:,Lnew))
!
!  Normalize Hessian eigenvectors:
!
!    ad_var(Lnew) = fac * ad_var(Lnew)
!
        fac=1.0_r8/SQRT(dot(0))

        CALL state_scale (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Lnew, Lnew, fac,                              &
#ifdef MASKING
     &                    rmask, umask, vmask,                          &
#endif
#ifdef ADJUST_WSTRESS
     &                    ad_ustr, ad_vstr,                             &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                    ad_tflux,                                     &
# endif
     &                    ad_t, ad_u, ad_v,                             &
#else
     &                    ad_ubar, ad_vbar,                             &
#endif
     &                    ad_zeta)
!
!  Write out converged Ritz eigenvalues and is associated accuracy.
!
        IF (OutThread) THEN
          status=nf90_inq_varid(ncHSSid(ng),'Ritz',varid)
          start(1)=nvec
          total(1)=1
          status=nf90_put_var(ncHSSid(ng), varid, Ritz(nvec:),          &
     &                        start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'Ritz', TRIM(HSSname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
          status=nf90_inq_varid(ncHSSid(ng),'Ritz_error',varid)
          status=nf90_put_var(ncHSSid(ng), varid, RitzErr(nvec:),       &
     &                        start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'Ritz_error', TRIM(HSSname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
!  Replace record "nvec" of Hessian eigenvectors NetCDF with the
!  normalized value in adjoint state arrays at index Lnew.
!
        tHSSindx(ng)=nvec-1
        LwrtState2d(ng)=.TRUE.
        CALL wrt_hessian (ng, Lnew, Lnew)
        LwrtState2d(ng)=.FALSE.
        IF (exit_flag.ne.NoERRor) RETURN

      END DO

  10  FORMAT (/,' HESSIAN_EVECS - error while writing variable: ',a,/,  &
     &        17x,'into NetCDF file: ',a)
  20  FORMAT (a,'_',i3.3,'.nc')

      RETURN
      END SUBROUTINE hessian_evecs

      SUBROUTINE cg_write (ng, innLoop, outLoop)
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
      integer, intent(in) :: ng, innLoop, outLoop
!
!  Local variable declarations.
!
      logical, save :: First = .TRUE.

      integer :: i, status
      integer :: start(2), total(2)
      integer, save :: varid(18)
!
!-----------------------------------------------------------------------
!  Write out conjugate gradient vectors.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
        IF (First) THEN
          First=.FALSE.
          DO i=1,18
            varid(i)=0
          END DO
        END IF
!
!  Write out outer and inner iteration.
!
        IF (varid(1).eq.0) THEN
          status=nf90_inq_varid(ncMODid(ng), 'outer', varid(1))
        END IF        
        status=nf90_put_var(ncMODid(ng), varid(1), outer)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'outer', TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF

        IF (varid(2).eq.0) THEN
          status=nf90_inq_varid(ncMODid(ng), 'inner', varid(2))
        END IF
        status=nf90_put_var(ncMODid(ng), varid(2), inner)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'inner', TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Write out number of converged Ritz eigenvalues.
!
        IF (innLoop.eq.(Ninner-1)) THEN
          IF (varid(3).eq.0) THEN
            status=nf90_inq_varid(ncMODid(ng), 'nConvRitz', varid(3))
          END IF
          status=nf90_put_var(ncMODid(ng), varid(3), nConvRitz)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'nConvRitz', TRIM(MODname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
!  Write out converged Ritz eigenvalues.
!
        IF (innLoop.eq.(Ninner-1)) THEN
          IF (varid(4).eq.0) THEN
            status=nf90_inq_varid(ncMODid(ng), 'Ritz', varid(4))
          END IF
          start(1)=1
          total(1)=nConvRitz
          status=nf90_put_var(ncMODid(ng), varid(4), Ritz, start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'Ritz', TRIM(MODname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
!  Write out conjugate gradient norms.
!
        IF (innLoop.gt.0) THEN
          IF (varid(5).eq.0) THEN
            status=nf90_inq_varid(ncMODid(ng), 'cg_beta', varid(5))
          END IF
          start(1)=1
          total(1)=Ninner+1
          start(2)=1
          total(2)=Nouter
          status=nf90_put_var(ncMODid(ng), varid(5), cg_beta,           &
     &                        start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'cg_beta', TRIM(MODname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
        IF (varid(6).eq.0) THEN
          status=nf90_inq_varid(ncMODid(ng), 'cg_tau', varid(6))
        END IF
        start(1)=1
        total(1)=Ninner+1
        start(2)=1
        total(2)=Nouter
        status=nf90_put_var(ncMODid(ng), varid(6), cg_tau(0:,:),        &
     &                      start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'cg_tau', TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Write out Lanczos algorithm coefficients.
!
        IF (innLoop.gt.0) THEN
          IF (varid(7).eq.0) THEN
            status=nf90_inq_varid(ncMODid(ng), 'cg_delta', varid(7))
          END IF
          start(1)=1
          total(1)=Ninner
          start(2)=1
          total(2)=Nouter
          status=nf90_put_var(ncMODid(ng), varid(7), cg_delta,          &
     &                        start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'cg_delta', TRIM(MODname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
        IF (innLoop.gt.0) THEN
          IF (varid(8).eq.0) THEN
            status=nf90_inq_varid(ncMODid(ng), 'cg_gamma', varid(8))
          END IF
          start(1)=1
          total(1)=Ninner
          start(2)=1
          total(2)=Nouter
          status=nf90_put_var(ncMODid(ng), varid(8), cg_gamma,          &
     &                        start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'cg_gamma', TRIM(MODname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
!  Initial gradient normalization factor.
!
        IF (innLoop.eq.0) THEN
          IF (varid(9).eq.0) THEN
            status=nf90_inq_varid(ncMODid(ng), 'cg_Gnorm', varid(9))
          END IF
          start(1)=1
          total(1)=Nouter
          status=nf90_put_var(ncMODid(ng), varid(9), cg_Gnorm,          &
     &                        start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'cg_Gnorm', TRIM(MODname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
!  Lanczos vector normalization factor.
!
        IF (innLoop.gt.0) THEN
          IF (varid(10).eq.0) THEN
            status=nf90_inq_varid(ncMODid(ng), 'cg_QG', varid(10))
          END IF
          start(1)=1
          total(1)=Ninner
          start(2)=1
          total(2)=Nouter
          status=nf90_put_var(ncMODid(ng), varid(10), cg_QG,            &
     &                        start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'cg_QG', TRIM(MODname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
!  Reduction in the gradient norm.
!
        IF (innLoop.gt.0) THEN
          IF (varid(11).eq.0) THEN
            status=nf90_inq_varid(ncMODid(ng), 'cg_Greduc', varid(11))
          END IF
          start(1)=1
          total(1)=Ninner
          start(2)=1
          total(2)=Nouter
          status=nf90_put_var(ncMODid(ng), varid(11), cg_Greduc,        &
     &                        start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'cg_Greduc', TRIM(MODname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
!  Lanczos recurrence tridiagonal matrix.
!
        IF (innLoop.gt.0) THEN
          IF (varid(12).eq.0) THEN
            status=nf90_inq_varid(ncMODid(ng), 'cg_Tmatrix', varid(12))
          END IF
          start(1)=1
          total(1)=Ninner
          start(2)=1
          total(2)=3
          status=nf90_put_var(ncMODid(ng), varid(12), cg_Tmatrix,       &
     &                        start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'cg_Tmatrix', TRIM(MODname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
!  Lanczos tridiagonal matrix, upper diagonal elements.
!
        IF (innLoop.gt.0) THEN
          IF (varid(13).eq.0) THEN
            status=nf90_inq_varid(ncMODid(ng), 'cg_zu', varid(13))
          END IF
          start(1)=1
          total(1)=Ninner
          start(2)=1
          total(2)=Nouter
          status=nf90_put_var(ncMODid(ng), varid(13), cg_zu,            &
     &                        start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'cg_zu', TRIM(MODname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
!  Eigenvalues of Lanczos recurrence relationship.
!
        IF (innLoop.gt.0) THEN
          IF (varid(14).eq.0) THEN
            status=nf90_inq_varid(ncMODid(ng), 'cg_Ritz', varid(14))
          END IF
          start(1)=1
          total(1)=Ninner
          start(2)=1
          total(2)=Nouter
          status=nf90_put_var(ncMODid(ng), varid(14), cg_Ritz,          &
     &                        start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'cg_Ritz', TRIM(MODname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
!  Eigenvalues relative error.
!
        IF (innLoop.gt.0) THEN
          IF (varid(15).eq.0) THEN
            status=nf90_inq_varid(ncMODid(ng), 'cg_RitzErr', varid(15))
          END IF
          start(1)=1
          total(1)=Ninner
          start(2)=1
          total(2)=Nouter
          status=nf90_put_var(ncMODid(ng), varid(15), cg_RitzErr,       &
     &                        start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'cg_RitzErr', TRIM(MODname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
!  Eigenvectors of Lanczos recurrence relationship.
!
        IF (innLoop.gt.0) THEN
          IF (varid(16).eq.0) THEN
            status=nf90_inq_varid(ncMODid(ng), 'cg_zv', varid(16))
          END IF
          start(1)=1
          total(1)=Ninner
          start(2)=1
          total(2)=Ninner
          status=nf90_put_var(ncMODid(ng), varid(16), cg_zv,            &
     &                        start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'cg_zv', TRIM(MODname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
!  Write out Lanczos algorithm coefficients into the adjoint NetCDF
!  file.  These coefficients can be used to compute the sensitivity
!  of the observations to the 4DVAR data assimilation system.
!
        IF (innLoop.gt.0) THEN
          IF (varid(17).eq.0) THEN
            status=nf90_inq_varid(ncADJid(ng), 'cg_beta', varid(17))
          END IF
          start(1)=1
          total(1)=Ninner+1
          start(2)=1
          total(2)=Nouter
          status=nf90_put_var(ncADJid(ng), varid(17), cg_beta,          &
     &                        start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'cg_beta', TRIM(ADJname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
!
          IF (varid(18).eq.0) THEN
            status=nf90_inq_varid(ncADJid(ng), 'cg_delta', varid(18))
          END IF
          start(1)=1
          total(1)=Ninner
          start(2)=1
          total(2)=Nouter
          status=nf90_put_var(ncADJid(ng), varid(18), cg_delta,         &
     &                        start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'cg_delta', TRIM(ADJname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF

        END IF
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
