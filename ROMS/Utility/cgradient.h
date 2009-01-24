!
!svn $Id: cgradient.h 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module minimizes  incremental  4Dvar  quadratic cost function  !
!  using a preconditioned version of the conjugate gradient algorithm  !
!  proposed by Mike Fisher (ECMWF).                                    !
!                                                                      !
!  In the following,  M represents the preconditioner.  Specifically,  !
!                                                                      !
!    M = I + SUM_i [ (mu_i-1) h_i (h_i)'],                             !
!                                                                      !
!  where mu_i can take the following values:                           !
!                                                                      !
!    Lscale= 1:    mu_i = lambda_i                                     !
!    Lscale=-1:    mu_i = 1 / lambda_i                                 !
!    Lscale= 2:    mu_i = SQRT (lambda_i)                              !
!    Lscale=-2:    mu_i = 1 / SQRT(lambda_i)                           !
!                                                                      !
!  where lambda_i are the Hessian eigenvalues and h_i are the Hessian  !
!  eigenvectors. For Lscale=-1 the preconditioner is an approximation  !
!  of the inverse Hessian matrix constructed from the leading Hessian  !
!  eigenvectors.  A full description is given in  Fisher and Courtier  !
!  (1995), ECMWF Tech Memo 220.                                        !
!                                                                      !
!  Given an initial model state X(0), gradient G(0), descent direction !
!  d(0), and trial step size tau(1), the minimization algorithm at the !
!  k-iteration is :                                                    !
!                                                                      !
!  (1) Run tangent linear model initialized with trial step, Xhat(k):  !
!                                                                      !
!      Xhat(k) = X(k) + tau(k) * d(k)                          (Eq 5a) !
!                                                                      !
!  (2) Run adjoint model to compute gradient at trial point, Ghat(k):  !
!                                                                      !
!      Ghat(k) = GRAD[ f(Xhat(k)) ]                            (Eq 5b) !
!                                                                      !
!  (3) Compute optimum step size, alpha(k):                            !
!                                                                      !
!      alpha(k) = tau(k) * <d(k), G(k)> /                              !
!                                                                      !
!                 (<d(k),G(k)> - <d(k), Ghat(k)>)              (Eq 5c) !
!                                                                      !
!      here <...> denotes dot product                                  !
!                                                                      !
!  (4) Compute new starting point (TLM increments), X(k+1):            !
!                                                                      !
!      X(k+1) = X(k) + alpha(k) * d(k)                         (Eq 5d) !
!                                                                      !
!  (5) Compute gradient at new point, G(k+1):                          !
!                                                                      !
!      G(k+1) = G(k) + (alpha(k) / tau(k)) * (Ghat(k) - G(k))  (Eq 5e) !
!                                                                      !
!      overwrite G(k+1) in the NetCDF for latter use.                  !
!                                                                      !
!  (6) Orthogonalize new gradient, G(k+1), against all previous        !
!      gradients [G(k), ..., G(0)], in reverse order, using the        !
!      modified Gramm-Schmidt algorithm. Notice that we need to        !
!      save all inner loop gradient solutions.                         !
!                                                                      !
!      For the preconditioned case the appropriate inner-product       !
!      for the orthonormalizatio is <G,MG>.                            !
!                                                                      !
!  (7) Compute new descent direction, d(k+1):                          !
!                                                                      !
!      beta(k+1) = <G(k+1),M G(k+1)> / <G(k),M G(k)>           (Eq 5g) !
!                                                                      !
!      d(k+1) = - MG(k+1) + beta(k+1) * d(k)                   (Eq 5f) !
!                                                                      !
!  After the first iteration, the trial step size is:                  !
!                                                                      !
!      tau(k) = alpha(k-1)                                             !
!                                                                      !
!  NOTE: In all of the following computations we are using the NLM     !
!        state variable arrays as temporary arrays.                    !
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
     &                     FORCES(ng) % ustr,                           &
     &                     FORCES(ng) % vstr,                           &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                     FORCES(ng) % tflux,                          &
# endif
     &                     OCEAN(ng) % t,                               &
     &                     OCEAN(ng) % u,                               &
     &                     OCEAN(ng) % v,                               &
#else
     &                     OCEAN(ng) % ubar,                            &
     &                     OCEAN(ng) % vbar,                            &
#endif
     &                     OCEAN(ng) % zeta,                            &
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
     &                     OCEAN(ng) % ad_zeta,                         &
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
     &                     OCEAN(ng) % tl_zeta)
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
     &                           nl_ustr, nl_vstr,                      &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                           nl_tflux,                              &
# endif
     &                           nl_t, nl_u, nl_v,                      &
#else
     &                           nl_ubar, nl_vbar,                      &
#endif
     &                           nl_zeta,                               &
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
     &                           ad_zeta,                               &
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
     &                           tl_zeta)
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
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
#endif
      USE state_dotprod_mod, ONLY : state_dotprod
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
      real(r8), intent(inout) :: nl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: nl_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: nl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: nl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: nl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: nl_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: nl_zeta(LBi:,LBj:,:)
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
      real(r8), intent(inout) :: nl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: nl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: nl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: nl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: nl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: nl_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: nl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: nl_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: nl_zeta(LBi:UBi,LBj:UBj,3)
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
      integer :: Linp, Lout, Lwrk, L1, L2, i, Lscale
      integer :: status, varid
      integer :: start(4), total(4)

      real(r8) :: norm

      real(r8), dimension(0:NstateVar(ng)) :: Adjust
      real(r8), dimension(0:NstateVar(ng)) :: dot_old, dot_new
      real(r8), dimension(0:NstateVar(ng)) :: old_dot, new_dot

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize trial step size.
!-----------------------------------------------------------------------
!
      IF (innLoop.eq.0) THEN
        cg_tau(innLoop,outLoop)=CGstepI
        cg_alpha(innLoop,outLoop)=cg_tau(innLoop,outLoop)
        DO i=0,NstateVar(ng)
          dot_old(i)=0.0_r8
          dot_new(i)=0.0_r8
          old_dot(i)=0.0_r8
          new_dot(i)=0.0_r8
          FOURDVAR(ng)%CostGradDot(i)=0.0_r8
        END DO
      END IF
      IF (Master) THEN
        WRITE (stdout,10)
 10     FORMAT (/,' <<<< Descent Algorithm >>>>')
      END IF
!
!  If preconditioning, read in number of converged eigenvectors and their
!  associated eigenvalues.
!
      IF (Lprecond.and.((innLoop.eq.0).and.(outLoop.eq.1))) THEN
        IF (InpThread) THEN
          status=nf90_inq_varid(ncHSSid(ng), 'nConvRitz', varid)
          status=nf90_get_var(ncHSSid(ng), varid, nConvRitz)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,20) 'nConvRitz', TRIM(HSSname(ng))
 20         FORMAT (/,' CGRADIENT - error while reading variable: ',    &
     &              a,/,12x,'from NetCDF file: ',a)
            exit_flag=2
            ioerror=status
            RETURN
          END IF
          status=nf90_inq_varid(ncHSSid(ng), 'Ritz', varid)
          start(1)=1
          total(1)=nConvRitz
          status=nf90_get_var(ncHSSid(ng), varid, Ritz, start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,20) 'Ritz', TRIM(HSSname(ng))
            exit_flag=2
            ioerror=status
            RETURN
          END IF
        END IF
#ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iADM, nConvRitz, 1)
        CALL mp_bcastf (ng, iADM, Ritz, nConvRitz)
#endif
!
!  Reset number of eigenpairs to use to specified value.
!
        nConvRitz=NritzEV
      END IF
!
!-----------------------------------------------------------------------
!  Compute conjugate gradient optimum step size, alpha(k).
!-----------------------------------------------------------------------
!
      IF (innLoop.gt.0) THEN
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), dot_old(0:),                 &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef ADJUST_WSTRESS
     &                      d_sustr, ad_ustr(:,:,:,Lold),               &
     &                      d_svstr, ad_vstr(:,:,:,Lold),               &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      d_stflx, ad_tflux(:,:,:,Lold,:),            &
# endif
     &                      d_t, ad_t(:,:,:,Lold,:),                    &
     &                      d_u, ad_u(:,:,:,Lold),                      &
     &                      d_v, ad_v(:,:,:,Lold),                      &
#else
     &                      d_ubar, ad_ubar(:,:,Lold),                  &
     &                      d_vbar, ad_vbar(:,:,Lold),                  &
#endif
     &                      d_zeta, ad_zeta(:,:,Lold))
!
!  If preconditioning, compute new dot product, <d(k), H^-1 * Ghat(k)>.
!
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), dot_new(0:),                 &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef ADJUST_WSTRESS
     &                      d_sustr, ad_ustr(:,:,:,Lnew),               &
     &                      d_svstr, ad_vstr(:,:,:,Lnew),               &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      d_stflx, ad_tflux(:,:,:,Lnew,:),            &
# endif
     &                      d_t, ad_t(:,:,:,Lnew,:),                    &
     &                      d_u, ad_u(:,:,:,Lnew),                      &
     &                      d_v, ad_v(:,:,:,Lnew),                      &
#else
     &                      d_ubar, ad_ubar(:,:,Lnew),                  &
     &                      d_vbar, ad_vbar(:,:,Lnew),                  &
#endif
     &                      d_zeta, ad_zeta(:,:,Lnew))
!
!  Compute new optimal step size.
!
        cg_tau(innLoop,outLoop)=cg_alpha(innLoop-1,outLoop)
        cg_alpha(innLoop,outLoop)=cg_tau(innLoop,outLoop)*              &
     &                            dot_old(0)/(dot_old(0)-dot_new(0))
      END IF
!
!  Adjust the cost function for the previous inner-loop iteration.
!  This is based on a first-order Taylor expansion of the cost function.
!  Let vhat=v+tau*d. During each inner-loop the tangent linear
!  model provides J(vhat). What we require is J(v). Using a 1st-order
!  Taylor expansion we have: J(vhat)=J(v)+tau*<d,grad> where grad is
!  the cost function gradient computed during the last inner-loop
!  immediately prior to the orthogonalization. Rearranging this
!  equation we have: J(v)=J(vhat)-tau*<d,grad>. In the code
!  J(vhat)=CostFun(:) and <d,grad>=CostFunDot(:). Remember though
!  that J(v) is the cost function associated with v from the previous
!  inner-loop.
!
      DO i=0,NstateVar(ng)
        Adjust(i)=cg_tau(innLoop,outLoop)*FOURDVAR(ng)%CostGradDot(i)
        FOURDVAR(ng)%CostFun(i)=FOURDVAR(ng)%CostFun(i)-Adjust(i)
      END DO
!
!-----------------------------------------------------------------------
!  Estimate the gradient for the new state vector, G(k+1).
!-----------------------------------------------------------------------
!
!  If preconditioning, compute old dot product, <G(k), H^-1 * G(k)>.
!  The ADM arrays, index Lold, will be used a as temporary storage
!  after this.
!
      IF (Lprecond) THEN
        Lscale=-1
        Lwrk=2
        CALL precond (ng, tile, model,                                  &
     &                LBi, UBi, LBj, UBj,                               &
     &                NstateVar(ng), Lold, Lwrk, Lscale,                &
     &                nConvRitz, Ritz,                                  &
#ifdef MASKING
     &                rmask, umask, vmask,                              &
#endif
#ifdef ADJUST_WSTRESS
     &                ad_ustr, nl_ustr, tl_ustr,                        &
     &                ad_vstr, nl_vstr, tl_vstr,                        &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                ad_tflux, nl_tflux, tl_tflux,                     &
# endif
     &                ad_t, nl_t, tl_t,                                 &
     &                ad_u, nl_u, tl_u,                                 &
     &                ad_v, nl_v, tl_v,                                 &
#else
     &                ad_ubar, nl_ubar, tl_ubar,                        &
     &                ad_vbar, nl_vbar, tl_vbar,                        &
#endif
     &                ad_zeta, nl_zeta, tl_zeta)
!
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), old_dot(0:),                 &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef ADJUST_WSTRESS
     &                      ad_ustr(:,:,:,Lold), tl_ustr(:,:,:,Lwrk),   &
     &                      ad_vstr(:,:,:,Lold), tl_vstr(:,:,:,Lwrk),   &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      ad_tflux(:,:,:,Lold,:),                     &
     &                      tl_tflux(:,:,:,Lwrk,:),                     &
# endif
     &                      ad_t(:,:,:,Lold,:), tl_t(:,:,:,Lwrk,:),     &
     &                      ad_u(:,:,:,Lold), tl_u(:,:,:,Lwrk),         &
     &                      ad_v(:,:,:,Lold), tl_v(:,:,:,Lwrk),         &
#else
     &                      ad_ubar(:,:,Lold), tl_ubar(:,:,Lwrk),       &
     &                      ad_vbar(:,:,Lold), tl_vbar(:,:,Lwrk),       &
#endif
     &                      ad_zeta(:,:,Lold), tl_zeta(:,:,Lwrk))
!
!  If not preconditioning, compute old dot product, <G(k), G(k)>.
!  The ADM arrays, index Lold, will be used a as temporary storage
!  after this.
!
      ELSE
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), old_dot(0:),                 &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef ADJUST_WSTRESS
     &                      ad_ustr(:,:,:,Lold), ad_ustr(:,:,:,Lold),   &
     &                      ad_vstr(:,:,:,Lold), ad_vstr(:,:,:,Lold),   &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      ad_tflux(:,:,:,Lold,:),                     &
     &                      ad_tflux(:,:,:,Lold,:),                     &
# endif
     &                      ad_t(:,:,:,Lold,:), ad_t(:,:,:,Lold,:),     &
     &                      ad_u(:,:,:,Lold), ad_u(:,:,:,Lold),         &
     &                      ad_v(:,:,:,Lold), ad_v(:,:,:,Lold),         &
#else
     &                      ad_ubar(:,:,Lold), ad_ubar(:,:,Lold),       &
     &                      ad_vbar(:,:,Lold), ad_vbar(:,:,Lold),       &
#endif
     &                      ad_zeta(:,:,Lold), ad_zeta(:,:,Lold))
      END IF
!
!  Notice that the current gradient Ghat(k) in time index Lnew is
!  overwritten with the new gradient G(k+1).
!
!    G(k+1) = G(k) + (alpha(k) / tau(k)) * (Ghat(k) - G(k))
!    Lnew     Lold                          Lnew      Lold      index
!
!  Also save G(k+1) in time index Lold as a non-orthogonalized new
!  gradient.
!
      CALL ad_new_state (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Lold, Lnew,                                    &
     &                   cg_alpha(innLoop,outLoop),                     &
     &                   cg_tau(innLoop,outLoop),                       &
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

#ifdef ORTHOGONALIZATION
!
!  Orthogonalize new gradient, G(k+1), against all previous gradients
!  G(0) to G(k). Use TLM state arrays at time index Lwrk=2, to load
!  each of the previous gradients.
!
      IF (innLoop.gt.0) THEN
        Lwrk=2
        CALL orthogonalize (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Lold, Lnew, Lwrk,                           &
     &                      innLoop, outLoop,                           &
# ifdef MASKING
     &                      rmask, umask, vmask,                        &
# endif
# ifdef ADJUST_WSTRESS
     &                      nl_ustr, nl_vstr,                           &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                      nl_tflux,                                   &
#  endif
     &                      nl_t, nl_u, nl_v,                           &
# else
     &                      nl_ubar, nl_vbar,                           &
# endif
     &                      nl_zeta,                                    &
# ifdef ADJUST_WSTRESS
     &                      tl_ustr, tl_vstr,                           &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                      tl_tflux,                                   &
#  endif
     &                      tl_t, tl_u, tl_v,                           &
# else
     &                      tl_ubar, tl_vbar,                           &
# endif
     &                      tl_zeta,                                    &
# ifdef ADJUST_WSTRESS
     &                      ad_ustr, ad_vstr,                           &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                      ad_tflux,                                   &
#  endif
     &                      ad_t, ad_u, ad_v,                           &
# else
     &                      ad_ubar, ad_vbar,                           &
# endif
     &                      ad_zeta)
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Compute new starting tangent linear state vector, X(k+1).
!-----------------------------------------------------------------------
!
!  Here we are doing step (4), equation 5d, the new TLM increment for
!  the initial conditions are always saved at time level Lout=1.
!
!    X(k+1) = X(k) + alpha(k) * d(k)
!    Lout     Linp                      index
!
      IF (innLoop.gt.0) THEN
        Linp=1
        Lout=1
        CALL tl_new_state (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Linp, Lout,                                  &
     &                     cg_alpha(innLoop,outLoop),                   &
#ifdef MASKING
     &                     rmask, umask, vmask,                         &
#endif
#ifdef ADJUST_WSTRESS
     &                     d_sustr, d_svstr,                            &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                     d_stflx,                                     &
# endif
     &                     d_t, d_u, d_v,                               &
#else
     &                     d_ubar, d_vbar,                              &
#endif
     &                     d_zeta,                                      &
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
!  If last iteration of inner loop, skip remaining computations. The
!  TLM increments computed here are the ones that are needed update
!  the NLM model initial conditions.
!
!!      IF (innLoop.eq.Ninner) RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Compute new conjugate descent direction, d(k+1).
!-----------------------------------------------------------------------
!
!  If preconditioning, multiply the new gradient by H^-1 and save in
!  nl_var(Lwrk).
!
      IF (Lprecond) THEN
        Lscale=-1
        Lwrk=2
        CALL precond (ng, tile, model,                                  &
     &                LBi, UBi, LBj, UBj,                               &
     &                NstateVar(ng), Lnew, Lwrk, Lscale,                &
     &                nConvRitz, Ritz,                                  &
#ifdef MASKING
     &                rmask, umask, vmask,                              &
#endif
#ifdef ADJUST_WSTRESS
     &                ad_ustr, nl_ustr, nl_ustr,                        &
     &                ad_vstr, nl_vstr, nl_vstr,                        &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                ad_tflux, nl_tflux, nl_tflux,                     &
# endif
     &                ad_t, nl_t, nl_t,                                 &
     &                ad_u, nl_u, nl_u,                                 &
     &                ad_v, nl_v, nl_v,                                 &
#else
     &                ad_ubar, nl_ubar, nl_ubar,                        &
     &                ad_vbar, nl_vbar, nl_vbar,                        &
#endif
     &                ad_zeta, nl_zeta, nl_zeta)
       END IF
!
!  If preconditioning, compute new dot product, <G(k+1), H^-1 * G(k+1)>.
!
      IF (innLoop.gt.0) THEN
        IF (Lprecond) THEN
          CALL state_dotprod (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NstateVar(ng), new_dot(0:),               &
#ifdef MASKING
     &                        rmask, umask, vmask,                      &
#endif
#ifdef ADJUST_WSTRESS
     &                        ad_ustr(:,:,:,Lnew), nl_ustr(:,:,:,Lwrk), &
     &                        ad_vstr(:,:,:,Lnew), nl_vstr(:,:,:,Lwrk), &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                        ad_tflux(:,:,:,Lnew,:),                   &
     &                        nl_tflux(:,:,:,Lwrk,:),                   &
# endif
     &                        ad_t(:,:,:,Lnew,:), nl_t(:,:,:,Lwrk,:),   &
     &                        ad_u(:,:,:,Lnew), nl_u(:,:,:,Lwrk),       &
     &                        ad_v(:,:,:,Lnew), nl_v(:,:,:,Lwrk),       &
#else
     &                        ad_ubar(:,:,Lnew), nl_ubar(:,:,Lwrk),     &
     &                        ad_vbar(:,:,Lnew), nl_vbar(:,:,Lwrk),     &
#endif
     &                        ad_zeta(:,:,Lnew), nl_zeta(:,:,Lwrk))
        ELSE
!
!  If not preconditioning, compute new dot product, <G(k+1), G(k+1)>.
!
          CALL state_dotprod (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NstateVar(ng), new_dot(0:),               &
#ifdef MASKING
     &                        rmask, umask, vmask,                      &
#endif
#ifdef ADJUST_WSTRESS
     &                        ad_ustr(:,:,:,Lnew), ad_ustr(:,:,:,Lnew), &
     &                        ad_vstr(:,:,:,Lnew), ad_vstr(:,:,:,Lnew), &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                        ad_tflux(:,:,:,Lnew,:),                   &
     &                        ad_tflux(:,:,:,Lnew,:),                   &
# endif
     &                        ad_t(:,:,:,Lnew,:), ad_t(:,:,:,Lnew,:),   &
     &                        ad_u(:,:,:,Lnew), ad_u(:,:,:,Lnew),       &
     &                        ad_v(:,:,:,Lnew), ad_v(:,:,:,Lnew),       &
#else
     &                        ad_ubar(:,:,Lnew), ad_ubar(:,:,Lnew),     &
     &                        ad_vbar(:,:,Lnew), ad_vbar(:,:,Lnew),     &
#endif
     &                        ad_zeta(:,:,Lnew), ad_zeta(:,:,Lnew))
        END IF
!
!  Compute conjugate direction coefficient, beta(k+1).
!
        cg_beta(innLoop,outLoop)=new_dot(0)/old_dot(0)
      ELSE
        cg_beta(innLoop,outLoop)=0.0_r8
      END IF
!
!  If preconditioning, compute new conjugate direction, d(k+1).  Notice
!  that the preconditined gradient is in NLM (index Lwrk) state arrays.
!
      IF (Lprecond) THEN
        CALL new_direction (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Lold, Lwrk,                                 &
     &                      cg_beta(innLoop,outLoop),                   &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef ADJUST_WSTRESS
     &                      nl_ustr, nl_vstr,                           &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      nl_tflux,                                   &
# endif
     &                      nl_t, nl_u, nl_v,                           &
#else
     &                      nl_ubar, nl_vbar,                           &
#endif
     &                      nl_zeta,                                    &
#ifdef ADJUST_WSTRESS
     &                      d_sustr, d_svstr,                           &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      d_stflx,                                    &
# endif
     &                      d_t, d_u, d_v,                              &
#else
     &                      d_ubar, d_vbar,                             &
#endif
     &                      d_zeta)
!
!  If not preconditioning, compute new conjugate direction, d(k+1).
!
      ELSE
        CALL new_direction (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Lold, Lnew,                                 &
     &                      cg_beta(innLoop,outLoop),                   &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef ADJUST_WSTRESS
     &                      ad_ustr, ad_vstr,                           &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      ad_tflux,                                   &
# endif
     &                      ad_t, ad_u, ad_v,                           &
#else
     &                      ad_ubar, ad_vbar,                           &
#endif
     &                      ad_zeta,                                    &
#ifdef ADJUST_WSTRESS
     &                      d_sustr, d_svstr,                           &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      d_stflx,                                    &
# endif
     &                      d_t, d_u, d_v,                              &
#else
     &                      d_ubar, d_vbar,                             &
#endif
     &                      d_zeta)
      END IF
!
!  Compute next iteration dot product, <d(k), G(k)>, using new d(k+1)
!  and non-orthogonalized G(k+1) used to adjust cost function.
!
      CALL state_dotprod (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NstateVar(ng),                                &
     &                    FOURDVAR(ng)%CostGradDot(0:),                 &
#ifdef MASKING
     &                    rmask, umask, vmask,                          &
#endif
#ifdef ADJUST_WSTRESS
     &                    d_sustr, ad_ustr(:,:,:,Lold),                 &
     &                    d_svstr, ad_vstr(:,:,:,Lold),                 &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                    d_stflx, ad_tflux(:,:,:,Lold,:),              &
# endif
     &                    d_t, ad_t(:,:,:,Lold,:),                      &
     &                    d_u, ad_u(:,:,:,Lold),                        &
     &                    d_v, ad_v(:,:,:,Lold),                        &
#else
     &                    d_ubar, ad_ubar(:,:,Lold),                    &
     &                    d_vbar, ad_vbar(:,:,Lold),                    &
#endif
     &                    d_zeta, ad_zeta(:,:,Lold))
!
!-----------------------------------------------------------------------
!  Set TLM initial conditions for next inner loop, Xhat(k+1).
!-----------------------------------------------------------------------
!
!  Here we are doing step (1), equation 5a, the new TLM initial
!  conditions for the next inner loop are always saved at Lout=2.
!
!    Xhat(k+1) = X(k+1) + tau(k+1) * d(k+1),  where  tau(k+1)=alpha(k)
!    Lout        Linp                         index
!
      Linp=1
      Lout=2
      CALL tl_new_state (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Linp, Lout,                                    &
     &                   cg_alpha(innLoop,outLoop),                     &
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
     &                   tl_zeta)
!
!-----------------------------------------------------------------------
!  Write out conjugate gradient information into NetCDF file.
!-----------------------------------------------------------------------
!
      CALL cg_write (ng, innLoop, outLoop)
!
!  Report  algorithm parameters.
!
      IF (Master) THEN
        WRITE (stdout,30) outLoop, innLoop,                             &
     &                    cg_tau(innLoop,outLoop),                      &
     &                    cg_alpha(innLoop,outLoop),                    &
     &                    cg_beta(innLoop,outLoop),                     &
     &                    outLoop, MAX(0,innLoop-1), Adjust(0),         &
     &                    outLoop, innLoop,                             &
     &                    'dot product', innLoop, innLoop,              &
     &                    dot_old(0), 'alpha',                          &
     &                    'dot product', innLoop, innLoop,              &
     &                    dot_new(0), 'alpha',                          &
     &                    'dot product', innLoop, innLoop,              &
     &                    old_dot(0), 'beta',                           &
     &                    'dot product', innLoop+1, innLoop+1,          &
     &                     new_dot(0), 'beta'
 30     FORMAT (/,1x,'(',i3.3,',',i3.3,'): ',                           &
     &          'tau = ',1p,e14.7,                                      &
     &          ', alpha = ',1p,e14.7,                                  &
     &          ', Beta = ',1p,e14.7,                                   &
     &          /,1x,'(',i3.3,',',i3.3,'): ',                           &
     &          'Total COST Function Adjustment = ',1p,e19.12,          &
     &          /,1x,'(',i3.3,',',i3.3,'): ',                           &
     &          a,' <d(',i3.3,'),G(',i3.3,')> = ',1p,e19.12,3x,a,/,12x, &
     &          a,' <d(',i3.3,'),g(',i3.3,')> = ',1p,e19.12,3x,a,/,12x, &
     &          a,' <G(',i3.3,'),G(',i3.3,')> = ',1p,e19.12,3x,a,/,12x, &
     &          a,' <G(',i3.3,'),G(',i3.3,')> = ',1p,e19.12,3x,a,/)
      END IF

      RETURN
      END SUBROUTINE cgradient_tile

!
!***********************************************************************
      SUBROUTINE tl_new_state (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Linp, Lout, alphaK,                      &
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
     &                         tl_zeta)
!***********************************************************************
!
      USE mod_param
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_scalars
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Linp, Lout

      real(r8), intent(in) :: alphaK
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
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
      integer :: i, j, k
#ifdef SOLVE3D
      integer :: itrc
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute new starting tangent linear state vector, X(k+1).
!-----------------------------------------------------------------------
!
!  Free-surface.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          tl_zeta(i,j,Lout)=tl_zeta(i,j,Linp)+                          &
     &                      alphaK*d_zeta(i,j)
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
          tl_ubar(i,j,Lout)=tl_ubar(i,j,Linp)+                          &
     &                      alphaK*d_ubar(i,j)
# ifdef MASKING
          tl_ubar(i,j,Lout)=tl_ubar(i,j,Lout)*umask(i,j)
# endif
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          tl_vbar(i,j,Lout)=tl_vbar(i,j,Linp)+                          &
     &                      alphaK*d_vbar(i,j)
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
            tl_ustr(i,j,k,Lout)=tl_ustr(i,j,k,Linp)+                    &
     &                          alphaK*d_sustr(i,j,k)
# ifdef MASKING
            tl_ustr(i,j,k,Lout)=tl_ustr(i,j,k,Lout)*umask(i,j)
# endif
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            tl_vstr(i,j,k,Lout)=tl_vstr(i,j,k,Linp)+                    &
     &                         alphaK*d_svstr(i,j,k)
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
            tl_u(i,j,k,Lout)=tl_u(i,j,k,Linp)+                          &
     &                       alphaK*d_u(i,j,k)
# ifdef MASKING
            tl_u(i,j,k,Lout)=tl_u(i,j,k,Lout)*umask(i,j)
# endif
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            tl_v(i,j,k,Lout)=tl_v(i,j,k,Linp)+                          &
     &                       alphaK*d_v(i,j,k)
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
              tl_t(i,j,k,Lout,itrc)=tl_t(i,j,k,Linp,itrc)+              &
     &                              alphaK*d_t(i,j,k,itrc)
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
              tl_tflux(i,j,k,Lout,itrc)=tl_tflux(i,j,k,Linp,itrc)+      &
     &                                  alphaK*d_stflx(i,j,k,itrc)
#  ifdef MASKING
              tl_tflux(i,j,k,Lout,itrc)=tl_tflux(i,j,k,Lout,itrc)*      &
     &                                  rmask(i,j)
#  endif
            END DO
          END DO
        END DO
      END DO
# endif
#endif

      RETURN
      END SUBROUTINE tl_new_state
!
!***********************************************************************
      SUBROUTINE ad_new_state (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Lold, Lnew, alphaK, tauK,                &
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
!***********************************************************************
!
      USE mod_param
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_scalars
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lold, Lnew

      real(r8), intent(in) :: alphaK, tauK
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
#endif
!
!  Local variable declarations.
!
      integer :: i, j, k
#ifdef SOLVE3D
      integer :: itrc
#endif
      real(r8) :: fac

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Estimate the gradient for the new state vector, G(k+1). Notice that
!  the Lnew record is overwritten.
!-----------------------------------------------------------------------
!
      fac=alphaK/tauK
!
!  Free-surface.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ad_zeta(i,j,Lnew)=ad_zeta(i,j,Lold)+                          &
     &                      fac*(ad_zeta(i,j,Lnew)-                     &
     &                           ad_zeta(i,j,Lold))
#ifdef MASKING
          ad_zeta(i,j,Lnew)=ad_zeta(i,j,Lnew)*rmask(i,j)
#endif
          ad_zeta(i,j,Lold)=ad_zeta(i,j,Lnew)
        END DO
      END DO

#ifndef SOLVE3D
!
!  2D momentum.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ad_ubar(i,j,Lnew)=ad_ubar(i,j,Lold)+                          &
     &                      fac*(ad_ubar(i,j,Lnew)-                     &
     &                           ad_ubar(i,j,Lold))
# ifdef MASKING
          ad_ubar(i,j,Lnew)=ad_ubar(i,j,Lnew)*umask(i,j)
# endif
          ad_ubar(i,j,Lold)=ad_ubar(i,j,Lnew)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          ad_vbar(i,j,Lnew)=ad_vbar(i,j,Lold)+                          &
     &                      fac*(ad_vbar(i,j,Lnew)-                     &
     &                           ad_vbar(i,j,Lold))
# ifdef MASKING
          ad_vbar(i,j,Lnew)=ad_vbar(i,j,Lnew)*vmask(i,j)
# endif
          ad_vbar(i,j,Lold)=ad_vbar(i,j,Lnew)
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
            ad_ustr(i,j,k,Lnew)=ad_ustr(i,j,k,Lold)+                    &
     &                          fac*(ad_ustr(i,j,k,Lnew)-               &
     &                               ad_ustr(i,j,k,Lold))
# ifdef MASKING
            ad_ustr(i,j,k,Lnew)=ad_ustr(i,j,k,Lnew)*umask(i,j)
# endif
            ad_ustr(i,j,k,Lold)=ad_ustr(i,j,k,Lnew)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_vstr(i,j,k,Lnew)=ad_vstr(i,j,k,Lold)+                    &
     &                          fac*(ad_vstr(i,j,k,Lnew)-               &
     &                               ad_vstr(i,j,k,Lold))
# ifdef MASKING
            ad_vstr(i,j,k,Lnew)=ad_vstr(i,j,k,Lnew)*vmask(i,j)
# endif
            ad_vstr(i,j,k,Lold)=ad_vstr(i,j,k,Lnew)
          END DO
        END DO
      END DO
#endif
#ifdef SOLVE3D
!
!  3D state variables.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_u(i,j,k,Lnew)=ad_u(i,j,k,Lold)+                          &
     &                       fac*(ad_u(i,j,k,Lnew)-                     &
     &                            ad_u(i,j,k,Lold))
# ifdef MASKING
            ad_u(i,j,k,Lnew)=ad_u(i,j,k,Lnew)*umask(i,j)
# endif
            ad_u(i,j,k,Lold)=ad_u(i,j,k,Lnew)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_v(i,j,k,Lnew)=ad_v(i,j,k,Lold)+                          &
     &                       fac*(ad_v(i,j,k,Lnew)-                     &
     &                            ad_v(i,j,k,Lold))
# ifdef MASKING
            ad_v(i,j,k,Lnew)=ad_v(i,j,k,Lnew)*vmask(i,j)
# endif
            ad_v(i,j,k,Lold)=ad_v(i,j,k,Lnew)
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
              ad_t(i,j,k,Lnew,itrc)=ad_t(i,j,k,Lold,itrc)+              &
     &                              fac*(ad_t(i,j,k,Lnew,itrc)-         &
     &                                   ad_t(i,j,k,Lold,itrc))
# ifdef MASKING
              ad_t(i,j,k,Lnew,itrc)=ad_t(i,j,k,Lnew,itrc)*rmask(i,j)
# endif
              ad_t(i,j,k,Lold,itrc)=ad_t(i,j,k,Lnew,itrc)
            END DO
          END DO
        END DO
      END DO

# ifdef ADJUST_STFLUX
!
!  Tracers.
!
      DO itrc=1,NT(ng)
        DO k=1,Nfrec(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              ad_tflux(i,j,k,Lnew,itrc)=ad_tflux(i,j,k,Lold,itrc)+      &
     &                                  fac*(ad_tflux(i,j,k,Lnew,itrc)- &
     &                                       ad_tflux(i,j,k,Lold,itrc))
#  ifdef MASKING
              ad_tflux(i,j,k,Lnew,itrc)=ad_tflux(i,j,k,Lnew,itrc)*      &
     &                                  rmask(i,j)
#  endif
              ad_tflux(i,j,k,Lold,itrc)=ad_tflux(i,j,k,Lnew,itrc)
            END DO
          END DO
        END DO
      END DO
# endif
#endif

      RETURN
      END SUBROUTINE ad_new_state
!
!***********************************************************************
      SUBROUTINE orthogonalize (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Lold, Lnew, Lwrk,                       &
     &                          innLoop, outLoop,                       &
#ifdef MASKING
     &                          rmask, umask, vmask,                    &
#endif
#ifdef ADJUST_WSTRESS
     &                          nl_ustr, nl_vstr,                       &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                          nl_tflux,                               &
# endif
     &                          nl_t, nl_u, nl_v,                       &
#else
     &                          nl_ubar, nl_vbar,                       &
#endif
     &                          nl_zeta,                                &
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
      real(r8), intent(inout) :: nl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: nl_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: nl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: nl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: nl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: nl_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: nl_zeta(LBi:,LBj:,:)
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
      real(r8), intent(inout) :: nl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: nl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: nl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: nl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: nl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: nl_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: nl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: nl_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: nl_zeta(LBi:UBi,LBj:UBj,3)
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
      integer :: i, j, k, lstr, rec, Lscale
#ifdef SOLVE3D
      integer :: itrc
#endif
      integer :: L1 = 1
      integer :: L2 = 2

      real(r8) :: fac, fac1, fac2

      real(r8), dimension(0:NstateVar(ng)) :: dot
      real(r8), dimension(0:Ninner) :: DotProd, dot_new, dot_old

      character (len=80) :: ncname

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Orthogonalize current gradient, G(k+1), against all previous
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
 10       FORMAT (a,'_',i3.3,'.nc')
        ELSE
          ncname=ADJname(ng)
        END IF
!
!  Read in each previous gradient state solutions, G(0) to G(k), and
!  compute its associated dot against current G(k+1). Each gradient
!  solution is loaded NLM (index L2, if preconditioning) or
!  TLM (index Lwrk, if not preconditioning) state arrays.
!
        IF (Lprecond) THEN
          CALL read_state (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     L2, rec,                                     &
     &                     ndefADJ(ng), ncADJid(ng), ncname,            &
#ifdef MASKING
     &                     rmask, umask, vmask,                         &
#endif
#ifdef ADJUST_WSTRESS
     &                     nl_ustr, nl_vstr,                            &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                     nl_tflux,                                    &
# endif
     &                     nl_t, nl_u, nl_v,                            &
#else
     &                     nl_ubar, nl_vbar,                            &
#endif
     &                     nl_zeta)
!
        ELSE
          CALL read_state (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Lwrk, rec,                                   &
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
        END IF
!
!  If preconditioning, compute  H^-1 * G(rec) and store it TLM state
!  arrays (index Lwrk).
!
        IF (Lprecond) THEN
          Lscale=-1
          CALL precond (ng, tile, model,                                &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  NstateVar(ng), L2, Lwrk, Lscale,                &
     &                  nConvRitz, Ritz,                                &
#ifdef MASKING
     &                  rmask, umask, vmask,                            &
#endif
#ifdef ADJUST_WSTRESS
     &                  nl_ustr, nl_ustr, tl_ustr,                      &
     &                  nl_vstr, nl_vstr, tl_vstr,                      &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                  nl_tflux, nl_tflux, tl_tflux,                   &
# endif
     &                  nl_t, nl_t, tl_t,                               &
     &                  nl_u, nl_u, tl_u,                               &
     &                  nl_v, nl_v, tl_v,                               &
#else
     &                  nl_ubar, nl_ubar, tl_ubar,                      &
     &                  nl_vbar, nl_vbar, tl_vbar,                      &
#endif
     &                  nl_zeta, nl_zeta, tl_zeta)
        END IF
!
!  If preconditioning, compute dot product <G(k+1), H^-1 G(rec)>.
!  Otherwise, compute <G(k+1), G(rec)>. Recall that the TLM
!  (index Lwrk) contains either H^-1 G(rec) or G(rec).
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
        dot_new(rec)=dot(0)
!
!  If preconditioning, compute dot product <G(rec), H^-1 * G(rec)>.
!
        IF (Lprecond) THEN
          CALL state_dotprod (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NstateVar(ng), dot(0:),                   &
#ifdef MASKING
     &                        rmask, umask, vmask,                      &
#endif
#ifdef ADJUST_WSTRESS
     &                        nl_ustr(:,:,:,L2), tl_ustr(:,:,:,Lwrk),   &
     &                        nl_vstr(:,:,:,L2), tl_vstr(:,:,:,Lwrk),   &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                        nl_tflux(:,:,:,L2,:),                     &
     &                        tl_tflux(:,:,:,Lwrk,:),                   &
# endif
     &                        nl_t(:,:,:,L2,:), tl_t(:,:,:,Lwrk,:),     &
     &                        nl_u(:,:,:,L2), tl_u(:,:,:,Lwrk),         &
     &                        nl_v(:,:,:,L2), tl_v(:,:,:,Lwrk),         &
#else
     &                        nl_ubar(:,:,L2), tl_ubar(:,:,Lwrk),       &
     &                        nl_vbar(:,:,L2), tl_vbar(:,:,Lwrk),       &
#endif
     &                        nl_zeta(:,:,L2), tl_zeta(:,:,Lwrk))
!
!  Otherwise, compute dot product <G(rec), G(rec)>.
!
        ELSE
          CALL state_dotprod (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NstateVar(ng), dot(0:),                   &
#ifdef MASKING
     &                        rmask, umask, vmask,                      &
#endif
#ifdef ADJUST_WSTRESS
     &                        tl_ustr(:,:,:,Lwrk), tl_ustr(:,:,:,Lwrk), &
     &                        tl_vstr(:,:,:,Lwrk), tl_vstr(:,:,:,Lwrk), &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                        tl_tflux(:,:,:,Lwrk,:),                   &
     &                        tl_tflux(:,:,:,Lwrk,:),                   &
# endif
     &                        tl_t(:,:,:,Lwrk,:), tl_t(:,:,:,Lwrk,:),   &
     &                        tl_u(:,:,:,Lwrk), tl_u(:,:,:,Lwrk),       &
     &                        tl_v(:,:,:,Lwrk), tl_v(:,:,:,Lwrk),       &
#else
     &                        tl_ubar(:,:,Lwrk), tl_ubar(:,:,Lwrk),     &
     &                        tl_vbar(:,:,Lwrk), tl_vbar(:,:,Lwrk),     &
#endif
     &                        tl_zeta(:,:,Lwrk), tl_zeta(:,:,Lwrk))
        END IF
        dot_old(rec)=dot(0)
!
!  Compute Gramm-Schmidt scaling coefficient.
!
        DotProd(rec)=dot_new(rec)/dot_old(rec)

        fac1=1.0_r8
        fac2=-DotProd(rec)
!
!  If preconditioning, perform Gramm-Schmidt orthonormalization as:
!
!    ad_var(Lnew) = fac1 * ad_var(Lnew) + fac2 * nl_var(L2)
!
        IF (Lprecond) THEN
          CALL state_addition (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Lnew, L2, Lnew, fac1, fac2,              &
#ifdef MASKING
     &                         rmask, umask, vmask,                     &
#endif
#ifdef ADJUST_WSTRESS
     &                         ad_ustr, nl_ustr,                        &
     &                         ad_vstr, nl_vstr,                        &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                         ad_tflux, nl_tflux,                      &
# endif
     &                         ad_t, nl_t,                              &
     &                         ad_u, nl_u,                              &
     &                         ad_v, nl_v,                              &
#else
     &                         ad_ubar, nl_ubar,                        &
     &                         ad_vbar, nl_vbar,                        &
#endif
     &                         ad_zeta, nl_zeta)
!
!  If not preconditioning, perform Gramm-Schmidt orthonormalization as:
!
!    ad_var(Lnew) = fac1 * ad_var(Lnew) + fac2 * tl_var(Lwrk)
!
        ELSE
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
        END IF
      END DO

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
!  Read in each previous gradient state solutions, G(0) to G(k), and
!  compute its associated dot angaint orthogonalized G(k+1). Again,
!  each gradient solution is loaded into TANGENT LINEAR STATE ARRAYS
!  at index Lwrk.
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
        dot_new(rec)=dot(0)
      END DO
!
!  Report dot products. If everything is working correctly, at the
!  end of the orthogonalization dot_new(rec) << dot_old(rec).
!
      IF (Master) THEN
        WRITE (stdout,20) outer, inner
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
      END SUBROUTINE orthogonalize

!
!***********************************************************************
      SUBROUTINE new_direction (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Lwrk, Lnew, betaK,                      &
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
      integer, intent(in) :: Lwrk, Lnew

      real(r8), intent(in) :: betaK
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
!  descent direction is overwritten. Also the initial value is just
!  d(0)=-G(0) since betaK=0 when inner=0.
!-----------------------------------------------------------------------
!
!  Free-surface.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          d_zeta(i,j)=-ad_zeta(i,j,Lnew)+betaK*d_zeta(i,j)
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
          d_ubar(i,j)=-ad_ubar(i,j,Lnew)+betaK*d_ubar(i,j)
# ifdef MASKING
          d_ubar(i,j)=d_ubar(i,j)*umask(i,j)
# endif
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          d_vbar(i,j)=-ad_vbar(i,j,Lnew)+betaK*d_vbar(i,j)
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
            d_sustr(i,j,k)=-ad_ustr(i,j,k,Lnew)+betaK*d_sustr(i,j,k)
# ifdef MASKING
            d_sustr(i,j,k)=d_sustr(i,j,k)*umask(i,j)
# endif
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            d_svstr(i,j,k)=-ad_vstr(i,j,k,Lnew)+betaK*d_svstr(i,j,k)
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
            d_u(i,j,k)=-ad_u(i,j,k,Lnew)+betaK*d_u(i,j,k)
# ifdef MASKING
            d_u(i,j,k)=d_u(i,j,k)*umask(i,j)
# endif
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            d_v(i,j,k)=-ad_v(i,j,k,Lnew)+betaK*d_v(i,j,k)
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
              d_t(i,j,k,itrc)=-ad_t(i,j,k,Lnew,itrc)+                   &
     &                        betaK*d_t(i,j,k,itrc)
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
              d_stflx(i,j,k,itrc)=-ad_tflux(i,j,k,Lnew,itrc)+           &
     &                            betaK*d_stflx(i,j,k,itrc)
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
!**********************************************************************
      SUBROUTINE precond (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NstateVars, Linp, Lwrk, Lscale,               &
     &                    nConvRitz, Ritz,                              &
#ifdef MASKING
     &                    rmask, umask, vmask,                          &
#endif
#ifdef ADJUST_WSTRESS
     &                    ad_ustr, nl_ustr, tl_ustr,                    &
     &                    ad_vstr, nl_vstr, tl_vstr,                    &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                    ad_tflux, nl_tflux, tl_tflux,                 &
# endif
     &                    ad_t, nl_t, tl_t,                             &
     &                    ad_u, nl_u, tl_u,                             &
     &                    ad_v, nl_v, tl_v,                             &
#else
     &                    ad_ubar, nl_ubar, tl_ubar,                    &
     &                    ad_vbar, nl_vbar, tl_vbar,                    &
#endif
     &                    ad_zeta, nl_zeta, tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_netcdf
      USE mod_iounits
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_scalars
#endif
!
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_reduce
#endif
      USE state_addition_mod, ONLY : state_addition
      USE state_copy_mod, ONLY : state_copy
      USE state_dotprod_mod, ONLY : state_dotprod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: NstateVars, Linp, Lwrk, Lscale
      integer, intent(in) :: nConvRitz
!
      real(r8), intent(in) :: Ritz(:)
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
#ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_vstr(LBi:,LBj:,:,:)
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
      real(r8), intent(in) :: ad_tflux(LBi:,LBj:,:,:,:)
# endif
      real(r8), intent(in) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_v(LBi:,LBj:,:,:)
#else
      real(r8), intent(in) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_vbar(LBi:,LBj:,:)
#endif
      real(r8), intent(in) :: ad_zeta(LBi:,LBj:,:)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: nl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: nl_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: nl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: nl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: nl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: nl_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: nl_zeta(LBi:,LBj:,:)
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
#ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(in) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
      real(r8), intent(in) :: ad_tflux(LBi:UBi,LBj:UBj,                 &
     &                                 Nfrec(ng),2,NT(ng))
# endif
      real(r8), intent(in) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#else
      real(r8), intent(in) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#endif
      real(r8), intent(in) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: nl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: nl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: nl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: nl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: nl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: nl_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: nl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: nl_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: nl_zeta(LBi:UBi,LBj:UBj,3)
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
      integer :: NSUB, i, j, k, L1, L2, nvec
#ifdef SOLVE3D
      integer :: itrc
#endif
      real(r8) :: cff, fac, fac1, fac2
      real(r8), dimension(0:NstateVars) :: Dotprod
#ifdef DISTRIBUTE
      character (len=3) :: op_handle
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Apply the preconditioner. The approximated Hessian matrix is computed
!  from the eigenvectors computed by the Lanczos algorithm which are
!  stored in HSSname NetCDF file.
!-----------------------------------------------------------------------
!
!  Copy ad_var(Linp) into tl_var(Lwrk)
!
      CALL state_copy (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 Linp, Lwrk,                                      &
#ifdef ADJUST_WSTRESS
     &                 tl_ustr, ad_ustr,                                &
     &                 tl_vstr, ad_vstr,                                &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                 tl_tflux, ad_tflux,                              &
# endif
     &                 tl_t, ad_t,                                      &
     &                 tl_u, ad_u,                                      &
     &                 tl_v, ad_v,                                      &
#else
     &                 tl_ubar, ad_ubar,                                &
     &                 tl_vbar, ad_vbar,                                &
#endif
     &                 tl_zeta, ad_zeta)
!
!  Read the converged Hessian eigenvectors into NLM state array,
!  index L1.
!
      DO nvec=1,nConvRitz
        L1=1
        L2=2
        CALL read_state (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   L1, nvec,                                      &
     &                   0, ncHSSid(ng), HSSname(ng),                   &
#ifdef MASKING
     &                   rmask, umask, vmask,                           &
#endif
#ifdef ADJUST_WSTRESS
     &                   nl_ustr, nl_vstr,                              &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                   nl_tflux,                                      &
# endif
     &                   nl_t, nl_u, nl_v,                              &
#else
     &                   nl_ubar, nl_vbar,                              &
#endif
     &                   nl_zeta)
!
!  Compute dot product between gradient and Hessian eigenvector.
!
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVars, Dotprod(0:),                    &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef ADJUST_WSTRESS
     &                      ad_ustr(:,:,:,Linp), nl_ustr(:,:,:,L1),     &
     &                      ad_vstr(:,:,:,Linp), nl_vstr(:,:,:,L1),     &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      ad_tflux(:,:,:,Linp,:),                     &
                            nl_tflux(:,:,:,L1,:),                       &
# endif
     &                      ad_t(:,:,:,Linp,:), nl_t(:,:,:,L1,:),       &
     &                      ad_u(:,:,:,Linp), nl_u(:,:,:,L1),           &
     &                      ad_v(:,:,:,Linp), nl_v(:,:,:,L1),           &
#else
     &                      ad_ubar(:,:,Linp), nl_ubar(:,:,L1),         &
     &                      ad_vbar(:,:,Linp), nl_vbar(:,:,L1),         &
#endif
     &                      ad_zeta(:,:,Linp), nl_zeta(:,:,L1))
!
!    Lscale determines the form of the preconditioner:
!
!       1= Hessian
!      -1= Inverse Hessian
!       2= Hessian square root
!      -2= Inverse Hessian square root
!
!    tl_var(Lwrk) = fac1 * tl_var(Lwrk) + fac2 * nl_var(L1)
!
        fac1=1.0_r8

        IF (Lscale.eq.1) THEN
          fac2=(Ritz(nvec)-1.0_r8)*Dotprod(0)
        ELSE IF (Lscale.eq.-1) THEN
          fac2=(1.0_r8/Ritz(nvec)-1.0_r8)*Dotprod(0)
        ELSE IF (Lscale.eq.2) THEN
          fac2=(SQRT(Ritz(nvec))-1.0_r8)*Dotprod(0)
        ELSE IF (Lscale.eq.-2) THEN
          fac2=(1.0_r8/SQRT(Ritz(nvec))-1.0_r8)*Dotprod(0)
        END IF

        CALL state_addition (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Lwrk, L1, Lwrk, fac1, fac2,                &
#ifdef MASKING
     &                       rmask, umask, vmask,                       &
#endif
#ifdef ADJUST_WSTRESS
     &                       tl_ustr, nl_ustr,                          &
     &                       tl_vstr, nl_vstr,                          &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                       tl_tflux, nl_tflux,                        &
# endif
     &                       tl_t, nl_t,                                &
     &                       tl_u, nl_u,                                &
     &                       tl_v, nl_v,                                &
#else
     &                       tl_ubar, nl_ubar,                          &
     &                       tl_vbar, nl_vbar,                          &
#endif
     &                       tl_zeta, nl_zeta)
      END DO

      RETURN
      END SUBROUTINE precond
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
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lwrk, rec, ndef, ncfileid

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
      integer :: i, j
#ifdef SOLVE3D
      integer :: itrc, k
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
            RETURN
          END IF
        END IF
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
      IF (ndef.gt.0) THEN
        status=nf90_close(ncid)
      END IF
!
 10   FORMAT (' READ_STATE - unable to open NetCDF file: ',a)
 20   FORMAT (' READ_STATE - error while reading variable: ',a,2x,      &
     &        'at time record = ',i3,/,16x,'in NetCDF file: ',a)

      RETURN
      END SUBROUTINE read_state

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
      integer, save :: varid(7)
!
!-----------------------------------------------------------------------
!  Write out conjugate gradient vectors.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
        IF (First) THEN
          First=.FALSE.
          DO i=1,7
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
        IF ((innLoop.eq.0).and.(outloop.eq.1)) THEN
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
        IF ((innLoop.eq.0).and.(outloop.eq.1)) THEN
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
        IF (varid(5).eq.0) THEN
          status=nf90_inq_varid(ncMODid(ng), 'cg_alpha', varid(5))
        END IF
        start(1)=1
        total(1)=Ninner+1
        start(2)=1
        total(2)=Nouter
        status=nf90_put_var(ncMODid(ng), varid(5), cg_alpha(0:,:),      &
     &                      start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'cg_alpha', TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        IF (varid(6).eq.0) THEN
          status=nf90_inq_varid(ncMODid(ng), 'cg_beta', varid(6))
        END IF
        start(1)=1
        total(1)=Ninner+1
        start(2)=1
        total(2)=Nouter
        status=nf90_put_var(ncMODid(ng), varid(6), cg_beta(0:,:),       &
     &                      start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'cg_beta', TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        IF (varid(7).eq.0) THEN
          status=nf90_inq_varid(ncMODid(ng), 'cg_tau', varid(7))
        END IF
        start(1)=1
        total(1)=Ninner+1
        start(2)=1
        total(2)=Nouter
        status=nf90_put_var(ncMODid(ng), varid(7), cg_tau(0:,:),        &
     &                      start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'cg_tau', TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
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
