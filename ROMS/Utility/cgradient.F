#include "cppdefs.h"

      MODULE cgradient_mod

#ifdef IS4DVAR
!
!svn $Id: cgradient.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module minimizes  incremental  4Dvar  quadratic cost function  !
!  using a preconditioned version of the conjugate gradient algorithm  !
!  proposed by Mike Fisher (ECMWF) and modified  by  Tshimanga et al.  !
!  (2008) using Limited-Memory Precondtioners (LMP).                   !
!                                                                      !
!  In the following,  M represents the preconditioner.  Specifically,  !
!                                                                      !
!    M = I + SUM_i [ (mu_i-1) h_i transpose(h_i)],                     !
!                                                                      !
!  where mu_i can take the following values:                           !
!                                                                      !
!    Lscale=-1:    mu_i = lambda_i                                     !
!    Lscale= 1:    mu_i = 1 / lambda_i                                 !
!    Lscale=-2:    mu_i = SQRT (lambda_i)                              !
!    Lscale= 2:    mu_i = 1 / SQRT(lambda_i)                           !
!                                                                      !
!  where lambda_i are the Hessian eigenvalues and h_i are the Hessian  !
!  eigenvectors.                                                       !
!                                                                      !
!  For Lscale= 1 spectral LMP is used as the preconditioner.           !
!  For Lscale=-1 inverse spectral LMP is used as the preconditioner.   !
!  For Lscale= 2 SQRT spectral LMP is used as the preconditioner.      !
!  For Lscale=-2 inverse SQRT spectral LMP is the preconditioner.      !
!                                                                      !
!  If Lritz=.TRUE. then Ritz LMP is used and the expressions for mu_i  !
!  are more complicated.                                               !
!                                                                      !
!  For some operations the tranpose of the preconditioner is required. !
!  For spectral LMP the preconditioner and its tranpose are identical. !
!  For Ritz LMP the preconditioner and its tranpose differ.            !
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
!    H Q(k) = Q(k) T(k) + Gamma(k) q(k+1) transpose[e(k)]              !
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
!    transpose[e(k)] = { 0, ...,0, 1 }                                 !
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
!    T(k) Z(k) = - transpose[Q(k)] g(0)                                !
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
!    e(k)        Tansposed unit vector                                 !
!    g(k)        Gradient vectors (adjoint solution: GRAD(J))          !
!    q(k)        Lanczos vectors                                       !
!    <...>       Dot product                                           !
!    ||...||     Euclidean norm, ||g(k)|| = SQRT( <g(k),g(k)> )        !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Fisher, M., 1997: Efficient Minimization of Quadratic Penalty     !
!      funtions, unpublish manuscript, 1-14.                           !
!                                                                      !
!    Fisher, M., 1998: Minimization algorithms for variational data    !
!      data assimilation. In Recent Developments in Numerical          !
!      Methods for Atmospheric Modelling, pp 364-385, ECMWF.           !
!                                                                      !
!    Tchimanga, J., S. Gratton, A.T. Weaver, and A. Sartenaer, 2008:   !
!      Limited-memory preconditioners, with application to incremental !
!      four-dimensional variational ocean data assimilation, Q.J.R.    !
!      Meteorol. Soc., 134, 753-771.                                   !
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
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
# ifdef SOLVE3D
      USE mod_coupling
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
# endif
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
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 36)
# endif
      CALL cgradient_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj, LBij, UBij,              &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     Lold(ng), Lnew(ng),                          &
     &                     innLoop, outLoop,                            &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                     BOUNDARY(ng) % t_obc,                        &
     &                     BOUNDARY(ng) % u_obc,                        &
     &                     BOUNDARY(ng) % v_obc,                        &
#  endif
     &                     BOUNDARY(ng) % ubar_obc,                     &
     &                     BOUNDARY(ng) % vbar_obc,                     &
     &                     BOUNDARY(ng) % zeta_obc,                     &
# endif
# ifdef ADJUST_WSTRESS
     &                     FORCES(ng) % ustr,                           &
     &                     FORCES(ng) % vstr,                           &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     FORCES(ng) % tflux,                          &
#  endif
     &                     OCEAN(ng) % t,                               &
     &                     OCEAN(ng) % u,                               &
     &                     OCEAN(ng) % v,                               &
# else
     &                     OCEAN(ng) % ubar,                            &
     &                     OCEAN(ng) % vbar,                            &
# endif
     &                     OCEAN(ng) % zeta,                            &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                     BOUNDARY(ng) % tl_t_obc,                     &
     &                     BOUNDARY(ng) % tl_u_obc,                     &
     &                     BOUNDARY(ng) % tl_v_obc,                     &
#  endif
     &                     BOUNDARY(ng) % tl_ubar_obc,                  &
     &                     BOUNDARY(ng) % tl_vbar_obc,                  &
     &                     BOUNDARY(ng) % tl_zeta_obc,                  &
# endif
# ifdef ADJUST_WSTRESS
     &                     FORCES(ng) % tl_ustr,                        &
     &                     FORCES(ng) % tl_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     FORCES(ng) % tl_tflux,                       &
#  endif
     &                     OCEAN(ng) % tl_t,                            &
     &                     OCEAN(ng) % tl_u,                            &
     &                     OCEAN(ng) % tl_v,                            &
# else
     &                     OCEAN(ng) % tl_ubar,                         &
     &                     OCEAN(ng) % tl_vbar,                         &
# endif
     &                     OCEAN(ng) % tl_zeta,                         &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                     BOUNDARY(ng) % d_t_obc,                      &
     &                     BOUNDARY(ng) % d_u_obc,                      &
     &                     BOUNDARY(ng) % d_v_obc,                      &
#  endif
     &                     BOUNDARY(ng) % d_ubar_obc,                   &
     &                     BOUNDARY(ng) % d_vbar_obc,                   &
     &                     BOUNDARY(ng) % d_zeta_obc,                   &
# endif
# ifdef ADJUST_WSTRESS
     &                     FORCES(ng) % d_sustr,                        &
     &                     FORCES(ng) % d_svstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     FORCES(ng) % d_stflx,                        &
#  endif
     &                     OCEAN(ng) % d_t,                             &
     &                     OCEAN(ng) % d_u,                             &
     &                     OCEAN(ng) % d_v,                             &
# else
     &                     OCEAN(ng) % d_ubar,                          &
     &                     OCEAN(ng) % d_vbar,                          &
# endif
     &                     OCEAN(ng) % d_zeta,                          &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                     BOUNDARY(ng) % ad_t_obc,                     &
     &                     BOUNDARY(ng) % ad_u_obc,                     &
     &                     BOUNDARY(ng) % ad_v_obc,                     &
#  endif
     &                     BOUNDARY(ng) % ad_ubar_obc,                  &
     &                     BOUNDARY(ng) % ad_vbar_obc,                  &
     &                     BOUNDARY(ng) % ad_zeta_obc,                  &
# endif
# ifdef ADJUST_WSTRESS
     &                     FORCES(ng) % ad_ustr,                        &
     &                     FORCES(ng) % ad_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     FORCES(ng) % ad_tflux,                       &
#  endif
     &                     OCEAN(ng) % ad_t,                            &
     &                     OCEAN(ng) % ad_u,                            &
     &                     OCEAN(ng) % ad_v,                            &
# else
     &                     OCEAN(ng) % ad_ubar,                         &
     &                     OCEAN(ng) % ad_vbar,                         &
# endif
     &                     OCEAN(ng) % ad_zeta)
# ifdef PROFILE
      CALL wclock_on (ng, model, 36)
# endif
      RETURN
      END SUBROUTINE cgradient
!
!***********************************************************************
      SUBROUTINE cgradient_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj, LBij, UBij,        &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           Lold, Lnew,                            &
     &                           innLoop, outLoop,                      &
# ifdef MASKING
     &                           rmask, umask, vmask,                   &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                           nl_t_obc, nl_u_obc, nl_v_obc,          &
#  endif
     &                           nl_ubar_obc, nl_vbar_obc,              &
     &                           nl_zeta_obc,                           &
# endif
# ifdef ADJUST_WSTRESS
     &                           nl_ustr, nl_vstr,                      &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                           nl_tflux,                              &
#  endif
     &                           nl_t, nl_u, nl_v,                      &
# else
     &                           nl_ubar, nl_vbar,                      &
# endif
     &                           nl_zeta,                               &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                           tl_t_obc, tl_u_obc, tl_v_obc,          &
#  endif
     &                           tl_ubar_obc, tl_vbar_obc,              &
     &                           tl_zeta_obc,                           &
# endif
# ifdef ADJUST_WSTRESS
     &                           tl_ustr, tl_vstr,                      &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                           tl_tflux,                              &
#  endif
     &                           tl_t, tl_u, tl_v,                      &
# else
     &                           tl_ubar, tl_vbar,                      &
# endif
     &                           tl_zeta,                               &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                           d_t_obc, d_u_obc, d_v_obc,             &
#  endif
     &                           d_ubar_obc, d_vbar_obc,                &
     &                           d_zeta_obc,                            &
# endif
# ifdef ADJUST_WSTRESS
     &                           d_sustr, d_svstr,                      &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                           d_stflx,                               &
#  endif
     &                           d_t, d_u, d_v,                         &
# else
     &                           d_ubar, d_vbar,                        &
# endif
     &                           d_zeta,                                &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                           ad_t_obc, ad_u_obc, ad_v_obc,          &
#  endif
     &                           ad_ubar_obc, ad_vbar_obc,              &
     &                           ad_zeta_obc,                           &
# endif
# ifdef ADJUST_WSTRESS
     &                           ad_ustr, ad_vstr,                      &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                           ad_tflux,                              &
#  endif
     &                           ad_t, ad_u, ad_v,                      &
# else
     &                           ad_ubar, ad_vbar,                      &
# endif
     &                           ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
# endif
      USE state_copy_mod, ONLY : state_copy
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lold, Lnew
      integer, intent(in) :: innLoop, outLoop
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: ad_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: ad_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: d_t_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: d_u_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: d_v_obc(LBij:,:,:,:)
#   endif
      real(r8), intent(inout) :: d_ubar_obc(LBij:,:,:)
      real(r8), intent(inout) :: d_vbar_obc(LBij:,:,:)
      real(r8), intent(inout) :: d_zeta_obc(LBij:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: d_sustr(LBi:,LBj:,:)
      real(r8), intent(inout) :: d_svstr(LBi:,LBj:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: d_stflx(LBi:,LBj:,:,:)
#   endif
      real(r8), intent(inout) :: d_t(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: d_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: d_v(LBi:,LBj:,:)
#  else
      real(r8), intent(inout) :: d_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: d_vbar(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: d_zeta(LBi:,LBj:)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: nl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: nl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: nl_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: nl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: nl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: nl_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: nl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: nl_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: nl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: nl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: nl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: nl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: nl_zeta(LBi:,LBj:,:)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: tl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: tl_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: ad_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: d_t_obc(LBij:UBij,N(ng),4,             &
     &                                   Nbrec(ng),NT(ng))
      real(r8), intent(inout) :: d_u_obc(LBij:UBij,N(ng),4,Nbrec(ng))
      real(r8), intent(inout) :: d_v_obc(LBij:UBij,N(ng),4,Nbrec(ng))
#   endif
      real(r8), intent(inout) :: d_ubar_obc(LBij:UBij,4,Nbrec(ng))
      real(r8), intent(inout) :: d_vbar_obc(LBij:UBij,4,Nbrec(ng))
      real(r8), intent(inout) :: d_zeta_obc(LBij:UBij,4,Nbrec(ng))
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: d_sustr(LBi:UBi,LBj:UBj,Nfrec(ng))
      real(r8), intent(inout) :: d_svstr(LBi:UBi,LBj:UBj,Nfrec(ng))
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: d_stflx(LBi:UBi,LBj:UBj,               &
     &                                   Nfrec(ng),NT(ng))
#   endif
      real(r8), intent(inout) :: d_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: d_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: d_v(LBi:UBi,LBj:UBj,N(ng))
#  else
      real(r8), intent(inout) :: d_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: d_vbar(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: d_zeta(LBi:UBi,LBj:UBj)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: nl_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: nl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: nl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: nl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: nl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: nl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: nl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: nl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: nl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: nl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: nl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: nl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: nl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: nl_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: nl_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: tl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      logical :: Ltrans

      integer :: L1 = 1
      integer :: L2 = 2

      integer :: Linp, Lout, Lscale, Lwrk, Lwrk1, i, j, ic
      integer :: info, itheta1

      real(r8) :: norm, zbeta, ztheta1

      real(r8), dimension(2*Ninner-2) :: work

      character (len=13) :: string

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize trial step size.
!-----------------------------------------------------------------------
!
      Ltrans=.FALSE.
!
      IF (Master) WRITE (stdout,10)
 10   FORMAT (/,' <<<< Descent Algorithm >>>>',/)
!
!  If preconditioning, convert the total gradient ad_var(L2)
!  from v-space to y-space.
!
      IF (Lprecond.and.(outLoop.gt.1)) THEN

        Lscale=2                 ! SQRT spectral LMP
        Ltrans=.TRUE.
!
!  Copy ad_var(L2) into nl_var(L1)
!
        CALL state_copy (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   L2, L1,                                        &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   nl_t_obc, ad_t_obc,                            &
     &                   nl_u_obc, ad_u_obc,                            &
     &                   nl_v_obc, ad_v_obc,                            &
#  endif
     &                   nl_ubar_obc, ad_ubar_obc,                      &
     &                   nl_vbar_obc, ad_vbar_obc,                      &
     &                   nl_zeta_obc, ad_zeta_obc,                      &
# endif
# ifdef ADJUST_WSTRESS
     &                   nl_ustr, ad_ustr,                              &
     &                   nl_vstr, ad_vstr,                              &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                   nl_tflux, ad_tflux,                            &
#  endif
     &                   nl_t, ad_t,                                    &
     &                   nl_u, ad_u,                                    &
     &                   nl_v, ad_v,                                    &
# else
     &                   nl_ubar, ad_ubar,                              &
     &                   nl_vbar, ad_vbar,                              &
# endif
     &                   nl_zeta, ad_zeta)
!
        CALL precond (ng, tile, model, 'convert gradient to y-space',   &
     &                LBi, UBi, LBj, UBj, LBij, UBij,                   &
     &                IminS, ImaxS, JminS, JmaxS,                       &
     &                NstateVar(ng), Lscale, Ltrans,                    &
     &                innLoop, outLoop,                                 &
# ifdef MASKING
     &                rmask, umask, vmask,                              &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                nl_t_obc, nl_u_obc, nl_v_obc,                     &
#  endif
     &                nl_ubar_obc, nl_vbar_obc,                         &
     &                nl_zeta_obc,                                      &
# endif
# ifdef ADJUST_WSTRESS
     &                nl_ustr, nl_vstr,                                 &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                nl_tflux,                                         &
#  endif
     &                nl_t, nl_u, nl_v,                                 &
# else
     &                nl_ubar, nl_vbar,                                 &
# endif
     &                nl_zeta)
        IF (exit_flag.ne.NoError) RETURN
!
!  Copy nl_var(L1) into ad_var(L2).
!
        CALL state_copy (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   L1, L2,                                        &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   ad_t_obc, nl_t_obc,                            &
     &                   ad_u_obc, nl_u_obc,                            &
     &                   ad_v_obc, nl_v_obc,                            &
#  endif
     &                   ad_ubar_obc, nl_ubar_obc,                      &
     &                   ad_vbar_obc, nl_vbar_obc,                      &
     &                   ad_zeta_obc, nl_zeta_obc,                      &
# endif
# ifdef ADJUST_WSTRESS
     &                   ad_ustr, nl_ustr,                              &
     &                   ad_vstr, nl_vstr,                              &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                   ad_tflux, nl_tflux,                            &
#  endif
     &                   ad_t, nl_t,                                    &
     &                   ad_u, nl_u,                                    &
     &                   ad_v, nl_v,                                    &
# else
     &                   ad_ubar, nl_ubar,                              &
     &                   ad_vbar, nl_vbar,                              &
# endif
     &                   ad_zeta, nl_zeta)
!
      END IF
!
!  Estimate the Hessian. Note that ad_var(Lold) will be in y-space
!  already if preconditioning since all of the Lanczos vectors saved
!  in the ADM(ng)%name file will be in y-space.
!
      IF (innLoop.gt.0) THEN
        Lwrk=2
        Linp=1
        Lout=2
        CALL hessian (ng, tile, model,                                  &
     &                LBi, UBi, LBj, UBj, LBij, UBij,                   &
     &                IminS, ImaxS, JminS, JmaxS,                       &
     &                Linp, Lout, Lwrk,                                 &
     &                innLoop, outLoop,                                 &
# ifdef MASKING
     &                rmask, umask, vmask,                              &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                ad_t_obc, ad_u_obc, ad_v_obc,                     &
#  endif
     &                ad_ubar_obc, ad_vbar_obc,                         &
     &                ad_zeta_obc,                                      &
# endif
# ifdef ADJUST_WSTRESS
     &                ad_ustr, ad_vstr,                                 &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                ad_tflux,                                         &
#  endif
     &                ad_t, ad_u, ad_v,                                 &
# else
     &                ad_ubar, ad_vbar,                                 &
# endif
     &                ad_zeta,                                          &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                tl_t_obc, tl_u_obc, tl_v_obc,                     &
#  endif
     &                tl_ubar_obc, tl_vbar_obc,                         &
     &                tl_zeta_obc,                                      &
# endif
# ifdef ADJUST_WSTRESS
     &                tl_ustr, tl_vstr,                                 &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                tl_tflux,                                         &
#  endif
     &                tl_t, tl_u, tl_v,                                 &
# else
     &                tl_ubar, tl_vbar,                                 &
# endif
     &                tl_zeta)
        IF (exit_flag.ne.NoError) RETURN
!
!  Check for positive Hessian of J.
!
        IF (cg_delta(innLoop,outLoop).le.0.0_r8) THEN
          PRINT *, ' CG_DELTA not positive.'
          PRINT *, ' CG_DELTA = ', cg_delta(innLoop,outLoop),           &
     &             ', outer = ', outLoop, ', inner = ', innLoop
          exit_flag=8
          RETURN
        END IF
      END IF
!
!  Apply the Lanczos recurrence and orthonormalize.
!  If preconditioning, the Lanczos recursion relation is identical
!  in v-space and y-space, and all ad_var are in y-space already.
!
      Linp=1
      Lout=2
      Lwrk=2
      CALL lanczos (ng, tile, model,                                    &
     &              LBi, UBi, LBj, UBj, LBij, UBij,                     &
     &              IminS, ImaxS, JminS, JmaxS,                         &
     &              Linp, Lout, Lwrk,                                   &
     &              innLoop, outLoop,                                   &
# ifdef MASKING
     &              rmask, umask, vmask,                                &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &              tl_t_obc, tl_u_obc, tl_v_obc,                       &
#  endif
     &              tl_ubar_obc, tl_vbar_obc,                           &
     &              tl_zeta_obc,                                        &
# endif
# ifdef ADJUST_WSTRESS
     &              tl_ustr, tl_vstr,                                   &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &              tl_tflux,                                           &
#  endif
     &              tl_t, tl_u, tl_v,                                   &
# else
     &              tl_ubar, tl_vbar,                                   &
# endif
     &              tl_zeta,                                            &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &              ad_t_obc, ad_u_obc, ad_v_obc,                       &
#  endif
     &              ad_ubar_obc, ad_vbar_obc,                           &
     &              ad_zeta_obc,                                        &
# endif
# ifdef ADJUST_WSTRESS
     &              ad_ustr, ad_vstr,                                   &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &              ad_tflux,                                           &
#  endif
     &              ad_t, ad_u, ad_v,                                   &
# else
     &              ad_ubar, ad_vbar,                                   &
# endif
     &              ad_zeta)
      IF (exit_flag.ne.NoError) RETURN
!
!  Compute new direction, d(k+1).
!
      CALL new_direction (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj, LBij, UBij,               &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    Linp, Lout,                                   &
# ifdef MASKING
     &                    rmask, umask, vmask,                          &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    ad_t_obc, ad_u_obc, ad_v_obc,                 &
#  endif
     &                    ad_ubar_obc, ad_vbar_obc,                     &
     &                    ad_zeta_obc,                                  &
# endif
# ifdef ADJUST_WSTRESS
     &                    ad_ustr, ad_vstr,                             &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    ad_tflux,                                     &
#  endif
     &                    ad_t, ad_u, ad_v,                             &
# else
     &                    ad_ubar, ad_vbar,                             &
# endif
     &                    ad_zeta,                                      &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    d_t_obc, d_u_obc, d_v_obc,                    &
#  endif
     &                    d_ubar_obc, d_vbar_obc,                       &
     &                    d_zeta_obc,                                   &
# endif
# ifdef ADJUST_WSTRESS
     &                    d_sustr, d_svstr,                             &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    d_stflx,                                      &
#  endif
     &                    d_t, d_u, d_v,                                &
# else
     &                    d_ubar, d_vbar,                               &
# endif
     &                    d_zeta)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Calculate the reduction in the gradient norm by solving a
!  tridiagonal system.
!-----------------------------------------------------------------------
!
!  Decomposition and forward substitution.
!
      IF (innLoop.gt.0) THEN
        zbeta=cg_delta(1,outLoop)
        cg_zu(1,outLoop)=-cg_QG(1,outLoop)/zbeta
      END IF
!
      IF (innLoop.gt.1) THEN
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
!  Compute gradient norm using ad_var(:,:,1) and tl_var(:,:,2) as
!  temporary storage.
!
        Linp=1
        Lout=2
        Lwrk=2
        CALL new_gradient (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj, LBij, UBij,              &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     Linp, Lout, Lwrk,                            &
     &                     innLoop, outLoop,                            &
# ifdef MASKING
     &                     rmask, umask, vmask,                         &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                     tl_t_obc, tl_u_obc, tl_v_obc,                &
#  endif
     &                     tl_ubar_obc, tl_vbar_obc,                    &
     &                     tl_zeta_obc,                                 &
# endif
# ifdef ADJUST_WSTRESS
     &                     tl_ustr, tl_vstr,                            &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     tl_tflux,                                    &
#  endif
     &                     tl_t, tl_u, tl_v,                            &
# else
     &                     tl_ubar, tl_vbar,                            &
# endif
     &                     tl_zeta,                                     &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                     ad_t_obc, ad_u_obc, ad_v_obc,                &
#  endif
     &                     ad_ubar_obc, ad_vbar_obc,                    &
     &                     ad_zeta_obc,                                 &
# endif
# ifdef ADJUST_WSTRESS
     &                     ad_ustr, ad_vstr,                            &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     ad_tflux,                                    &
#  endif
     &                     ad_t, ad_u, ad_v,                            &
# else
     &                     ad_ubar, ad_vbar,                            &
# endif
     &                     ad_zeta)
      END IF
!
!  Compute the new cost function.
!
      IF (innLoop.gt.0) THEN
        CALL new_cost (ng, tile, model,                                 &
     &                 LBi, UBi, LBj, UBj, LBij, UBij,                  &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 innLoop, outLoop,                                &
# ifdef MASKING
     &                 rmask, umask, vmask,                             &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                 nl_t_obc, nl_u_obc, nl_v_obc,                    &
#  endif
     &                 nl_ubar_obc, nl_vbar_obc,                        &
     &                 nl_zeta_obc,                                     &
# endif
# ifdef ADJUST_WSTRESS
     &                 nl_ustr, nl_vstr,                                &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                 nl_tflux,                                        &
#  endif
     &                 nl_t, nl_u, nl_v,                                &
# else
     &                 nl_ubar, nl_vbar,                                &
# endif
     &                 nl_zeta)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Determine the eigenvalues and eigenvectors of the tridiagonal matrix.
!  These will be used on the last inner-loop to compute the eigenvectors
!  of the Hessian.
!-----------------------------------------------------------------------
!
      IF (innLoop.gt.0) THEN
        IF (Lprecond.or.LhessianEV) THEN
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
            CALL DSTEQR ('I', innLoop, cg_Ritz(1,outLoop), cg_Tmatrix,  &
     &                   cg_zv, Ninner, work, info)
          END IF
# ifdef DISTRIBUTE
          CALL mp_bcasti (ng, model, info)
# endif
          IF (info.ne.0) THEN
            PRINT *, ' Error in DSTEQR: info=', info
            exit_flag=8
            RETURN
          END IF
# ifdef DISTRIBUTE
          CALL mp_bcastf (ng, model, cg_Ritz(:,outLoop))
          CALL mp_bcastf (ng, model, cg_zv)
# endif
!
!  Estimate the Ritz value error bounds.
!
          DO i=1,innLoop
            cg_RitzErr(i,outLoop)=ABS(cg_beta(innLoop+1,outLoop)*       &
     &                                cg_zv(innLoop,i))
          END DO
!
!  Check for exploding or negative Ritz values.
!
          DO i=1,innLoop
            IF (cg_Ritz(i,outLoop).lt.0.0_r8) THEN
              PRINT *, ' Negative Ritz value found.'
              exit_flag=8
              RETURN
            END IF
          END DO
!
!  Calculate the converged eigenvectors of the Hessian.
!
          IF (innLoop.eq.Ninner) THEN
            RitzMaxErr=HevecErr
            DO i=1,innLoop
              cg_RitzErr(i,outLoop)=cg_RitzErr(i,outLoop)/              &
     &                              cg_Ritz(Ninner,outLoop)
            END DO
            Lwrk=2
            Linp=1
            Lout=2
            CALL hessian_evecs (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj, LBij, UBij,         &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Linp, Lout, Lwrk,                       &
     &                          innLoop, outLoop,                       &
# ifdef MASKING
     &                          rmask, umask, vmask,                    &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                          nl_t_obc, nl_u_obc, nl_v_obc,           &
#  endif
     &                          nl_ubar_obc, nl_vbar_obc,               &
     &                          nl_zeta_obc,                            &
# endif
# ifdef ADJUST_WSTRESS
     &                          nl_ustr, nl_vstr,                       &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                          nl_tflux,                               &
#  endif
     &                          nl_t, nl_u, nl_v,                       &
# else
     &                          nl_ubar, nl_vbar,                       &
# endif
     &                          nl_zeta,                                &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                          tl_t_obc, tl_u_obc, tl_v_obc,           &
#  endif
     &                          tl_ubar_obc, tl_vbar_obc,               &

     &                          tl_zeta_obc,                            &
# endif
# ifdef ADJUST_WSTRESS
     &                          tl_ustr, tl_vstr,                       &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                          tl_tflux,                               &
#  endif
     &                          tl_t, tl_u, tl_v,                       &
# else
     &                          tl_ubar, tl_vbar,                       &
# endif
     &                          tl_zeta,                                &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                          ad_t_obc, ad_u_obc, ad_v_obc,           &
#  endif
     &                          ad_ubar_obc, ad_vbar_obc,               &
     &                          ad_zeta_obc,                            &
# endif
# ifdef ADJUST_WSTRESS
     &                          ad_ustr, ad_vstr,                       &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                          ad_tflux,                               &
#  endif
     &                          ad_t, ad_u, ad_v,                       &
# else
     &                          ad_ubar, ad_vbar,                       &
# endif
     &                          ad_zeta)
            IF (exit_flag.ne.NoError) RETURN

            IF (Master.and.(nConvRitz.eq.0)) THEN
              PRINT *, ' No converged Hesssian eigenvectors found.'
            END IF
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
!   For the Lanczos algorithm, X(Linp) is ALWAYS the starting TLM
!   initial condition which for incremental 4DVar is zero.
!
      Linp=1
      Lout=2
      CALL tl_new_state (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   Linp, Lout,                                    &
     &                   innLoop, outLoop,                              &
# ifdef MASKING
     &                   rmask, umask, vmask,                           &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   d_t_obc, d_u_obc, d_v_obc,                     &
#  endif
     &                   d_ubar_obc, d_vbar_obc,                        &
     &                   d_zeta_obc,                                    &
# endif
# ifdef ADJUST_WSTRESS
     &                   d_sustr, d_svstr,                              &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                   d_stflx,                                       &
#  endif
     &                   d_t, d_u, d_v,                                 &
# else
     &                   d_ubar, d_vbar,                                &
# endif
     &                   d_zeta,                                        &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   tl_t_obc, tl_u_obc, tl_v_obc,                  &
#  endif
     &                   tl_ubar_obc, tl_vbar_obc,                      &
     &                   tl_zeta_obc,                                   &
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
     &                   tl_zeta,                                       &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   ad_t_obc, ad_u_obc, ad_v_obc,                  &
#  endif
     &                   ad_ubar_obc, ad_vbar_obc,                      &
     &                   ad_zeta_obc,                                   &
# endif
# ifdef ADJUST_WSTRESS
     &                   ad_ustr, ad_vstr,                              &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                   ad_tflux,                                      &
#  endif
     &                   ad_t, ad_u, ad_v,                              &
# else
     &                   ad_ubar, ad_vbar,                              &
# endif
     &                   ad_zeta)
!
!  If preconditioning, convert tl_var(Lout) back into v-space.
!
      IF (Lprecond.and.(outLoop.gt.1)) THEN

        Lscale=2                 ! SQRT spectral LMP
        Ltrans=.FALSE.
!
!  Copy tl_var(Lout) into nl_var(L1).
!
        CALL state_copy (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Lout, L1,                                      &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   nl_t_obc, tl_t_obc,                            &
     &                   nl_u_obc, tl_u_obc,                            &
     &                   nl_v_obc, tl_v_obc,                            &
#  endif
     &                   nl_ubar_obc, tl_ubar_obc,                      &
     &                   nl_vbar_obc, tl_vbar_obc,                      &
     &                   nl_zeta_obc, tl_zeta_obc,                      &
# endif
# ifdef ADJUST_WSTRESS
     &                   nl_ustr, tl_ustr,                              &
     &                   nl_vstr, tl_vstr,                              &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                   nl_tflux, tl_tflux,                            &
#  endif
     &                   nl_t, tl_t,                                    &
     &                   nl_u, tl_u,                                    &
     &                   nl_v, tl_v,                                    &
# else
     &                   nl_ubar, tl_ubar,                              &
     &                   nl_vbar, tl_vbar,                              &
# endif
     &                   nl_zeta, tl_zeta)
!
        CALL precond (ng, tile, model, 'convert increment to v-space',  &
     &                LBi, UBi, LBj, UBj, LBij, UBij,                   &
     &                IminS, ImaxS, JminS, JmaxS,                       &
     &                NstateVar(ng), Lscale, Ltrans,                    &
     &                innLoop, outLoop,                                 &
# ifdef MASKING
     &                rmask, umask, vmask,                              &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                nl_t_obc, nl_u_obc, nl_v_obc,                     &
#  endif
     &                nl_ubar_obc, nl_vbar_obc,                         &
     &                nl_zeta_obc,                                      &
# endif
# ifdef ADJUST_WSTRESS
     &                nl_ustr, nl_vstr,                                 &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                nl_tflux,                                         &
#  endif
     &                nl_t, nl_u, nl_v,                                 &
# else
     &                nl_ubar, nl_vbar,                                 &
# endif
     &                nl_zeta)
        IF (exit_flag.ne.NoError) RETURN
!
!  Copy nl_var(L1) into tl_var(Lout)
!
        CALL state_copy (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   L1, Lout,                                      &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   tl_t_obc, nl_t_obc,                            &
     &                   tl_u_obc, nl_u_obc,                            &
     &                   tl_v_obc, nl_v_obc,                            &
#  endif
     &                   tl_ubar_obc, nl_ubar_obc,                      &
     &                   tl_vbar_obc, nl_vbar_obc,                      &
     &                   tl_zeta_obc, nl_zeta_obc,                      &
# endif
# ifdef ADJUST_WSTRESS
     &                   tl_ustr, nl_ustr,                              &
     &                   tl_vstr, nl_vstr,                              &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                   tl_tflux, nl_tflux,                            &
#  endif
     &                   tl_t, nl_t,                                    &
     &                   tl_u, nl_u,                                    &
     &                   tl_v, nl_v,                                    &
# else
     &                   tl_ubar, nl_ubar,                              &
     &                   tl_vbar, nl_vbar,                              &
# endif
     &                   tl_zeta, nl_zeta)
      END IF
!
!-----------------------------------------------------------------------
!  Write out conjugate gradient information into NetCDF file.
!-----------------------------------------------------------------------
!
      CALL cg_write (ng, model, innLoop, outLoop)
      IF (exit_flag.ne.NoError) RETURN
!
!  Report algorithm parameters.
!
      IF (Master) THEN
        IF (inner.eq.0) THEN
          WRITE (stdout,20) outLoop, innLoop,                           &
     &                      cg_Gnorm(outLoop)
 20       FORMAT (/,1x,'(',i3.3,',',i3.3,'): ',                         &
     &            'Initial gradient norm, Gnorm  = ',1p,e14.7)
        END IF
        IF (innLoop.gt.0) THEN
          WRITE (stdout,30) outLoop, innLoop,                           &
     &                      cg_Greduc(innLoop,outLoop),                 &
     &                      outLoop, innLoop,                           &
     &                      cg_delta(innLoop,outLoop)
 30       FORMAT (/,1x,'(',i3.3,',',i3.3,'): ',                         &
     &            'Reduction in the gradient norm,  Greduc = ',         &
     &            1p,e14.7,/,                                           &
     &            1x,'(',i3.3,',',i3.3,'): ',                           &
     &            'Lanczos algorithm coefficient,    delta = ',         &
     &            1p,e14.7)
          WRITE (stdout,40) RitzMaxErr
 40       FORMAT (/,' Ritz Eigenvalues and relative accuracy: ',        &
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
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         Linp, Lout,                              &
     &                         innLoop, outLoop,                        &
# ifdef MASKING
     &                         rmask, umask, vmask,                     &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                         d_t_obc, d_u_obc, d_v_obc,               &
#  endif
     &                         d_ubar_obc, d_vbar_obc,                  &
     &                         d_zeta_obc,                              &
# endif
# ifdef ADJUST_WSTRESS
     &                         d_sustr, d_svstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                         d_stflx,                                 &
#  endif
     &                         d_t, d_u, d_v,                           &
# else
     &                         d_ubar, d_vbar,                          &
# endif
     &                         d_zeta,                                  &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                         tl_t_obc, tl_u_obc, tl_v_obc,            &
#  endif
     &                         tl_ubar_obc, tl_vbar_obc,                &
     &                         tl_zeta_obc,                             &
# endif
# ifdef ADJUST_WSTRESS
     &                         tl_ustr, tl_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                         tl_tflux,                                &
#  endif
     &                         tl_t, tl_u, tl_v,                        &
# else
     &                         tl_ubar, tl_vbar,                        &
# endif
     &                         tl_zeta,                                 &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                         ad_t_obc, ad_u_obc, ad_v_obc,            &
#  endif
     &                         ad_ubar_obc, ad_vbar_obc,                &
     &                         ad_zeta_obc,                             &
# endif
# ifdef ADJUST_WSTRESS
     &                         ad_ustr, ad_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                         ad_tflux,                                &
#  endif
     &                         ad_t, ad_u, ad_v,                        &
# else
     &                         ad_ubar, ad_vbar,                        &
# endif
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
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Linp, Lout
      integer, intent(in) :: innLoop, outLoop
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: d_t_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: d_u_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: d_v_obc(LBij:,:,:,:)
#   endif
      real(r8), intent(inout) :: d_ubar_obc(LBij:,:,:)
      real(r8), intent(inout) :: d_vbar_obc(LBij:,:,:)
      real(r8), intent(inout) :: d_zeta_obc(LBij:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: d_sustr(LBi:,LBj:,:)
      real(r8), intent(in) :: d_svstr(LBi:,LBj:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(in) :: d_stflx(LBi:,LBj:,:,:)
#   endif
      real(r8), intent(in) :: d_t(LBi:,LBj:,:,:)
      real(r8), intent(in) :: d_u(LBi:,LBj:,:)
      real(r8), intent(in) :: d_v(LBi:,LBj:,:)
#  else
      real(r8), intent(in) :: d_ubar(LBi:,LBj:)
      real(r8), intent(in) :: d_vbar(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: d_zeta(LBi:,LBj:)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: ad_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: ad_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: tl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: tl_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(in) :: d_t_obc(LBij:UBij,N(ng),4,                &
     &                                Nbrec(ng),NT(ng))
      real(r8), intent(in) :: d_u_obc(LBij:UBij,N(ng),4,Nbrec(ng))
      real(r8), intent(in) :: d_v_obc(LBij:UBij,N(ng),4,Nbrec(ng))
#   endif
      real(r8), intent(in) :: d_ubar_obc(LBij:UBij,4,Nbrec(ng))
      real(r8), intent(in) :: d_vbar_obc(LBij:UBij,4,Nbrec(ng))
      real(r8), intent(in) :: d_zeta_obc(LBij:UBij,4,Nbrec(ng))
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: d_sustr(LBi:UBi,LBj:UBj,Nfrec(ng))
      real(r8), intent(in) :: d_svstr(LBi:UBi,LBj:UBj,Nfrec(ng))
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(in) :: d_stflx(LBi:UBi,LBj:UBj,                  &
     &                                Nfrec(ng),NT(ng))
#   endif
      real(r8), intent(in) :: d_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(in) :: d_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: d_v(LBi:UBi,LBj:UBj,N(ng))
#  else
      real(r8), intent(in) :: d_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: d_vbar(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: d_zeta(LBi:UBi,LBj:UBj)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: ad_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: tl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k, rec
      integer :: ib, ir, it

      real(r8) :: fac, fac1, fac2

      character (len=256) :: ncname

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute new starting tangent linear state vector, X(k+1).
!-----------------------------------------------------------------------
!
      IF (innLoop.ne.Ninner) THEN
!
!  Free-surface.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            tl_zeta(i,j,Lout)=d_zeta(i,j)
# ifdef MASKING
            tl_zeta(i,j,Lout)=tl_zeta(i,j,Lout)*rmask(i,j)
# endif
          END DO
        END DO

# ifdef ADJUST_BOUNDARY
!
!  Free-surface open boundaries.
!
        IF (ANY(Lobc(:,isFsur,ng))) THEN
          DO ir=1,Nbrec(ng)
            IF ((Lobc(iwest,isFsur,ng)).and.                            &
     &          DOMAIN(ng)%Western_Edge(tile)) THEN
              ib=iwest
              DO j=Jstr,Jend
                tl_zeta_obc(j,ib,ir,Lout)=d_zeta_obc(j,ib,ir)
#  ifdef MASKING
                tl_zeta_obc(j,ib,ir,Lout)=tl_zeta_obc(j,ib,ir,Lout)*    &
     &                                    rmask(Istr-1,j)
#  endif
              END DO
            END IF
            IF ((Lobc(ieast,isFsur,ng)).and.                            &
     &          DOMAIN(ng)%Eastern_Edge(tile)) THEN
              ib=ieast
              DO j=Jstr,Jend
                tl_zeta_obc(j,ib,ir,Lout)=d_zeta_obc(j,ib,ir)
#  ifdef MASKING
                tl_zeta_obc(j,ib,ir,Lout)=tl_zeta_obc(j,ib,ir,Lout)*    &
     &                                    rmask(Iend+1,j)
#  endif
              END DO
            END IF
            IF ((Lobc(isouth,isFsur,ng)).and.                           &
     &          DOMAIN(ng)%Southern_Edge(tile)) THEN
              ib=isouth
              DO i=Istr,Iend
                tl_zeta_obc(i,ib,ir,Lout)=d_zeta_obc(i,ib,ir)
#  ifdef MASKING
                tl_zeta_obc(i,ib,ir,Lout)=tl_zeta_obc(i,ib,ir,Lout)*    &
     &                                    rmask(i,Jstr-1)
#  endif
              END DO
            END IF
            IF ((Lobc(inorth,isFsur,ng)).and.                           &
     &          DOMAIN(ng)%Northern_Edge(tile)) THEN
              ib=inorth
              DO i=Istr,Iend
                tl_zeta_obc(i,ib,ir,Lout)=d_zeta_obc(i,ib,ir)
#  ifdef MASKING
                tl_zeta_obc(i,ib,ir,Lout)=tl_zeta_obc(i,ib,ir,Lout)*    &
     &                                    rmask(i,Jend+1)
#  endif
              END DO
            END IF
          END DO
        END IF
# endif

# ifndef SOLVE3D
!
!  2D U-momentum.
!
        DO j=JstrR,JendR
          DO i=Istr,IendR
            tl_ubar(i,j,Lout)=d_ubar(i,j)
#  ifdef MASKING
            tl_ubar(i,j,Lout)=tl_ubar(i,j,Lout)*umask(i,j)
#  endif
          END DO
        END DO
# endif

# ifdef ADJUST_BOUNDARY
!
!  2D U-momentum open boundaries.
!
        IF (ANY(Lobc(:,isUbar,ng))) THEN
          DO ir=1,Nbrec(ng)
            IF ((Lobc(iwest,isUbar,ng)).and.                            &
     &          DOMAIN(ng)%Western_Edge(tile)) THEN
              ib=iwest
              DO j=Jstr,Jend
                tl_ubar_obc(j,ib,ir,Lout)=d_ubar_obc(j,ib,ir)
#  ifdef MASKING
                tl_ubar_obc(j,ib,ir,Lout)=tl_ubar_obc(j,ib,ir,Lout)*    &
     &                                    umask(Istr,j)
#  endif
              END DO
            END IF
            IF ((Lobc(ieast,isUbar,ng)).and.                            &
     &          DOMAIN(ng)%Eastern_Edge(tile)) THEN
              ib=ieast
              DO j=Jstr,Jend
                tl_ubar_obc(j,ib,ir,Lout)=d_ubar_obc(j,ib,ir)
#  ifdef MASKING
                tl_ubar_obc(j,ib,ir,Lout)=tl_ubar_obc(j,ib,ir,Lout)*    &
     &                                    umask(Iend+1,j)
#  endif
              END DO
            END IF
            IF ((Lobc(isouth,isUbar,ng)).and.                           &
     &          DOMAIN(ng)%Southern_Edge(tile)) THEN
              ib=isouth
              DO i=IstrU,Iend
                tl_ubar_obc(i,ib,ir,Lout)=d_ubar_obc(i,ib,ir)
#  ifdef MASKING
                tl_ubar_obc(i,ib,ir,Lout)=tl_ubar_obc(i,ib,ir,Lout)*    &
     &                                    umask(i,Jstr-1)
#  endif
              END DO
            END IF
            IF ((Lobc(inorth,isUbar,ng)).and.                           &
     &          DOMAIN(ng)%Northern_Edge(tile)) THEN
              ib=inorth
              DO i=IstrU,Iend
                tl_ubar_obc(i,ib,ir,Lout)=d_ubar_obc(i,ib,ir)
#  ifdef MASKING
                tl_ubar_obc(i,ib,ir,Lout)=tl_ubar_obc(i,ib,ir,Lout)*    &
     &                                    umask(i,Jend+1)
#  endif
              END DO
            END IF
          END DO
        END IF
# endif

# ifndef SOLVE3D
!
!  2D V-momentum.
!
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            tl_vbar(i,j,Lout)=d_vbar(i,j)
#  ifdef MASKING
            tl_vbar(i,j,Lout)=tl_vbar(i,j,Lout)*vmask(i,j)
#  endif
          END DO
        END DO
# endif

# ifdef ADJUST_BOUNDARY
!
!  2D V-momentum open boundaries.
!
        IF (ANY(Lobc(:,isVbar,ng))) THEN
          DO ir=1,Nbrec(ng)
            IF ((Lobc(iwest,isVbar,ng)).and.                            &
     &          DOMAIN(ng)%Western_Edge(tile)) THEN
              ib=iwest
              DO j=JstrV,Jend
                tl_vbar_obc(j,ib,ir,Lout)=d_vbar_obc(j,ib,ir)
#  ifdef MASKING
                tl_vbar_obc(j,ib,ir,Lout)=tl_vbar_obc(j,ib,ir,Lout)*    &
     &                                    vmask(Istr-1,j)
#  endif
              END DO
            END IF
            IF ((Lobc(ieast,isVbar,ng)).and.                            &
     &          DOMAIN(ng)%Eastern_Edge(tile)) THEN
              ib=ieast
              DO j=JstrV,Jend
                tl_vbar_obc(j,ib,ir,Lout)=d_vbar_obc(j,ib,ir)
#  ifdef MASKING
                tl_vbar_obc(j,ib,ir,Lout)=tl_vbar_obc(j,ib,ir,Lout)*    &
     &                                    vmask(Iend+1,j)
#  endif
              END DO
            END IF
            IF ((Lobc(isouth,isVbar,ng)).and.                           &
     &          DOMAIN(ng)%Southern_Edge(tile)) THEN
              ib=isouth
              DO i=Istr,Iend
                tl_vbar_obc(i,ib,ir,Lout)=d_vbar_obc(i,ib,ir)
#  ifdef MASKING
                tl_vbar_obc(i,ib,ir,Lout)=tl_vbar_obc(i,ib,ir,Lout)*    &
     &                                    vmask(i,Jstr)
#  endif
              END DO
            END IF
            IF ((Lobc(inorth,isVbar,ng)).and.                           &
     &          DOMAIN(ng)%Northern_Edge(tile)) THEN
              ib=inorth
              DO i=Istr,Iend
                tl_vbar_obc(i,ib,ir,Lout)=d_vbar_obc(i,ib,ir)
#  ifdef MASKING
                tl_vbar_obc(i,ib,ir,Lout)=tl_vbar_obc(i,ib,ir,Lout)*    &
     &                                    vmask(i,Jend+1)
#  endif
              END DO
            END IF
          END DO
        END IF
# endif

# ifdef ADJUST_WSTRESS
!
!  Surface momentum stress.
!
        DO ir=1,Nfrec(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              tl_ustr(i,j,ir,Lout)=d_sustr(i,j,ir)
#  ifdef MASKING
              tl_ustr(i,j,ir,Lout)=tl_ustr(i,j,ir,Lout)*umask(i,j)
#  endif
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              tl_vstr(i,j,ir,Lout)=d_svstr(i,j,ir)
#  ifdef MASKING
              tl_vstr(i,j,ir,Lout)=tl_vstr(i,j,ir,Lout)*vmask(i,j)
#  endif
            END DO
          END DO
        END DO
# endif

# ifdef SOLVE3D
!
!  3D U-momentum.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              tl_u(i,j,k,Lout)=d_u(i,j,k)
#  ifdef MASKING
              tl_u(i,j,k,Lout)=tl_u(i,j,k,Lout)*umask(i,j)
#  endif
            END DO
          END DO
        END DO

#  ifdef ADJUST_BOUNDARY
!
!  3D U-momentum open boundaries.
!
        IF (ANY(Lobc(:,isUvel,ng))) THEN
          DO ir=1,Nbrec(ng)
            IF ((Lobc(iwest,isUvel,ng)).and.                            &
     &          DOMAIN(ng)%Western_Edge(tile)) THEN
              ib=iwest
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  tl_u_obc(j,k,ib,ir,Lout)=d_u_obc(j,k,ib,ir)
#   ifdef MASKING
                  tl_u_obc(j,k,ib,ir,Lout)=tl_u_obc(j,k,ib,ir,Lout)*    &
     &                                     umask(Istr,j)
#   endif
                END DO
              END DO
            END IF
            IF ((Lobc(ieast,isUvel,ng)).and.                            &
     &          DOMAIN(ng)%Eastern_Edge(tile)) THEN
              ib=ieast
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  tl_u_obc(j,k,ib,ir,Lout)=d_u_obc(j,k,ib,ir)
#   ifdef MASKING
                  tl_u_obc(j,k,ib,ir,Lout)=tl_u_obc(j,k,ib,ir,Lout)*    &
     &                                     umask(Iend+1,j)
#   endif
                END DO
              END DO
            END IF
            IF ((Lobc(isouth,isUvel,ng)).and.                           &
     &          DOMAIN(ng)%Southern_Edge(tile)) THEN
              ib=isouth
              DO k=1,N(ng)
                DO i=IstrU,Iend
                  tl_u_obc(i,k,ib,ir,Lout)=d_u_obc(i,k,ib,ir)
#   ifdef MASKING
                  tl_u_obc(i,k,ib,ir,Lout)=tl_u_obc(i,k,ib,ir,Lout)*    &
     &                                     umask(i,Jstr-1)
#   endif
                END DO
              END DO
            END IF
            IF ((Lobc(inorth,isUvel,ng)).and.                           &
     &          DOMAIN(ng)%Northern_Edge(tile)) THEN
              ib=inorth
              DO k=1,N(ng)
                DO i=IstrU,Iend
                  tl_u_obc(i,k,ib,ir,Lout)=d_u_obc(i,k,ib,ir)
#   ifdef MASKING
                  tl_u_obc(i,k,ib,ir,Lout)=tl_u_obc(i,k,ib,ir,Lout)*    &
     &                                     umask(i,Jend+1)
#   endif
                END DO
              END DO
            END IF
          END DO
        END IF
#  endif
!
!  3D V-momentum.
!
        DO k=1,N(ng)
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              tl_v(i,j,k,Lout)=d_v(i,j,k)
#  ifdef MASKING
              tl_v(i,j,k,Lout)=tl_v(i,j,k,Lout)*vmask(i,j)
#  endif
            END DO
          END DO
        END DO

#  ifdef ADJUST_BOUNDARY
!
!  3D V-momentum open boundaries.
!
        IF (ANY(Lobc(:,isVvel,ng))) THEN
          DO ir=1,Nbrec(ng)
            IF ((Lobc(iwest,isVvel,ng)).and.                            &
     &          DOMAIN(ng)%Western_Edge(tile)) THEN
              ib=iwest
              DO k=1,N(ng)
                DO j=JstrV,Jend
                  tl_v_obc(j,k,ib,ir,Lout)=d_v_obc(j,k,ib,ir)
#   ifdef MASKING
                  tl_v_obc(j,k,ib,ir,Lout)=tl_v_obc(j,k,ib,ir,Lout)*    &
     &                                     vmask(Istr-1,j)
#   endif
                END DO
              END DO
            END IF
            IF ((Lobc(ieast,isVvel,ng)).and.                            &
     &          DOMAIN(ng)%Eastern_Edge(tile)) THEN
              ib=ieast
              DO k=1,N(ng)
                DO j=JstrV,Jend
                  tl_v_obc(j,k,ib,ir,Lout)=d_v_obc(j,k,ib,ir)
#   ifdef MASKING
                  tl_v_obc(j,k,ib,ir,Lout)=tl_v_obc(j,k,ib,ir,Lout)*    &
     &                                     vmask(Iend+1,j)
#   endif
                END DO
              END DO
            END IF
            IF ((Lobc(isouth,isVvel,ng)).and.                           &
     &          DOMAIN(ng)%Southern_Edge(tile)) THEN
              ib=isouth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  tl_v_obc(i,k,ib,ir,Lout)=d_v_obc(i,k,ib,ir)
#   ifdef MASKING
                  tl_v_obc(i,k,ib,ir,Lout)=tl_v_obc(i,k,ib,ir,Lout)*    &
     &                                     vmask(i,Jstr)
#   endif
                END DO
              END DO
            END IF
            IF ((Lobc(inorth,isVvel,ng)).and.                           &
     &          DOMAIN(ng)%Northern_Edge(tile)) THEN
              ib=inorth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  tl_v_obc(i,k,ib,ir,Lout)=d_v_obc(i,k,ib,ir)
#   ifdef MASKING
                  tl_v_obc(i,k,ib,ir,Lout)=tl_v_obc(i,k,ib,ir,Lout)*    &
     &                                     vmask(i,Jend+1)
#   endif
                END DO
              END DO
            END IF
          END DO
        END IF
#  endif
!
!  Tracers.
!
        DO it=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                tl_t(i,j,k,Lout,it)=d_t(i,j,k,it)
#  ifdef MASKING
                tl_t(i,j,k,Lout,it)=tl_t(i,j,k,Lout,it)*rmask(i,j)
#  endif
              END DO
            END DO
          END DO
        END DO

#  ifdef ADJUST_BOUNDARY
!
!  Tracers open boundaries.
!
        DO it=1,NT(ng)
          IF (ANY(Lobc(:,isTvar(it),ng))) THEN
            DO ir=1,Nbrec(ng)
              IF ((Lobc(iwest,isTvar(it),ng)).and.                      &
     &            DOMAIN(ng)%Western_Edge(tile)) THEN
                ib=iwest
                DO k=1,N(ng)
                  DO j=Jstr,Jend
                    tl_t_obc(j,k,ib,ir,Lout,it)=d_t_obc(j,k,ib,ir,it)
#   ifdef MASKING
                    tl_t_obc(j,k,ib,ir,Lout,it)=                        &
     &                      tl_t_obc(j,k,ib,ir,Lout,it)*rmask(Istr-1,j)
#   endif
                  END DO
                END DO
              END IF
              IF ((Lobc(ieast,isTvar(it),ng)).and.                      &
     &            DOMAIN(ng)%Eastern_Edge(tile)) THEN
                ib=ieast
                DO k=1,N(ng)
                  DO j=Jstr,Jend
                    tl_t_obc(j,k,ib,ir,Lout,it)=d_t_obc(j,k,ib,ir,it)
#   ifdef MASKING
                    tl_t_obc(j,k,ib,ir,Lout,it)=                        &
     &                      tl_t_obc(j,k,ib,ir,Lout,it)*rmask(Iend+1,j)
#   endif
                  END DO
                END DO
              END IF
              IF ((Lobc(isouth,isTvar(it),ng)).and.                     &
     &            DOMAIN(ng)%Southern_Edge(tile)) THEN
                ib=isouth
                DO k=1,N(ng)
                  DO i=Istr,Iend
                    tl_t_obc(i,k,ib,ir,Lout,it)=d_t_obc(i,k,ib,ir,it)
#   ifdef MASKING
                    tl_t_obc(i,k,ib,ir,Lout,it)=                        &
     &                      tl_t_obc(i,k,ib,ir,Lout,it)*rmask(i,Jstr-1)
#   endif
                  END DO
                END DO
              END IF
              IF ((Lobc(inorth,isTvar(it),ng)).and.                     &
     &            DOMAIN(ng)%Northern_Edge(tile)) THEN
                ib=inorth
                DO k=1,N(ng)
                  DO i=Istr,Iend
                    tl_t_obc(i,k,ib,ir,Lout,it)=d_t_obc(i,k,ib,ir,it)
#   ifdef MASKING
                    tl_t_obc(i,k,ib,ir,Lout,it)=                        &
     &                      tl_t_obc(i,k,ib,ir,Lout,it)*rmask(i,Jend+1)
#   endif
                  END DO
                END DO
              END IF
            END DO
          END IF
        END DO
#  endif

#  ifdef ADJUST_STFLUX
!
!  Surface tracers flux.
!
        DO it=1,NT(ng)
          IF (Lstflux(it,ng)) THEN
            DO ir=1,Nfrec(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  tl_tflux(i,j,ir,Lout,it)=d_stflx(i,j,ir,it)
#   ifdef MASKING
                  tl_tflux(i,j,ir,Lout,it)=tl_tflux(i,j,ir,Lout,it)*    &
     &                                     rmask(i,j)
#   endif
                END DO
              END DO
            END DO
          END IF
        END DO
#  endif

# endif
!
!-----------------------------------------------------------------------
!  If last inner-loop, compute the tangent linear model initial
!  conditions from the Lanczos algorithm. Use adjoint state arrays,
!  index Linp, as temporary storage.
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
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         Linp, fac,                               &
# ifdef MASKING
     &                         rmask, umask, vmask,                     &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                         ad_t_obc, ad_u_obc, ad_v_obc,            &
#  endif
     &                         ad_ubar_obc, ad_vbar_obc,                &
     &                         ad_zeta_obc,                             &
# endif
# ifdef ADJUST_WSTRESS
     &                         ad_ustr, ad_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                         ad_tflux,                                &
#  endif
     &                         ad_t, ad_u, ad_v,                        &
# else
     &                         ad_ubar, ad_vbar,                        &
# endif
     &                         ad_zeta)
!
!  Read in each previous gradient state solutions, g(0) to g(k).
!
        IF (ndefADJ(ng).gt.0) THEN
          WRITE (ncname,10) TRIM(ADM(ng)%base), outLoop
 10       FORMAT (a,'_',i3.3,'.nc')
        ELSE
          ncname=ADM(ng)%name
        END IF
!
        DO rec=1,innLoop
!
!  Read gradient solution and load it into TANGENT LINEAR STATE ARRAYS
!  at index Lout.
!
          CALL read_state (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj, LBij, UBij,              &
     &                     Lout, rec,                                   &
     &                     ndefADJ(ng), ADM(ng)%ncid, ncname,           &
# ifdef MASKING
     &                     rmask, umask, vmask,                         &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                     tl_t_obc, tl_u_obc, tl_v_obc,                &
#  endif
     &                     tl_ubar_obc, tl_vbar_obc,                    &
     &                     tl_zeta_obc,                                 &
# endif
# ifdef ADJUST_WSTRESS
     &                     tl_ustr, tl_vstr,                            &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     tl_tflux,                                    &
#  endif
     &                     tl_t, tl_u, tl_v,                            &
# else
     &                     tl_ubar, tl_vbar,                            &
# endif
     &                     tl_zeta)
          IF (exit_flag.ne.NoError) RETURN
!
!  Sum all previous normalized gradients:
!
!    ad_var(Linp) = fac1 * ad_var(Linp) + fac2 * tl_var(Lout)
!
          fac1=1.0_r8
          fac2=cg_zu(rec,outLoop)

          CALL state_addition (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         Linp, Lout, Linp, fac1, fac2,            &
# ifdef MASKING
     &                         rmask, umask, vmask,                     &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                         ad_t_obc, tl_t_obc,                      &
     &                         ad_u_obc, tl_u_obc,                      &
     &                         ad_v_obc, tl_v_obc,                      &
#  endif
     &                         ad_ubar_obc, tl_ubar_obc,                &
     &                         ad_vbar_obc, tl_vbar_obc,                &
     &                         ad_zeta_obc, tl_zeta_obc,                &
# endif
# ifdef ADJUST_WSTRESS
     &                         ad_ustr, tl_ustr,                        &
     &                         ad_vstr, tl_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                         ad_tflux, tl_tflux,                      &
#  endif
     &                         ad_t, tl_t,                              &
     &                         ad_u, tl_u,                              &
     &                         ad_v, tl_v,                              &
# else
     &                         ad_ubar, tl_ubar,                        &
     &                         ad_vbar, tl_vbar,                        &
# endif
     &                         ad_zeta, tl_zeta)
        END DO
!
!  Load new tangent linear model initial conditions to respective state
!  arrays, index Lout:
!
!    tl_var(Lout) = ad_var(Linp)
!
        CALL state_copy (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Linp, Lout,                                    &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   tl_t_obc, ad_t_obc,                            &
     &                   tl_u_obc, ad_u_obc,                            &
     &                   tl_v_obc, ad_v_obc,                            &
#  endif
     &                   tl_ubar_obc, ad_ubar_obc,                      &
     &                   tl_vbar_obc, ad_vbar_obc,                      &
     &                   tl_zeta_obc, ad_zeta_obc,                      &
# endif
# ifdef ADJUST_WSTRESS
     &                   tl_ustr, ad_ustr,                              &
     &                   tl_vstr, ad_vstr,                              &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                   tl_tflux, ad_tflux,                            &
#  endif
     &                   tl_t, ad_t,                                    &
     &                   tl_u, ad_u,                                    &
     &                   tl_v, ad_v,                                    &
# else
     &                   tl_ubar, ad_ubar,                              &
     &                   tl_vbar, ad_vbar,                              &
# endif
     &                   tl_zeta, ad_zeta)
      END IF

      RETURN
      END SUBROUTINE tl_new_state
!
!***********************************************************************
      SUBROUTINE read_state (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       Lwrk, rec,                                 &
     &                       ndef, ncfileid, ncname,                    &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       s_t_obc, s_u_obc, s_v_obc,                 &
#  endif
     &                       s_ubar_obc, s_vbar_obc,                    &
     &                       s_zeta_obc,                                &
# endif
# ifdef ADJUST_WSTRESS
     &                       s_ustr, s_vstr,                            &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                       s_tflux,                                   &
#  endif
     &                       s_t, s_u, s_v,                             &
# else
     &                       s_ubar, s_vbar,                            &
# endif
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
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
# endif
# ifdef ADJUST_BOUNDARY
      USE nf_fread2d_bry_mod, ONLY : nf_fread2d_bry
#  ifdef SOLVE3D
      USE nf_fread3d_bry_mod, ONLY : nf_fread3d_bry
#  endif
# endif
      USE nf_fread2d_mod, ONLY : nf_fread2d
# ifdef SOLVE3D
      USE nf_fread3d_mod, ONLY : nf_fread3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: Lwrk, rec, ndef

      integer, intent(inout) :: ncfileid

      character (len=*), intent(in) :: ncname
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: s_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: s_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: s_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: s_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: s_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: s_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: s_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: s_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: s_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: s_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: s_zeta(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: s_t_obc(LBij:UBij,N(ng),4,             &
     &                                   Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: s_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: s_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: s_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: s_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: s_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: s_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s_tflux(LBi:UBi,LBj:UBj,               &
     &                                   Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: s_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: s_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: s_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: s_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: s_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: s_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k
      integer :: ifield, it
      integer :: gtype, ncid, status, varid

      integer, dimension(4) :: Vsize

      real(r8) :: Fmin, Fmax, scale

# include "set_bounds.h"
!
      SourceFile='cgradient_lanczos.h, read_state'
!
!-----------------------------------------------------------------------
!  Read in requested model state record. Load data into state array
!  index Lwrk.
!-----------------------------------------------------------------------
!
!  Determine file and variables ids.
!
      IF (ndef.gt.0) THEN
        CALL netcdf_open (ng, model, ncname, 0, ncid)
        IF (exit_flag.ne.NoError) THEN
          WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
        ncfileid=ncid
      ELSE
        ncid=ncfileid
      END IF

      DO i=1,4
        Vsize(i)=0
      END DO
!
!  Read in free-surface.
!
      gtype=r2dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idFsur),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread2d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idFsur), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  rmask,                                          &
# endif
     &                  s_zeta(:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idFsur)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Read in free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        ifield=idSbry(isFsur)
        gtype=r2dvar
        scale=1.0_r8
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),      &
     &                         ncid, varid)
        IF (exit_flag.ne.NoError) RETURN

        status=nf_fread2d_bry (ng, model, ncname, ncid,                 &
     &                         Vname(1,ifield), varid,                  &
     &                         rec, gtype,                              &
     &                         LBij, UBij, Nbrec(ng),                   &
     &                         scale, Fmin, Fmax,                       &
     &                         s_zeta_obc(:,:,:,Lwrk))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifndef SOLVE3D
!
!  Read in 2D U-momentum component.
!
      gtype=u2dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idUbar),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread2d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idUbar), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  umask,                                          &
#  endif
     &                  s_ubar(:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUbar)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Read in 2D V-momentum component.
!
      gtype=v2dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idVbar),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread2d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idVbar), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  vmask,                                          &
#  endif
     &                  s_vbar(:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVbar)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# ifdef ADJUST_BOUNDARY
!
!  Read in 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        ifield=idSbry(isUbar)
        gtype=u2dvar
        scale=1.0_r8
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),      &
     &                         ncid, varid)
        IF (exit_flag.ne.NoError) RETURN

        status=nf_fread2d_bry (ng, model, ncname, ncid,                 &
     &                         Vname(1,ifield), varid,                  &
     &                         rec, gtype,                              &
     &                         LBij, UBij, Nbrec(ng),                   &
     &                         scale, Fmin, Fmax,                       &
     &                         s_ubar_obc(:,:,:,Lwrk))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Read in 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        ifield=idSbry(isVbar)
        gtype=u2dvar
        scale=1.0_r8
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),      &
     &                         ncid, varid)
        IF (exit_flag.ne.NoError) RETURN

        status=nf_fread2d_bry (ng, model, ncname, ncid,                 &
     &                         Vname(1,ifield), varid,                  &
     &                         rec, gtype,                              &
     &                         LBij, UBij, Nbrec(ng),                   &
     &                         scale, Fmin, Fmax,                       &
     &                         s_vbar_obc(:,:,:,Lwrk))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef ADJUST_WSTRESS
!
!  Read surface momentum stress.
!
      gtype=u3dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idUsms),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread3d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idUsms), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj, 1, Nfrec(ng),               &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  umask,                                          &
#  endif
     &                  s_ustr(:,:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUsms)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

      gtype=v3dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idVsms),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread3d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idVsms), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj, 1, Nfrec(ng),               &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  vmask,                                          &
#  endif
     &                  s_vstr(:,:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVsms)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif

# ifdef SOLVE3D
!
!  Read in 3D U-momentum component.
!
      gtype=u3dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idUvel),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread3d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idUvel), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  umask,                                          &
#  endif
     &                  s_u(:,:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUvel)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Read in 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        ifield=idSbry(isUvel)
        gtype=u3dvar
        scale=1.0_r8
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),      &
     &                         ncid, varid)
        IF (exit_flag.ne.NoError) RETURN

        status=nf_fread3d_bry (ng, model, ncname, ncid,                 &
     &                         Vname(1,ifield), varid,                  &
     &                         rec, gtype,                              &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng),         &
     &                         scale, Fmin, Fmax,                       &
     &                         s_u_obc(:,:,:,:,Lwrk))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Read in 3D U-momentum component.
!
      gtype=v3dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idVvel),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread3d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idVvel), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  vmask,                                          &
#  endif
     &                  s_v(:,:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVvel)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Read in 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        ifield=idSbry(isVvel)
        gtype=u3dvar
        scale=1.0_r8
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),      &
     &                         ncid, varid)
        IF (exit_flag.ne.NoError) RETURN

        status=nf_fread3d_bry (ng, model, ncname, ncid,                 &
     &                         Vname(1,ifield), varid,                  &
     &                         rec, gtype,                              &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng),         &
     &                         scale, Fmin, Fmax,                       &
     &                         s_v_obc(:,:,:,:,Lwrk))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Read in tracers.
!
      gtype=r3dvar
      scale=1.0_r8
      DO it=1,NT(ng)
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idTvar(it)),  &
     &                         ncid, varid)
        IF (exit_flag.ne.NoError) RETURN

        status=nf_fread3d(ng, model, ncname, ncid,                      &
     &                    Vname(1,idTvar(it)), varid,                   &
     &                    rec, gtype, Vsize,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    scale, Fmin, Fmax,                            &
#  ifdef MASKING
     &                    rmask,                                        &
#  endif
     &                    s_t(:,:,:,Lwrk,it))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTvar(it))), rec,           &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  Read in tracers open boundaries.
!
      DO it=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(it),ng))) THEN
          ifield=idSbry(isTvar(it))
          gtype=r3dvar
          scale=1.0_r8
          CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),    &
     &                           ncid, varid)
          IF (exit_flag.ne.NoError) RETURN

          status=nf_fread3d_bry (ng, model, ncname, ncid,               &
     &                           Vname(1,ifield), varid,                &
     &                           rec, gtype,                            &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng),       &
     &                           scale, Fmin, Fmax,                     &
     &                           s_t_obc(:,:,:,:,Lwrk,it))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
#  ifdef ADJUST_STFLUX
!
!  Read in surface tracers flux.
!
      gtype=r3dvar
      scale=1.0_r8
      DO it=1,NT(ng)
        IF (Lstflux(it,ng)) THEN
          CALL netcdf_inq_varid (ng, model, ncname,                     &
     &                           Vname(1,idTsur(it)), ncid, varid)
          IF (exit_flag.ne.NoError) RETURN

          status=nf_fread3d(ng, model, ncname, ncid,                    &
     &                      Vname(1,idTsur(it)), varid,                 &
     &                      rec, gtype, Vsize,                          &
     &                      LBi, UBi, LBj, UBj, 1, Nfrec(ng),           &
     &                      scale, Fmin, Fmax,                          &
#   ifdef MASKING
     &                      rmask,                                      &
#   endif
     &                      s_tflux(:,:,:,Lwrk,it))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTsur(it))), rec,         &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
# endif
!
!  If multiple files, close current file.
!
      IF (ndef.gt.0) THEN
        CALL netcdf_close (ng, model, ncid, ncname, .FALSE.)
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
     &                          LBi, UBi, LBj, UBj, LBij, UBij,         &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Lold, Lnew,                             &
# ifdef MASKING
     &                          rmask, umask, vmask,                    &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                          ad_t_obc, ad_u_obc, ad_v_obc,           &
#  endif
     &                          ad_ubar_obc, ad_vbar_obc,               &
     &                          ad_zeta_obc,                            &
# endif
# ifdef ADJUST_WSTRESS
     &                          ad_ustr, ad_vstr,                       &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                          ad_tflux,                               &
#  endif
     &                          ad_t, ad_u, ad_v,                       &
# else
     &                          ad_ubar, ad_vbar,                       &
# endif
     &                          ad_zeta,                                &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                          d_t_obc, d_u_obc, d_v_obc,              &
#  endif
     &                          d_ubar_obc, d_vbar_obc,                 &
     &                          d_zeta_obc,                             &
# endif
# ifdef ADJUST_WSTRESS
     &                          d_sustr, d_svstr,                       &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                          d_stflx,                                &
#  endif
     &                          d_t, d_u, d_v,                          &
# else
     &                          d_ubar, d_vbar,                         &
# endif
     &                          d_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_parallel
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS || \
    defined ADJUST_BOUNDARY
      USE mod_scalars
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lold, Lnew
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(in) :: ad_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(in) :: ad_u_obc(LBij:,:,:,:,:)
      real(r8), intent(in) :: ad_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(in) :: ad_ubar_obc(LBij:,:,:,:)
      real(r8), intent(in) :: ad_vbar_obc(LBij:,:,:,:)
      real(r8), intent(in) :: ad_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(in) :: ad_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(in) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(in) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: ad_zeta(LBi:,LBj:,:)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: d_t_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: d_u_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: d_v_obc(LBij:,:,:,:)
#   endif
      real(r8), intent(inout) :: d_ubar_obc(LBij:,:,:)
      real(r8), intent(inout) :: d_vbar_obc(LBij:,:,:)
      real(r8), intent(inout) :: d_zeta_obc(LBij:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: d_sustr(LBi:,LBj:,:)
      real(r8), intent(inout) :: d_svstr(LBi:,LBj:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: d_stflx(LBi:,LBj:,:,:)
#   endif
      real(r8), intent(inout) :: d_t(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: d_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: d_v(LBi:,LBj:,:)
#  else
      real(r8), intent(inout) :: d_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: d_vbar(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: d_zeta(LBi:,LBj:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(in) :: ad_t_obc(LBij:UBij,N(ng),4,               &
     &                                 Nbrec(ng),2,NT(ng))
      real(r8), intent(in) :: ad_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(in) :: ad_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(in) :: ad_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(in) :: ad_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(in) :: ad_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(in) :: ad_vstr(LBi:UBI,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(in) :: ad_tflux(LBi:UBi,LBj:UBj,                 &
     &                                 Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(in) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(in) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(in) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: d_t_obc(LBij:UBij,N(ng),4,             &
     &                                   Nbrec(ng),NT(ng))
      real(r8), intent(inout) :: d_u_obc(LBij:UBij,N(ng),4,Nbrec(ng))
      real(r8), intent(inout) :: d_v_obc(LBij:UBij,N(ng),4,Nbrec(ng))
#   endif
      real(r8), intent(inout) :: d_ubar_obc(LBij:UBij,4,Nbrec(ng))
      real(r8), intent(inout) :: d_vbar_obc(LBij:UBij,4,Nbrec(ng))
      real(r8), intent(inout) :: d_zeta_obc(LBij:UBij,4,Nbrec(ng))
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: d_sustr(LBi:UBi,LBj:UBj,Nfrec(ng))
      real(r8), intent(inout) :: d_svstr(LBi:UBI,LBj:UBj,Nfrec(ng))
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: d_stflx(LBi:UBi,LBj:UBj,               &
     &                                   Nfrec(ng),NT(ng))
#   endif
      real(r8), intent(inout) :: d_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: d_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: d_v(LBi:UBi,LBj:UBj,N(ng))
#  else
      real(r8), intent(inout) :: d_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: d_vbar(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: d_zeta(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k
      integer :: ib, ir, it

# include "set_bounds.h"
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
# ifdef MASKING
          d_zeta(i,j)=d_zeta(i,j)*rmask(i,j)
# endif
        END DO
      END DO

# ifdef ADJUST_BOUNDARY
!
!  Free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isFsur,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO j=Jstr,Jend
              d_zeta_obc(j,ib,ir)=ad_zeta_obc(j,ib,ir,Lnew)
#  ifdef MASKING
              d_zeta_obc(j,ib,ir)=d_zeta_obc(j,ib,ir)*                  &
     &                            rmask(Istr-1,j)
#  endif
            END DO
          END IF
          IF ((Lobc(ieast,isFsur,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=Jstr,Jend
              d_zeta_obc(j,ib,ir)=ad_zeta_obc(j,ib,ir,Lnew)
#  ifdef MASKING
              d_zeta_obc(j,ib,ir)=d_zeta_obc(j,ib,ir)*                  &
     &                            rmask(Iend+1,j)
#  endif
            END DO
          END IF
          IF ((Lobc(isouth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=Istr,Iend
              d_zeta_obc(i,ib,ir)=ad_zeta_obc(i,ib,ir,Lnew)
#  ifdef MASKING
              d_zeta_obc(i,ib,ir)=d_zeta_obc(i,ib,ir)*                  &
     &                            rmask(i,Jstr-1)
#  endif
            END DO
          END IF
          IF ((Lobc(inorth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=Istr,Iend
              d_zeta_obc(i,ib,ir)=ad_zeta_obc(i,ib,ir,Lnew)
#  ifdef MASKING
              d_zeta_obc(i,ib,ir)=d_zeta_obc(i,ib,ir)*                  &
     &                            rmask(i,Jend+1)
#  endif
            END DO
          END IF
        END DO
      END IF
# endif

# ifndef SOLVE3D
!
!  2D U-momentum.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          d_ubar(i,j)=ad_ubar(i,j,Lnew)
#  ifdef MASKING
          d_ubar(i,j)=d_ubar(i,j)*umask(i,j)
#  endif
        END DO
      END DO
# endif

# ifdef ADJUST_BOUNDARY
!
!  2D U-momentum open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isUbar,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO j=Jstr,Jend
              d_ubar_obc(j,ib,ir)=ad_ubar_obc(j,ib,ir,Lnew)
#  ifdef MASKING
              d_ubar_obc(j,ib,ir)=d_ubar_obc(j,ib,ir)*                  &
     &                            umask(Istr,j)
#  endif
            END DO
          END IF
          IF ((Lobc(ieast,isUbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=Jstr,Jend
              d_ubar_obc(j,ib,ir)=ad_ubar_obc(j,ib,ir,Lnew)
#  ifdef MASKING
              d_ubar_obc(j,ib,ir)=d_ubar_obc(j,ib,ir)*                  &
     &                            umask(Iend+1,j)
#  endif
            END DO
          END IF
          IF ((Lobc(isouth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=IstrU,Iend
              d_ubar_obc(i,ib,ir)=ad_ubar_obc(i,ib,ir,Lnew)
#  ifdef MASKING
              d_ubar_obc(i,ib,ir)=d_ubar_obc(i,ib,ir)*                  &
     &                            umask(i,Jstr-1)
#  endif
            END DO
          END IF
          IF ((Lobc(inorth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=IstrU,Iend
              d_ubar_obc(i,ib,ir)=ad_ubar_obc(i,ib,ir,Lnew)
#  ifdef MASKING
              d_ubar_obc(i,ib,ir)=d_ubar_obc(i,ib,ir)*                  &
     &                            umask(i,Jend+1)
#  endif
            END DO
          END IF
        END DO
      END IF
# endif

# ifndef SOLVE3D
!
!  2D V-momentum.
!
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          d_vbar(i,j)=ad_vbar(i,j,Lnew)
#  ifdef MASKING
          d_vbar(i,j)=d_vbar(i,j)*vmask(i,j)
#  endif
        END DO
      END DO
# endif

# ifdef ADJUST_BOUNDARY
!
!  2D V-momentum open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isVbar,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO j=JstrV,Jend
              d_vbar_obc(j,ib,ir)=ad_vbar_obc(j,ib,ir,Lnew)
#  ifdef MASKING
              d_vbar_obc(j,ib,ir)=d_vbar_obc(j,ib,ir)*                  &
     &                            vmask(Istr-1,j)
#  endif
            END DO
          END IF
          IF ((Lobc(ieast,isVbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=JstrV,Jend
              d_vbar_obc(j,ib,ir)=ad_vbar_obc(j,ib,ir,Lnew)
#  ifdef MASKING
              d_vbar_obc(j,ib,ir)=d_vbar_obc(j,ib,ir)*                  &
     &                            vmask(Iend+1,j)
#  endif
            END DO
          END IF
          IF ((Lobc(isouth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=Istr,Iend
              d_vbar_obc(i,ib,ir)=ad_vbar_obc(i,ib,ir,Lnew)
#  ifdef MASKING
              d_vbar_obc(i,ib,ir)=d_vbar_obc(i,ib,ir)*                  &
     &                            vmask(i,Jstr)
#  endif
            END DO
          END IF
          IF ((Lobc(inorth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=Istr,Iend
              d_vbar_obc(i,ib,ir)=ad_vbar_obc(i,ib,ir,Lnew)
#  ifdef MASKING
              d_vbar_obc(i,ib,ir)=d_vbar_obc(i,ib,ir)*                  &
     &                            vmask(i,Jend+1)
#  endif
            END DO
          END IF
        END DO
      END IF
# endif

# ifdef ADJUST_WSTRESS
!
!  Surface momentum stress.
!
      DO ir=1,Nfrec(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            d_sustr(i,j,ir)=ad_ustr(i,j,ir,Lnew)
#  ifdef MASKING
            d_sustr(i,j,ir)=d_sustr(i,j,ir)*umask(i,j)
#  endif
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            d_svstr(i,j,ir)=ad_vstr(i,j,ir,Lnew)
#  ifdef MASKING
            d_svstr(i,j,ir)=d_svstr(i,j,ir)*vmask(i,j)
#  endif
          END DO
        END DO
      END DO
# endif

# ifdef SOLVE3D
!
!  3D U-momentum.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            d_u(i,j,k)=ad_u(i,j,k,Lnew)
#  ifdef MASKING
            d_u(i,j,k)=d_u(i,j,k)*umask(i,j)
#  endif
          END DO
        END DO
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  3D U-momentum open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isUvel,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO k=1,N(ng)
              DO j=Jstr,Jend
                d_u_obc(j,k,ib,ir)=ad_u_obc(j,k,ib,ir,Lnew)
#   ifdef MASKING
                d_u_obc(j,k,ib,ir)=d_u_obc(j,k,ib,ir)*                  &
     &                             umask(Istr,j)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isUvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=Jstr,Jend
                d_u_obc(j,k,ib,ir)=ad_u_obc(j,k,ib,ir,Lnew)
#   ifdef MASKING
                d_u_obc(j,k,ib,ir)=d_u_obc(j,k,ib,ir)*                  &
     &                             umask(Iend+1,j)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                d_u_obc(i,k,ib,ir)=ad_u_obc(i,k,ib,ir,Lnew)
#   ifdef MASKING
                d_u_obc(i,k,ib,ir)=d_u_obc(i,k,ib,ir)*                  &
     &                             umask(i,Jstr-1)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                d_u_obc(i,k,ib,ir)=ad_u_obc(i,k,ib,ir,Lnew)
#   ifdef MASKING
                d_u_obc(i,k,ib,ir)=d_u_obc(i,k,ib,ir)*                  &
     &                             umask(i,Jend+1)
#   endif
              END DO
            END DO
          END IF
        END DO
      END IF
#  endif
!
!  3D V-momentum.
!
      DO k=1,N(ng)
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            d_v(i,j,k)=ad_v(i,j,k,Lnew)
#  ifdef MASKING
            d_v(i,j,k)=d_v(i,j,k)*vmask(i,j)
#  endif
          END DO
        END DO
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  3D V-momentum open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isVvel,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO k=1,N(ng)
              DO j=JstrV,Jend
                d_v_obc(j,k,ib,ir)=ad_v_obc(j,k,ib,ir,Lnew)
#   ifdef MASKING
                d_v_obc(j,k,ib,ir)=d_v_obc(j,k,ib,ir)*                  &
     &                             vmask(Istr-1,j)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isVvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=JstrV,Jend
                d_v_obc(j,k,ib,ir)=ad_v_obc(j,k,ib,ir,Lnew)
#   ifdef MASKING
                d_v_obc(j,k,ib,ir)=d_v_obc(j,k,ib,ir)*                  &
     &                             vmask(Iend+1,j)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=Istr,Iend
                d_v_obc(i,k,ib,ir)=ad_v_obc(i,k,ib,ir,Lnew)
#   ifdef MASKING
                d_v_obc(i,k,ib,ir)=d_v_obc(i,k,ib,ir)*                  &
     &                             vmask(i,Jstr)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=Istr,Iend
                d_v_obc(i,k,ib,ir)=ad_v_obc(i,k,ib,ir,Lnew)
#   ifdef MASKING
                d_v_obc(i,k,ib,ir)=d_v_obc(i,k,ib,ir)*                  &
     &                             vmask(i,Jend+1)
#   endif
              END DO
            END DO
          END IF
        END DO
      END IF
#  endif
!
!  Tracers.
!
      DO it=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              d_t(i,j,k,it)=ad_t(i,j,k,Lnew,it)
#  ifdef MASKING
              d_t(i,j,k,it)=d_t(i,j,k,it)*rmask(i,j)
#  endif
            END DO
          END DO
        END DO
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  Tracers open boundaries.
!
      DO it=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(it),ng))) THEN
          DO ir=1,Nbrec(ng)
            IF ((Lobc(iwest,isTvar(it),ng)).and.                        &
     &          DOMAIN(ng)%Western_Edge(tile)) THEN
              ib=iwest
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  d_t_obc(j,k,ib,ir,it)=ad_t_obc(j,k,ib,ir,Lnew,it)
#   ifdef MASKING
                  d_t_obc(j,k,ib,ir,it)=d_t_obc(j,k,ib,ir,it)*          &
     &                                  rmask(Istr-1,j)
#   endif
                END DO
              END DO
            END IF
            IF ((Lobc(ieast,isTvar(it),ng)).and.                        &
     &          DOMAIN(ng)%Eastern_Edge(tile)) THEN
              ib=ieast
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  d_t_obc(j,k,ib,ir,it)=ad_t_obc(j,k,ib,ir,Lnew,it)
#   ifdef MASKING
                  d_t_obc(j,k,ib,ir,it)=d_t_obc(j,k,ib,ir,it)*          &
     &                                  rmask(Iend+1,j)
#   endif
                END DO
              END DO
            END IF
            IF ((Lobc(isouth,isTvar(it),ng)).and.                       &
     &          DOMAIN(ng)%Southern_Edge(tile)) THEN
              ib=isouth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  d_t_obc(i,k,ib,ir,it)=ad_t_obc(i,k,ib,ir,Lnew,it)
#   ifdef MASKING
                  d_t_obc(i,k,ib,ir,it)=d_t_obc(i,k,ib,ir,it)*          &
     &                                  rmask(i,Jstr-1)
#   endif
                END DO
              END DO
            END IF
            IF ((Lobc(inorth,isTvar(it),ng)).and.                       &
     &          DOMAIN(ng)%Northern_Edge(tile)) THEN
              ib=inorth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  d_t_obc(i,k,ib,ir,it)=ad_t_obc(i,k,ib,ir,Lnew,it)
#   ifdef MASKING
                  d_t_obc(i,k,ib,ir,it)=d_t_obc(i,k,ib,ir,it)*          &
     &                                  rmask(i,Jend+1)
#   endif
                END DO
              END DO
            END IF
          END DO
        END IF
      END DO
#  endif

#  ifdef ADJUST_STFLUX
!
!  Surface tracers flux.
!
      DO it=1,NT(ng)
        IF (Lstflux(it,ng)) THEN
          DO ir=1,Nfrec(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                d_stflx(i,j,ir,it)=ad_tflux(i,j,ir,Lnew,it)
#   ifdef MASKING
                d_stflx(i,j,ir,it)=d_stflx(i,j,ir,it)*rmask(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
      END DO
#  endif
# endif

      RETURN
      END SUBROUTINE new_direction
!
!***********************************************************************
      SUBROUTINE hessian (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj, LBij, UBij,               &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    Lold, Lnew, Lwrk,                             &
     &                    innLoop, outLoop,                             &
# ifdef MASKING
     &                    rmask, umask, vmask,                          &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    ad_t_obc, ad_u_obc, ad_v_obc,                 &
#  endif
     &                    ad_ubar_obc, ad_vbar_obc,                     &
     &                    ad_zeta_obc,                                  &
# endif
# ifdef ADJUST_WSTRESS
     &                    ad_ustr, ad_vstr,                             &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    ad_tflux,                                     &
#  endif
     &                    ad_t, ad_u, ad_v,                             &
# else
     &                    ad_ubar, ad_vbar,                             &
# endif
     &                    ad_zeta,                                      &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    tl_t_obc, tl_u_obc, tl_v_obc,                 &
#  endif
     &                    tl_ubar_obc, tl_vbar_obc,                     &
     &                    tl_zeta_obc,                                  &
# endif
# ifdef ADJUST_WSTRESS
     &                    tl_ustr, tl_vstr,                             &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    tl_tflux,                                     &
#  endif
     &                    tl_t, tl_u, tl_v,                             &
# else
     &                    tl_ubar, tl_vbar,                             &
# endif
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
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lold, Lnew, Lwrk
      integer, intent(in) :: innLoop, outLoop
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: ad_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: ad_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: tl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: tl_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: ad_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: tl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k
      integer :: ib, ir, it

      real(r8) :: fac

      real(r8), dimension(0:NstateVar(ng)) :: dot

      character (len=256) :: ncname

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Estimate the action of the Hessian according to:
!
!    grad(v) = Hv + grad(0)
!
!  where grad(v) is the gradient for current value of v, and grad(0)
!  is the gradient on the first inner-loop when v=0. Therefore,
!
!    Hv = grad(v) - grad(0).
!-----------------------------------------------------------------------
!
!  Need to multiply the adjoint state arrays (index Lold) by zgnorm to
!  convert back to the non-normalized gradient.
!
!  Free-surface.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ad_zeta(i,j,Lnew)=ad_zeta(i,j,Lnew)-                          &
     &                      ad_zeta(i,j,Lold)*                          &
     &                      cg_Gnorm(outLoop)
# ifdef MASKING
          ad_zeta(i,j,Lnew)=ad_zeta(i,j,Lnew)*rmask(i,j)
# endif
        END DO
      END DO

# ifdef ADJUST_BOUNDARY
!
!  Free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isFsur,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO j=Jstr,Jend
              ad_zeta_obc(j,ib,ir,Lnew)=ad_zeta_obc(j,ib,ir,Lnew)-      &
     &                                  ad_zeta_obc(j,ib,ir,Lold)*      &
     &                                  cg_Gnorm(outLoop)
#  ifdef MASKING
              ad_zeta_obc(j,ib,ir,Lnew)=ad_zeta_obc(j,ib,ir,Lnew)*      &
     &                                  rmask(Istr-1,j)
#  endif
            END DO
          END IF
          IF ((Lobc(ieast,isFsur,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=Jstr,Jend
              ad_zeta_obc(j,ib,ir,Lnew)=ad_zeta_obc(j,ib,ir,Lnew)-      &
     &                                  ad_zeta_obc(j,ib,ir,Lold)*      &
     &                                  cg_Gnorm(outLoop)
#  ifdef MASKING
              ad_zeta_obc(j,ib,ir,Lnew)=ad_zeta_obc(j,ib,ir,Lnew)*      &
     &                                  rmask(Iend+1,j)
#  endif
            END DO
          END IF
          IF ((Lobc(isouth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=Istr,Iend
              ad_zeta_obc(i,ib,ir,Lnew)=ad_zeta_obc(i,ib,ir,Lnew)-      &
     &                                  ad_zeta_obc(i,ib,ir,Lold)*      &
     &                                  cg_Gnorm(outLoop)
#  ifdef MASKING
              ad_zeta_obc(i,ib,ir,Lnew)=ad_zeta_obc(i,ib,ir,Lnew)*      &
     &                                  rmask(i,Jstr-1)
#  endif
            END DO
          END IF
          IF ((Lobc(inorth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=Istr,Iend
              ad_zeta_obc(i,ib,ir,Lnew)=ad_zeta_obc(i,ib,ir,Lnew)-      &
     &                                  ad_zeta_obc(i,ib,ir,Lold)*      &
     &                                  cg_Gnorm(outLoop)
#  ifdef MASKING
              ad_zeta_obc(i,ib,ir,Lnew)=ad_zeta_obc(i,ib,ir,Lnew)*      &
     &                                  rmask(i,Jend+1)
#  endif
            END DO
          END IF
        END DO
      END IF
# endif

# ifndef SOLVE3D
!
!  2D U-momentum.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ad_ubar(i,j,Lnew)=ad_ubar(i,j,Lnew)-                          &
     &                      ad_ubar(i,j,Lold)*                          &
     &                      cg_Gnorm(outLoop)
#  ifdef MASKING
          ad_ubar(i,j,Lnew)=ad_ubar(i,j,Lnew)*umask(i,j)
#  endif
        END DO
      END DO
# endif

# ifdef ADJUST_BOUNDARY
!
!  2D U-momentum open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isUbar,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO j=Jstr,Jend
              ad_ubar_obc(j,ib,ir,Lnew)=ad_ubar_obc(j,ib,ir,Lnew)-      &
     &                                  ad_ubar_obc(j,ib,ir,Lold)*      &
     &                                  cg_Gnorm(outLoop)
#  ifdef MASKING
              ad_ubar_obc(j,ib,ir,Lnew)=ad_ubar_obc(j,ib,ir,Lnew)*      &
     &                                  umask(Istr,j)
#  endif
            END DO
          END IF
          IF ((Lobc(ieast,isUbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=Jstr,Jend
              ad_ubar_obc(j,ib,ir,Lnew)=ad_ubar_obc(j,ib,ir,Lnew)-      &
     &                                  ad_ubar_obc(j,ib,ir,Lold)*      &
     &                                  cg_Gnorm(outLoop)
#  ifdef MASKING
              ad_ubar_obc(j,ib,ir,Lnew)=ad_ubar_obc(j,ib,ir,Lnew)*      &
     &                                  umask(Iend+1,j)
#  endif
            END DO
          END IF
          IF ((Lobc(isouth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=IstrU,Iend
              ad_ubar_obc(i,ib,ir,Lnew)=ad_ubar_obc(i,ib,ir,Lnew)-      &
     &                                  ad_ubar_obc(i,ib,ir,Lold)*      &
     &                                  cg_Gnorm(outLoop)
#  ifdef MASKING
              ad_ubar_obc(i,ib,ir,Lnew)=ad_ubar_obc(i,ib,ir,Lnew)*      &
     &                                  umask(i,Jstr-1)
#  endif
            END DO
          END IF
          IF ((Lobc(inorth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=IstrU,Iend
              ad_ubar_obc(i,ib,ir,Lnew)=ad_ubar_obc(i,ib,ir,Lnew)-      &
     &                                  ad_ubar_obc(i,ib,ir,Lold)*      &
     &                                  cg_Gnorm(outLoop)
#  ifdef MASKING
              ad_ubar_obc(i,ib,ir,Lnew)=ad_ubar_obc(i,ib,ir,Lnew)*      &
     &                                  umask(i,Jend+1)
#  endif
            END DO
          END IF
        END DO
      END IF
# endif

# ifndef SOLVE3D
!
!  2D V-momentum.
!
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          ad_vbar(i,j,Lnew)=ad_vbar(i,j,Lnew)-                          &
     &                      ad_vbar(i,j,Lold)*                          &
     &                      cg_Gnorm(outLoop)
#  ifdef MASKING
          ad_vbar(i,j,Lnew)=ad_vbar(i,j,Lnew)*vmask(i,j)
#  endif
        END DO
      END DO
# endif

# ifdef ADJUST_BOUNDARY
!
!  2D V-momentum open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isVbar,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO j=JstrV,Jend
              ad_vbar_obc(j,ib,ir,Lnew)=ad_vbar_obc(j,ib,ir,Lnew)-      &
     &                                  ad_vbar_obc(j,ib,ir,Lold)*      &
     &                                  cg_Gnorm(outLoop)
#  ifdef MASKING
              ad_vbar_obc(j,ib,ir,Lnew)=ad_vbar_obc(j,ib,ir,Lnew)*      &
     &                                  vmask(Istr-1,j)
#  endif
            END DO
          END IF
          IF ((Lobc(ieast,isVbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=JstrV,Jend
              ad_vbar_obc(j,ib,ir,Lnew)=ad_vbar_obc(j,ib,ir,Lnew)-      &
     &                                  ad_vbar_obc(j,ib,ir,Lold)*      &
     &                                  cg_Gnorm(outLoop)
#  ifdef MASKING
              ad_vbar_obc(j,ib,ir,Lnew)=ad_vbar_obc(j,ib,ir,Lnew)*      &
     &                                  vmask(Iend+1,j)
#  endif
            END DO
          END IF
          IF ((Lobc(isouth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=Istr,Iend
              ad_vbar_obc(i,ib,ir,Lnew)=ad_vbar_obc(i,ib,ir,Lnew)-      &
     &                                  ad_vbar_obc(i,ib,ir,Lold)*      &
     &                                  cg_Gnorm(outLoop)
#  ifdef MASKING
              ad_vbar_obc(i,ib,ir,Lnew)=ad_vbar_obc(i,ib,ir,Lnew)*      &
     &                                  vmask(i,Jstr)
#  endif
            END DO
          END IF
          IF ((Lobc(inorth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=Istr,Iend
              ad_vbar_obc(i,ib,ir,Lnew)=ad_vbar_obc(i,ib,ir,Lnew)-      &
     &                                  ad_vbar_obc(i,ib,ir,Lold)*      &
     &                                  cg_Gnorm(outLoop)
#  ifdef MASKING
              ad_vbar_obc(i,ib,ir,Lnew)=ad_vbar_obc(i,ib,ir,Lnew)*      &
     &                                  vmask(i,Jend+1)
#  endif
            END DO
          END IF
        END DO
      END IF
# endif

# ifdef ADJUST_WSTRESS
!
!  Surface momentum stress.
!
      DO ir=1,Nfrec(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_ustr(i,j,ir,Lnew)=ad_ustr(i,j,ir,Lnew)-                  &
     &                           ad_ustr(i,j,ir,Lold)*                  &
     &                           cg_Gnorm(outLoop)
#  ifdef MASKING
            ad_ustr(i,j,ir,Lnew)=ad_ustr(i,j,ir,Lnew)*umask(i,j)
#  endif
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_vstr(i,j,ir,Lnew)=ad_vstr(i,j,ir,Lnew)-                  &
     &                           ad_vstr(i,j,ir,Lold)*                  &
     &                           cg_Gnorm(outLoop)
#  ifdef MASKING
            ad_vstr(i,j,ir,Lnew)=ad_vstr(i,j,ir,Lnew)*vmask(i,j)
#  endif
          END DO
        END DO
      END DO
# endif

# ifdef SOLVE3D
!
!  3D U-momentum.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_u(i,j,k,Lnew)=ad_u(i,j,k,Lnew)-                          &
     &                       ad_u(i,j,k,Lold)*                          &
     &                       cg_Gnorm(outLoop)
#  ifdef MASKING
            ad_u(i,j,k,Lnew)=ad_u(i,j,k,Lnew)*umask(i,j)
#  endif
          END DO
        END DO
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  3D U-momentum open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isUvel,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO k=1,N(ng)
              DO j=Jstr,Jend
                ad_u_obc(j,k,ib,ir,Lnew)=ad_u_obc(j,k,ib,ir,Lnew)-      &
     &                                   ad_u_obc(j,k,ib,ir,Lold)*      &
     &                                   cg_Gnorm(outLoop)
#   ifdef MASKING
                ad_u_obc(j,k,ib,ir,Lnew)=ad_u_obc(j,k,ib,ir,Lnew)*      &
     &                                   umask(Istr,j)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isUvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=Jstr,Jend
                ad_u_obc(j,k,ib,ir,Lnew)=ad_u_obc(j,k,ib,ir,Lnew)-      &
     &                                   ad_u_obc(j,k,ib,ir,Lold)*      &
     &                                   cg_Gnorm(outLoop)
#   ifdef MASKING
                ad_u_obc(j,k,ib,ir,Lnew)=ad_u_obc(j,k,ib,ir,Lnew)*      &
     &                                   umask(Iend+1,j)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                ad_u_obc(i,k,ib,ir,Lnew)=ad_u_obc(i,k,ib,ir,Lnew)-      &
     &                                   ad_u_obc(i,k,ib,ir,Lold)*      &
     &                                   cg_Gnorm(outLoop)
#   ifdef MASKING
                ad_u_obc(i,k,ib,ir,Lnew)=ad_u_obc(i,k,ib,ir,Lnew)*      &
     &                                   umask(i,Jstr-1)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                ad_u_obc(i,k,ib,ir,Lnew)=ad_u_obc(i,k,ib,ir,Lnew)-      &
     &                                   ad_u_obc(i,k,ib,ir,Lold)*      &
     &                                   cg_Gnorm(outLoop)
#   ifdef MASKING
                ad_u_obc(i,k,ib,ir,Lnew)=ad_u_obc(i,k,ib,ir,Lnew)*      &
     &                                   umask(i,Jend+1)
#   endif
              END DO
            END DO
          END IF
        END DO
      END IF
#  endif
!
!  3D V-momentum.
!
      DO k=1,N(ng)
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_v(i,j,k,Lnew)=ad_v(i,j,k,Lnew)-                          &
     &                       ad_v(i,j,k,Lold)*                          &
     &                       cg_Gnorm(outLoop)
#  ifdef MASKING
            ad_v(i,j,k,Lnew)=ad_v(i,j,k,Lnew)*vmask(i,j)
#  endif
          END DO
        END DO
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  3D V-momentum open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isVvel,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO k=1,N(ng)
              DO j=JstrV,Jend
                ad_v_obc(j,k,ib,ir,Lnew)=ad_v_obc(j,k,ib,ir,Lnew)-      &
     &                                   ad_v_obc(j,k,ib,ir,Lold)*      &
     &                                   cg_Gnorm(outLoop)
#   ifdef MASKING
                ad_v_obc(j,k,ib,ir,Lnew)=ad_v_obc(j,k,ib,ir,Lnew)*      &
     &                                   vmask(Istr-1,j)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isVvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=JstrV,Jend
                ad_v_obc(j,k,ib,ir,Lnew)=ad_v_obc(j,k,ib,ir,Lnew)-      &
     &                                   ad_v_obc(j,k,ib,ir,Lold)*      &
     &                                   cg_Gnorm(outLoop)
#   ifdef MASKING
                ad_v_obc(j,k,ib,ir,Lnew)=ad_v_obc(j,k,ib,ir,Lnew)*      &
     &                                   vmask(Iend+1,j)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=Istr,Iend
                ad_v_obc(i,k,ib,ir,Lnew)=ad_v_obc(i,k,ib,ir,Lnew)-      &
     &                                   ad_v_obc(i,k,ib,ir,Lold)*      &
     &                                   cg_Gnorm(outLoop)
#   ifdef MASKING
                ad_v_obc(i,k,ib,ir,Lnew)=ad_v_obc(i,k,ib,ir,Lnew)*      &
     &                                   vmask(i,Jstr)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=Istr,Iend
                ad_v_obc(i,k,ib,ir,Lnew)=ad_v_obc(i,k,ib,ir,Lnew)-      &
     &                                   ad_v_obc(i,k,ib,ir,Lold)*      &
     &                                   cg_Gnorm(outLoop)
#   ifdef MASKING
                ad_v_obc(i,k,ib,ir,Lnew)=ad_v_obc(i,k,ib,ir,Lnew)*      &
     &                                   vmask(i,Jend+1)
#   endif
              END DO
            END DO
          END IF
        END DO
      END IF
#  endif
!
!  Tracers.
!
      DO it=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              ad_t(i,j,k,Lnew,it)=ad_t(i,j,k,Lnew,it)-                  &
     &                            ad_t(i,j,k,Lold,it)*                  &
     &                            cg_Gnorm(outLoop)
#  ifdef MASKING
              ad_t(i,j,k,Lnew,it)=ad_t(i,j,k,Lnew,it)*rmask(i,j)
#  endif
            END DO
          END DO
        END DO
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  Tracers open boundaries.
!
      DO it=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(it),ng))) THEN
          DO ir=1,Nbrec(ng)
            IF ((Lobc(iwest,isTvar(it),ng)).and.                        &
     &          DOMAIN(ng)%Western_Edge(tile)) THEN
              ib=iwest
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  ad_t_obc(j,k,ib,ir,Lnew,it)=                          &
     &                                    ad_t_obc(j,k,ib,ir,Lnew,it)-  &
     &                                    ad_t_obc(j,k,ib,ir,Lold,it)*  &
     &                                    cg_Gnorm(outLoop)
#   ifdef MASKING
                  ad_t_obc(j,k,ib,ir,Lnew,it)=                          &
     &                    ad_t_obc(j,k,ib,ir,Lnew,it)*rmask(Istr-1,j)
#   endif
                END DO
              END DO
            END IF
            IF ((Lobc(ieast,isTvar(it),ng)).and.                        &
     &          DOMAIN(ng)%Eastern_Edge(tile)) THEN
              ib=ieast
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  ad_t_obc(j,k,ib,ir,Lnew,it)=                          &
     &                                    ad_t_obc(j,k,ib,ir,Lnew,it)-  &
     &                                    ad_t_obc(j,k,ib,ir,Lold,it)*  &
     &                                    cg_Gnorm(outLoop)
#   ifdef MASKING
                  ad_t_obc(j,k,ib,ir,Lnew,it)=                          &
     &                    ad_t_obc(j,k,ib,ir,Lnew,it)*rmask(Iend+1,j)
#   endif
                END DO
              END DO
            END IF
            IF ((Lobc(isouth,isTvar(it),ng)).and.                       &
     &          DOMAIN(ng)%Southern_Edge(tile)) THEN
              ib=isouth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  ad_t_obc(i,k,ib,ir,Lnew,it)=                          &
     &                                    ad_t_obc(i,k,ib,ir,Lnew,it)-  &
     &                                    ad_t_obc(i,k,ib,ir,Lold,it)*  &
     &                                    cg_Gnorm(outLoop)
#   ifdef MASKING
                  ad_t_obc(i,k,ib,ir,Lnew,it)=                          &
     &                    ad_t_obc(i,k,ib,ir,Lnew,it)*rmask(i,Jstr-1)
#   endif
                END DO
              END DO
            END IF
            IF ((Lobc(inorth,isTvar(it),ng)).and.                       &
     &          DOMAIN(ng)%Northern_Edge(tile)) THEN
              ib=inorth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  ad_t_obc(i,k,ib,ir,Lnew,it)=                          &
     &                                    ad_t_obc(i,k,ib,ir,Lnew,it)-  &
     &                                    ad_t_obc(i,k,ib,ir,Lold,it)*  &
     &                                    cg_Gnorm(outLoop)
#   ifdef MASKING
                  ad_t_obc(i,k,ib,ir,Lnew,it)=                          &
     &                    ad_t_obc(i,k,ib,ir,Lnew,it)*rmask(i,Jend+1)
#   endif
                END DO
              END DO
            END IF
          END DO
        END IF
      END DO
#  endif

#  ifdef ADJUST_STFLUX
!
!  Surface tracers flux.
!
      DO it=1,NT(ng)
        IF (Lstflux(it,ng)) THEN
          DO ir=1,Nfrec(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                ad_tflux(i,j,ir,Lnew,it)=ad_tflux(i,j,ir,Lnew,it)-      &
     &                                   ad_tflux(i,j,ir,Lold,it)*      &
     &                                   cg_Gnorm(outLoop)
#   ifdef MASKING
                ad_tflux(i,j,ir,Lnew,it)=ad_tflux(i,j,ir,Lnew,it)*      &
     &                                   rmask(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
      END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Compute norm Delta(k) as the dot product between the new gradient
!  and current iteration gradient solution.
!-----------------------------------------------------------------------
!
!  Determine gradient file to process.
!
      IF (ndefADJ(ng).gt.0) THEN
        WRITE (ncname,10) TRIM(ADM(ng)%base), outLoop
 10     FORMAT (a,'_',i3.3,'.nc')
      ELSE
        ncname=ADM(ng)%name
      END IF
!
!  Read Lanczos vector on which the Hessian matrix is operating
!  into tangent linear state array, index Lwrk.
!
      CALL read_state (ng, tile, model,                                 &
     &                 LBi, UBi, LBj, UBj, LBij, UBij,                  &
     &                 Lwrk, innLoop,                                   &
     &                 ndefADJ(ng), ADM(ng)%ncid, ncname,               &
# ifdef MASKING
     &                 rmask, umask, vmask,                             &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                 tl_t_obc, tl_u_obc, tl_v_obc,                    &
#  endif
     &                 tl_ubar_obc, tl_vbar_obc,                        &
     &                 tl_zeta_obc,                                     &
# endif
# ifdef ADJUST_WSTRESS
     &                 tl_ustr, tl_vstr,                                &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                 tl_tflux,                                        &
#  endif
     &                 tl_t, tl_u, tl_v,                                &
# else
     &                 tl_ubar, tl_vbar,                                &
# endif
     &                 tl_zeta)
      IF (exit_flag.ne.NoError) RETURN
!
!  Compute current iteration norm Delta(k) used to compute tri-diagonal
!  matrix T(k) in the Lanczos recurrence.
!
      CALL state_dotprod (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj, LBij, UBij,               &
     &                    NstateVar(ng), dot(0:),                       &
# ifdef MASKING
     &                    rmask, umask, vmask,                          &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    ad_t_obc(:,:,:,:,Lnew,:),                     &
     &                    tl_t_obc(:,:,:,:,Lwrk,:),                     &
     &                    ad_u_obc(:,:,:,:,Lnew),                       &
     &                    tl_u_obc(:,:,:,:,Lwrk),                       &
     &                    ad_v_obc(:,:,:,:,Lnew),                       &
     &                    tl_v_obc(:,:,:,:,Lwrk),                       &
#  endif
     &                    ad_ubar_obc(:,:,:,Lnew),                      &
     &                    tl_ubar_obc(:,:,:,Lwrk),                      &
     &                    ad_vbar_obc(:,:,:,Lnew),                      &
     &                    tl_vbar_obc(:,:,:,Lwrk),                      &
     &                    ad_zeta_obc(:,:,:,Lnew),                      &
     &                    tl_zeta_obc(:,:,:,Lwrk),                      &
# endif
# ifdef ADJUST_WSTRESS
     &                    ad_ustr(:,:,:,Lnew), tl_ustr(:,:,:,Lwrk),     &
     &                    ad_vstr(:,:,:,Lnew), tl_vstr(:,:,:,Lwrk),     &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    ad_tflux(:,:,:,Lnew,:),                       &
     &                    tl_tflux(:,:,:,Lwrk,:),                       &
#  endif
     &                    ad_t(:,:,:,Lnew,:), tl_t(:,:,:,Lwrk,:),       &
     &                    ad_u(:,:,:,Lnew), tl_u(:,:,:,Lwrk),           &
     &                    ad_v(:,:,:,Lnew), tl_v(:,:,:,Lwrk),           &
# else
     &                    ad_ubar(:,:,Lnew), tl_ubar(:,:,Lwrk),         &
     &                    ad_vbar(:,:,Lnew), tl_vbar(:,:,Lwrk),         &
# endif
     &                    ad_zeta(:,:,Lnew), tl_zeta(:,:,Lwrk))

      cg_delta(innLoop,outLoop)=dot(0)

      RETURN
      END SUBROUTINE hessian
!
!***********************************************************************
      SUBROUTINE lanczos (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj, LBij, UBij,               &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    Lold, Lnew, Lwrk,                             &
     &                    innLoop, outLoop,                             &
# ifdef MASKING
     &                    rmask, umask, vmask,                          &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    tl_t_obc, tl_u_obc, tl_v_obc,                 &
#  endif
     &                    tl_ubar_obc, tl_vbar_obc,                     &
     &                    tl_zeta_obc,                                  &
# endif
# ifdef ADJUST_WSTRESS
     &                    tl_ustr, tl_vstr,                             &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    tl_tflux,                                     &
#  endif
     &                    tl_t, tl_u, tl_v,                             &
# else
     &                    tl_ubar, tl_vbar,                             &
# endif
     &                    tl_zeta,                                      &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    ad_t_obc, ad_u_obc, ad_v_obc,                 &
#  endif
     &                    ad_ubar_obc, ad_vbar_obc,                     &
     &                    ad_zeta_obc,                                  &
# endif
# ifdef ADJUST_WSTRESS
     &                    ad_ustr, ad_vstr,                             &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    ad_tflux,                                     &
#  endif
     &                    ad_t, ad_u, ad_v,                             &
# else
     &                    ad_ubar, ad_vbar,                             &
# endif
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
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lold, Lnew, Lwrk
      integer, intent(in) :: innLoop, outLoop
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: ad_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: ad_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: tl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: tl_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: ad_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: tl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, rec

      real(r8) :: fac, fac1, fac2

      real(r8), dimension(0:NstateVar(ng)) :: dot
      real(r8), dimension(0:Ninner) :: DotProd, dot_new, dot_old

      character (len=256) :: ncname

# include "set_bounds.h"
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
!  tangent linear state arrays (index Lwrk) - it was read in the
!  routine hessian.
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
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       Lnew, Lwrk, Lnew, fac1, fac2,              &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       ad_t_obc, tl_t_obc,                        &
     &                       ad_u_obc, tl_u_obc,                        &
     &                       ad_v_obc, tl_v_obc,                        &
#  endif
     &                       ad_ubar_obc, tl_ubar_obc,                  &
     &                       ad_vbar_obc, tl_vbar_obc,                  &
     &                       ad_zeta_obc, tl_zeta_obc,                  &
# endif
# ifdef ADJUST_WSTRESS
     &                       ad_ustr, tl_ustr,                          &
     &                       ad_vstr, tl_vstr,                          &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                       ad_tflux, tl_tflux,                        &
#  endif
     &                       ad_t, tl_t,                                &
     &                       ad_u, tl_u,                                &
     &                       ad_v, tl_v,                                &
# else
     &                       ad_ubar, tl_ubar,                          &
     &                       ad_vbar, tl_vbar,                          &
# endif
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
          WRITE (ncname,10) TRIM(ADM(ng)%base), outLoop
 10       FORMAT (a,'_',i3.3,'.nc')
        ELSE
          ncname=ADM(ng)%name
        END IF
!
!  Read in the previous (innLoop-1) orthonormal Lanczos vector.
!
        CALL read_state (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Lwrk, innLoop-1,                               &
     &                   ndefADJ(ng), ADM(ng)%ncid, ncname,             &
# ifdef MASKING
     &                   rmask, umask, vmask,                           &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   tl_t_obc, tl_u_obc, tl_v_obc,                  &
#  endif
     &                   tl_ubar_obc, tl_vbar_obc,                      &
     &                   tl_zeta_obc,                                   &
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
        IF (exit_flag.ne.NoError) RETURN
!
!  Substract previous orthonormal Lanczos vector:
!
!    ad_var(Lnew) = fac1 * ad_var(Lnew) + fac2 * tl_var(Lwrk)
!
        fac1=1.0_r8
        fac2=-cg_beta(innLoop,outLoop)

        CALL state_addition (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       Lnew, Lwrk, Lnew, fac1, fac2,              &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       ad_t_obc, tl_t_obc,                        &
     &                       ad_u_obc, tl_u_obc,                        &
     &                       ad_v_obc, tl_v_obc,                        &
#  endif
     &                       ad_ubar_obc, tl_ubar_obc,                  &
     &                       ad_vbar_obc, tl_vbar_obc,                  &
     &                       ad_zeta_obc, tl_zeta_obc,                  &
# endif
# ifdef ADJUST_WSTRESS
     &                       ad_ustr, tl_ustr,                          &
     &                       ad_vstr, tl_vstr,                          &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                       ad_tflux, tl_tflux,                        &
#  endif
     &                       ad_t, tl_t,                                &
     &                       ad_u, tl_u,                                &
     &                       ad_v, tl_v,                                &
# else
     &                       ad_ubar, tl_ubar,                          &
     &                       ad_vbar, tl_vbar,                          &
# endif
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
      IF (ndefADJ(ng).gt.0) THEN
        WRITE (ncname,10) TRIM(ADM(ng)%base), outLoop
      ELSE
        ncname=ADM(ng)%name
      END IF
!
      DO rec=innLoop,1,-1
!
!  Read in each previous gradient state solutions, G(0) to G(k), and
!  compute its associated dot angaint curret G(k+1). Each gradient
!  solution is loaded into TANGENT LINEAR STATE ARRAYS at index Lwrk.
!
        CALL read_state (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Lwrk, rec,                                     &
     &                   ndefADJ(ng), ADM(ng)%ncid, ncname,             &
# ifdef MASKING
     &                   rmask, umask, vmask,                           &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   tl_t_obc, tl_u_obc, tl_v_obc,                  &
#  endif
     &                   tl_ubar_obc, tl_vbar_obc,                      &
     &                   tl_zeta_obc,                                   &
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
        IF (exit_flag.ne.NoError) RETURN
!
!  Compute dot product <q(k+1), q(rec)>.
!
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj, LBij, UBij,             &
     &                      NstateVar(ng), dot(0:),                     &
# ifdef MASKING
     &                      rmask, umask, vmask,                        &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                      ad_t_obc(:,:,:,:,Lnew,:),                   &
     &                      tl_t_obc(:,:,:,:,Lwrk,:),                   &
     &                      ad_u_obc(:,:,:,:,Lnew),                     &
     &                      tl_u_obc(:,:,:,:,Lwrk),                     &
     &                      ad_v_obc(:,:,:,:,Lnew),                     &
     &                      tl_v_obc(:,:,:,:,Lwrk),                     &
#  endif
     &                      ad_ubar_obc(:,:,:,Lnew),                    &
     &                      tl_ubar_obc(:,:,:,Lwrk),                    &
     &                      ad_vbar_obc(:,:,:,Lnew),                    &
     &                      tl_vbar_obc(:,:,:,Lwrk),                    &
     &                      ad_zeta_obc(:,:,:,Lnew),                    &
     &                      tl_zeta_obc(:,:,:,Lwrk),                    &
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
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       Lnew, Lwrk, Lnew, fac1, fac2,              &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       ad_t_obc, tl_t_obc,                        &
     &                       ad_u_obc, tl_u_obc,                        &
     &                       ad_v_obc, tl_v_obc,                        &
#  endif
     &                       ad_ubar_obc, tl_ubar_obc,                  &
     &                       ad_vbar_obc, tl_vbar_obc,                  &
     &                       ad_zeta_obc, tl_zeta_obc,                  &
# endif
# ifdef ADJUST_WSTRESS
     &                       ad_ustr, tl_ustr,                          &
     &                       ad_vstr, tl_vstr,                          &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                       ad_tflux, tl_tflux,                        &
#  endif
     &                       ad_t, tl_t,                                &
     &                       ad_u, tl_u,                                &
     &                       ad_v, tl_v,                                &
# else
     &                       ad_ubar, tl_ubar,                          &
     &                       ad_vbar, tl_vbar,                          &
# endif
     &                       ad_zeta, tl_zeta)
      END DO
!
!-----------------------------------------------------------------------
!  Normalize current orthogonal gradient vector.
!-----------------------------------------------------------------------
!
      CALL state_dotprod (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj, LBij, UBij,               &
     &                    NstateVar(ng), dot(0:),                       &
# ifdef MASKING
     &                    rmask, umask, vmask,                          &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    ad_t_obc(:,:,:,:,Lnew,:),                     &
     &                    ad_t_obc(:,:,:,:,Lnew,:),                     &
     &                    ad_u_obc(:,:,:,:,Lnew),                       &
     &                    ad_u_obc(:,:,:,:,Lnew),                       &
     &                    ad_v_obc(:,:,:,:,Lnew),                       &
     &                    ad_v_obc(:,:,:,:,Lnew),                       &
#  endif
     &                    ad_ubar_obc(:,:,:,Lnew),                      &
     &                    ad_ubar_obc(:,:,:,Lnew),                      &
     &                    ad_vbar_obc(:,:,:,Lnew),                      &
     &                    ad_vbar_obc(:,:,:,Lnew),                      &
     &                    ad_zeta_obc(:,:,:,Lnew),                      &
     &                    ad_zeta_obc(:,:,:,Lnew),                      &
# endif
# ifdef ADJUST_WSTRESS
     &                    ad_ustr(:,:,:,Lnew), ad_ustr(:,:,:,Lnew),     &
     &                    ad_vstr(:,:,:,Lnew), ad_vstr(:,:,:,Lnew),     &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    ad_tflux(:,:,:,Lnew,:),                       &
     &                    ad_tflux(:,:,:,Lnew,:),                       &
#  endif
     &                    ad_t(:,:,:,Lnew,:), ad_t(:,:,:,Lnew,:),       &
     &                    ad_u(:,:,:,Lnew), ad_u(:,:,:,Lnew),           &
     &                    ad_v(:,:,:,Lnew), ad_v(:,:,:,Lnew),           &
# else
     &                    ad_ubar(:,:,Lnew), ad_ubar(:,:,Lnew),         &
     &                    ad_vbar(:,:,Lnew), ad_vbar(:,:,Lnew),         &
# endif
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
     &                  LBi, UBi, LBj, UBj, LBij, UBij,                 &
     &                  Lnew, Lnew, fac,                                &
# ifdef MASKING
     &                  rmask, umask, vmask,                            &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                  ad_t_obc, ad_u_obc, ad_v_obc,                   &
#  endif
     &                  ad_ubar_obc, ad_vbar_obc,                       &
     &                  ad_zeta_obc,                                    &
# endif
# ifdef ADJUST_WSTRESS
     &                  ad_ustr, ad_vstr,                               &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                  ad_tflux,                                       &
#  endif
     &                  ad_t, ad_u, ad_v,                               &
# else
     &                  ad_ubar, ad_vbar,                               &
# endif
     &                  ad_zeta)
!
!-----------------------------------------------------------------------
!  Compute dot product of new Lanczos vector with gradient.
!-----------------------------------------------------------------------
!
      IF (innLoop.eq.0) THEN
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj, LBij, UBij,             &
     &                      NstateVar(ng), dot(0:),                     &
# ifdef MASKING
     &                      rmask, umask, vmask,                        &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                      ad_t_obc(:,:,:,:,Lnew,:),                   &
     &                      ad_t_obc(:,:,:,:,Lnew,:),                   &
     &                      ad_u_obc(:,:,:,:,Lnew),                     &
     &                      ad_u_obc(:,:,:,:,Lnew),                     &
     &                      ad_v_obc(:,:,:,:,Lnew),                     &
     &                      ad_v_obc(:,:,:,:,Lnew),                     &
#  endif
     &                      ad_ubar_obc(:,:,:,Lnew),                    &
     &                      ad_ubar_obc(:,:,:,Lnew),                    &
     &                      ad_vbar_obc(:,:,:,Lnew),                    &
     &                      ad_vbar_obc(:,:,:,Lnew),                    &
     &                      ad_zeta_obc(:,:,:,Lnew),                    &
     &                      ad_zeta_obc(:,:,:,Lnew),                    &
# endif
# ifdef ADJUST_WSTRESS
     &                      ad_ustr(:,:,:,Lnew), ad_ustr(:,:,:,Lnew),   &
     &                      ad_vstr(:,:,:,Lnew), ad_vstr(:,:,:,Lnew),   &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                      ad_tflux(:,:,:,Lnew,:),                     &
     &                      ad_tflux(:,:,:,Lnew,:),                     &
#  endif
     &                      ad_t(:,:,:,Lnew,:), ad_t(:,:,:,Lnew,:),     &
     &                      ad_u(:,:,:,Lnew), ad_u(:,:,:,Lnew),         &
     &                      ad_v(:,:,:,Lnew), ad_v(:,:,:,Lnew),         &
# else
     &                      ad_ubar(:,:,Lnew), ad_ubar(:,:,Lnew),       &
     &                      ad_vbar(:,:,Lnew), ad_vbar(:,:,Lnew),       &
# endif
     &                      ad_zeta(:,:,Lnew), ad_zeta(:,:,Lnew))
      ELSE
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj, LBij, UBij,             &
     &                      NstateVar(ng), dot(0:),                     &
# ifdef MASKING
     &                      rmask, umask, vmask,                        &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                      ad_t_obc(:,:,:,:,Lold,:),                   &
     &                      ad_t_obc(:,:,:,:,Lnew,:),                   &
     &                      ad_u_obc(:,:,:,:,Lold),                     &
     &                      ad_u_obc(:,:,:,:,Lnew),                     &
     &                      ad_v_obc(:,:,:,:,Lold),                     &
     &                      ad_v_obc(:,:,:,:,Lnew),                     &
#  endif
     &                      ad_ubar_obc(:,:,:,Lold),                    &
     &                      ad_ubar_obc(:,:,:,Lnew),                    &
     &                      ad_vbar_obc(:,:,:,Lold),                    &
     &                      ad_vbar_obc(:,:,:,Lnew),                    &
     &                      ad_zeta_obc(:,:,:,Lold),                    &
     &                      ad_zeta_obc(:,:,:,Lnew),                    &
# endif
# ifdef ADJUST_WSTRESS
     &                      ad_ustr(:,:,:,Lold), ad_ustr(:,:,:,Lnew),   &
     &                      ad_vstr(:,:,:,Lold), ad_vstr(:,:,:,Lnew),   &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                      ad_tflux(:,:,:,Lold,:),                     &
     &                      ad_tflux(:,:,:,Lnew,:),                     &
#  endif
     &                      ad_t(:,:,:,Lold,:), ad_t(:,:,:,Lnew,:),     &
     &                      ad_u(:,:,:,Lold), ad_u(:,:,:,Lnew),         &
     &                      ad_v(:,:,:,Lold), ad_v(:,:,:,Lnew),         &
# else
     &                      ad_ubar(:,:,Lold), ad_ubar(:,:,Lnew),       &
     &                      ad_vbar(:,:,Lold), ad_vbar(:,:,Lnew),       &
# endif
     &                      ad_zeta(:,:,Lold), ad_zeta(:,:,Lnew))
      ENDIF
!
!  Need to multiply dot(0) by zgnorm because the gradient (index Lold)
!  has been normalized.
!
      cg_QG(innLoop+1,outLoop)=cg_Gnorm(outLoop)*dot(0)

# ifdef TEST_ORTHOGONALIZATION
!
!-----------------------------------------------------------------------
!  Test orthogonal properties of the new gradient.
!-----------------------------------------------------------------------
!
!  Determine adjoint file to process.
!
      IF (ndefADJ(ng).gt.0) THEN
        WRITE (ncname,10) TRIM(ADM(ng)%base), outLoop
      ELSE
        ncname=ADM(ng)%name
      END IF
!
      DO rec=innLoop,1,-1
!
!  Read in each previous gradient state solutions, q(0) to q(k), and
!  compute its associated dot angaint orthogonalized q(k+1). Again,
!  each gradient solution is loaded into TANGENT LINEAR STATE ARRAYS
!  at index Lwrk.
!
        CALL read_state (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Lwrk, rec,                                     &
     &                   ndefADJ(ng), ADM(ng)%ncid, ncname,             &
#  ifdef MASKING
     &                   rmask, umask, vmask,                           &
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
     &                   tl_t_obc, tl_u_obc, tl_v_obc,                  &
#   endif
     &                   tl_ubar_obc, tl_vbar_obc,                      &
     &                   tl_zeta_obc,                                   &
#  endif
#  ifdef ADJUST_WSTRESS
     &                   tl_ustr, tl_vstr,                              &
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
     &                   tl_tflux,                                      &
#   endif
     &                   tl_t, tl_u, tl_v,                              &
#  else
     &                   tl_ubar, tl_vbar,                              &
#  endif
     &                   tl_zeta)
        IF (exit_flag.ne.NoError) RETURN
!
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj, LBij, UBij,             &
     &                      NstateVar(ng), dot(0:),                     &
#  ifdef MASKING
     &                      rmask, umask, vmask,                        &
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
     &                      ad_t_obc(:,:,:,:,Lnew,:),                   &
     &                      tl_t_obc(:,:,:,:,Lwrk,:),                   &
     &                      ad_u_obc(:,:,:,:,Lnew),                     &
     &                      tl_u_obc(:,:,:,:,Lwrk),                     &
     &                      ad_v_obc(:,:,:,:,Lnew),                     &
     &                      tl_v_obc(:,:,:,:,Lwrk),                     &
#   endif
     &                      ad_ubar_obc(:,:,:,Lnew),                    &
     &                      tl_ubar_obc(:,:,:,Lwrk),                    &
     &                      ad_vbar_obc(:,:,:,Lnew),                    &
     &                      tl_vbar_obc(:,:,:,Lwrk),                    &
     &                      ad_zeta_obc(:,:,:,Lnew),                    &
     &                      tl_zeta_obc(:,:,:,Lwrk),                    &
#  endif
#  ifdef ADJUST_WSTRESS
     &                      ad_ustr(:,:,:,Lnew), tl_ustr(:,:,:,Lwrk),   &
     &                      ad_vstr(:,:,:,Lnew), tl_vstr(:,:,:,Lwrk),   &
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
     &                      ad_tflux(:,:,:,Lnew,:),                     &
     &                      tl_tflux(:,:,:,Lwrk,:),                     &
#   endif
     &                      ad_t(:,:,:,Lnew,:), tl_t(:,:,:,Lwrk,:),     &
     &                      ad_u(:,:,:,Lnew), tl_u(:,:,:,Lwrk),         &
     &                      ad_v(:,:,:,Lnew), tl_v(:,:,:,Lwrk),         &
#  else
     &                      ad_ubar(:,:,Lnew), tl_ubar(:,:,Lwrk),       &
     &                      ad_vbar(:,:,Lnew), tl_vbar(:,:,Lwrk),       &
#  endif
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
# endif

      RETURN
      END SUBROUTINE lanczos
!
!***********************************************************************
      SUBROUTINE new_gradient (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         Lold, Lnew, Lwrk,                        &
     &                         innLoop, outLoop,                        &
# ifdef MASKING
     &                         rmask, umask, vmask,                     &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                         tl_t_obc, tl_u_obc, tl_v_obc,            &
#  endif
     &                         tl_ubar_obc, tl_vbar_obc,                &
     &                         tl_zeta_obc,                             &
# endif
# ifdef ADJUST_WSTRESS
     &                         tl_ustr, tl_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                         tl_tflux,                                &
#  endif
     &                         tl_t, tl_u, tl_v,                        &
# else
     &                         tl_ubar, tl_vbar,                        &
# endif
     &                         tl_zeta,                                 &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                         ad_t_obc, ad_u_obc, ad_v_obc,            &
#  endif
     &                         ad_ubar_obc, ad_vbar_obc,                &
     &                         ad_zeta_obc,                             &
# endif
# ifdef ADJUST_WSTRESS
     &                         ad_ustr, ad_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                         ad_tflux,                                &
#  endif
     &                         ad_t, ad_u, ad_v,                        &
# else
     &                         ad_ubar, ad_vbar,                        &
# endif
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
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lold, Lnew, Lwrk
      integer, intent(in) :: innLoop, outLoop
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: ad_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: ad_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:,:,:,:)
#  endif
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: tl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: tl_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: ad_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: tl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: rec

      real(r8) :: fac1, fac2

      real(r8), dimension(0:NstateVar(ng)) :: dot
      real(r8), dimension(0:Ninner) :: DotProd, dot_new, dot_old

      character (len=256) :: ncname

# include "set_bounds.h"
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
     &                     LBi, UBi, LBj, UBj, LBij, UBij,              &
     &                     Lold, Lnew, Lold, fac1, fac2,                &
# ifdef MASKING
     &                     rmask, umask, vmask,                         &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                     ad_t_obc, ad_t_obc,                          &
     &                     ad_u_obc, ad_u_obc,                          &
     &                     ad_v_obc, ad_v_obc,                          &
#  endif
     &                     ad_ubar_obc, ad_ubar_obc,                    &
     &                     ad_vbar_obc, ad_vbar_obc,                    &
     &                     ad_zeta_obc, ad_zeta_obc,                    &
# endif
# ifdef ADJUST_WSTRESS
     &                     ad_ustr, ad_ustr,                            &
     &                     ad_vstr, ad_vstr,                            &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     ad_tflux, ad_tflux,                          &
#  endif
     &                     ad_t, ad_t,                                  &
     &                     ad_u, ad_u,                                  &
     &                     ad_v, ad_v,                                  &
# else
     &                     ad_ubar, ad_ubar,                            &
     &                     ad_vbar, ad_vbar,                            &
# endif
     &                     ad_zeta, ad_zeta)
!
!  Adjust gradient against all previous gradients.
!
      IF (ndefADJ(ng).gt.0) THEN
        WRITE (ncname,10) TRIM(ADM(ng)%base), outLoop
 10     FORMAT (a,'_',i3.3,'.nc')
      ELSE
        ncname=ADM(ng)%name
      END IF
!
      DO rec=1,innLoop
!
!  Read in each previous gradient state solutions, G(0) to G(k), and
!  compute its associated dot angaint curret G(k+1). Each gradient
!  solution is loaded into TANGENT LINEAR STATE ARRAYS at index Lwrk.
!
        CALL read_state (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Lwrk, rec,                                     &
     &                   ndefADJ(ng), ADM(ng)%ncid, ncname,             &
# ifdef MASKING
     &                   rmask, umask, vmask,                           &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   tl_t_obc, tl_u_obc, tl_v_obc,                  &
#  endif
     &                   tl_ubar_obc, tl_vbar_obc,                      &
     &                   tl_zeta_obc,                                   &
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
        IF (exit_flag.ne.NoError) RETURN
!
!  In this expression for FAC2, the term cg_QG gives the contribution
!  to the gradient of Jo, and the term cg_Tmatrix gives the contribution
!  of Jb:
!
!    ad_var(Lold) = fac1 * ad_var(Lold) + fac2 * tl_var(Lwrk)
!
!  AMM: I do not think the second term in fac2 is needed now since we
!       are always working in terms of the total gradient of J=Jb+Jo.
!       CHECK:
!
!       fac2=-(cg_Tmatrix(rec,3)+cg_QG(rec,outLoop))
!
        fac1=1.0_r8
        fac2=-cg_QG(rec,outLoop)

        CALL state_addition (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       Lold, Lwrk, Lold, fac1, fac2,              &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       ad_t_obc, tl_t_obc,                        &
     &                       ad_u_obc, tl_u_obc,                        &
     &                       ad_v_obc, tl_v_obc,                        &
#  endif
     &                       ad_ubar_obc, tl_ubar_obc,                  &
     &                       ad_vbar_obc, tl_vbar_obc,                  &
     &                       ad_zeta_obc, tl_zeta_obc,                  &
# endif
# ifdef ADJUST_WSTRESS
     &                       ad_ustr, tl_ustr,                          &
     &                       ad_vstr, tl_vstr,                          &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                       ad_tflux, tl_tflux,                        &
#  endif
     &                       ad_t, tl_t,                                &
     &                       ad_u, tl_u,                                &
     &                       ad_v, tl_v,                                &
# else
     &                       ad_ubar, tl_ubar,                          &
     &                       ad_vbar, tl_vbar,                          &
# endif
     &                       ad_zeta, tl_zeta)
      END DO
!
!  Compute cost function gradient reduction.
!
      CALL state_dotprod (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj, LBij, UBij,               &
     &                    NstateVar(ng), dot(0:),                       &
# ifdef MASKING
     &                    rmask, umask, vmask,                          &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    ad_t_obc(:,:,:,:,Lold,:),                     &
     &                    ad_t_obc(:,:,:,:,Lold,:),                     &
     &                    ad_u_obc(:,:,:,:,Lold),                       &
     &                    ad_u_obc(:,:,:,:,Lold),                       &
     &                    ad_v_obc(:,:,:,:,Lold),                       &
     &                    ad_v_obc(:,:,:,:,Lold),                       &
#  endif
     &                    ad_ubar_obc(:,:,:,Lold),                      &
     &                    ad_ubar_obc(:,:,:,Lold),                      &
     &                    ad_vbar_obc(:,:,:,Lold),                      &
     &                    ad_vbar_obc(:,:,:,Lold),                      &
     &                    ad_zeta_obc(:,:,:,Lold),                      &
     &                    ad_zeta_obc(:,:,:,Lold),                      &
# endif
# ifdef ADJUST_WSTRESS
     &                    ad_ustr(:,:,:,Lold), ad_ustr(:,:,:,Lold),     &
     &                    ad_vstr(:,:,:,Lold), ad_vstr(:,:,:,Lold),     &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    ad_tflux(:,:,:,Lold,:),                       &
     &                    ad_tflux(:,:,:,Lold,:),                       &
#  endif
     &                    ad_t(:,:,:,Lold,:), ad_t(:,:,:,Lold,:),       &
     &                    ad_u(:,:,:,Lold), ad_u(:,:,:,Lold),           &
     &                    ad_v(:,:,:,Lold), ad_v(:,:,:,Lold),           &
# else
     &                    ad_ubar(:,:,Lold), ad_ubar(:,:,Lold),         &
     &                    ad_vbar(:,:,Lold), ad_vbar(:,:,Lold),         &
# endif
     &                    ad_zeta(:,:,Lold), ad_zeta(:,:,Lold))

      cg_Greduc(innLoop,outLoop)=SQRT(dot(0))/cg_Gnorm(outLoop)

      RETURN
      END SUBROUTINE new_gradient
!
!***********************************************************************
      SUBROUTINE hessian_evecs (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj, LBij, UBij,         &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Lold, Lnew, Lwrk,                       &
     &                          innLoop, outLoop,                       &
# ifdef MASKING
     &                          rmask, umask, vmask,                    &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                          nl_t_obc, nl_u_obc, nl_v_obc,           &
#  endif
     &                          nl_ubar_obc, nl_vbar_obc,               &
     &                          nl_zeta_obc,                            &
# endif
# ifdef ADJUST_WSTRESS
     &                          nl_ustr, nl_vstr,                       &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                          nl_tflux,                               &
#  endif
     &                          nl_t, nl_u, nl_v,                       &
# else
     &                          nl_ubar, nl_vbar,                       &
# endif
     &                          nl_zeta,                                &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                          tl_t_obc, tl_u_obc, tl_v_obc,           &
#  endif
     &                          tl_ubar_obc, tl_vbar_obc,               &
     &                          tl_zeta_obc,                            &
# endif
# ifdef ADJUST_WSTRESS
     &                          tl_ustr, tl_vstr,                       &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                          tl_tflux,                               &
#  endif
     &                          tl_t, tl_u, tl_v,                       &
# else
     &                          tl_ubar, tl_vbar,                       &
# endif
     &                          tl_zeta,                                &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                          ad_t_obc, ad_u_obc, ad_v_obc,           &
#  endif
     &                          ad_ubar_obc, ad_vbar_obc,               &
     &                          ad_zeta_obc,                            &
# endif
# ifdef ADJUST_WSTRESS
     &                          ad_ustr, ad_vstr,                       &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                          ad_tflux,                               &
#  endif
     &                          ad_t, ad_u, ad_v,                       &
# else
     &                          ad_ubar, ad_vbar,                       &
# endif
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
      USE state_copy_mod, ONLY : state_copy
      USE state_dotprod_mod, ONLY : state_dotprod
      USE state_initialize_mod, ONLY : state_initialize
      USE state_scale_mod, ONLY : state_scale
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lold, Lnew, Lwrk
      integer, intent(in) :: innLoop, outLoop
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: ad_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: ad_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: tl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: tl_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: nl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: nl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: nl_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: nl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: nl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: nl_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: nl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: nl_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: nl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: nl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: nl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: nl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: nl_zeta(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: ad_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: tl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: nl_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: nl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: nl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: nl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: nl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: nl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: nl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: nl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: nl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: nl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: nl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: nl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: nl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: nl_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: nl_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, ingood, rec, nvec, status, varid
      integer :: L1
      integer :: start(4), total(4)

      real(r8) :: fac, fac1, fac2

      real(r8), dimension(Ninner) :: RitzErr

      real(r8), dimension(0:NstateVar(ng)) :: dot
      real(r8), dimension(0:Ninner) :: DotProd, dot_new, dot_old

      character (len=256) :: ncname

# include "set_bounds.h"
!
      SourceFile='cgradient_lanczos.h, hessian_evecs'
!
!-----------------------------------------------------------------------
!  Calculate converged eigenvectors of the Hessian.
!-----------------------------------------------------------------------
!
!  Count and collect the converged eigenvalues.
!
      ingood=0
      DO i=innLoop,1,-1
        ingood=ingood+1
        Ritz(ingood)=cg_Ritz(i,outLoop)
        RitzErr(ingood)=cg_RitzErr(i,outLoop)
      END DO
      nConvRitz=ingood
!
!  Write out number of converged eigenvalues.
!
      CALL netcdf_put_ivar (ng, model, HSS(ng)%name, 'nConvRitz',       &
     &                      nConvRitz, (/0/), (/0/),                    &
     &                      ncid = HSS(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  First, premultiply the converged eigenvectors of the tridiagonal
!  matrix T(k) by the matrix of Lanczos vectors Q(k).  Use tangent
!  linear (index Lwrk) and adjoint (index Lold) state arrays as
!  temporary storage.
!-----------------------------------------------------------------------
!
      IF (Master) WRITE (stdout,10)
!
      COLUMNS : DO nvec=innLoop,1,-1
!
!  Initialize adjoint state arrays: ad_var(Lold) = fac
!
        fac=0.0_r8

        CALL state_initialize (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         Lold, fac,                               &
# ifdef MASKING
     &                         rmask, umask, vmask,                     &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                         ad_t_obc, ad_u_obc, ad_v_obc,            &
#  endif
     &                         ad_ubar_obc, ad_vbar_obc,                &
     &                         ad_zeta_obc,                             &
# endif
# ifdef ADJUST_WSTRESS
     &                         ad_ustr, ad_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                         ad_tflux,                                &
#  endif
     &                         ad_t, ad_u, ad_v,                        &
# else
     &                         ad_ubar, ad_vbar,                        &
# endif
     &                         ad_zeta)
!
!  Compute Hessian eigenvectors.
!
        IF (ndefADJ(ng).gt.0) THEN
          WRITE (ncname,20) TRIM(ADM(ng)%base), outLoop
        ELSE
          ncname=ADM(ng)%name
        END IF
!
        ROWS : DO rec=1,innLoop
!
!  Read gradient solution and load it into TANGENT LINEAR STATE ARRAYS
!  at index Lwrk.
!
          CALL read_state (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj, LBij, UBij,              &
     &                     Lwrk, rec,                                   &
     &                     ndefADJ(ng), ADM(ng)%ncid, ncname,           &
# ifdef MASKING
     &                     rmask, umask, vmask,                         &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                     tl_t_obc, tl_u_obc, tl_v_obc,                &
#  endif
     &                     tl_ubar_obc, tl_vbar_obc,                    &
     &                     tl_zeta_obc,                                 &
# endif
# ifdef ADJUST_WSTRESS
     &                     tl_ustr, tl_vstr,                            &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     tl_tflux,                                    &
#  endif
     &                     tl_t, tl_u, tl_v,                            &
# else
     &                     tl_ubar, tl_vbar,                            &
# endif
     &                     tl_zeta)
          IF (exit_flag.ne.NoError) RETURN
!
!  Compute Hessian eigenvectors:
!
!    ad_var(Lold) = fac1 * ad_var(Lold) + fac2 * tl_var(Lwrk)
!
          fac1=1.0_r8
          fac2=cg_zv(rec,nvec)

          CALL state_addition (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         Lold, Lwrk, Lold, fac1, fac2,            &
# ifdef MASKING
     &                         rmask, umask, vmask,                     &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                         ad_t_obc, tl_t_obc,                      &
     &                         ad_u_obc, tl_u_obc,                      &
     &                         ad_v_obc, tl_v_obc,                      &
#  endif
     &                         ad_ubar_obc, tl_ubar_obc,                &
     &                         ad_vbar_obc, tl_vbar_obc,                &
     &                         ad_zeta_obc, tl_zeta_obc,                &
# endif
# ifdef ADJUST_WSTRESS
     &                         ad_ustr, tl_ustr,                        &
     &                         ad_vstr, tl_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                         ad_tflux, tl_tflux,                      &
#  endif
     &                         ad_t, tl_t,                              &
     &                         ad_u, tl_u,                              &
     &                         ad_v, tl_v,                              &
# else
     &                         ad_ubar, tl_ubar,                        &
     &                         ad_vbar, tl_vbar,                        &
# endif
     &                         ad_zeta, tl_zeta)
        END DO ROWS
!
!  Write eigenvectors into Hessian NetCDF.
!
        LwrtState2d(ng)=.TRUE.
        CALL wrt_hessian (ng, Lold, Lold)
        LwrtState2d(ng)=.FALSE.
        IF (exit_flag.ne.NoERRor) RETURN

      END DO COLUMNS
!
!-----------------------------------------------------------------------
!  Second, orthonormalize the converged Hessian vectors against each
!  other. Use tangent linear state arrays (index Lwrk) as temporary
!  storage.
!-----------------------------------------------------------------------
!
!  Use nl_var(1) as temporary storage since we need to preserve
!  ad_var(Lnew).
!
      IF (ndefADJ(ng).gt.0) THEN
        WRITE (ncname,20) TRIM(HSS(ng)%base), outLoop
      ELSE
        ncname=HSS(ng)%name
      END IF
      IF (Master) WRITE (stdout,30)
!
      DO nvec=1,innLoop
!
!  Read in just computed Hessian eigenvectors into adjoint state array
!  index Lold.
!
        CALL read_state (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Lold, nvec,                                    &
     &                   0, HSS(ng)%ncid, ncname,                       &
# ifdef MASKING
     &                   rmask, umask, vmask,                           &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   ad_t_obc, ad_u_obc, ad_v_obc,                  &
#  endif
     &                   ad_ubar_obc, ad_vbar_obc,                      &
     &                   ad_zeta_obc,                                   &
# endif
# ifdef ADJUST_WSTRESS
     &                   ad_ustr, ad_vstr,                              &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                   ad_tflux,                                      &
#  endif
     &                   ad_t, ad_u, ad_v,                              &
# else
     &                   ad_ubar, ad_vbar,                              &
# endif
     &                   ad_zeta)
        IF (exit_flag.ne.NoError) RETURN
!
!  Initialize nonlinear state arrays index L1 with just read Hessian
!  vector in index Lold (initialize the summation):
!
!  Copy ad_var(Lold) into nl_var(L1).
!
        L1=1

        CALL state_copy (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Lold, L1,                                      &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   nl_t_obc, ad_t_obc,                            &
     &                   nl_u_obc, ad_u_obc,                            &
     &                   nl_v_obc, ad_v_obc,                            &
#  endif
     &                   nl_ubar_obc, ad_ubar_obc,                      &
     &                   nl_vbar_obc, ad_vbar_obc,                      &
     &                   nl_zeta_obc, ad_zeta_obc,                      &
# endif
# ifdef ADJUST_WSTRESS
     &                   nl_ustr, ad_ustr,                              &
     &                   nl_vstr, ad_vstr,                              &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                   nl_tflux, ad_tflux,                            &
#  endif
     &                   nl_t, ad_t,                                    &
     &                   nl_u, ad_u,                                    &
     &                   nl_v, ad_v,                                    &
# else
     &                   nl_ubar, ad_ubar,                              &
     &                   nl_vbar, ad_vbar,                              &
# endif
     &                   nl_zeta, ad_zeta)
!
!  Orthogonalize Hessian eigenvectors against each other.
!
        DO rec=1,nvec-1
!
!  Read in gradient just computed Hessian eigenvectors into tangent
!  linear state array index Lwrk.
!
          CALL read_state (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj, LBij, UBij,              &
     &                     Lwrk, rec,                                   &
     &                     0, HSS(ng)%ncid, ncname,                     &
# ifdef MASKING
     &                     rmask, umask, vmask,                         &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                     tl_t_obc, tl_u_obc, tl_v_obc,                &
#  endif
     &                     tl_ubar_obc, tl_vbar_obc,                    &
     &                     tl_zeta_obc,                                 &
# endif
# ifdef ADJUST_WSTRESS
     &                     tl_ustr, tl_vstr,                            &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     tl_tflux,                                    &
#  endif
     &                     tl_t, tl_u, tl_v,                            &
# else
     &                     tl_ubar, tl_vbar,                            &
# endif
     &                     tl_zeta)
          IF (exit_flag.ne.NoError) RETURN
!
!  Compute dot product.
!
          CALL state_dotprod (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj, LBij, UBij,           &
     &                        NstateVar(ng), dot(0:),                   &
# ifdef MASKING
     &                        rmask, umask, vmask,                      &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                        ad_t_obc(:,:,:,:,Lold,:),                 &
     &                        tl_t_obc(:,:,:,:,Lwrk,:),                 &
     &                        ad_u_obc(:,:,:,:,Lold),                   &
     &                        tl_u_obc(:,:,:,:,Lwrk),                   &
     &                        ad_v_obc(:,:,:,:,Lold),                   &
     &                        tl_v_obc(:,:,:,:,Lwrk),                   &
#  endif
     &                        ad_ubar_obc(:,:,:,Lold),                  &
     &                        tl_ubar_obc(:,:,:,Lwrk),                  &
     &                        ad_vbar_obc(:,:,:,Lold),                  &
     &                        tl_vbar_obc(:,:,:,Lwrk),                  &
     &                        ad_zeta_obc(:,:,:,Lold),                  &
     &                        tl_zeta_obc(:,:,:,Lwrk),                  &
# endif
# ifdef ADJUST_WSTRESS
     &                        ad_ustr(:,:,:,Lold), tl_ustr(:,:,:,Lwrk), &
     &                        ad_vstr(:,:,:,Lold), tl_vstr(:,:,:,Lwrk), &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                        ad_tflux(:,:,:,Lold,:),                   &
     &                        tl_tflux(:,:,:,Lwrk,:),                   &
#  endif
     &                        ad_t(:,:,:,Lold,:), tl_t(:,:,:,Lwrk,:),   &
     &                        ad_u(:,:,:,Lold), tl_u(:,:,:,Lwrk),       &
     &                        ad_v(:,:,:,Lold), tl_v(:,:,:,Lwrk),       &
# else
     &                        ad_ubar(:,:,Lold), tl_ubar(:,:,Lwrk),     &
     &                        ad_vbar(:,:,Lold), tl_vbar(:,:,Lwrk),     &
# endif
     &                        ad_zeta(:,:,Lold), tl_zeta(:,:,Lwrk))
!
!  Orthogonalize Hessian eigenvectors:
!
!    nl_var(L1) = fac1 * nl_var(L1) + fac2 * tl_var(Lwrk)
!
          fac1=1.0_r8
          fac2=-dot(0)

          CALL state_addition (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         L1, Lwrk, L1, fac1, fac2,                &
# ifdef MASKING
     &                         rmask, umask, vmask,                     &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                         nl_t_obc, tl_t_obc,                      &
     &                         nl_u_obc, tl_u_obc,                      &
     &                         nl_v_obc, tl_v_obc,                      &
#  endif
     &                         nl_ubar_obc, tl_ubar_obc,                &
     &                         nl_vbar_obc, tl_vbar_obc,                &
     &                         nl_zeta_obc, tl_zeta_obc,                &
# endif
# ifdef ADJUST_WSTRESS
     &                         nl_ustr, tl_ustr,                        &
     &                         nl_vstr, tl_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                         nl_tflux, tl_tflux,                      &
#  endif
     &                         nl_t, tl_t,                              &
     &                         nl_u, tl_u,                              &
     &                         nl_v, tl_v,                              &
# else
     &                         nl_ubar, tl_ubar,                        &
     &                         nl_vbar, tl_vbar,                        &
# endif
     &                         nl_zeta, tl_zeta)
        END DO
!
!  Compute normalization factor.
!
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj, LBij, UBij,             &
     &                      NstateVar(ng), dot(0:),                     &
# ifdef MASKING
     &                      rmask, umask, vmask,                        &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                      nl_t_obc(:,:,:,:,L1,:),                     &
     &                      nl_t_obc(:,:,:,:,L1,:),                     &
     &                      nl_u_obc(:,:,:,:,L1),                       &
     &                      nl_u_obc(:,:,:,:,L1),                       &
     &                      nl_v_obc(:,:,:,:,L1),                       &
     &                      nl_v_obc(:,:,:,:,L1),                       &
#  endif
     &                      nl_ubar_obc(:,:,:,L1),                      &
     &                      nl_ubar_obc(:,:,:,L1),                      &
     &                      nl_vbar_obc(:,:,:,L1),                      &
     &                      nl_vbar_obc(:,:,:,L1),                      &
     &                      nl_zeta_obc(:,:,:,L1),                      &
     &                      nl_zeta_obc(:,:,:,L1),                      &
# endif
# ifdef ADJUST_WSTRESS
     &                      nl_ustr(:,:,:,L1), nl_ustr(:,:,:,L1),       &
     &                      nl_vstr(:,:,:,L1), nl_vstr(:,:,:,L1),       &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                      nl_tflux(:,:,:,L1,:),                       &
     &                      nl_tflux(:,:,:,L1,:),                       &
#  endif
     &                      nl_t(:,:,:,L1,:), nl_t(:,:,:,L1,:),         &
     &                      nl_u(:,:,:,L1), nl_u(:,:,:,L1),             &
     &                      nl_v(:,:,:,L1), nl_v(:,:,:,L1),             &
# else
     &                      nl_ubar(:,:,L1), nl_ubar(:,:,L1),           &
     &                      nl_vbar(:,:,L1), nl_vbar(:,:,L1),           &
# endif
     &                      nl_zeta(:,:,L1), nl_zeta(:,:,L1))
!
!  Normalize Hessian eigenvectors:
!
!    nl_var(L1) = fac * nl_var(L1)
!
        fac=1.0_r8/SQRT(dot(0))

        CALL state_scale (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, LBij, UBij,               &
     &                    L1, L1, fac,                                  &
# ifdef MASKING
     &                    rmask, umask, vmask,                          &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    nl_t_obc, nl_u_obc, nl_v_obc,                 &
#  endif
     &                    nl_ubar_obc, nl_vbar_obc,                     &
     &                    nl_zeta_obc,                                  &
# endif
# ifdef ADJUST_WSTRESS
     &                    nl_ustr, nl_vstr,                             &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    nl_tflux,                                     &
#  endif
     &                    nl_t, nl_u, nl_v,                             &
# else
     &                    nl_ubar, nl_vbar,                             &
# endif
     &                    nl_zeta)
!
!  Copy nl_var(L1) into  ad_var(Lold).
!
        CALL state_copy (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   L1, Lold,                                      &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   ad_t_obc, nl_t_obc,                            &
     &                   ad_u_obc, nl_u_obc,                            &
     &                   ad_v_obc, nl_v_obc,                            &
#  endif
     &                   ad_ubar_obc, nl_ubar_obc,                      &
     &                   ad_vbar_obc, nl_vbar_obc,                      &
     &                   ad_zeta_obc, nl_zeta_obc,                      &
# endif
# ifdef ADJUST_WSTRESS
     &                   ad_ustr, nl_ustr,                              &
     &                   ad_vstr, nl_vstr,                              &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                   ad_tflux, nl_tflux,                            &
#  endif
     &                   ad_t, nl_t,                                    &
     &                   ad_u, nl_u,                                    &
     &                   ad_v, nl_v,                                    &
# else
     &                   ad_ubar, nl_ubar,                              &
     &                   ad_vbar, nl_vbar,                              &
# endif
     &                   ad_zeta, nl_zeta)
!
!  Write out converged Ritz eigenvalues and is associated accuracy.
!
        CALL netcdf_put_fvar (ng, model, HSS(ng)%name, 'Ritz',          &
     &                        Ritz(nvec:), (/nvec/), (/1/),             &
     &                        ncid = HSS(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN

        CALL netcdf_put_fvar (ng, model, HSS(ng)%name, 'Ritz_error',    &
     &                        RitzErr(nvec:), (/nvec/), (/1/),          &
     &                        ncid = HSS(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
!  Replace record "nvec" of Hessian eigenvectors NetCDF with the
!  normalized value in adjoint state arrays at index Lold.
!
        HSS(ng)%Rindex=nvec-1
        LwrtState2d(ng)=.TRUE.
        CALL wrt_hessian (ng, Lold, Lold)
        LwrtState2d(ng)=.FALSE.
        IF (exit_flag.ne.NoERRor) RETURN

      END DO

  10  FORMAT (/,' Computing converged Hessian eigenvectors...',/)
  20  FORMAT (a,'_',i3.3,'.nc')
  30  FORMAT (/,' Orthonormalizing converged Hessian eigenvectors...',/)

      RETURN
      END SUBROUTINE hessian_evecs

      SUBROUTINE cg_write (ng, model, innLoop, outLoop)
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
!
!  Local variable declarations.
!
      integer :: status
!
      SourceFile='cgradient_lanczos.h, cg_write'
!
!-----------------------------------------------------------------------
!  Write out conjugate gradient vectors.
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
     &                        nConvRitz, (/0/), (/0/),                  &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out converged Ritz eigenvalues.
!
      IF (innLoop.eq.Ninner) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'Ritz',          &
     &                        Ritz, (/1/), (/nConvRitz/),               &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out conjugate gradient norms.
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
      IF (innLoop.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_gamma',      &
     &                        cg_gamma, (/1,1/), (/Ninner,Nouter/),     &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Initial gradient normalization factor.
!
      IF (innLoop.eq.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_Gnorm',      &
     &                        cg_Gnorm, (/1/), (/Nouter/),              &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Lanczos vector normalization factor.
!
      IF (innLoop.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_QG',         &
     &                        cg_QG, (/1,1/), (/Ninner,Nouter/),        &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Reduction in the gradient norm.
!
      IF (innLoop.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_Greduc',     &
     &                        cg_Greduc, (/1,1/), (/Ninner,Nouter/),    &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Lanczos recurrence tridiagonal matrix.
!
      IF (innLoop.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_Tmatrix',    &
     &                        cg_Tmatrix, (/1,1/), (/Ninner,3/),        &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Lanczos tridiagonal matrix, upper diagonal elements.
!
      IF (innLoop.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_zu',         &
     &                        cg_zu, (/1,1/), (/Ninner,Nouter/),        &
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
     &                        cg_zv, (/1,1/), (/Ninner,Ninner/),        &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out Lanczos algorithm coefficients into Lanczos vectors
!  output file (for now adjoint history file). These coefficients
!  can be used for preconditioning or to compute the sensitivity
!  of the observations to the 4DVAR data assimilation system.
!
      IF (innLoop.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, ADM(ng)%name, 'cg_beta',       &
     &                        cg_beta, (/1,1/), (/Ninner+1,Nouter/),    &
     &                        ncid = ADM(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
        CALL netcdf_put_fvar (ng, model, ADM(ng)%name, 'cg_delta',      &
     &                        cg_delta, (/1,1/), (/Ninner,Nouter/),     &
     &                        ncid = ADM(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
        CALL netcdf_put_fvar (ng, model, ADM(ng)%name, 'cg_zv',         &
     &                        cg_zv, (/1,1/), (/Ninner,Ninner/),        &
     &                        ncid = ADM(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
# ifdef EVOLVED_LCZ
        CALL netcdf_put_fvar (ng, model, HSS(ng)%name, 'cg_beta',       &
     &                        cg_beta, (/1,1/), (/Ninner+1,Nouter/),    &
     &                        ncid = HSS(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
        CALL netcdf_put_fvar (ng, model, HSS(ng)%name, 'cg_delta',      &
     &                        cg_delta, (/1,1/), (/Ninner,Nouter/),     &
     &                        ncid = HSS(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
        CALL netcdf_put_fvar (ng, model, HSS(ng)%name, 'cg_zv',         &
     &                        cg_zv, (/1,1/), (/Ninner,Ninner/),        &
     &                        ncid = HSS(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize model/observation NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, model, DAV(ng)%name, DAV(ng)%ncid)

      RETURN
      END SUBROUTINE cg_write
!
!=======================================================================
      SUBROUTINE new_cost (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj, LBij, UBij,              &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     innLoop, outLoop,                            &
# ifdef MASKING
     &                     rmask, umask, vmask,                         &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                     nl_t_obc, nl_u_obc, nl_v_obc,                &
#  endif
     &                     nl_ubar_obc, nl_vbar_obc,                    &
     &                     nl_zeta_obc,                                 &
# endif
# ifdef ADJUST_WSTRESS
     &                     nl_ustr, nl_vstr,                            &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     nl_tflux,                                    &
#  endif
     &                     nl_t, nl_u, nl_v,                            &
# else
     &                     nl_ubar, nl_vbar,                            &
# endif
     &                     nl_zeta)
!=======================================================================
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
      USE state_initialize_mod, ONLY : state_initialize
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: innLoop, outLoop
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: nl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: nl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: nl_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: nl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: nl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: nl_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: nl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: nl_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: nl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: nl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: nl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: nl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: nl_zeta(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: nl_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: nl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: nl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: nl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: nl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: nl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: nl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: nl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: nl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: nl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: nl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: nl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: nl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: nl_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: nl_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, rec, Lscale
      integer :: L1 = 1
      integer :: L2 = 2

      real(r8) :: fac, fac1, fac2

      real(r8), dimension(0:NstateVar(ng)) :: dot

      logical :: Ltrans

      character (len=256) :: ncname

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute the cost function based on the formula of Tshimanga
!  (PhD thesis, p 154, eqn A.15):
!
!    J = J_initial + 0.5 * transpose(r) Q z
!
!  where J_initial is the value of the cost function when inner=0
!  (i.e. Cost0), r is the initial cost function gradient when inner=0,
!  Q is the matrix of Lanczos vectors, and z is the solution of
!  Tz=-transpose(Q) r, T being the associated tridiagonal matrix of
!  the Lanczos recursion. Note that even when r and x are in y-space,
!  as in the preconditioned case, their dot-product is equal to that
!  of the same variables transformed to v-space.
!-----------------------------------------------------------------------
!
!  Compute the current increment and save in nl_var(L1).
!
!  Clear the adjoint working arrays (index Linp) since the tangent
!  linear model initial condition on the first inner-loop is zero:
!
!    nl_var(L1) = fac
!
      fac=0.0_r8

      CALL state_initialize (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       L1, fac,                                   &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       nl_t_obc, nl_u_obc, nl_v_obc,              &
#  endif
     &                       nl_ubar_obc, nl_vbar_obc,                  &
     &                       nl_zeta_obc,                               &
# endif
# ifdef ADJUST_WSTRESS
     &                       nl_ustr, nl_vstr,                          &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                       nl_tflux,                                  &
#  endif
     &                       nl_t, nl_u, nl_v,                          &
# else
     &                       nl_ubar, nl_vbar,                          &
# endif
     &                       nl_zeta)
!
!  Now compute nl_var(L1) = Q * cg_zu.
!
      IF (ndefADJ(ng).gt.0) THEN
        WRITE (ncname,10) TRIM(ADM(ng)%base), outLoop
 10     FORMAT (a,'_',i3.3,'.nc')
      ELSE
        ncname=ADM(ng)%name
      END IF
!
      DO rec=1,innLoop
!
!  Read gradient solution and load it into nl_var(L2).
!
        CALL read_state (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   L2, rec,                                       &
     &                   ndefADJ(ng), ADM(ng)%ncid, ncname,             &
# ifdef MASKING
     &                   rmask, umask, vmask,                           &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   nl_t_obc, nl_u_obc, nl_v_obc,                  &
#  endif
     &                   nl_ubar_obc, nl_vbar_obc,                      &
     &                   nl_zeta_obc,                                   &
# endif
# ifdef ADJUST_WSTRESS
     &                   nl_ustr, nl_vstr,                              &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                   nl_tflux,                                      &
#  endif
     &                   nl_t, nl_u, nl_v,                              &
# else
     &                   nl_ubar, nl_vbar,                              &
# endif
     &                   nl_zeta)
        IF (exit_flag.ne.NoError) RETURN
!
!  Perform the summation:
!
!    nl_var(L1) = fac1 * nl_var(L1) + fac2 * nl_var(L2)
!
        fac1=1.0_r8
        fac2=cg_zu(rec,outLoop)

        CALL state_addition (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       L1, L2, L1, fac1, fac2,                    &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       nl_t_obc, nl_t_obc,                        &
     &                       nl_u_obc, nl_u_obc,                        &
     &                       nl_v_obc, nl_v_obc,                        &
#  endif
     &                       nl_ubar_obc, nl_ubar_obc,                  &
     &                       nl_vbar_obc, nl_vbar_obc,                  &
     &                       nl_zeta_obc, nl_zeta_obc,                  &
# endif
# ifdef ADJUST_WSTRESS
     &                       nl_ustr, nl_ustr,                          &
     &                       nl_vstr, nl_vstr,                          &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                       nl_tflux, nl_tflux,                        &
#  endif
     &                       nl_t, nl_t,                                &
     &                       nl_u, nl_u,                                &
     &                       nl_v, nl_v,                                &
# else
     &                       nl_ubar, nl_ubar,                          &
     &                       nl_vbar, nl_vbar,                          &
# endif
     &                       nl_zeta, nl_zeta)
      END DO
!
!  Now read the initial Lanczos vector again into nl_var(L2).
!
      rec=1
      IF (ndefADJ(ng).gt.0) THEN
        WRITE (ncname,10) TRIM(ADM(ng)%base), outLoop
      ELSE
        ncname=ADM(ng)%name
      END IF
!
      CALL read_state (ng, tile, model,                                 &
     &                 LBi, UBi, LBj, UBj, LBij, UBij,                  &
     &                 L2, rec,                                         &
     &                 ndefADJ(ng), ADM(ng)%ncid, ncname,               &
# ifdef MASKING
     &                 rmask, umask, vmask,                             &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                 nl_t_obc, nl_u_obc, nl_v_obc,                    &
#  endif
     &                 nl_ubar_obc, nl_vbar_obc,                        &
     &                 nl_zeta_obc,                                     &
# endif
# ifdef ADJUST_WSTRESS
     &                 nl_ustr, nl_vstr,                                &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                 nl_tflux,                                        &
#  endif
     &                 nl_t, nl_u, nl_v,                                &
# else
     &                 nl_ubar, nl_vbar,                                &
# endif
     &                 nl_zeta)
      IF (exit_flag.ne.NoError) RETURN
!
!  Compute the dot-product of the initial Lanczos vector with the
!  current increment.
!
      CALL state_dotprod (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj, LBij, UBij,               &
     &                    NstateVar(ng), dot(0:),                       &
# ifdef MASKING
     &                    rmask, umask, vmask,                          &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    nl_t_obc(:,:,:,:,L1,:),                       &
     &                    nl_t_obc(:,:,:,:,L2,:),                       &
     &                    nl_u_obc(:,:,:,:,L1),                         &
     &                    nl_u_obc(:,:,:,:,L2),                         &
     &                    nl_v_obc(:,:,:,:,L1),                         &
     &                    nl_v_obc(:,:,:,:,L2),                         &
#  endif
     &                    nl_ubar_obc(:,:,:,L1),                        &
     &                    nl_ubar_obc(:,:,:,L2),                        &
     &                    nl_vbar_obc(:,:,:,L1),                        &
     &                    nl_vbar_obc(:,:,:,L2),                        &
     &                    nl_zeta_obc(:,:,:,L1),                        &
     &                    nl_zeta_obc(:,:,:,L2),                        &
# endif
# ifdef ADJUST_WSTRESS
     &                    nl_ustr(:,:,:,L1), nl_ustr(:,:,:,L2),         &
     &                    nl_vstr(:,:,:,L1), nl_vstr(:,:,:,L2),         &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    nl_tflux(:,:,:,L1,:),                         &
     &                    nl_tflux(:,:,:,L2,:),                         &
#  endif
     &                    nl_t(:,:,:,L1,:), nl_t(:,:,:,L2,:),           &
     &                    nl_u(:,:,:,L1), nl_u(:,:,:,L2),               &
     &                    nl_v(:,:,:,L1), nl_v(:,:,:,L2),               &
# else
     &                    nl_ubar(:,:,L1), nl_ubar(:,:,L2),             &
     &                    nl_vbar(:,:,L1), nl_vbar(:,:,L2),             &
# endif
     &                    nl_zeta(:,:,L1), nl_zeta(:,:,L2))
!
!  Compute the new cost function. Only the total value is meaningful.
!  Note that we need to multiply dot(0) by cg_Gnorm(outLoop) because
!  the initial gradient is cg_Gnorm*q(0).
!
      FOURDVAR(ng)%CostFun(0)=FOURDVAR(ng)%Cost0(outLoop)+              &
     &                        0.5_r8*dot(0)*cg_Gnorm(outLoop)
      DO i=1,NstateVar(ng)
        FOURDVAR(ng)%CostFun(i)=0.0_r8
      END DO
!
!  Compute the background cost function.
!
!    If preconditioning, convert nl_var(L1) from y-space into v-space.
!    This must be called before reading the sum of previous v-space
!    gradient since nl_var(L2) is used a temporary storage in precond.
!
      IF (Lprecond.and.(outLoop.gt.1)) THEN

        Lscale=2                 ! SQRT spectral LMP
        Ltrans=.FALSE.
!
        CALL precond (ng, tile, model, 'new cost function',             &
     &                LBi, UBi, LBj, UBj, LBij, UBij,                   &
     &                IminS, ImaxS, JminS, JmaxS,                       &
     &                NstateVar(ng), Lscale, Ltrans,                    &
     &                innLoop, outLoop,                                 &
# ifdef MASKING
     &                rmask, umask, vmask,                              &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                nl_t_obc, nl_u_obc, nl_v_obc,                     &
#  endif
     &                nl_ubar_obc, nl_vbar_obc,                         &
     &                nl_zeta_obc,                                      &
# endif
# ifdef ADJUST_WSTRESS
     &                nl_ustr, nl_vstr,                                 &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                nl_tflux,                                         &
#  endif
     &                nl_t, nl_u, nl_v,                                 &
# else
     &                nl_ubar, nl_vbar,                                 &
# endif
     &                nl_zeta)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Read the sum of previous v-space gradients from record 4 of ITL
!  tile into nl_var(L2). Note that all fields in the ITL file
!  are in v-space so there is no need to apply the preconditioner
!  to nl_var(L2).
!
      CALL read_state (ng, tile, model,                                 &
     &                 LBi, UBi, LBj, UBj, LBij, UBij,                  &
     &                 L2, 4,                                           &
     &                 ndefTLM(ng), ITL(ng)%ncid, ITL(ng)%name,         &
# ifdef MASKING
     &                 rmask, umask, vmask,                             &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                 nl_t_obc, nl_u_obc, nl_v_obc,                    &
#  endif
     &                 nl_ubar_obc, nl_vbar_obc,                        &
     &                 nl_zeta_obc,                                     &
# endif
# ifdef ADJUST_WSTRESS
     &                 nl_ustr, nl_vstr,                                &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                 nl_tflux,                                        &
#  endif
     &                 nl_t, nl_u, nl_v,                                &
# else
     &                 nl_ubar, nl_vbar,                                &
# endif
     &                 nl_zeta)
      IF (exit_flag.ne.NoError) RETURN
!
!  Perform the summation:
!
!  nl_var(L1) = fac1 * nl_var(L1) + fac2 * nl_var(L2)
!
      fac1=1.0_r8
      fac2=1.0_r8

      CALL state_addition (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj, LBij, UBij,              &
     &                     L1, L2, L1, fac1, fac2,                      &
# ifdef MASKING
     &                     rmask, umask, vmask,                         &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                     nl_t_obc, nl_t_obc,                          &
     &                     nl_u_obc, nl_u_obc,                          &
     &                     nl_v_obc, nl_v_obc,                          &
#  endif
     &                     nl_ubar_obc, nl_ubar_obc,                    &
     &                     nl_vbar_obc, nl_vbar_obc,                    &
     &                     nl_zeta_obc, nl_zeta_obc,                    &
# endif
# ifdef ADJUST_WSTRESS
     &                     nl_ustr, nl_ustr,                            &
     &                     nl_vstr, nl_vstr,                            &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     nl_tflux, nl_tflux,                          &
#  endif
     &                     nl_t, nl_t,                                  &
     &                     nl_u, nl_u,                                  &
     &                     nl_v, nl_v,                                  &
# else
     &                     nl_ubar, nl_ubar,                            &
     &                     nl_vbar, nl_vbar,                            &
# endif
     &                     nl_zeta, nl_zeta)
!
      CALL state_dotprod (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj, LBij, UBij,               &
     &                    NstateVar(ng), dot(0:),                       &
# ifdef MASKING
     &                    rmask, umask, vmask,                          &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    nl_t_obc(:,:,:,:,L1,:),                       &
     &                    nl_t_obc(:,:,:,:,L1,:),                       &
     &                    nl_u_obc(:,:,:,:,L1),                         &
     &                    nl_u_obc(:,:,:,:,L1),                         &
     &                    nl_v_obc(:,:,:,:,L1),                         &
     &                    nl_v_obc(:,:,:,:,L1),                         &
#  endif
     &                    nl_ubar_obc(:,:,:,L1),                        &
     &                    nl_ubar_obc(:,:,:,L1),                        &
     &                    nl_vbar_obc(:,:,:,L1),                        &
     &                    nl_vbar_obc(:,:,:,L1),                        &
     &                    nl_zeta_obc(:,:,:,L1),                        &
     &                    nl_zeta_obc(:,:,:,L1),                        &
# endif
# ifdef ADJUST_WSTRESS
     &                    nl_ustr(:,:,:,L1), nl_ustr(:,:,:,L1),         &
     &                    nl_vstr(:,:,:,L1), nl_vstr(:,:,:,L1),         &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    nl_tflux(:,:,:,L1,:),                         &
     &                    nl_tflux(:,:,:,L1,:),                         &
#  endif
     &                    nl_t(:,:,:,L1,:), nl_t(:,:,:,L1,:),           &
     &                    nl_u(:,:,:,L1), nl_u(:,:,:,L1),               &
     &                    nl_v(:,:,:,L1), nl_v(:,:,:,L1),               &
# else
     &                    nl_ubar(:,:,L1), nl_ubar(:,:,L1),             &
     &                    nl_vbar(:,:,L1), nl_vbar(:,:,L1),             &
# endif
     &                    nl_zeta(:,:,L1), nl_zeta(:,:,L1))
!
      FOURDVAR(ng)%BackCost(0)=0.5_r8*dot(0)
      FOURDVAR(ng)%ObsCost(0)=FOURDVAR(ng)%CostFun(0)-                  &
     &                        FOURDVAR(ng)%BackCost(0)
      DO i=1,NstateVar(ng)
        FOURDVAR(ng)%BackCost(i)=0.0_r8
        FOURDVAR(ng)%ObsCost(i)=0.0_r8
      END DO

      RETURN
      END SUBROUTINE new_cost
!
!***********************************************************************
      SUBROUTINE precond (ng, tile, model, message,                     &
     &                    LBi, UBi, LBj, UBj, LBij, UBij,               &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    NstateVars, Lscale, Ltrans,                   &
     &                    innLoop, outLoop,                             &
# ifdef MASKING
     &                    rmask, umask, vmask,                          &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    nl_t_obc, nl_u_obc, nl_v_obc,                 &
#  endif
     &                    nl_ubar_obc, nl_vbar_obc,                     &
     &                    nl_zeta_obc,                                  &
# endif
# ifdef ADJUST_WSTRESS
     &                    nl_ustr, nl_vstr,                             &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    nl_tflux,                                     &
#  endif
     &                    nl_t, nl_u, nl_v,                             &
# else
     &                    nl_ubar, nl_vbar,                             &
# endif
     &                    nl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_ncparam
      USE mod_netcdf
      USE mod_iounits
      USE mod_scalars
!
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_reduce
# endif
      USE state_addition_mod, ONLY : state_addition
      USE state_copy_mod, ONLY : state_copy
      USE state_dotprod_mod, ONLY : state_dotprod
!
!  Imported variable declarations.
!
      logical, intent(in) :: Ltrans

      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: NstateVars, Lscale
      integer, intent(in) :: innLoop, outLoop

      character (len=*), intent(in) :: message
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: nl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: nl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: nl_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: nl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: nl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: nl_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: nl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: nl_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: nl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: nl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: nl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: nl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: nl_zeta(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: nl_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: nl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: nl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: nl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: nl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: nl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: nl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: nl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: nl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: nl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: nl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: nl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: nl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: nl_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: nl_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: NSUB, i, j, k, L1, L2, nvec, rec, ndefLCZ
      integer :: is, ie, inc, iss, ncid
      integer :: nol, nols, nole, ninc
      integer :: ingood
      integer :: namm
# ifdef SOLVE3D
      integer :: it
# endif
      integer, parameter :: ndef = 1

      real(r8) :: cff, fac, fac1, fac2, facritz
      real(r8), dimension(0:NstateVars) :: Dotprod
      real(r8), dimension(1:Ninner+1,Nouter) :: beta_lcz
      real(r8), dimension(Ninner,Ninner) :: zv_lcz

      character (len=256) :: ncname

# ifdef DISTRIBUTE
      character (len=3) :: op_handle
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  THIS PRECONDITIONER IS WRITTEN IN PRODUCT FORM AS DESCRIBED BY
!  TSHIMANGA - PhD thesis, page 75, proof of proposition 2.3.1.
!  IT IS THEREFORE IMPORTANT THAT THE EIGENVECTORS/RITZ VECTORS THAT
!  ARE COMPUTED BY is4dvar_lanczos.h ARE ORTHONORMALIZED.
!
!  Apply the preconditioner. The approximated Hessian matrix is computed
!  from the eigenvectors computed by the Lanczos algorithm which are
!  stored in HSS(ng)%name NetCDF file.
!-----------------------------------------------------------------------
!
      L1=1
      L2=2
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

      IF (Master) THEN
        IF (Lritz) THEN
          WRITE (stdout,10) outLoop, innLoop, 'Ritz', TRIM(message)
        ELSE
          WRITE (stdout,10) outLoop, innLoop, 'Spectral', TRIM(message)
        END IF
      END IF
!
!  Apply the preconditioners derived from all previous outer-loops
!  sequentially.
!
      DO nol=nols,nole,ninc
!
!  Read the primitive Ritz vectors cg_v and cg_beta.
!
        WRITE (ncname,20) TRIM(ADM(ng)%base), nol
        IF (Master) THEN
          WRITE (stdout,30) outLoop, innLoop, TRIM(ncname)
        END IF
!
        CALL netcdf_get_fvar (ng, model, ncname, 'cg_beta',             &
     &                        beta_lcz)
        IF (exit_flag.ne. NoError) RETURN

        CALL netcdf_get_fvar (ng, model, ncname, 'cg_zv',               &
     &                        zv_lcz)
        IF (exit_flag.ne. NoError) RETURN
!
!  Determine the number of Ritz vectors to use.
!  For Lritz=.TRUE., choose HvecErr to be larger.
!
        ingood=0
        DO i=1,Ninner
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
        IF (Master) THEN
          WRITE (stdout,40) outLoop, innLoop, ingood
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
        DO nvec=is,ie,inc
!
          fac2=0.0_r8
!
!  If using the Ritz preconditioner, read information from the
!  Lanczos vector file.
!
          IF (Lritz) THEN
!
            IF (.not.Ltrans)THEN
!
!  Determine adjoint file to process.
!
              WRITE (ncname,20) TRIM(ADM(ng)%base), nol
              IF (Master.and.(nvec.eq.is)) THEN
                WRITE (stdout,50) outLoop, innLoop, TRIM(ncname)
              END IF
!
!  Read in the Lanczos vector q_k+1 computed from the incremental
!  4DVar algorithm previous outers loop, where k=Ninner+1. Load
!  Lanczos vectors into NL state arrays at index L2.
!
              rec=Ninner+1
              CALL read_state (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         L2, rec,                                 &
     &                         ndef, ncid, ncname,                      &
# ifdef MASKING
     &                         rmask, umask, vmask,                     &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                         nl_t_obc, nl_u_obc, nl_v_obc,            &
#  endif
     &                         nl_ubar_obc, nl_vbar_obc,                &
     &                         nl_zeta_obc,                             &
# endif
# ifdef ADJUST_WSTRESS
     &                         nl_ustr, nl_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                         nl_tflux,                                &
#  endif
     &                         nl_t, nl_u, nl_v,                        &
# else
     &                         nl_ubar, nl_vbar,                        &
# endif
     &                         nl_zeta)
              IF (exit_flag.ne.NoError) RETURN
!
!  Compute the dot-product between the input vector and the Ninner+1
!  Lanczos vector.
!
              CALL state_dotprod (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj, LBij, UBij,       &
     &                            NstateVars, Dotprod(0:),              &
# ifdef MASKING
     &                            rmask, umask, vmask,                  &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                            nl_t_obc(:,:,:,:,L1,:),               &
     &                            nl_t_obc(:,:,:,:,L2,:),               &
     &                            nl_u_obc(:,:,:,:,L1),                 &
     &                            nl_u_obc(:,:,:,:,L2),                 &
     &                            nl_v_obc(:,:,:,:,L1),                 &
     &                            nl_v_obc(:,:,:,:,L2),                 &
#  endif
     &                            nl_ubar_obc(:,:,:,L1),                &
     &                            nl_ubar_obc(:,:,:,L2),                &
     &                            nl_vbar_obc(:,:,:,L1),                &
     &                            nl_vbar_obc(:,:,:,L2),                &
     &                            nl_zeta_obc(:,:,:,L1),                &
     &                            nl_zeta_obc(:,:,:,L2),                &
# endif
# ifdef ADJUST_WSTRESS
     &                            nl_ustr(:,:,:,L1), nl_ustr(:,:,:,L2), &
     &                            nl_vstr(:,:,:,L1), nl_vstr(:,:,:,L2), &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                            nl_tflux(:,:,:,L1,:),                 &
     &                            nl_tflux(:,:,:,L2,:),                 &
#  endif
     &                            nl_t(:,:,:,L1,:), nl_t(:,:,:,L2,:),   &
     &                            nl_u(:,:,:,L1), nl_u(:,:,:,L2),       &
     &                            nl_v(:,:,:,L1), nl_v(:,:,:,L2),       &
# else
     &                            nl_ubar(:,:,L1), nl_ubar(:,:,L2),     &
     &                            nl_vbar(:,:,L1), nl_vbar(:,:,L2),     &
# endif
     &                            nl_zeta(:,:,L1), nl_zeta(:,:,L2))

            END IF
!
!  Note: the primitive Ritz vectors zv_lcz are arranged in order of
!        ASCENDING eigenvalue while the Hessian eigenvectors are
!        arranged in DESCENDING order.
!
            facritz=beta_lcz(Ninner+1,nol)*zv_lcz(Ninner,Ninner+1-nvec)

            IF (.not.Ltrans) THEN
              facritz=facritz*Dotprod(0)
            END IF

          END IF
!
!  Read the converged Hessian eigenvectors into NLM state array,
!  index L2.
!
          WRITE (ncname,20) TRIM(HSS(ng)%base), nol
          IF (Master.and.(nvec.eq.is)) THEN
            WRITE (stdout,60) outLoop, innLoop, TRIM(ncname)
          END IF
!
          CALL read_state (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj, LBij, UBij,              &
     &                     L2, nvec,                                    &
     &                     ndef, ncid, ncname,                          &
# ifdef MASKING
     &                     rmask, umask, vmask,                         &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                     nl_t_obc, nl_u_obc, nl_v_obc,                &
#  endif
     &                     nl_ubar_obc, nl_vbar_obc,                    &
     &                     nl_zeta_obc,                                 &
# endif
# ifdef ADJUST_WSTRESS
     &                     nl_ustr, nl_vstr,                            &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     nl_tflux,                                    &
#  endif
     &                     nl_t, nl_u, nl_v,                            &
# else
     &                     nl_ubar, nl_vbar,                            &
# endif
     &                     nl_zeta)
          IF (exit_flag.ne.NoError) RETURN
!
!  Compute dot product between input vector and Hessian eigenvector.
!  The input vector is in nl_var(L1) and the Hessian vector in
!  nl_var(L2)
!
          CALL state_dotprod (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj, LBij, UBij,           &
     &                        NstateVars, Dotprod(0:),                  &
# ifdef MASKING
     &                        rmask, umask, vmask,                      &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                        nl_t_obc(:,:,:,:,L1,:),                   &
     &                        nl_t_obc(:,:,:,:,L2,:),                   &
     &                        nl_u_obc(:,:,:,:,L1),                     &
     &                        nl_u_obc(:,:,:,:,L2),                     &
     &                        nl_v_obc(:,:,:,:,L1),                     &
     &                        nl_v_obc(:,:,:,:,L2),                     &
#  endif
     &                        nl_ubar_obc(:,:,:,L1),                    &
     &                        nl_ubar_obc(:,:,:,L2),                    &
     &                        nl_vbar_obc(:,:,:,L1),                    &
     &                        nl_vbar_obc(:,:,:,L2),                    &
     &                        nl_zeta_obc(:,:,:,L1),                    &
     &                        nl_zeta_obc(:,:,:,L2),                    &
# endif
# ifdef ADJUST_WSTRESS
     &                        nl_ustr(:,:,:,L1), nl_ustr(:,:,:,L2),     &
     &                        nl_vstr(:,:,:,L1), nl_vstr(:,:,:,L2),     &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                        nl_tflux(:,:,:,L1,:),                     &
     &                        nl_tflux(:,:,:,L2,:),                     &
#  endif
     &                        nl_t(:,:,:,L1,:), nl_t(:,:,:,L2,:),       &
     &                        nl_u(:,:,:,L1), nl_u(:,:,:,L2),           &
     &                        nl_v(:,:,:,L1), nl_v(:,:,:,L2),           &
# else
     &                        nl_ubar(:,:,L1), nl_ubar(:,:,L2),         &
     &                        nl_vbar(:,:,L1), nl_vbar(:,:,L2),         &
# endif
     &                        nl_zeta(:,:,L1), nl_zeta(:,:,L2))
!
!  Lscale determines the form of the preconditioner:
!
!     1= spectral LMP
!    -1= Inverse spectral LMP
!     2= Square root spectral LMP
!    -2= Inverse square spectral root LMP
!
!    nl_var(L1) = fac1 * nl_var(L1) + fac2 * nl_var(L2)
!
!  Note: cg_Ritz contains the Ritz values written in ASCENDING order.
!
          fac1=1.0_r8

          IF (Lscale.eq.-1) THEN
            fac2=(cg_Ritz(Ninner+1-nvec,nol)-1.0_r8)*Dotprod(0)
          ELSE IF (Lscale.eq.1) THEN
            fac2=(1.0_r8/cg_Ritz(Ninner+1-nvec,nol)-1.0_r8)*Dotprod(0)
          ELSE IF (Lscale.eq.-2) THEN
            fac2=(SQRT(cg_Ritz(Ninner+1-nvec,nol))-1.0_r8)*Dotprod(0)
          ELSE IF (Lscale.eq.2) THEN
            fac2=(1.0_r8/SQRT(cg_Ritz(Ninner+1-nvec,nol))-1.0_r8)*      &
     &           Dotprod(0)
          END IF
!
          IF (.not.Ltrans) THEN
            IF (Lritz.and.(Lscale.eq.-2)) THEN
              fac2=fac2+facritz/SQRT(cg_Ritz(Ninner+1-nvec,nol))
            END IF
            IF (Lritz.and.(Lscale.eq.2)) THEN
              fac2=fac2-facritz/cg_Ritz(Ninner+1-nvec,nol)
            END IF
          END IF
!
          CALL state_addition (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         L1, L2, L1, fac1, fac2,                  &
# ifdef MASKING
     &                         rmask, umask, vmask,                     &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                         nl_t_obc, nl_t_obc,                      &
     &                         nl_u_obc, nl_u_obc,                      &
     &                         nl_v_obc, nl_v_obc,                      &
#  endif
     &                         nl_ubar_obc, nl_ubar_obc,                &
     &                         nl_vbar_obc, nl_vbar_obc,                &
     &                         nl_zeta_obc, nl_zeta_obc,                &
# endif
# ifdef ADJUST_WSTRESS
     &                         nl_ustr, nl_ustr,                        &
     &                         nl_vstr, nl_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                         nl_tflux, nl_tflux,                      &
#  endif
     &                         nl_t, nl_t,                              &
     &                         nl_u, nl_u,                              &
     &                         nl_v, nl_v,                              &
# else
     &                         nl_ubar, nl_ubar,                        &
     &                         nl_vbar, nl_vbar,                        &
# endif
     &                         nl_zeta, nl_zeta)
!
          IF (Lritz.and.Ltrans) THEN

            WRITE (ncname,20) TRIM(ADM(ng)%base), nol
            IF (Master.and.(nvec.eq.is)) THEN
              WRITE (stdout,50) outLoop, innLoop, TRIM(ncname)
            END IF
!
!  Read in the Lanczos vector q_k+1 computed from the incremental
!  4DVar algorithm first outer loop, where k=Ninner+1. Load Lanczos
!  vectors into NL state arrays at index L2.
!
            rec=Ninner+1
            CALL read_state (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       L2, rec,                                   &
     &                       ndef, ncid, ncname,                        &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       nl_t_obc, nl_u_obc, nl_v_obc,              &
#  endif
     &                       nl_ubar_obc, nl_vbar_obc,                  &
     &                       nl_zeta_obc,                               &
# endif
# ifdef ADJUST_WSTRESS
     &                       nl_ustr, nl_vstr,                          &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                       nl_tflux,                                  &
#  endif
     &                       nl_t, nl_u, nl_v,                          &
# else
     &                       nl_ubar, nl_vbar,                          &
# endif
     &                       nl_zeta)
            IF (exit_flag.ne.NoError) RETURN
!
            IF (Lscale.eq.2) THEN
              fac2=-facritz*Dotprod(0)/cg_Ritz(Ninner+1-nvec,nol)
            END IF
            IF (Lscale.eq.-2) THEN
              fac2=facritz*Dotprod(0)/SQRT(cg_Ritz(Ninner+1-nvec,nol))
            END IF
!
            CALL state_addition (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, LBij, UBij,        &
     &                           L1, L2, L1, fac1, fac2,                &
# ifdef MASKING
     &                           rmask, umask, vmask,                   &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                           nl_t_obc, nl_t_obc,                    &
     &                           nl_u_obc, nl_u_obc,                    &
     &                           nl_v_obc, nl_v_obc,                    &
#  endif
     &                           nl_ubar_obc, nl_ubar_obc,              &
     &                           nl_vbar_obc, nl_vbar_obc,              &
     &                           nl_zeta_obc, nl_zeta_obc,              &
# endif
# ifdef ADJUST_WSTRESS
     &                           nl_ustr, nl_ustr,                      &
     &                           nl_vstr, nl_vstr,                      &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                           nl_tflux, nl_tflux,                    &
#  endif
     &                           nl_t, nl_t,                            &
     &                           nl_u, nl_u,                            &
     &                           nl_v, nl_v,                            &
# else
     &                           nl_ubar, nl_ubar,                      &
     &                           nl_vbar, nl_vbar,                      &
# endif
     &                           nl_zeta, nl_zeta)

          END IF

        END DO
      END DO

 10   FORMAT (/,1x,'(',i3.3,',',i3.3,'): PRECOND -',1x,                 &
     &        a,1x,'preconditioning:',1x,a/)
 20   FORMAT (a,'_',i3.3,'.nc')
 30   FORMAT (1x,'(',i3.3,',',i3.3,'): PRECOND -',1x,                   &
     &        'Reading Lanczos eigenpairs from:',t58,a)
 40   FORMAT (1x,'(',i3.3,',',i3.3,'): PRECOND -',1x,                   &
     &        'Number of good Ritz eigenvalues,',t58,'ingood = ',i3)
 50   FORMAT (1x,'(',i3.3,',',i3.3,'): PRECOND -',1x,                   &
     &        'Processing Lanczos vectors from:',t58,a)
 60   FORMAT (1x,'(',i3.3,',',i3.3,'): PRECOND -',1x,                   &
     &        'Processing Hessian vectors from:',t58,a)

      RETURN
      END SUBROUTINE precond
#endif
      END MODULE cgradient_mod
