#include "cppdefs.h"
      MODULE ad_convolution_mod
#if defined ADJOINT && defined FOUR_DVAR
!
!svn $Id: ad_convolution.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine performs a spatial convolution of the adjoint state    !
!  solution to model the  background error correlations,  C,  using    !
!  a generalized diffusion operator.  This allows the observational    !
!  information to spread spatially in 4DVAR data assimilation.         !
!                                                                      !
!  The background error covariance is defined as:                      !
!                                                                      !
!    B = S C S                                                         !
!                                                                      !
!    C = C^(1/2) C^(T/2)                                               !
!                                                                      !
!    C^(1/2) = G L^(1/2) W^(-1/2)                               TLM    !
!    C^(T/2) = W^(-1/2) L^(T/2) G                               ADM    !
!                                                                      !
!  where                                                               !
!                                                                      !
!    B : background-error covariance matrix                            !
!    S : diagonal matrix of background-error standard deviations       !
!    C : symmetric matrix of background-error correlations             !
!    G : normalization coefficients matrix used to ensure that the     !
!          diagonal variances of C are equal to unity.                 !
!    L : tangent linear and adjoint diffusion operators                !
!    W : diagonal matrix of local area or volume metrics used to       !
!          convert L into a symmetric matrix: LW^(-1).                 !
!                                                                      !
!  Here, T/2 denote the transpose if a squared-root factor.            !
!                                                                      !
!  This routine is used to provide a  better preconditioning of the    !
!  minimization problem,  which is expressed as a function of a new    !
!  state vector, v, given by:                                          !
!                                                                      !
!                    v = B^(-1/2) delta_x      (v-space)               !
!  or                                                                  !
!              delta_x = B^(1/2) v                                     !
!                                                                      !
!  where                                                               !
!                                                                      !
!                    B = tranpose{B^(1/2)} B^(1/2)                     !
!                                                                      !
!  Therefore, the cost function, J, gradient becomes:                  !
!                                                                      !
!            GRAD_v(J) = v + transpose{B^(1/2)} GRAD_x(J)              !
!                                                                      !
!  In incremental  4DVAR,  these spatial convolutions constitutes a    !
!  smoothing action on the  correlation operator  and they are used    !
!  to transform between model space to minimization space and  vice    !
!  versa:                                                              !
!                                                                      !
!    ad_convolution     compute GRAD_v(J) from GRAD_x(J)               !
!    tl_convolution     compute x from v                               !
!                                                                      !
!  The minimization of of J in the descent algorithm is in v-space.    !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Weaver, A. and P. Courtier, 2001: Correlation modeling on the     !
!      sphere using a generalized diffusion equation, Q.J.R. Meteo.    !
!      Soc, 127, 1815-1846.                                            !
!                                                                      !
!======================================================================!
!
      USE mod_kinds

      implicit none

      PRIVATE
      PUBLIC :: ad_convolution

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_convolution (ng, tile, Linp, Lweak, ifac)
!***********************************************************************
!
      USE mod_param
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
# endif
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
# if defined SEDIMENT && defined SED_MORPH && defined SOLVE3D
      USE mod_sedbed
# endif
      USE mod_stepping
!
!  Imported variable declarations.
!
      logical, intent(in) :: Lweak

      integer, intent(in) :: ng, tile, Linp, ifac
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL ad_convolution_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, LBij, UBij,         &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nstp(ng), nnew(ng), Linp, Lweak, ifac,  &
     &                          GRID(ng) % pm,                          &
     &                          GRID(ng) % om_p,                        &
     &                          GRID(ng) % om_r,                        &
     &                          GRID(ng) % om_u,                        &
     &                          GRID(ng) % om_v,                        &
     &                          GRID(ng) % pn,                          &
     &                          GRID(ng) % on_p,                        &
     &                          GRID(ng) % on_r,                        &
     &                          GRID(ng) % on_u,                        &
     &                          GRID(ng) % on_v,                        &
     &                          GRID(ng) % pmon_p,                      &
     &                          GRID(ng) % pmon_r,                      &
     &                          GRID(ng) % pmon_u,                      &
     &                          GRID(ng) % pnom_p,                      &
     &                          GRID(ng) % pnom_r,                      &
     &                          GRID(ng) % pnom_v,                      &
# ifdef MASKING
     &                          GRID(ng) % rmask,                       &
     &                          GRID(ng) % pmask,                       &
     &                          GRID(ng) % umask,                       &
     &                          GRID(ng) % vmask,                       &
# endif
# ifdef SOLVE3D
     &                          GRID(ng) % h,                           &
#  ifdef ICESHELF
     &                          GRID(ng) % zice,                        &
#  endif
#  if defined SEDIMENT && defined SED_MORPH
     &                          SEDBED(ng) % bed_thick,                 &
#  endif
     &                          GRID(ng) % Hz,                          &
     &                          GRID(ng) % z_r,                         &
     &                          GRID(ng) % z_w,                         &
# endif
     &                          MIXING(ng) % Kh,                        &
# ifdef SOLVE3D
     &                          MIXING(ng) % Kv,                        &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                          BOUNDARY(ng) % b_t_obc,                 &
     &                          BOUNDARY(ng) % b_u_obc,                 &
     &                          BOUNDARY(ng) % b_v_obc,                 &
#  endif
     &                          BOUNDARY(ng) % b_ubar_obc,              &
     &                          BOUNDARY(ng) % b_vbar_obc,              &
     &                          BOUNDARY(ng) % b_zeta_obc,              &
# endif
# ifdef ADJUST_WSTRESS
     &                          FORCES(ng) % b_sustr,                   &
     &                          FORCES(ng) % b_svstr,                   &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                          FORCES(ng) % b_stflx,                   &
# endif
# ifdef SOLVE3D
     &                          OCEAN(ng) % b_t,                        &
     &                          OCEAN(ng) % b_u,                        &
     &                          OCEAN(ng) % b_v,                        &
# endif
     &                          OCEAN(ng) % b_zeta,                     &
     &                          OCEAN(ng) % b_ubar,                     &
     &                          OCEAN(ng) % b_vbar,                     &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                          BOUNDARY(ng) % ad_t_obc,                &
     &                          BOUNDARY(ng) % ad_u_obc,                &
     &                          BOUNDARY(ng) % ad_v_obc,                &
#  endif
     &                          BOUNDARY(ng) % ad_ubar_obc,             &
     &                          BOUNDARY(ng) % ad_vbar_obc,             &
     &                          BOUNDARY(ng) % ad_zeta_obc,             &
# endif
# ifdef ADJUST_WSTRESS
     &                          FORCES(ng) % ad_ustr,                   &
     &                          FORCES(ng) % ad_vstr,                   &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                          FORCES(ng) % ad_tflux,                  &
# endif
# ifdef SOLVE3D
     &                          OCEAN(ng) % ad_t,                       &
     &                          OCEAN(ng) % ad_u,                       &
     &                          OCEAN(ng) % ad_v,                       &
# endif
     &                          OCEAN(ng) % ad_ubar,                    &
     &                          OCEAN(ng) % ad_vbar,                    &
     &                          OCEAN(ng) % ad_zeta)
      RETURN
      END SUBROUTINE ad_convolution
!
!***********************************************************************
      SUBROUTINE ad_convolution_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, LBij, UBij,   &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                nstp, nnew, Linp, Lweak, ifac,    &
     &                                pm, om_p, om_r, om_u, om_v,       &
     &                                pn, on_p, on_r, on_u, on_v,       &
     &                                pmon_p, pmon_r, pmon_u,           &
     &                                pnom_p, pnom_r, pnom_v,           &
# ifdef MASKING
     &                                rmask, pmask, umask, vmask,       &
# endif
# ifdef SOLVE3D
     &                                h,                                &
#  ifdef ICESHELF
     &                                zice,                             &
#  endif
#  if defined SEDIMENT && defined SED_MORPH
     &                                bed_thick,                        &
#  endif
     &                                Hz, z_r, z_w,                     &
# endif
     &                                Kh,                               &
# ifdef SOLVE3D
     &                                Kv,                               &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                                VnormRobc, VnormUobc, VnormVobc,  &
#  endif
     &                                HnormRobc, HnormUobc, HnormVobc,  &
# endif
# ifdef ADJUST_WSTRESS
     &                                HnormSUS, HnormSVS,               &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                                HnormSTF,                         &
# endif
# ifdef SOLVE3D
     &                                VnormR, VnormU, VnormV,           &
# endif
     &                                HnormR, HnormU, HnormV,           &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                                ad_t_obc, ad_u_obc, ad_v_obc,     &
#  endif
     &                                ad_ubar_obc, ad_vbar_obc,         &
     &                                ad_zeta_obc,                      &
# endif
# ifdef ADJUST_WSTRESS
     &                                ad_ustr, ad_vstr,                 &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                                ad_tflux,                         &
# endif
# ifdef SOLVE3D
     &                                ad_t, ad_u, ad_v,                 &
# endif
     &                                ad_ubar, ad_vbar, ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_fourdvar
      USE mod_ncparam
      USE mod_scalars
!
      USE ad_conv_2d_mod
# ifdef SOLVE3D
      USE ad_conv_3d_mod
# endif
# ifdef ADJUST_BOUNDARY
      USE ad_conv_bry2d_mod
#  ifdef SOLVE3D
      USE ad_conv_bry3d_mod
#  endif
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod
# endif
      USE set_depth_mod
!
!  Imported variable declarations.
!
      logical, intent(in) :: Lweak

      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew, Linp, ifac
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: om_p(LBi:,LBj:)
      real(r8), intent(in) :: om_r(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: on_p(LBi:,LBj:)
      real(r8), intent(in) :: on_r(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: pmon_p(LBi:,LBj:)
      real(r8), intent(in) :: pmon_r(LBi:,LBj:)
      real(r8), intent(in) :: pmon_u(LBi:,LBj:)
      real(r8), intent(in) :: pnom_p(LBi:,LBj:)
      real(r8), intent(in) :: pnom_r(LBi:,LBj:)
      real(r8), intent(in) :: pnom_v(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: pmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Kh(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(in) :: Kv(LBi:,LBj:,0:)
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#   endif
#   if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(inout):: bed_thick(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: h(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent (in) :: VnormRobc(LBij:,:,:,:)
      real(r8), intent (in) :: VnormUobc(LBij:,:,:)
      real(r8), intent (in) :: VnormVobc(LBij:,:,:)
#   endif
      real(r8), intent (in) :: HnormRobc(LBij:,:)
      real(r8), intent (in) :: HnormUobc(LBij:,:)
      real(r8), intent (in) :: HnormVobc(LBij:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: HnormSUS(LBi:,LBj:)
      real(r8), intent(in) :: HnormSVS(LBi:,LBj:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(in) :: HnormSTF(LBi:,LBj:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: VnormR(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: VnormU(LBi:,LBj:,:,:)
      real(r8), intent(in) :: VnormV(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(in) :: HnormR(LBi:,LBj:,:)
      real(r8), intent(in) :: HnormU(LBi:,LBj:,:)
      real(r8), intent(in) :: HnormV(LBi:,LBj:,:)
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
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
#  ifdef SOLVE3D
      real(r8), intent(out) :: Hz(LBi:,LBj:,:)
      real(r8), intent(out) :: z_r(LBi:,LBj:,:)
      real(r8), intent(out) :: z_w(LBi:,LBj:,0:)
#  endif
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_v(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Kh(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(in) :: Kv(LBi:UBi,LBj:UBj,0:N(ng))
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#   endif
#   if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(inout):: bed_thick(LBi:UBi,LBj:UBj,2)
#   endif
      real(r8), intent(inout) :: h(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent (in) :: VnormRobc(LBij:UBij,N(ng),4,NT(ng))
      real(r8), intent (in) :: VnormUobc(LBij:UBij,N(ng),4)
      real(r8), intent (in) :: VnormVobc(LBij:UBij,N(ng),4)
#   endif
      real(r8), intent (in) :: HnormRobc(LBij:UBij,4)
      real(r8), intent (in) :: HnormUobc(LBij:UBij,4)
      real(r8), intent (in) :: HnormVobc(LBij:UBij,4)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: HnormSUS(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: HnormSVS(LBi:UBi,LBj:UBj)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(in) :: HnormSTF(LBi:UBi,LBj:UBj,NT(ng))
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: VnormR(LBi:UBi,LBj:UBj,N(ng),NSA,NT(ng))
      real(r8), intent(in) :: VnormU(LBi:UBi,LBj:UBj,NSA,N(ng))
      real(r8), intent(in) :: VnormV(LBi:UBi,LBj:UBj,NSA,N(ng))
#  endif
      real(r8), intent(in) :: HnormR(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(in) :: HnormU(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(in) :: HnormV(LBi:UBi,LBj:UBj,NSA)
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
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#  endif
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef SOLVE3D
      real(r8), intent(out) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
# endif
!
!  Local variable declarations.
!
# ifdef ADJUST_BOUNDARY
      logical, dimension(4) :: Lconvolve
# endif
      integer :: i, ib, ir, is, it, j, k, rec
      real(r8) :: cff
# ifdef SOLVE3D
      real(r8) :: fac
# endif
# ifdef SOLVE3D
      real(r8), dimension(LBi:UBi,LBj:UBj) :: work
# endif
!
# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Determine error covariance normalization factors to use.
!-----------------------------------------------------------------------
!
      IF (Lweak) THEN
        rec=2                        ! weak constraint
      ELSE
        rec=1                        ! strong constraint
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Set switch to convolve boundary segments by the appropriate
!  tiles.
!
      Lconvolve(iwest )=DOMAIN(ng)%Western_Edge (tile)
      Lconvolve(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      Lconvolve(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      Lconvolve(inorth)=DOMAIN(ng)%Northern_Edge(tile)
# endif

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute time invariant depths (use zero free-surface).
!-----------------------------------------------------------------------
!
      DO i=LBi,UBi
        DO j=LBj,UBj
          work(i,j)=0.0_r8
        END DO
      END DO

      CALL set_depth_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp, nnew,                                  &
     &                     h,                                           &
#  ifdef ICESHELF
     &                     zice,                                        &
#  endif
#  if defined SEDIMENT && defined SED_MORPH
     &                     bed_thick,                                   &
#  endif
     &                     work,                                        &
     &                     Hz, z_r, z_w)
# endif
!
!-----------------------------------------------------------------------
!  Multiply adjoint state by its corresponding normalization factor.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
      CALL ad_mp_exchange2d (ng, tile, iADM, 3,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_zeta(:,:,Linp),                         &
     &                       ad_ubar(:,:,Linp),                         &
     &                       ad_vbar(:,:,Linp))
# endif
!
!  Adjoint free-surface.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ad_zeta(i,j,Linp)=ad_zeta(i,j,Linp)*HnormR(i,j,rec)
        END DO
      END DO
!
!  Adjoint 2D momentum.
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          ad_ubar(i,j,Linp)=ad_ubar(i,j,Linp)*HnormU(i,j,rec)
        END DO
      END DO
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          ad_vbar(i,j,Linp)=ad_vbar(i,j,Linp)*HnormV(i,j,rec)
        END DO
      END DO
# ifdef SOLVE3D
!
!  Adjoint 3D momentum.
!
#  ifdef DISTRIBUTE
      CALL ad_mp_exchange3d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_u(:,:,:,Linp),                          &
     &                       ad_v(:,:,:,Linp))
#  endif
      DO k=1,N(ng)
        DO j=JstrT,JendT
          DO i=IstrP,IendT
            ad_u(i,j,k,Linp)=ad_u(i,j,k,Linp)*VnormU(i,j,k,rec)
          END DO
        END DO
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            ad_v(i,j,k,Linp)=ad_v(i,j,k,Linp)*VnormV(i,j,k,rec)
          END DO
        END DO
      END DO
!
!  Adjoint tracers.
!
#  ifdef DISTRIBUTE
      CALL ad_mp_exchange4d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),   &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_t(:,:,:,Linp,:))
#  endif
      DO it=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              ad_t(i,j,k,Linp,it)=ad_t(i,j,k,Linp,it)*                  &
     &                            VnormR(i,j,k,rec,it)
            END DO
          END DO
        END DO
      END DO
# endif

# ifdef ADJUST_BOUNDARY
!
!  Adjoint free-surface open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isFsur,ng)) THEN
#  ifdef DISTRIBUTE
            CALL ad_mp_exchange2d_bry (ng, tile, iADM, 1, ib,           &
     &                                 LBij, UBij,                      &
     &                                 NghostPoints,                    &
     &                                 EWperiodic(ng), NSperiodic(ng),  &
     &                                 ad_zeta_obc(:,ib,ir,Linp))
#  endif
            IF (Lconvolve(ib)) THEN
              IF ((ib.eq.iwest).or.(ib.eq.ieast)) THEN
                DO j=JstrT,JendT
                  ad_zeta_obc(j,ib,ir,Linp)=ad_zeta_obc(j,ib,ir,Linp)*  &
     &                                      HnormRobc(j,ib)
                END DO
              ELSE IF ((ib.eq.isouth).or.(ib.eq.inorth)) THEN
                DO i=IstrT,IendT
                  ad_zeta_obc(i,ib,ir,Linp)=ad_zeta_obc(i,ib,ir,Linp)*  &
     &                                      HnormRobc(i,ib)
                END DO
              END IF
            END IF
          END IF
        END DO
      END DO
!
!  Tangent linear 2D U-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isUbar,ng)) THEN
#  ifdef DISTRIBUTE
            CALL ad_mp_exchange2d_bry (ng, tile, iADM, 1, ib,           &
     &                                 LBij, UBij,                      &
     &                                 NghostPoints,                    &
     &                                 EWperiodic(ng), NSperiodic(ng),  &
     &                                 ad_ubar_obc(:,ib,ir,Linp))
#  endif
            IF (Lconvolve(ib)) THEN
              IF ((ib.eq.iwest).or.(ib.eq.ieast)) THEN
                DO j=JstrT,JendT
                  ad_ubar_obc(j,ib,ir,Linp)=ad_ubar_obc(j,ib,ir,Linp)*  &
     &                                      HnormUobc(j,ib)
                END DO
              ELSE IF ((ib.eq.isouth).or.(ib.eq.inorth)) THEN
                DO i=IstrP,IendT
                  ad_ubar_obc(i,ib,ir,Linp)=ad_ubar_obc(i,ib,ir,Linp)*  &
     &                                      HnormUobc(i,ib)
                END DO
              END IF
            END IF
          END IF
        END DO
      END DO
!
!  Tangent linear 2D V-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isVbar,ng)) THEN
#  ifdef DISTRIBUTE
            CALL ad_mp_exchange2d_bry (ng, tile, iADM, 1, ib,           &
     &                                 LBij, UBij,                      &
     &                                 NghostPoints,                    &
     &                                 EWperiodic(ng), NSperiodic(ng),  &
     &                                 ad_vbar_obc(:,ib,ir,Linp))
#  endif
            IF (Lconvolve(ib)) THEN
              IF ((ib.eq.iwest).or.(ib.eq.ieast)) THEN
                DO j=JstrP,JendT
                  ad_vbar_obc(j,ib,ir,Linp)=ad_vbar_obc(j,ib,ir,Linp)*  &
     &                                      HnormVobc(j,ib)
                END DO
              ELSE IF ((ib.eq.isouth).or.(ib.eq.inorth)) THEN
                DO i=IstrT,IendT
                  ad_vbar_obc(i,ib,ir,Linp)=ad_vbar_obc(i,ib,ir,Linp)*  &
     &                                      HnormVobc(i,ib)
                END DO
              END IF
            END IF
          END IF
        END DO
      END DO

#  ifdef SOLVE3D
!
!  Tangent linear 3D U-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isUvel,ng)) THEN
#   ifdef DISTRIBUTE
            CALL ad_mp_exchange3d_bry (ng, tile, iADM, 1, ib,           &
     &                                 LBij, UBij, 1, N(ng),            &
     &                                 NghostPoints,                    &
     &                                 EWperiodic(ng), NSperiodic(ng),  &
     &                                 ad_u_obc(:,:,ib,ir,Linp))
#   endif
            IF (Lconvolve(ib)) THEN
              IF ((ib.eq.iwest).or.(ib.eq.ieast)) THEN
                DO k=1,N(ng)
                  DO j=JstrT,JendT
                    ad_u_obc(j,k,ib,ir,Linp)=ad_u_obc(j,k,ib,ir,Linp)*  &
     &                                       VnormUobc(j,k,ib)
                  END DO
                END DO
              ELSE IF ((ib.eq.isouth).or.(ib.eq.inorth)) THEN
                DO k=1,N(ng)
                  DO i=IstrP,IendT
                    ad_u_obc(i,k,ib,ir,Linp)=ad_u_obc(i,k,ib,ir,Linp)*  &
     &                                       VnormUobc(i,k,ib)
                  END DO
                END DO
              END IF
            END IF
          END IF
        END DO
      END DO
!
!  Tangent linear 3D V-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isVvel,ng)) THEN
#   ifdef DISTRIBUTE
            CALL ad_mp_exchange3d_bry (ng, tile, iADM, 1, ib,           &
     &                                 LBij, UBij, 1, N(ng),            &
     &                                 NghostPoints,                    &
     &                                 EWperiodic(ng), NSperiodic(ng),  &
     &                                 ad_v_obc(:,:,ib,ir,Linp))
#   endif
            IF (Lconvolve(ib)) THEN
              IF ((ib.eq.iwest).or.(ib.eq.ieast)) THEN
                DO k=1,N(ng)
                  DO j=JstrP,JendT
                    ad_v_obc(j,k,ib,ir,Linp)=ad_v_obc(j,k,ib,ir,Linp)*  &
     &                                       VnormVobc(j,k,ib)
                  END DO
                END DO
              ELSE IF ((ib.eq.isouth).or.(ib.eq.inorth)) THEN
                DO k=1,N(ng)
                  DO i=IstrT,IendT
                    ad_v_obc(i,k,ib,ir,Linp)=ad_v_obc(i,k,ib,ir,Linp)*  &
     &                                       VnormVobc(i,k,ib)
                  END DO
                END DO
              END IF
            END IF
          END IF
        END DO
      END DO
!
!  Tangent linear tracers open boundaries.
!
      DO it=1,NT(ng)
        DO ir=1,Nbrec(ng)
          DO ib=1,4
            IF (.not.Lweak.and.Lobc(ib,isTvar(it),ng)) THEN
#   ifdef DISTRIBUTE
              CALL ad_mp_exchange3d_bry (ng, tile, iADM, 1, ib,         &
     &                                   LBij, UBij, 1, N(ng),          &
     &                                   NghostPoints,                  &
     &                                   EWperiodic(ng), NSperiodic(ng),&
     &                                   ad_t_obc(:,:,ib,ir,Linp,it))
#   endif
              IF (Lconvolve(ib)) THEN
                IF ((ib.eq.iwest).or.(ib.eq.ieast)) THEN
                  DO k=1,N(ng)
                    DO j=JstrT,JendT
                      ad_t_obc(j,k,ib,ir,Linp,it)=                      &
     &                                   ad_t_obc(j,k,ib,ir,Linp,it)*   &
     &                                   VnormRobc(j,k,ib,it)
                    END DO
                  END DO
                ELSE IF ((ib.eq.isouth).or.(ib.eq.inorth)) THEN
                  DO k=1,N(ng)
                    DO i=IstrT,IendT
                      ad_t_obc(i,k,ib,ir,Linp,it)=                      &
     &                                   ad_t_obc(i,k,ib,ir,Linp,it)*   &
     &                                   VnormRobc(i,k,ib,it)
                    END DO
                  END DO
                END IF
              END IF
            END IF
          END DO
        END DO
      END DO
#  endif
# endif

# ifdef ADJUST_WSTRESS
!
!  Adjoint surface momentum stress.
!
      IF (.not.Lweak) THEN
#  ifdef DISTRIBUTE
        CALL ad_mp_exchange3d (ng, tile, iADM, 2,                       &
     &                         LBi, UBi, LBj, UBj, 1, Nfrec(ng),        &
     &                         NghostPoints,                            &
     &                         EWperiodic(ng), NSperiodic(ng),          &
     &                         ad_ustr(:,:,:,Linp),                     &
     &                         ad_vstr(:,:,:,Linp))
#  endif
        DO k=1,Nfrec(ng)
          DO j=JstrT,JendT
            DO i=IstrP,IendT
              ad_ustr(i,j,k,Linp)=ad_ustr(i,j,k,Linp)*HnormSUS(i,j)
            END DO
          END DO
          DO j=JstrP,JendT
            DO i=IstrT,IendT
              ad_vstr(i,j,k,Linp)=ad_vstr(i,j,k,Linp)*HnormSVS(i,j)
            END DO
          END DO
        END DO
      END IF
# endif
# ifdef ADJUST_STFLUX
!
!  Adjoint surface tracers flux.
!
      IF (.not.Lweak) THEN
        DO it=1,NT(ng)
          IF (Lstflux(it,ng)) THEN
#  ifdef DISTRIBUTE
            CALL ad_mp_exchange3d (ng, tile, iADM, 1,                   &
     &                             LBi, UBi, LBj, UBj, 1, Nfrec(ng),    &
     &                             NghostPoints,                        &
     &                             EWperiodic(ng), NSperiodic(ng),      &
     &                             ad_tflux(:,:,:,Linp,it))
#  endif
            DO k=1,Nfrec(ng)
              DO j=JstrT,JendT
                DO i=IstrT,IendT
                  ad_tflux(i,j,k,Linp,it)=ad_tflux(i,j,k,Linp,it)*      &
     &                                    HnormSTF(i,j,it)
                END DO
              END DO
            END DO
          END IF
        END DO
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Initial conditions or model error covariance: Convolve adjoint state
!  vector with a generalized adjoint diffusion equation to filter
!  solution with specified horizontal scales. Convert from model space
!  to minimization space (v-space).
!-----------------------------------------------------------------------
!
!  Adjoint free-surface.
!
      CALL ad_conv_r2d_tile (ng, tile, iADM,                            &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       NghostPoints,                              &
     &                       NHsteps(rec,isFsur)/ifac,                  &
     &                       DTsizeH(rec,isFsur),                       &
     &                       Kh,                                        &
     &                       pm, pn, pmon_u, pnom_v,                    &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
     &                       ad_zeta(:,:,Linp))
!
!  Adjoint 2D momentum.
!
      CALL ad_conv_u2d_tile (ng, tile, iADM,                            &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       NghostPoints,                              &
     &                       NHsteps(rec,isUbar)/ifac,                  &
     &                       DTsizeH(rec,isUbar),                       &
     &                       Kh,                                        &
     &                       pm, pn, pmon_r, pnom_p,                    &
# ifdef MASKING
     &                       umask, pmask,                              &
# endif
     &                       ad_ubar(:,:,Linp))

      CALL ad_conv_v2d_tile (ng, tile, iADM,                            &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       NghostPoints,                              &
     &                       NHsteps(rec,isVbar)/ifac,                  &
     &                       DTsizeH(rec,isVbar),                       &
     &                       Kh,                                        &
     &                       pm, pn, pmon_p, pnom_r,                    &
# ifdef MASKING
     &                       vmask, pmask,                              &
# endif
     &                       ad_vbar(:,:,Linp))
# ifdef SOLVE3D
!
!  Adjoint 3D momentum.
!
      CALL ad_conv_u3d_tile (ng, tile, iADM,                            &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       NghostPoints,                              &
     &                       NHsteps(rec,isUvel)/ifac,                  &
     &                       NVsteps(rec,isUvel)/ifac,                  &
     &                       DTsizeH(rec,isUvel),                       &
     &                       DTsizeV(rec,isUvel),                       &
     &                       Kh, Kv,                                    &
     &                       pm, pn,                                    &
#  ifdef GEOPOTENTIAL_HCONV
     &                       on_r, om_p,                                &
#  else
     &                       pmon_r, pnom_p,                            &
#  endif
#  ifdef MASKING
#   ifdef GEOPOTENTIAL_HCONV
     &                       pmask, rmask, umask, vmask,                &
#   else
     &                       umask, pmask,                              &
#   endif
#  endif
     &                       Hz, z_r,                                   &
     &                       ad_u(:,:,:,Linp))

      CALL ad_conv_v3d_tile (ng, tile, iADM,                            &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       NghostPoints,                              &
     &                       NHsteps(rec,isUvel)/ifac,                  &
     &                       NVsteps(rec,isUvel)/ifac,                  &
     &                       DTsizeH(rec,isUvel),                       &
     &                       DTsizeV(rec,isUvel),                       &
     &                       Kh, Kv,                                    &
     &                       pm, pn,                                    &
#  ifdef GEOPOTENTIAL_HCONV
     &                       on_p, om_r,                                &
#  else
     &                       pmon_p, pnom_r,                            &
#  endif
#  ifdef MASKING
#   ifdef GEOPOTENTIAL_HCONV
     &                       pmask, rmask, umask, vmask,                &
#   else
     &                       vmask, pmask,                              &
#   endif
#  endif
     &                       Hz, z_r,                                   &
     &                       ad_v(:,:,:,Linp))
!
!  Adjoint tracers.
!
      DO it=1,NT(ng)
        is=isTvar(it)
        CALL ad_conv_r3d_tile (ng, tile, iADM,                          &
     &                         LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         NghostPoints,                            &
     &                         NHsteps(rec,is)/ifac,                    &
     &                         NVsteps(rec,is)/ifac,                    &
     &                         DTsizeH(rec,is),                         &
     &                         DTsizeV(rec,is),                         &
     &                         Kh, Kv,                                  &
     &                         pm, pn,                                  &
#  ifdef GEOPOTENTIAL_HCONV
     &                         on_u, om_v,                              &
#  else
     &                         pmon_u, pnom_v,                          &
#  endif
#  ifdef MASKING
     &                         rmask, umask, vmask,                     &
#  endif
     &                         Hz, z_r,                                 &
     &                         ad_t(:,:,:,Linp,it))
      END DO
# endif

# ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  Open boundaries error convariance:  Convolve adjoint state boundary
!  edges with a generalized adjoint diffusion equation to filter
!  solution with specified horizontal scales. Convert from model space
!  to minimization space (v-space).
!-----------------------------------------------------------------------
!
!  Adjoint free-surface open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isFsur,ng)) THEN
            CALL ad_conv_r2d_bry_tile (ng, tile, iADM, ib,              &
     &                                 BOUNDS(ng)%edge(:,r2dvar),       &
     &                                 LBij, UBij,                      &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 NghostPoints,                    &
     &                                 NHstepsB(ib,isFsur)/ifac,        &
     &                                 DTsizeHB(ib,isFsur),             &
     &                                 Kh,                              &
     &                                 pm, pn, pmon_u, pnom_v,          &
#  ifdef MASKING
     &                                 rmask, umask, vmask,             &
#  endif
     &                                 ad_zeta_obc(:,ib,ir,Linp))
          END IF
        END DO
      END DO
!
!  Tangent linear 2D U-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isUbar,ng)) THEN
            CALL ad_conv_u2d_bry_tile (ng, tile, iADM, ib,              &
     &                                 BOUNDS(ng)%edge(:,u2dvar),       &
     &                                 LBij, UBij,                      &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 NghostPoints,                    &
     &                                 NHstepsB(ib,isUbar)/ifac,        &
     &                                 DTsizeHB(ib,isUbar),             &
     &                                 Kh,                              &
     &                                 pm, pn, pmon_r, pnom_p,          &
#  ifdef MASKING
     &                                 umask, pmask,                    &
#  endif
     &                                 ad_ubar_obc(:,ib,ir,Linp))
          END IF
        END DO
      END DO
!
!  Tangent linear 2D V-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isVbar,ng)) THEN
            CALL ad_conv_v2d_bry_tile (ng, tile, iADM, ib,              &
     &                                 BOUNDS(ng)%edge(:,v2dvar),       &
     &                                 LBij, UBij,                      &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 NghostPoints,                    &
     &                                 NHstepsB(ib,isVbar)/ifac,        &
     &                                 DTsizeHB(ib,isVbar),             &
     &                                 Kh,                              &
     &                                 pm, pn, pmon_p, pnom_r,          &
#  ifdef MASKING
     &                                 vmask, pmask,                    &
#  endif
     &                                 ad_vbar_obc(:,ib,ir,Linp))
          END IF
        END DO
      END DO

#  ifdef SOLVE3D
!
!  Tangent linear 3D U-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isUvel,ng)) THEN
            CALL ad_conv_u3d_bry_tile (ng, tile, iADM, ib,              &
     &                                 BOUNDS(ng)%edge(:,u2dvar),       &
     &                                 LBij, UBij,                      &
     &                                 LBi, UBi, LBj, UBj, 1, N(ng),    &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 NghostPoints,                    &
     &                                 NHstepsB(ib,isUvel)/ifac,        &
     &                                 NVstepsB(ib,isUvel)/ifac,        &
     &                                 DTsizeHB(ib,isUvel),             &
     &                                 DTsizeVB(ib,isUvel),             &
     &                                 Kh, Kv,                          &
     &                                 pm, pn, pmon_r, pnom_p,          &
#   ifdef MASKING
     &                                 umask, pmask,                    &
#   endif
     &                                 Hz, z_r,                         &
     &                                 ad_u_obc(:,:,ib,ir,Linp))
          END IF
        END DO
      END DO
!
!  Tangent linear 3D V-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isVvel,ng)) THEN
            CALL ad_conv_v3d_bry_tile (ng, tile, iADM, ib,              &
     &                                 BOUNDS(ng)%edge(:,v2dvar),       &
     &                                 LBij, UBij,                      &
     &                                 LBi, UBi, LBj, UBj, 1, N(ng),    &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 NghostPoints,                    &
     &                                 NHstepsB(ib,isVvel)/ifac,        &
     &                                 NVstepsB(ib,isVvel)/ifac,        &
     &                                 DTsizeHB(ib,isVvel),             &
     &                                 DTsizeVB(ib,isVvel),             &
     &                                 Kh, Kv,                          &
     &                                 pm, pn, pmon_p, pnom_r,          &
#   ifdef MASKING
     &                                 vmask, pmask,                    &
#   endif
     &                                 Hz, z_r,                         &
     &                                 ad_v_obc(:,:,ib,ir,Linp))
          END IF
        END DO
      END DO
!
!  Tangent linear tracers open boundaries.
!
      DO it=1,NT(ng)
        is=isTvar(it)
        DO ir=1,Nbrec(ng)
          DO ib=1,4
            IF (.not.Lweak.and.Lobc(ib,is,ng)) THEN
              CALL ad_conv_r3d_bry_tile (ng, tile, iADM, ib,            &
     &                                   BOUNDS(ng)%edge(:,r2dvar),     &
     &                                   LBij, UBij,                    &
     &                                   LBi, UBi, LBj, UBj, 1, N(ng),  &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   NghostPoints,                  &
     &                                   NHstepsB(ib,is)/ifac,          &
     &                                   NVstepsB(ib,is)/ifac,          &
     &                                   DTsizeHB(ib,is),               &
     &                                   DTsizeVB(ib,is),               &
     &                                   Kh, Kv,                        &
     &                                   pm, pn, pmon_u, pnom_v,        &
#   ifdef MASKING
     &                                   rmask, umask, vmask,           &
#   endif
     &                                   Hz, z_r,                       &
     &                                   ad_t_obc(:,:,ib,ir,Linp,it))
            END IF
          END DO
        END DO
      END DO
#  endif
# endif

# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!-----------------------------------------------------------------------
!  Surface forcing error covariance: Convolve adjoint state vector with
!  a generalized adjoint diffusion equation to filter solution with
!  specified horizontal scales. Convert from model spaceto minimization
!  space (v-space).
!-----------------------------------------------------------------------

#  ifdef ADJUST_WSTRESS
!
!  Adjoint surface momentum stress.
!
      IF (.not.Lweak) THEN
        DO k=1,Nfrec(ng)
          CALL ad_conv_u2d_tile (ng, tile, iADM,                        &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           NghostPoints,                          &
     &                           NHsteps(rec,isUstr)/ifac,              &
     &                           DTsizeH(rec,isUstr),                   &
     &                           Kh,                                    &
     &                           pm, pn, pmon_r, pnom_p,                &
#   ifdef MASKING
     &                           umask, pmask,                          &
#   endif
     &                           ad_ustr(:,:,k,Linp))

          CALL ad_conv_v2d_tile (ng, tile, iADM,                        &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           NghostPoints,                          &
     &                           NHsteps(rec,isVstr)/ifac,              &
     &                           DTsizeH(rec,isVstr),                   &
     &                           Kh,                                    &
     &                           pm, pn, pmon_p, pnom_r,                &
#   ifdef MASKING
     &                           vmask, pmask,                          &
#   endif
     &                           ad_vstr(:,:,k,Linp))
        END DO
      END IF
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Adjoint surface tracers flux.
!
      IF (.not.Lweak) THEN
        DO it=1,NT(ng)
          IF (Lstflux(it,ng)) THEN
            is=isTsur(it)
            DO k=1,Nfrec(ng)
              CALL ad_conv_r2d_tile (ng, tile, iADM,                    &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               NghostPoints,                      &
     &                               NHsteps(rec,is)/ifac,              &
     &                               DTsizeH(rec,is),                   &
     &                               Kh,                                &
     &                               pm, pn, pmon_u, pnom_v,            &
#   ifdef MASKING
     &                               rmask, umask, vmask,               &
#   endif
     &                               ad_tflux(:,:,k,Linp,it))
            END DO
          END IF
        END DO
      END IF
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Multiply convolved adjoint state by the inverse squared root of its
!  associated area (2D) or volume (3D).
!-----------------------------------------------------------------------
!
!  Adjoint free-surface.
!
# ifdef DISTRIBUTE
      CALL ad_mp_exchange2d (ng, tile, iADM, 3,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_zeta(:,:,Linp),                         &
     &                       ad_ubar(:,:,Linp),                         &
     &                       ad_vbar(:,:,Linp))
# endif
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ad_zeta(i,j,Linp)=ad_zeta(i,j,Linp)/                          &
     &                      SQRT(om_r(i,j)*on_r(i,j))
        END DO
      END DO
!
!  Adjoint 2D momentum.
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          ad_ubar(i,j,Linp)=ad_ubar(i,j,Linp)/                          &
     &                      SQRT(om_u(i,j)*on_u(i,j))
        END DO
      END DO
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          ad_vbar(i,j,Linp)=ad_vbar(i,j,Linp)/                          &
     &                      SQRT(om_v(i,j)*on_v(i,j))
        END DO
      END DO
# ifdef SOLVE3D
!
!  Adjoint 3D momentum.
!
#  ifdef DISTRIBUTE
      CALL ad_mp_exchange3d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_u(:,:,:,Linp),                          &
     &                       ad_v(:,:,:,Linp))
#  endif
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          cff=om_u(i,j)*on_u(i,j)*0.5_r8
          DO k=1,N(ng)
            ad_u(i,j,k,Linp)=ad_u(i,j,k,Linp)/                          &
     &                       SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
          END DO
        END DO
      END DO
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          cff=om_v(i,j)*on_v(i,j)*0.5_r8
          DO k=1,N(ng)
            ad_v(i,j,k,Linp)=ad_v(i,j,k,Linp)/                          &
     &                       SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
          END DO
        END DO
      END DO
!
!  Adjoint tracers.
!
#  ifdef DISTRIBUTE
      CALL ad_mp_exchange4d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),   &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_t(:,:,:,Linp,:))
#  endif
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          cff=om_r(i,j)*on_r(i,j)
          DO k=1,N(ng)
            fac=1.0_r8/SQRT(cff*Hz(i,j,k))
            DO it=1,NT(ng)
              ad_t(i,j,k,Linp,it)=fac*ad_t(i,j,k,Linp,it)
            END DO
          END DO
        END DO
      END DO
# endif

# ifdef ADJUST_BOUNDARY
!
!  Adjoint free-surface open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isFsur,ng)) THEN
#  ifdef DISTRIBUTE
            CALL ad_mp_exchange2d_bry (ng, tile, iADM, 1, ib,           &
     &                                 LBij, UBij,                      &
     &                                 NghostPoints,                    &
     &                                 EWperiodic(ng), NSperiodic(ng),  &
     &                                 ad_zeta_obc(:,ib,ir,Linp))
#  endif
            IF (Lconvolve(ib)) THEN
              IF ((ib.eq.iwest).or.(ib.eq.ieast)) THEN
                i=BOUNDS(ng)%edge(ib,r2dvar)
                DO j=JstrT,JendT
                  ad_zeta_obc(j,ib,ir,Linp)=ad_zeta_obc(j,ib,ir,Linp)/  &
     &                                      SQRT(on_r(i,j))
                END DO
              ELSE IF ((ib.eq.isouth).or.(ib.eq.inorth)) THEN
                j=BOUNDS(ng)%edge(ib,r2dvar)
                DO i=IstrT,IendT
                  ad_zeta_obc(i,ib,ir,Linp)=ad_zeta_obc(i,ib,ir,Linp)/  &
     &                                      SQRT(om_r(i,j))
                END DO
              END IF
            END IF
          END IF
        END DO
      END DO
!
!  Tangent linear 2D U-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isUbar,ng)) THEN
#  ifdef DISTRIBUTE
            CALL ad_mp_exchange2d_bry (ng, tile, iADM, 1, ib,           &
     &                                 LBij, UBij,                      &
     &                                 NghostPoints,                    &
     &                                 EWperiodic(ng), NSperiodic(ng),  &
     &                                 ad_ubar_obc(:,ib,ir,Linp))
#  endif
            IF (Lconvolve(ib)) THEN
              IF ((ib.eq.iwest).or.(ib.eq.ieast)) THEN
                i=BOUNDS(ng)%edge(ib,u2dvar)
                DO j=JstrT,JendT
                  ad_ubar_obc(j,ib,ir,Linp)=ad_ubar_obc(j,ib,ir,Linp)/  &
     &                                      SQRT(on_u(i,j))
                END DO
              ELSE IF ((ib.eq.isouth).or.(ib.eq.inorth)) THEN
                j=BOUNDS(ng)%edge(ib,u2dvar)
                DO i=IstrP,IendT
                  ad_ubar_obc(i,ib,ir,Linp)=ad_ubar_obc(i,ib,ir,Linp)/  &
     &                                      SQRT(om_u(i,j))
                END DO
              END IF
            END IF
          END IF
        END DO
      END DO
!
!  Tangent linear 2D V-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isVbar,ng)) THEN
#  ifdef DISTRIBUTE
            CALL ad_mp_exchange2d_bry (ng, tile, iADM, 1, ib,           &
     &                                 LBij, UBij,                      &
     &                                 NghostPoints,                    &
     &                                 EWperiodic(ng), NSperiodic(ng),  &
     &                                 ad_vbar_obc(:,ib,ir,Linp))
#  endif
            IF (Lconvolve(ib)) THEN
              IF ((ib.eq.iwest).or.(ib.eq.ieast)) THEN
                i=BOUNDS(ng)%edge(ib,v2dvar)
                DO j=JstrP,JendT
                  ad_vbar_obc(j,ib,ir,Linp)=ad_vbar_obc(j,ib,ir,Linp)/  &
     &                                      SQRT(on_v(i,j))
                END DO
              ELSE IF ((ib.eq.isouth).or.(ib.eq.inorth)) THEN
                j=BOUNDS(ng)%edge(ib,v2dvar)
                DO i=IstrT,IendT
                  ad_vbar_obc(i,ib,ir,Linp)=ad_vbar_obc(i,ib,ir,Linp)/  &
     &                                      SQRT(om_v(i,j))
                END DO
              END IF
            END IF
          END IF
        END DO
      END DO

#  ifdef SOLVE3D
!
!  Tangent linear 3D U-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isUvel,ng)) THEN
#   ifdef DISTRIBUTE
            CALL ad_mp_exchange3d_bry (ng, tile, iADM, 1, ib,           &
     &                                 LBij, UBij, 1, N(ng),            &
     &                                 NghostPoints,                    &
     &                                 EWperiodic(ng), NSperiodic(ng),  &
     &                                 ad_u_obc(:,:,ib,ir,Linp))
#   endif
            IF (Lconvolve(ib)) THEN
              IF ((ib.eq.iwest).or.(ib.eq.ieast)) THEN
                i=BOUNDS(ng)%edge(ib,u2dvar)
                DO j=JstrT,JendT
                  cff=on_u(i,j)*0.5_r8
                  DO k=1,N(ng)
                    ad_u_obc(j,k,ib,ir,Linp)=ad_u_obc(j,k,ib,ir,Linp)/  &
     &                                       SQRT(cff*                  &
     &                                            (Hz(i-1,j,k)+         &
     &                                             Hz(i  ,j,k)))
                  END DO
                END DO
              ELSE IF ((ib.eq.isouth).or.(ib.eq.inorth)) THEN
                j=BOUNDS(ng)%edge(ib,u2dvar)
                DO i=IstrP,IendT
                  cff=om_u(i,j)*0.5_r8
                  DO k=1,N(ng)
                    ad_u_obc(i,k,ib,ir,Linp)=ad_u_obc(i,k,ib,ir,Linp)/  &
     &                                       SQRT(cff*                  &
     &                                            (Hz(i-1,j,k)+         &
     &                                             Hz(i  ,j,k)))
                  END DO
                END DO
              END IF
            END IF
          END IF
        END DO
      END DO
!
!  Tangent linear 3D V-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isVvel,ng)) THEN
#   ifdef DISTRIBUTE
            CALL ad_mp_exchange3d_bry (ng, tile, iADM, 1, ib,           &
     &                                 LBij, UBij, 1, N(ng),            &
     &                                 NghostPoints,                    &
     &                                 EWperiodic(ng), NSperiodic(ng),  &
     &                                 ad_v_obc(:,:,ib,ir,Linp))
#   endif
            IF (Lconvolve(ib)) THEN
              IF ((ib.eq.iwest).or.(ib.eq.ieast)) THEN
                i=BOUNDS(ng)%edge(ib,v2dvar)
                DO j=JstrP,JendT
                  cff=on_v(i,j)*0.5_r8
                  DO k=1,N(ng)
                    ad_v_obc(j,k,ib,ir,Linp)=ad_v_obc(j,k,ib,ir,Linp)/  &
     &                                       SQRT(cff*                  &
     &                                            (Hz(i,j-1,k)+         &
     &                                             Hz(i,j  ,k)))
                  END DO
                END DO
              ELSE IF ((ib.eq.isouth).or.(ib.eq.inorth)) THEN
                j=BOUNDS(ng)%edge(ib,v2dvar)
                DO i=IstrT,IendT
                  cff=om_v(i,j)*0.5_r8
                  DO k=1,N(ng)
                    ad_v_obc(i,k,ib,ir,Linp)=ad_v_obc(i,k,ib,ir,Linp)/  &
     &                                       SQRT(cff*                  &
     &                                            (Hz(i,j-1,k)+         &
     &                                             Hz(i,j  ,k)))
                  END DO
                END DO
              END IF
            END IF
          END IF
        END DO
      END DO
!
!  Tangent linear tracers open boundaries.
!
      DO it=1,NT(ng)
        DO ir=1,Nbrec(ng)
          DO ib=1,4
            IF (.not.Lweak.and.Lobc(ib,isTvar(it),ng)) THEN
#   ifdef DISTRIBUTE
              CALL ad_mp_exchange3d_bry (ng, tile, iADM, 1, ib,         &
     &                                   LBij, UBij, 1, N(ng),          &
     &                                   NghostPoints,                  &
     &                                   EWperiodic(ng), NSperiodic(ng),&
     &                                   ad_t_obc(:,:,ib,ir,Linp,it))
#   endif
              IF (Lconvolve(ib)) THEN
                IF ((ib.eq.iwest).or.(ib.eq.ieast)) THEN
                  i=BOUNDS(ng)%edge(ib,r2dvar)
                  DO j=JstrT,JendT
                    cff=on_r(i,j)
                    DO k=1,N(ng)
                      ad_t_obc(j,k,ib,ir,Linp,it)=                      &
     &                                   ad_t_obc(j,k,ib,ir,Linp,it)/   &
     &                                   SQRT(cff*Hz(i,j,k))
                    END DO
                  END DO
                ELSE IF ((ib.eq.isouth).or.(ib.eq.inorth)) THEN
                  j=BOUNDS(ng)%edge(ib,r2dvar)
                  DO i=IstrT,IendT
                    cff=om_r(i,j)
                    DO k=1,N(ng)
                      ad_t_obc(i,k,ib,ir,Linp,it)=                      &
     &                                   ad_t_obc(i,k,ib,ir,Linp,it)/   &
     &                                   SQRT(cff*Hz(i,j,k))
                    END DO
                  END DO
                END IF
              END IF
            END IF
          END DO
        END DO
      END DO
#  endif
# endif

# ifdef ADJUST_WSTRESS
!
!  Adjoint surface momentum stress.
!
      IF (.not.Lweak) THEN
#  ifdef DISTRIBUTE
        CALL ad_mp_exchange3d (ng, tile, iADM, 2,                       &
     &                         LBi, UBi, LBj, UBj, 1, Nfrec(ng),        &
     &                         NghostPoints,                            &
     &                         EWperiodic(ng), NSperiodic(ng),          &
     &                         ad_ustr(:,:,:,Linp),                     &
     &                         ad_vstr(:,:,:,Linp))
#  endif
        DO k=1,Nfrec(ng)
          DO j=JstrT,JendT
            DO i=IstrP,IendT
              ad_ustr(i,j,k,Linp)=ad_ustr(i,j,k,Linp)/                  &
     &                            SQRT(om_u(i,j)*on_u(i,j))
            END DO
          END DO
          DO j=JstrP,JendT
            DO i=IstrT,IendT
              ad_vstr(i,j,k,Linp)=ad_vstr(i,j,k,Linp)/                  &
     &                            SQRT(om_v(i,j)*on_v(i,j))
            END DO
          END DO
        END DO
      END IF
# endif

# if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Adjoint surface tracers flux.
!
      IF (.not.Lweak) THEN
#  ifdef DISTRIBUTE
        DO it=1,NT(ng)
          IF (Lstflux(it,ng)) THEN
            CALL ad_mp_exchange3d (ng, tile, iADM, 1,                   &
     &                             LBi, UBi, LBj, UBj, 1, Nfrec(ng),    &
     &                             NghostPoints,                        &
     &                             EWperiodic(ng), NSperiodic(ng),      &
     &                             ad_tflux(:,:,:,Linp,it))
          END IF
        END DO
#  endif
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            fac=1.0_r8/SQRT(om_r(i,j)*on_r(i,j))
            DO it=1,NT(ng)
              IF (Lstflux(it,ng)) THEN
                DO k=1,Nfrec(ng)
                  ad_tflux(i,j,k,Linp,it)=fac*ad_tflux(i,j,k,Linp,it)
                END DO
              END IF
            END DO
          END DO
        END DO
      END IF
# endif

      RETURN
      END SUBROUTINE ad_convolution_tile
#endif
      END MODULE ad_convolution_mod
