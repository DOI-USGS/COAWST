#include "cppdefs.h"
      MODULE ad_convolution_mod

#if defined ADJOINT && defined FOUR_DVAR
# ifdef EW_PERIODIC
#  define IR_RANGE Istr-1,Iend+1
#  define IU_RANGE Istr,Iend
#  define IV_RANGE Istr,Iend
# else
#  define IR_RANGE IstrR,IendR
#  define IU_RANGE Istr,IendR
#  define IV_RANGE IstrR,IendR
# endif
# ifdef NS_PERIODIC
#  define JR_RANGE Jstr-1,Jend+1
#  define JU_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JR_RANGE JstrR,JendR
#  define JU_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
# endif
!
!svn $Id: ad_convolution.F 694 2008-08-08 18:33:05Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
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
!  smoothing action on the  correlation operator  and they are use!    !
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
      SUBROUTINE ad_convolution (ng, tile, Linp, ifac)
!***********************************************************************
!
      USE mod_param
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
# endif
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Linp, ifac
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL ad_convolution_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nstp(ng), nnew(ng), Linp, ifac,         &
     &                          GRID(ng) % pm,                          &
     &                          GRID(ng) % om_r,                        &
     &                          GRID(ng) % om_u,                        &
     &                          GRID(ng) % om_v,                        &
     &                          GRID(ng) % pn,                          &
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
     &                          GRID(ng) % bed_thick,                   &
#  endif
     &                          GRID(ng) % Hz,                          &
     &                          GRID(ng) % z_r,                         &
     &                          GRID(ng) % z_w,                         &
# endif
     &                          MIXING(ng) % Kh,                        &
# ifdef SOLVE3D
     &                          MIXING(ng) % Kv,                        &
     &                          OCEAN(ng) % b_t,                        &
     &                          OCEAN(ng) % b_u,                        &
     &                          OCEAN(ng) % b_v,                        &
#  ifdef ADJUST_STFLUX
     &                          FORCES(ng) % b_stflx,                   &
#  endif
# endif
# ifdef ADJUST_WSTRESS
     &                          FORCES(ng) % b_sustr,                   &
     &                          FORCES(ng) % b_svstr,                   &
# endif
     &                          OCEAN(ng) % b_zeta,                     &
     &                          OCEAN(ng) % b_ubar,                     &
     &                          OCEAN(ng) % b_vbar,                     &
# ifdef SOLVE3D
     &                          OCEAN(ng) % ad_t,                       &
     &                          OCEAN(ng) % ad_u,                       &
     &                          OCEAN(ng) % ad_v,                       &
#  ifdef ADJUST_STFLUX
     &                          FORCES(ng) % ad_tflux,                  &
#  endif
# endif
# ifdef ADJUST_WSTRESS
     &                          FORCES(ng) % ad_ustr,                   &
     &                          FORCES(ng) % ad_vstr,                   &
# endif
     &                          OCEAN(ng) % ad_ubar,                    &
     &                          OCEAN(ng) % ad_vbar,                    &
     &                          OCEAN(ng) % ad_zeta)
      RETURN
      END SUBROUTINE ad_convolution
!
!***********************************************************************
      SUBROUTINE ad_convolution_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                nstp, nnew, Linp, ifac,           &
     &                                pm, om_r, om_u, om_v,             &
     &                                pn, on_r, on_u, on_v,             &
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
     &                                VnormR, VnormU, VnormV,           &
#  ifdef ADJUST_STFLUX
     &                                HnormSTF,                         &
#  endif
# endif
# ifdef ADJUST_WSTRESS
     &                                HnormSUS, HnormSVS,               &
# endif
     &                                HnormR, HnormU, HnormV,           &
# ifdef SOLVE3D
     &                                ad_t, ad_u, ad_v,                 &
#  ifdef ADJUST_STFLUX
     &                                ad_tflux,                         &
#  endif
# endif
# ifdef ADJUST_WSTRESS
     &                                ad_ustr, ad_vstr,                 &
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
# ifdef DISTRIBUTE
      USE mp_exchange_mod
# endif
      USE set_depth_mod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew, Linp, ifac
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: om_r(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
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
#  ifdef SOLVE3D
      real(r8), intent(in) :: VnormR(LBi:,LBj:,:,:)
      real(r8), intent(in) :: VnormU(LBi:,LBj:,:)
      real(r8), intent(in) :: VnormV(LBi:,LBj:,:)
#   ifdef ADJUST_STFLUX
      real(r8), intent(in) :: HnormSTF(LBi:,LBj:,:)
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: HnormSUS(LBi:,LBj:)
      real(r8), intent(in) :: HnormSVS(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: HnormR(LBi:,LBj:)
      real(r8), intent(in) :: HnormU(LBi:,LBj:)
      real(r8), intent(in) :: HnormV(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
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
      real(r8), intent(in) :: om_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
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
#  ifdef SOLVE3D
      real(r8), intent(in) :: VnormR(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(in) :: VnormU(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: VnormV(LBi:UBi,LBj:UBj,N(ng))
#   ifdef ADJUST_STFLUX
      real(r8), intent(in) :: HnormSTF(LBi:UBi,LBj:UBj,NT(ng))
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: HnormSUS(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: HnormSVS(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: HnormR(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: HnormU(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: HnormV(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
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
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      integer :: i, is, j, k
# ifdef SOLVE3D
      integer :: itrc
# endif
      real(r8) :: cff
# ifdef SOLVE3D
      real(r8) :: fac
# endif
# ifdef SOLVE3D
      real(r8), dimension(LBi:UBi,LBj:UBj) :: work
# endif
!
# include "set_bounds.h"

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
!  Adjoint free-surface.
!
      DO j=JR_RANGE
        DO i=IR_RANGE
          ad_zeta(i,j,Linp)=ad_zeta(i,j,Linp)*HnormR(i,j)
        END DO
      END DO
!
!  Adjoint 2D momentum.
!
      DO j=JU_RANGE
        DO i=IU_RANGE
          ad_ubar(i,j,Linp)=ad_ubar(i,j,Linp)*HnormU(i,j)
        END DO
      END DO
      DO j=JV_RANGE
        DO i=IV_RANGE
          ad_vbar(i,j,Linp)=ad_vbar(i,j,Linp)*HnormV(i,j)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL ad_mp_exchange2d (ng, tile, iADM, 3,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_zeta(:,:,Linp),                         &
     &                       ad_ubar(:,:,Linp),                         &
     &                       ad_vbar(:,:,Linp))
# endif
# ifdef ADJUST_WSTRESS
!
!  Adjoint surface momentum stress.
!
      DO k=1,Nfrec(ng)
        DO j=JU_RANGE
          DO i=IU_RANGE
            ad_ustr(i,j,k,Linp)=ad_ustr(i,j,k,Linp)*HnormSUS(i,j)
          END DO
        END DO
        DO j=JV_RANGE
          DO i=IV_RANGE
            ad_vstr(i,j,k,Linp)=ad_vstr(i,j,k,Linp)*HnormSVS(i,j)
          END DO
        END DO
      END DO
#  ifdef DISTRIBUTE
      CALL ad_mp_exchange3d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj, 1, Nfrec(ng),          &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_ustr(:,:,:,Linp),                       &
     &                       ad_vstr(:,:,:,Linp))
#  endif
# endif
# ifdef SOLVE3D
!
!  Adjoint 3D momentum.
!
      DO k=1,N(ng)
        DO j=JU_RANGE
          DO i=IU_RANGE
            ad_u(i,j,k,Linp)=ad_u(i,j,k,Linp)*VnormU(i,j,k)
          END DO
        END DO
        DO j=JV_RANGE
          DO i=IV_RANGE
            ad_v(i,j,k,Linp)=ad_v(i,j,k,Linp)*VnormV(i,j,k)
          END DO
        END DO
      END DO
#  ifdef DISTRIBUTE
      CALL ad_mp_exchange3d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_u(:,:,:,Linp),                          &
     &                       ad_v(:,:,:,Linp))
#  endif
!
!  Adjoint tracers.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JR_RANGE
            DO i=IR_RANGE
              ad_t(i,j,k,Linp,itrc)=ad_t(i,j,k,Linp,itrc)*              &
     &                              VnormR(i,j,k,itrc)
            END DO
          END DO
        END DO
      END DO
#  ifdef DISTRIBUTE
      CALL ad_mp_exchange4d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),   &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_t(:,:,:,Linp,:))
#  endif
#  ifdef ADJUST_STFLUX
!
!  Adjoint surface tracers flux.
!
      DO itrc=1,NT(ng)
        DO k=1,Nfrec(ng)
          DO j=JR_RANGE
            DO i=IR_RANGE
              ad_tflux(i,j,k,Linp,itrc)=ad_tflux(i,j,k,Linp,itrc)*      &
     &                                  HnormSTF(i,j,itrc)
            END DO
          END DO
        END DO
      END DO
#   ifdef DISTRIBUTE
      CALL ad_mp_exchange4d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj, 1,Nfrec(ng), 1,NT(ng), &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_tflux(:,:,:,Linp,:))
#   endif
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Convolve adjoint state vector with a generalized adjoint diffusion
!  equation to filter solution with specified horizontal scales.
!  Convert from model space to minimization space (v-space).
!-----------------------------------------------------------------------
!
!  Adjoint free-surface.
!
      CALL ad_conv_r2d_tile (ng, tile, iADM,                            &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       NghostPoints,                              &
     &                       NHsteps(isFsur)/ifac,                      &
     &                       DTsizeH(isFsur),                           &
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
     &                       NHsteps(isUbar)/ifac,                      &
     &                       DTsizeH(isUbar),                           &
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
     &                       NHsteps(isVbar)/ifac,                      &
     &                       DTsizeH(isVbar),                           &
     &                       Kh,                                        &
     &                       pm, pn, pmon_p, pnom_r,                    &
# ifdef MASKING
     &                       vmask, pmask,                              &
# endif
     &                       ad_vbar(:,:,Linp))
# ifdef ADJUST_WSTRESS
!
!  Adjoint surface momentum stress.
!
      DO k=1,Nfrec(ng)
        CALL ad_conv_u2d_tile (ng, tile, iADM,                          &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         NghostPoints,                            &
     &                         NHsteps(isUstr)/ifac,                    &
     &                         DTsizeH(isUstr),                         &
     &                         Kh,                                      &
     &                         pm, pn, pmon_r, pnom_p,                  &
#  ifdef MASKING
     &                         umask, pmask,                            &
#  endif
     &                         ad_ustr(:,:,k,Linp))

        CALL ad_conv_v2d_tile (ng, tile, iADM,                          &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         NghostPoints,                            &
     &                         NHsteps(isVstr)/ifac,                    &
     &                         DTsizeH(isVstr),                         &
     &                         Kh,                                      &
     &                         pm, pn, pmon_p, pnom_r,                  &
#  ifdef MASKING
     &                         vmask, pmask,                            &
#  endif
     &                         ad_vstr(:,:,k,Linp))
      END DO
# endif
# ifdef SOLVE3D
!
!  Adjoint 3D momentum.
!
      CALL ad_conv_u3d_tile (ng, tile, iADM,                            &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       NghostPoints,                              &
     &                       NHsteps(isUvel)/ifac,                      &
     &                       NVsteps(isUvel)/ifac,                      &
     &                       DTsizeH(isUvel),                           &
     &                       DTsizeV(isUvel),                           &
     &                       Kh, Kv,                                    &
     &                       pm, pn, pmon_r, pnom_p,                    &
#  ifdef MASKING
     &                       umask, pmask,                              &
#  endif
     &                       Hz, z_r,                                   &
     &                       ad_u(:,:,:,Linp))

      CALL ad_conv_v3d_tile (ng, tile, iADM,                            &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       NghostPoints,                              &
     &                       NHsteps(isUvel)/ifac,                      &
     &                       NVsteps(isUvel)/ifac,                      &
     &                       DTsizeH(isUvel),                           &
     &                       DTsizeV(isUvel),                           &
     &                       Kh, Kv,                                    &
     &                       pm, pn, pmon_p, pnom_r,                    &
#  ifdef MASKING
     &                       vmask, pmask,                              &
#  endif
     &                       Hz, z_r,                                   &
     &                       ad_v(:,:,:,Linp))
!
!  Adjoint tracers.
!
      DO itrc=1,NT(ng)
        is=isTvar(itrc)
        CALL ad_conv_r3d_tile (ng, tile, iADM,                          &
     &                         LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         NghostPoints,                            &
     &                         NHsteps(is)/ifac,                        &
     &                         NVsteps(is)/ifac,                        &
     &                         DTsizeH(is),                             &
     &                         DTsizeV(is),                             &
     &                         Kh, Kv,                                  &
     &                         pm, pn, pmon_u, pnom_v,                  &
#  ifdef MASKING
     &                         rmask, umask, vmask,                     &
#  endif
     &                         Hz, z_r,                                 &
     &                         ad_t(:,:,:,Linp,itrc))
      END DO
#  ifdef ADJUST_STFLUX
!
!  Adjoint surface tracers flux.
!
      DO itrc=1,NT(ng)
        is=isTsur(itrc)
        DO k=1,Nfrec(ng)
          CALL ad_conv_r2d_tile (ng, tile, iADM,                        &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           NghostPoints,                          &
     &                           NHsteps(is)/ifac,                      &
     &                           DTsizeH(is),                           &
     &                           Kh,                                    &
     &                           pm, pn, pmon_u, pnom_v,                &
#   ifdef MASKING
     &                           rmask, umask, vmask,                   &
#   endif
     &                           ad_tflux(:,:,k,Linp,itrc))
        END DO
      END DO
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
      DO j=JR_RANGE
        DO i=IR_RANGE
          ad_zeta(i,j,Linp)=ad_zeta(i,j,Linp)/                          &
     &                      SQRT(om_r(i,j)*on_r(i,j))
        END DO
      END DO
!
!  Adjoint 2D momentum.
!
      DO j=JU_RANGE
        DO i=IU_RANGE
          ad_ubar(i,j,Linp)=ad_ubar(i,j,Linp)/                          &
     &                      SQRT(om_u(i,j)*on_u(i,j))
        END DO
      END DO
      DO j=JV_RANGE
        DO i=IV_RANGE
          ad_vbar(i,j,Linp)=ad_vbar(i,j,Linp)/                          &
     &                      SQRT(om_v(i,j)*on_v(i,j))
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL ad_mp_exchange2d (ng, tile, iADM, 3,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_zeta(:,:,Linp),                         &
     &                       ad_ubar(:,:,Linp),                         &
     &                       ad_vbar(:,:,Linp))
# endif
# ifdef ADJUST_WSTRESS
!
!  Adjoint surface momentum stress.
!
      DO k=1,Nfrec(ng)
        DO j=JU_RANGE
          DO i=IU_RANGE
            ad_ustr(i,j,k,Linp)=ad_ustr(i,j,k,Linp)/                    &
     &                          SQRT(om_u(i,j)*on_u(i,j))
          END DO
        END DO
        DO j=JV_RANGE
          DO i=IV_RANGE
            ad_vstr(i,j,k,Linp)=ad_vstr(i,j,k,Linp)/                    &
     &                          SQRT(om_v(i,j)*on_v(i,j))
          END DO
        END DO
      END DO
#  ifdef DISTRIBUTE
      CALL ad_mp_exchange3d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj, 1, Nfrec(ng),          &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_ustr(:,:,:,Linp),                       &
     &                       ad_vstr(:,:,:,Linp))
#  endif
# endif
# ifdef SOLVE3D
!
!  Adjoint 3D momentum.
!
      DO j=JU_RANGE
        DO i=IU_RANGE
          cff=om_u(i,j)*on_u(i,j)*0.5_r8
          DO k=1,N(ng)
            ad_u(i,j,k,Linp)=ad_u(i,j,k,Linp)/                          &
     &                       SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
          END DO
        END DO
      END DO
      DO j=JV_RANGE
        DO i=IV_RANGE
          cff=om_v(i,j)*on_v(i,j)*0.5_r8
          DO k=1,N(ng)
            ad_v(i,j,k,Linp)=ad_v(i,j,k,Linp)/                          &
     &                       SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
          END DO
        END DO
      END DO
#  ifdef DISTRIBUTE
      CALL ad_mp_exchange3d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_u(:,:,:,Linp),                          &
     &                       ad_v(:,:,:,Linp))
#  endif
!
!  Adjoint tracers.
!
      DO j=JR_RANGE
        DO i=IR_RANGE
          cff=om_r(i,j)*on_r(i,j)
          DO k=1,N(ng)
            fac=1.0_r8/SQRT(cff*Hz(i,j,k))
            DO itrc=1,NT(ng)
              ad_t(i,j,k,Linp,itrc)=fac*ad_t(i,j,k,Linp,itrc)
            END DO
          END DO
        END DO
      END DO
#  ifdef DISTRIBUTE
      CALL ad_mp_exchange4d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),   &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_t(:,:,:,Linp,:))
#  endif
#  ifdef ADJUST_STFLUX
!
!  Adjoint surface tracers flux.
!
      DO j=JR_RANGE
        DO i=IR_RANGE
          fac=1.0_r8/SQRT(om_r(i,j)*on_r(i,j))
          DO itrc=1,NT(ng)
            DO k=1,Nfrec(ng)
              ad_tflux(i,j,k,Linp,itrc)=fac*ad_tflux(i,j,k,Linp,itrc)
            END DO
          END DO
        END DO
      END DO
#   ifdef DISTRIBUTE
      CALL ad_mp_exchange4d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj, 1,Nfrec(ng), 1,NT(ng), &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_tflux(:,:,:,Linp,:))
#   endif
#  endif
# endif
# undef IR_RANGE
# undef IU_RANGE
# undef IV_RANGE
# undef JR_RANGE
# undef JU_RANGE
# undef JV_RANGE

      RETURN
      END SUBROUTINE ad_convolution_tile
#endif
      END MODULE ad_convolution_mod
