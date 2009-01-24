#include "cppdefs.h"
      MODULE ad_variability_mod

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
!svn $Id: ad_variability.F 694 2008-08-08 18:33:05Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine converts adjoint error correlations, C, to background  !
!  error covariances,  B,  by multiplying by the standard deviations,  !
!  S.                                                                  !
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
!    G : normalization coefficients matrix (convert to correlations)   !
!    L : tangent linear and adjoint diffusion operators                !
!    W : diagonal matrix of local area or volume metrics               !
!                                                                      !
!  Here, T/2 denote the transpose of a squared-root factor.            !
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
      PUBLIC :: ad_variability

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_variability (ng, tile, Linp, Lweak)
!***********************************************************************
!
      USE mod_param
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
# endif
      USE mod_ocean
!
!  Imported variable declarations.
!
      logical, intent(in) :: Lweak

      integer, intent(in) :: ng, tile, Linp
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL ad_variability_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Linp, Lweak,                            &
# ifdef SOLVE3D
     &                          OCEAN(ng) % e_t,                        &
     &                          OCEAN(ng) % e_u,                        &
     &                          OCEAN(ng) % e_v,                        &
#  ifdef ADJUST_STFLUX
     &                          FORCES(ng) % e_stflx,                   &
#  endif
# endif
# ifdef ADJUST_WSTRESS
     &                          FORCES(ng) % e_sustr,                   &
     &                          FORCES(ng) % e_svstr,                   &
# endif
     &                          OCEAN(ng) % e_ubar,                     &
     &                          OCEAN(ng) % e_vbar,                     &
     &                          OCEAN(ng) % e_zeta,                     &
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
      END SUBROUTINE ad_variability
!
!***********************************************************************
      SUBROUTINE ad_variability_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                Linp, Lweak,                      &
# ifdef SOLVE3D
     &                                t_std, u_std, v_std,              &
#  ifdef ADJUST_STFLUX
     &                                stflx_std,                        &
#  endif
# endif
# ifdef ADJUST_WSTRESS
     &                                sustr_std, svstr_std,             &
# endif
     &                                ubar_std, vbar_std, zeta_std,     &
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
      USE mod_ncparam
      USE mod_scalars
!
# ifdef DISTRIBUTE
      USE mp_exchange_mod
# endif
!
!  Imported variable declarations.
!
      logical, intent(in) :: Lweak

      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Linp
!
# ifdef ASSUMED_SHAPE
#  ifdef SOLVE3D
      real(r8), intent(in) :: t_std(LBi:,LBj:,:,:)
      real(r8), intent(in) :: u_std(LBi:,LBj:,:)
      real(r8), intent(in) :: v_std(LBi:,LBj:,:)
#   ifdef ADJUST_STFLUX
      real(r8), intent(in) :: stflx_std(LBi:,LBj:,:)
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: sustr_std(LBi:,LBj:)
      real(r8), intent(in) :: svstr_std(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: ubar_std(LBi:,LBj:)
      real(r8), intent(in) :: vbar_std(LBi:,LBj:)
      real(r8), intent(in) :: zeta_std(LBi:,LBj:)
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
# else
#  ifdef SOLVE3D
      real(r8), intent(in) :: t_std(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(in) :: u_std(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: v_std(LBi:UBi,LBj:UBj,N(ng))
#   ifdef ADJUST_STFLUX
      real(r8), intent(in) :: stflx_std(LBi:UBi,LBj:UBj,NT(ng))
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: sustr_std(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr_std(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: ubar_std(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vbar_std(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: zeta_std(LBi:UBi,LBj:UBj)
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
      integer :: i, j, k
# ifdef SOLVE3D
      integer :: itrc
# endif
      real(r8), dimension(MstateVar) :: fac
!
# include "set_bounds.h"

!
!-----------------------------------------------------------------------
!  Multiply adjoint state by its corresponding background-error
!  standard deviation.
!-----------------------------------------------------------------------
!
      IF (Lweak) THEN
        DO i=1,MstateVar
          fac(i)=SQRT(Cfscale(i,ng))
        END DO
      ELSE
        DO i=1,MstateVar
          fac(i)=1.0_r8
        END DO
      END IF
!
!  Adjoint free-surface.
!
      DO j=JR_RANGE
        DO i=IR_RANGE
          ad_zeta(i,j,Linp)=ad_zeta(i,j,Linp)*zeta_std(i,j)*            &
     &                      fac(isFsur)
        END DO
      END DO
!
!  Adjoint 2D momentum.
!
      DO j=JU_RANGE
        DO i=IU_RANGE
          ad_ubar(i,j,Linp)=ad_ubar(i,j,Linp)*ubar_std(i,j)*            &
     &                      fac(isUbar)
        END DO
      END DO
      DO j=JV_RANGE
        DO i=IV_RANGE
          ad_vbar(i,j,Linp)=ad_vbar(i,j,Linp)*vbar_std(i,j)*            &
     &                      fac(isVbar)
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
            ad_ustr(i,j,k,Linp)=ad_ustr(i,j,k,Linp)*sustr_std(i,j)*     &
     &                          fac(isUstr)
          END DO
        END DO
        DO j=JV_RANGE
          DO i=IV_RANGE
            ad_vstr(i,j,k,Linp)=ad_vstr(i,j,k,Linp)*svstr_std(i,j)*     &
     &                          fac(isVstr)
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
            ad_u(i,j,k,Linp)=ad_u(i,j,k,Linp)*u_std(i,j,k)*             &
     &                       fac(isUvel)
          END DO
        END DO
        DO j=JV_RANGE
          DO i=IV_RANGE
            ad_v(i,j,k,Linp)=ad_v(i,j,k,Linp)*v_std(i,j,k)*             &
     &                       fac(isVvel)
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
     &                              t_std(i,j,k,itrc)*                  &
     &                              fac(isTvar(itrc))
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
     &                                  stflx_std(i,j,itrc)*            &
     &                                  fac(isTsur(itrc))
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
      END SUBROUTINE ad_variability_tile
#endif
      END MODULE ad_variability_mod
