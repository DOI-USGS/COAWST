#include "cppdefs.h"
      MODULE tl_variability_mod

#if defined TANGENT && defined FOUR_DVAR
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
!svn $Id: tl_variability.F 694 2008-08-08 18:33:05Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!======================================================================= 
!                                                                      !
!  This routine converts tangent error correlations, C, to background  !
!  error covariances,  B,  by multiplying by the standard deviations,  !
!  S.                                                                  !
!                                                                      !
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
      PUBLIC :: tl_variability

      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_variability (ng, tile, Linp, Lweak)
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
      CALL tl_variability_tile (ng, tile,                               &
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
     &                          OCEAN(ng) % tl_t,                       &
     &                          OCEAN(ng) % tl_u,                       &
     &                          OCEAN(ng) % tl_v,                       &
#  ifdef ADJUST_STFLUX
     &                          FORCES(ng) % tl_tflux,                  &
#  endif
# endif
# ifdef ADJUST_WSTRESS
     &                          FORCES(ng) % tl_ustr,                   &
     &                          FORCES(ng) % tl_vstr,                   &
# endif
     &                          OCEAN(ng) % tl_ubar,                    &
     &                          OCEAN(ng) % tl_vbar,                    &
     &                          OCEAN(ng) % tl_zeta)

      RETURN
      END SUBROUTINE tl_variability
!
!***********************************************************************
      SUBROUTINE tl_variability_tile (ng, tile,                         &
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
     &                                tl_t, tl_u, tl_v,                 &
#  ifdef ADJUST_STFLUX
     &                                tl_tflux,                         &
#  endif
# endif
# ifdef ADJUST_WSTRESS
     &                                tl_ustr, tl_vstr,               &
# endif
     &                                tl_ubar, tl_vbar, tl_zeta)
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
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
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
      real(r8), intent(in) :: sustr_std(LBi:,LBj:)
      real(r8), intent(in) :: svstr_std(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: ubar_std(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vbar_std(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: zeta_std(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,2,N(ng))
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,2,N(ng))
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
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
!  Multiply tangent linear state by its corresponding background-error
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
!  Tangent linear free-surface.
!
      DO j=JR_RANGE
        DO i=IR_RANGE
          tl_zeta(i,j,Linp)=tl_zeta(i,j,Linp)*zeta_std(i,j)*            &
     &                      fac(isFsur)
        END DO
      END DO
!
!  Tangent linear 2D momentum.
!
      DO j=JU_RANGE
        DO i=IU_RANGE
          tl_ubar(i,j,Linp)=tl_ubar(i,j,Linp)*ubar_std(i,j)*            &
     &                      fac(isUbar)
        END DO                                                          
      END DO                                                            
      DO j=JV_RANGE                                                     
        DO i=IV_RANGE                                                   
          tl_vbar(i,j,Linp)=tl_vbar(i,j,Linp)*vbar_std(i,j)*            &
     &                      fac(isVbar)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iTLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_zeta(:,:,Linp),                            &
     &                    tl_ubar(:,:,Linp),                            &
     &                    tl_vbar(:,:,Linp))
# endif
# ifdef ADJUST_WSTRESS
!
!  Tangent linear surface momentum stress.
!
      DO k=1,Nfrec(ng)
        DO j=JU_RANGE
          DO i=IU_RANGE
            tl_ustr(i,j,k,Linp)=tl_ustr(i,j,k,Linp)*sustr_std(i,j)*     &
     &                          fac(isUstr)
          END DO
        END DO
        DO j=JV_RANGE
          DO i=IV_RANGE
            tl_vstr(i,j,k,Linp)=tl_vstr(i,j,k,Linp)*svstr_std(i,j)*     &
     &                          fac(isVstr)
          END DO
        END DO
      END DO
#  ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iTLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, Nfrec(ng),             &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_ustr(:,:,:,Linp),                          &
     &                    tl_vstr(:,:,:,Linp))
#  endif
# endif
# ifdef SOLVE3D
!
!  Tangent linear 3D momentum.
!
      DO k=1,N(ng)
        DO j=JU_RANGE
          DO i=IU_RANGE
            tl_u(i,j,k,Linp)=tl_u(i,j,k,Linp)*u_std(i,j,k)*             &
     &                       fac(isUvel)
          END DO
        END DO
        DO j=JV_RANGE
          DO i=IV_RANGE
            tl_v(i,j,k,Linp)=tl_v(i,j,k,Linp)*v_std(i,j,k)*             &
     &                       fac(isVvel)
          END DO
        END DO
      END DO
#  ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iTLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_u(:,:,:,Linp),                             &
     &                    tl_v(:,:,:,Linp))
#  endif
!
!  Tangent linear tracers.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JR_RANGE
            DO i=IR_RANGE
              tl_t(i,j,k,Linp,itrc)=tl_t(i,j,k,Linp,itrc)*              &
     &                              t_std(i,j,k,itrc)*                  &
     &                              fac(isTvar(itrc))
            END DO
          END DO
        END DO
      END DO
#  ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_t(:,:,:,Linp,:))
#  endif
#  ifdef ADJUST_STFLUX
!
!  Tangent linear surface tracers flux.
!
      DO itrc=1,NT(ng)
        DO k=1,Nfrec(ng)
          DO j=JR_RANGE
            DO i=IR_RANGE
              tl_tflux(i,j,k,Linp,itrc)=tl_tflux(i,j,k,Linp,itrc)*      &
     &                                  stflx_std(i,j,itrc)*            &
     &                                  fac(isTsur(itrc))
            END DO
          END DO
        END DO
      END DO
#   ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, Nfrec(ng), 1, NT(ng),  &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_tflux(:,:,:,Linp,:))
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
      END SUBROUTINE tl_variability_tile
#endif
      END MODULE tl_variability_mod

