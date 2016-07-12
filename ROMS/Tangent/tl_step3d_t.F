#include "cppdefs.h"
      MODULE tl_step3d_t_mod
#if !defined TS_FIXED && (defined TANGENT && defined SOLVE3D)
!
!svn $Id: tl_step3d_t.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine time-steps tangent linear tracer equations.  Notice    !
!  that advective and diffusive terms are time-stepped differently.    !
!  It applies the corrector time-step for horizontal and vertical      !
!  advection, vertical diffusion, nudging if necessary, and lateral    !
!  boundary conditions.                                                !
!                                                                      !
!  Notice that at input the tracer arrays have:                        !
!                                                                      !
!    t(:,:,:,nnew,:)   m Tunits  n+1     horizontal/vertical diffusion !
!                                        terms plus source/sink terms  !
!                                        (biology, sediment), if any   !
!                                                                      !
!    t(:,:,:,3   ,:)   Tunits    n+1/2   advective terms and vertical  !
!                                        diffusion predictor step      !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: tl_step3d_t
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_step3d_t (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_clima
# ifdef DIAGNOSTICS_TS
!!    USE mod_diags
# endif
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iTLM, 35)
# endif
      CALL tl_step3d_t_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nrhs(ng), nstp(ng), nnew(ng),              &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
# endif
# ifdef TS_MPDATA_NOT_YET
#  ifdef WET_DRY
     &                       GRID(ng) % rmask_wet,                      &
     &                       GRID(ng) % umask_wet,                      &
     &                       GRID(ng) % vmask_wet,                      &
#  endif
     &                       GRID(ng) % omn,                            &
     &                       GRID(ng) % om_u,                           &
     &                       GRID(ng) % om_v,                           &
     &                       GRID(ng) % on_u,                           &
     &                       GRID(ng) % on_v,                           &
# endif
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pn,                             &
     &                       GRID(ng) % Hz,                             &
     &                       GRID(ng) % tl_Hz,                          &
     &                       GRID(ng) % Huon,                           &
     &                       GRID(ng) % tl_Huon,                        &
     &                       GRID(ng) % Hvom,                           &
     &                       GRID(ng) % tl_Hvom,                        &
     &                       GRID(ng) % z_r,                            &
     &                       GRID(ng) % tl_z_r,                         &
     &                       MIXING(ng) % Akt,                          &
     &                       MIXING(ng) % tl_Akt,                       &
     &                       OCEAN(ng) % W,                             &
     &                       OCEAN(ng) % tl_W,                          &
# if defined FLOATS_NOT_YET && defined FLOAT_VWALK
     &                       MIXING(ng) % dAktdz,                       &
# endif
# ifdef DIAGNOSTICS_TS
!!   &                       DIAGS(ng) % DiaTwrk,                       &
# endif
     &                       OCEAN(ng) % t,                             &
     &                       OCEAN(ng) % tl_t)
# ifdef PROFILE
      CALL wclock_off (ng, iTLM, 35)
# endif

      RETURN
      END SUBROUTINE tl_step3d_t
!
!***********************************************************************
      SUBROUTINE tl_step3d_t_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nrhs, nstp, nnew,                    &
# ifdef MASKING
     &                             rmask, umask, vmask,                 &
# endif
# ifdef TS_MPDATA_NOT_YET
#  ifdef WET_DRY
     &                             rmask_wet, umask_wet, vmask_wet,     &
#  endif
     &                             omn, om_u, om_v, on_u, on_v,         &
# endif
     &                             pm, pn,                              &
     &                             Hz, tl_Hz,                           &
     &                             Huon, tl_Huon,                       &
     &                             Hvom, tl_Hvom,                       &
     &                             z_r, tl_z_r,                         &
     &                             Akt, tl_Akt,                         &
     &                             W, tl_W,                             &
# if defined FLOATS_NOT_YET && defined FLOAT_VWALK
     &                             dAktdz,                              &
# endif
# ifdef DIAGNOSTICS_TS
!!   &                             DiaTwrk,                             &
# endif
     &                             t, tl_t)
!***********************************************************************
!
      USE mod_param
      USE mod_clima
      USE mod_scalars
      USE mod_sources
!
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
# ifdef DISTRIBUTE
#  if defined FLOATS_NOT_YET && defined FLOAT_VWALK
      USE mp_exchange_mod, ONLY : mp_exchange3d
#  endif
      USE mp_exchange_mod, ONLY : mp_exchange4d
# endif
# ifdef TS_MPDATA_NOT_YET
!!    USE tl_mpdata_adiff_mod
# endif
      USE tl_t3dbc_mod, ONLY : tl_t3dbc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, nstp, nnew
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef TS_MPDATA_NOT_YET
#   ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
#  ifdef SUN
      real(r8), intent(in) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#  else
      real(r8), intent(in) :: Akt(LBi:,LBj:,0:,:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(in) :: W(LBi:,LBj:,0:)

      real(r8), intent(in) :: tl_Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_z_r(LBi:,LBj:,:)
#  ifdef SUN
      real(r8), intent(in) :: tl_Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
#  else
      real(r8), intent(in) :: tl_Akt(LBi:,LBj:,0:,:)
#  endif
      real(r8), intent(in) :: tl_W(LBi:,LBj:,0:)
#  ifdef DIAGNOSTICS_TS
!!    real(r8), intent(inout) :: DiaTwrk(LBi:,LBj:,:,:,:)
#  endif
#  ifdef SUN
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#  else
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
#  endif
#  if defined FLOATS_NOT_YET && defined FLOAT_VWALK
      real(r8), intent(out) :: dAktdz(LBi:,LBj:,:)
#  endif

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef TS_MPDATA_NOT_YET
#   ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: W(LBi:UBi,LBj:UBj,0:N(ng))

      real(r8), intent(in) :: tl_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_Hvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
      real(r8), intent(in) :: tl_W(LBi:UBi,LBj:UBj,0:N(ng))
#  ifdef DIAGNOSTICS_TS
!!    real(r8), intent(inout) :: DiaTwrk(LBi:UBi,LBj:UBj,N(ng),NT(ng),  &
!!   &                                   NDT)
#  endif
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#  if defined FLOATS_NOT_YET && defined FLOAT_VWALK
      real(r8), intent(out) :: dAktdz(LBi:UBi,LBj:UBj,N(ng))
#  endif
# endif
!
!  Local variable declarations.
!
      integer :: Isrc, Jsrc
      integer :: i, ic, is, itrc, j, k, ltrc
# if defined AGE_MEAN && defined T_PASSIVE
      integer :: iage
# endif
# ifdef DIAGNOSTICS_TS
      integer :: idiag
# endif
      real(r8), parameter :: eps = 1.0E-16_r8

      real(r8) :: cff, cff1, cff2, cff3
      real(r8) :: tl_cff, tl_cff1, tl_cff2, tl_cff3

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: BC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: tl_CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: tl_BC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: tl_DC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: tl_FC

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: curv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_curv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_grad

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: oHz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: tl_oHz
# ifdef TS_MPDATA_NOT_YET
#  ifdef DIAGNOSTICS_TS
!!    real(r8), dimension(IminS:ImaxS,JminS:JmaxS,3) :: Dhadv
!!    real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng),NT(ng)) :: Dvadv
#  endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,                      &
     &                    N(ng),NT(ng)) :: Ta
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,                      &
     &                    N(ng),NT(ng)) :: tl_Ta

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: Ua
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: Va
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: Wa

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: tl_Ua
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: tl_Va
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: tl_Wa
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Time-step horizontal advection term.
!-----------------------------------------------------------------------
!
!  Compute inverse thickness.
!
# ifdef TS_MPDATA_NOT_YET
#  define I_RANGE Istrm2,Iendp2
#  define J_RANGE Jstrm2,Jendp2
# else
#  define I_RANGE Istr,Iend
#  define J_RANGE Jstr,Jend
# endif

      DO k=1,N(ng)
        DO j=J_RANGE
          DO i=I_RANGE
            oHz(i,j,k)=1.0_r8/Hz(i,j,k)
            tl_oHz(i,j,k)=-oHz(i,j,k)*oHz(i,j,k)*tl_Hz(i,j,k)
          END DO
        END DO
      END DO
# undef I_RANGE
# undef J_RANGE
# ifdef TS_MPDATA_NOT_YET
!
!  The MPDATA algorithm requires a three-point footprint, so exchange
!  boundary data on t(:,:,:,nnew,:) so other processes computed earlier
!  (horizontal diffusion, biology, or sediment) are accounted.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        DO itrc=1,NT(ng)
!>        CALL exchange_r3d_tile (ng, tile,                             &
!>   &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
!>   &                            t(:,:,:,nnew,itrc))
!>
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            tl_t(:,:,:,nnew,itrc))
        END DO
      END IF

#  ifdef DISTRIBUTE
!>    CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    t(:,:,:,nnew,:))
!>
      CALL mp_exchange4d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_t(:,:,:,nnew,:))
#  endif
# endif
!
!  Compute tangent linear horizontal tracer advection fluxes.
!
      T_LOOP : DO itrc=1,NT(ng)
        K_LOOP : DO k=1,N(ng)

# if defined TS_C2HADVECTION_TL
!
!  Second-order, centered differences horizontal advective fluxes.
!
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
!>            FX(i,j)=Huon(i,j,k)*                                      &
!>   &                0.5_r8*(t(i-1,j,k,3,itrc)+                        &
!>   &                        t(i  ,j,k,3,itrc))
!>
              tl_FX(i,j)=0.5_r8*                                        &
     &                   (tl_Huon(i,j,k)*(t(i-1,j,k,3,itrc)+            &
     &                                    t(i  ,j,k,3,itrc))+           &
     &                    Huon(i,j,k)*(tl_t(i-1,j,k,3,itrc)+            &
     &                                 tl_t(i  ,j,k,3,itrc)))
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
!>            FE(i,j)=Hvom(i,j,k)*                                      &
!>   &                0.5_r8*(t(i,j-1,k,3,itrc)+                        &
!>   &                        t(i,j  ,k,3,itrc))
!>
              tl_FE(i,j)=0.5_r8*                                        &
     &                   (tl_Hvom(i,j,k)*(t(i,j-1,k,3,itrc)+            &
     &                                    t(i,j  ,k,3,itrc))+           &
     &                    Hvom(i,j,k)*(tl_t(i,j-1,k,3,itrc)+            &
     &                                 tl_t(i,j  ,k,3,itrc)))
            END DO
          END DO

# elif defined TS_MPDATA_NOT_YET
!
!  First-order, upstream differences horizontal advective fluxes.
!
          DO j=JstrVm2,Jendp2i
            DO i=IstrUm2,Iendp3
              cff1=MAX(Huon(i,j,k),0.0_r8)
              cff2=MIN(Huon(i,j,k),0.0_r8)
              tl_cff1=(0.5_r8+SIGN(0.5_r8, Huon(i,j,k)))*tl_Huon(i,j,k)
              tl_cff2=(0.5_r8+SIGN(0.5_r8,-Huon(i,j,k)))*tl_Huon(i,j,k)
!>            FX(i,j)=cff1*t(i-1,j,k,3,itrc)+                           &
!>   &                cff2*t(i  ,j,k,3,itrc)
!>
              tl_FX(i,j)=tl_cff1*t(i-1,j,k,3,itrc)+                     &
     &                   cff1*tl_t(i-1,j,k,3,itrc)+                     &
     &                   tl_cff2*t(i  ,j,k,3,itrc)+                     &
     &                   cff2*tl_t(i  ,j,k,3,itrc)
            END DO
          END DO
          DO j=JstrVm2,Jendp3
            DO i=IstrUm2,Iendp2i
              cff1=MAX(Hvom(i,j,k),0.0_r8)
              cff2=MIN(Hvom(i,j,k),0.0_r8)
              tl_cff1=(0.5_r8+SIGN(0.5_r8, Hvom(i,j,k)))*tl_Hvom(i,j,k)
              tl_cff2=(0.5_r8+SIGN(0.5_r8,-Hvom(i,j,k)))*tl_Hvom(i,j,k)
!>            FE(i,j)=cff1*t(i,j-1,k,3,itrc)+                           &
!>   &                cff2*t(i,j  ,k,3,itrc)
!>
              tl_FE(i,j)=tl_cff1*t(i,j-1,k,3,itrc)+                     &
     &                   cff1*tl_t(i,j-1,k,3,itrc)+                     &
     &                   tl_cff2*t(i,j  ,k,3,itrc)+                     &
     &                   cff2*tl_t(i,j  ,k,3,itrc)
            END DO
          END DO

# else
!
#  if defined TS_U3HADVECTION_TL
!  Third-order, uptream-biased horizontal advective fluxes.
#  elif defined TS_A4HADVECTION_TL
!  Fourth-order, Akima horizontal advective fluxes.
#  else
!  Fourth-order, centered differences horizontal advective fluxes.
#  endif
!
          DO j=Jstr,Jend
            DO i=Istrm1,Iendp2
              FX(i,j)=t(i  ,j,k,3,itrc)-                                &
     &                t(i-1,j,k,3,itrc)
              tl_FX(i,j)=tl_t(i  ,j,k,3,itrc)-                          &
     &                   tl_t(i-1,j,k,3,itrc)
#  ifdef MASKING
              FX(i,j)=FX(i,j)*umask(i,j)
              tl_FX(i,j)=tl_FX(i,j)*umask(i,j)
#  endif
            END DO
          END DO
          IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Western_Edge(tile)) THEN
              DO j=Jstr,Jend
                FX(Istr-1,j)=FX(Istr,j)
                tl_FX(Istr-1,j)=tl_FX(Istr,j)
              END DO
            END IF
          END IF
          IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
              DO j=Jstr,Jend
                FX(Iend+2,j)=FX(Iend+1,j)
                tl_FX(Iend+2,j)=tl_FX(Iend+1,j)
              END DO
            END IF
          END IF
!
          DO j=Jstr,Jend
            DO i=Istr-1,Iend+1
#  if defined TS_U3HADVECTION_TL
              curv(i,j)=FX(i+1,j)-FX(i,j)
              tl_curv(i,j)=tl_FX(i+1,j)-tl_FX(i,j)
#  elif defined TS_A4HADVECTION_TL
              cff=2.0_r8*FX(i+1,j)*FX(i,j)
              tl_cff=2.0_r8*(tl_FX(i+1,j)*FX(i,j)+                      &
                             FX(i+1,j)*tl_FX(i,j))
              IF (cff.gt.eps) THEN
                grad(i,j)=cff/(FX(i+1,j)+FX(i,j))
                tl_grad(i,j)=((FX(i+1,j)+FX(i,j))*tl_cff-               &
     &                        cff*(tl_FX(i+1,j)+tl_FX(i,j)))/           &
     &                       ((FX(i+1,j)+FX(i,j))*(FX(i+1,j)+FX(i,j)))
              ELSE
                grad(i,j)=0.0_r8
                tl_grad(i,j)=0.0_r8
              END IF
#  else
              grad(i,j)=0.5_r8*(FX(i+1,j)+FX(i,j))
              tl_grad(i,j)=0.5_r8*(tl_FX(i+1,j)+tl_FX(i,j))
#  endif
            END DO
          END DO
!
          cff1=1.0_r8/6.0_r8
          cff2=1.0_r8/3.0_r8
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
#  ifdef TS_U3HADVECTION_TL
!>            FX(i,j)=Huon(i,j,k)*0.5_r8*                               &
!>   &                (t(i-1,j,k,3,itrc)+                               &
!>   &                 t(i  ,j,k,3,itrc))-                              &
!>   &                cff1*(curv(i-1,j)*MAX(Huon(i,j,k),0.0_r8)+        &
!>   &                      curv(i  ,j)*MIN(Huon(i,j,k),0.0_r8))
!>
              tl_FX(i,j)=0.5_r8*                                        &
     &                   (tl_Huon(i,j,k)*                               &
     &                    (t(i-1,j,k,3,itrc)+                           &
     &                     t(i  ,j,k,3,itrc))+                          &
     &                    Huon(i,j,k)*                                  &
     &                    (tl_t(i-1,j,k,3,itrc)+                        &
     &                     tl_t(i  ,j,k,3,itrc)))-                      &
     &                   cff1*                                          &
     &                   (tl_curv(i-1,j)*MAX(Huon(i,j,k),0.0_r8)+       &
     &                    curv(i-1,j)*                                  &
     &                    (0.5_r8+SIGN(0.5_r8, Huon(i,j,k)))*           &
     &                    tl_Huon(i,j,k)+                               &
     &                    tl_curv(i  ,j)*MIN(Huon(i,j,k),0.0_r8)+       &
     &                    curv(i  ,j)*                                  &
     &                    (0.5_r8+SIGN(0.5_r8,-Huon(i,j,k)))*           &
     &                    tl_Huon(i,j,k))
#  else
!>            FX(i,j)=Huon(i,j,k)*0.5_r8*                               &
!>   &                (t(i-1,j,k,3,itrc)+                               &
!>   &                 t(i  ,j,k,3,itrc)-                               &
!>   &                 cff2*(grad(i  ,j)-                               &
!>   &                       grad(i-1,j)))
!>
              tl_FX(i,j)=0.5_r8*                                        &
     &                   (tl_Huon(i,j,k)*                               &
     &                    (t(i-1,j,k,3,itrc)+                           &
     &                     t(i  ,j,k,3,itrc)-                           &
     &                     cff2*(grad(i  ,j)-                           &
     &                           grad(i-1,j)))+                         &
     &                    Huon(i,j,k)*                                  &
     &                    (tl_t(i-1,j,k,3,itrc)+                        &
     &                     tl_t(i  ,j,k,3,itrc)-                        &
     &                     cff2*(tl_grad(i  ,j)-                        &
     &                           tl_grad(i-1,j))))
#  endif
            END DO
          END DO
!
          DO j=Jstrm1,Jendp2
            DO i=Istr,Iend
              FE(i,j)=t(i,j  ,k,3,itrc)-                                &
     &                t(i,j-1,k,3,itrc)
              tl_FE(i,j)=tl_t(i,j  ,k,3,itrc)-                          &
     &                   tl_t(i,j-1,k,3,itrc)
#  ifdef MASKING
              FE(i,j)=FE(i,j)*vmask(i,j)
              tl_FE(i,j)=tl_FE(i,j)*vmask(i,j)
#  endif
            END DO
          END DO
          IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
            IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
              DO i=Istr,Iend
                FE(i,Jstr-1)=FE(i,Jstr)
                tl_FE(i,Jstr-1)=tl_FE(i,Jstr)
              END DO
            END IF
          END IF
          IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
            IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
              DO i=Istr,Iend
                FE(i,Jend+2)=FE(i,Jend+1)
                tl_FE(i,Jend+2)=tl_FE(i,Jend+1)
              END DO
            END IF
          END IF
!
          DO j=Jstr-1,Jend+1
            DO i=Istr,Iend
#  if defined TS_U3HADVECTION_TL
              curv(i,j)=FE(i,j+1)-FE(i,j)
              tl_curv(i,j)=tl_FE(i,j+1)-tl_FE(i,j)
#  elif defined TS_A4HADVECTION_TL
              cff=2.0_r8*FE(i,j+1)*FE(i,j)
              tl_cff=2.0_r8*(tl_FE(i,j+1)*FE(i,j)+                      &
     &                       FE(i,j+1)*tl_FE(i,j))
              IF (cff.gt.eps) THEN
                grad(i,j)=cff/(FE(i,j+1)+FE(i,j))
                tl_grad(i,j)=((FE(i,j+1)+FE(i,j))*tl_cff-               &
     &                        cff*(tl_FE(i,j+1)+tl_FE(i,j)))/           &
     &                       ((FE(i,j+1)+FE(i,j))*(FE(i,j+1)+FE(i,j)))
              ELSE
                grad(i,j)=0.0_r8
                tl_grad(i,j)=0.0_r8
              END IF
#  else
              grad(i,j)=0.5_r8*(FE(i,j+1)+FE(i,j))
              tl_grad(i,j)=0.5_r8*(tl_FE(i,j+1)+tl_FE(i,j))
#  endif
            END DO
          END DO
!
          cff1=1.0_r8/6.0_r8
          cff2=1.0_r8/3.0_r8
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
#  ifdef TS_U3HADVECTION_TL
!>            FE(i,j)=Hvom(i,j,k)*0.5_r8*                               &
!>   &                (t(i,j-1,k,3,itrc)+                               &
!>   &                 t(i,j  ,k,3,itrc))-                              &
!>   &                cff1*(curv(i,j-1)*MAX(Hvom(i,j,k),0.0_r8)+        &
!>   &                      curv(i,j  )*MIN(Hvom(i,j,k),0.0_r8))
!>
              tl_FE(i,j)=0.5_r8*                                        &
     &                   (tl_Hvom(i,j,k)*                               &
     &                    (t(i,j-1,k,3,itrc)+                           &
     &                     t(i,j  ,k,3,itrc))+                          &
     &                    Hvom(i,j,k)*                                  &
     &                    (tl_t(i,j-1,k,3,itrc)+                        &
     &                     tl_t(i,j  ,k,3,itrc)))-                      &
     &                    cff1*                                         &
     &                    (tl_curv(i,j-1)*MAX(Hvom(i,j,k),0.0_r8)+      &
     &                     curv(i,j-1)*                                 &
     &                     (0.5_r8+SIGN(0.5_r8, Hvom(i,j,k)))*          &
     &                     tl_Hvom(i,j,k)+                              &
     &                     tl_curv(i,j  )*MIN(Hvom(i,j,k),0.0_r8)+      &
     &                     curv(i,j  )*                                 &
     &                     (0.5_r8+SIGN(0.5_r8,-Hvom(i,j,k)))*          &
     &                     tl_Hvom(i,j,k))
#  else
!>            FE(i,j)=Hvom(i,j,k)*0.5_r8*                               &
!>   &                (t(i,j-1,k,3,itrc)+                               &
!>   &                 t(i,j  ,k,3,itrc)-                               &
!>   &                 cff2*(grad(i,j  )-                               &
!>   &                       grad(i,j-1)))
!>
              tl_FE(i,j)=0.5_r8*                                        &
     &                   (tl_Hvom(i,j,k)*                               &
     &                    (t(i,j-1,k,3,itrc)+                           &
     &                     t(i,j  ,k,3,itrc)-                           &
     &                     cff2*(grad(i,j  )-                           &
     &                           grad(i,j-1)))+                         &
     &                    Hvom(i,j,k)*                                  &
     &                    (tl_t(i,j-1,k,3,itrc)+                        &
     &                     tl_t(i,j  ,k,3,itrc)-                        &
     &                     cff2*(tl_grad(i,j  )-                        &
     &                           tl_grad(i,j-1))))
#  endif
            END DO
          END DO
# endif
!
!  Apply tracers point sources to the horizontal advection terms,
!  if any.
!
          IF (LuvSrc(ng).and.ANY(LtracerSrc(:,ng))) THEN
            DO is=1,Nsrc(ng)
              Isrc=SOURCES(ng)%Isrc(is)
              Jsrc=SOURCES(ng)%Jsrc(is)
              IF (INT(SOURCES(ng)%Dsrc(is)).eq.0) THEN
# ifdef TS_MPDATA_NOT_YET
                IF (((IstrUm2.le.Isrc).and.(Isrc.le.Iendp3)).and.       &
     &              ((JstrVm2.le.Jsrc).and.(Jsrc.le.Jendp2i))) THEN
# else
                IF (((Istr.le.Isrc).and.(Isrc.le.Iend+1)).and.          &
     &              ((Jstr.le.Jsrc).and.(Jsrc.le.Jend))) THEN
# endif
                  IF (LtracerSrc(itrc,ng)) THEN
!>                  FX(Isrc,Jsrc)=Huon(Isrc,Jsrc,k)*                    &
!>   &                            SOURCES(ng)%Tsrc(is,k,itrc)
!>
                    tl_FX(Isrc,Jsrc)=tl_Huon(Isrc,Jsrc,k)*              &
     &                               SOURCES(ng)%Tsrc(is,k,itrc)
# ifdef MASKING
                  ELSE
                    IF ((rmask(Isrc  ,Jsrc).eq.0.0_r8).and.             &
     &                  (rmask(Isrc-1,Jsrc).eq.1.0_r8)) THEN
!>                    FX(Isrc,Jsrc)=Huon(Isrc,Jsrc,k)*                  &
!>   &                              t(Isrc-1,Jsrc,k,3,itrc)
!>
                      tl_FX(Isrc,Jsrc)=tl_Huon(Isrc,Jsrc,k)*            &
     &                                 t(Isrc-1,Jsrc,k,3,itrc)+         &
     &                                 Huon(Isrc,Jsrc,k)*               &
     &                                 tl_t(Isrc-1,Jsrc,k,3,itrc)
                    ELSE IF ((rmask(Isrc  ,Jsrc).eq.1.0_r8).and.        &
     &                       (rmask(Isrc-1,Jsrc).eq.0.0_r8)) THEN
!>                    FX(Isrc,Jsrc)=Huon(Isrc,Jsrc,k)*                  &
!>   &                              t(Isrc  ,Jsrc,k,3,itrc)
!>
                      tl_FX(Isrc,Jsrc)=tl_Huon(Isrc,Jsrc,k)*            &
     &                                 t(Isrc  ,Jsrc,k,3,itrc)+         &
     &                                 Huon(Isrc,Jsrc,k)*               &
     &                                 tl_t(Isrc  ,Jsrc,k,3,itrc)
                    END IF
# endif
                  END IF
                END IF
              ELSE IF (INT(SOURCES(ng)%Dsrc(is)).eq.1) THEN
# ifdef TS_MPDATA_NOT_YET
                IF (((IstrUm2.le.Isrc).and.(Isrc.le.Iendp2i)).and.      &
     &              ((JstrVm2.le.Jsrc).and.(Jsrc.le.Jendp3))) THEN
# else
                IF (((Istr.le.Isrc).and.(Isrc.le.Iend)).and.            &
     &              ((Jstr.le.Jsrc).and.(Jsrc.le.Jend+1))) THEN
# endif
                  IF (LtracerSrc(itrc,ng)) THEN
!>                  FE(Isrc,Jsrc)=Hvom(Isrc,Jsrc,k)*                    &
!>   &                            SOURCES(ng)%Tsrc(is,k,itrc)
!>
                    tl_FE(Isrc,Jsrc)=tl_Hvom(Isrc,Jsrc,k)*              &
     &                               SOURCES(ng)%Tsrc(is,k,itrc)
# ifdef MASKING
                  ELSE
                    IF ((rmask(Isrc,Jsrc  ).eq.0.0_r8).and.             &
     &                  (rmask(Isrc,Jsrc-1).eq.1.0_r8)) THEN
!>                    FE(Isrc,Jsrc)=Hvom(Isrc,Jsrc,k)*                  &
!>   &                              t(Isrc,Jsrc-1,k,3,itrc)
!>
                      tl_FE(Isrc,Jsrc)=tl_Hvom(Isrc,Jsrc,k)*            &
     &                                 t(Isrc,Jsrc-1,k,3,itrc)+         &
     &                                 Hvom(Isrc,Jsrc,k)*               &
     &                                 tl_t(Isrc,Jsrc-1,k,3,itrc)
                    ELSE IF ((rmask(Isrc,Jsrc  ).eq.1.0_r8).and.        &
     &                       (rmask(Isrc,Jsrc-1).eq.0.0_r8)) THEN
!>                    FE(Isrc,Jsrc)=Hvom(Isrc,Jsrc,k)*                  &
!>   &                              t(Isrc,Jsrc  ,k,3,itrc)
!>
                      tl_FE(Isrc,Jsrc)=tl_Hvom(Isrc,Jsrc,k)*            &
     &                                 t(Isrc,Jsrc  ,k,3,itrc)+         &
     &                                 Hvom(Isrc,Jsrc,k)*               &
     &                                 tl_t(Isrc,Jsrc  ,k,3,itrc)
                    END IF
# endif
                  END IF
                END IF
              END IF
            END DO
          END IF
!
# ifdef TS_MPDATA_NOT_YET
!  Time-step horizontal advection for intermediate diffusive tracer, Ta.
!  Advective fluxes have units of Tunits m3/s.  The new tracer has
!  units of m Tunits.
# else
!  Time-step horizontal advection term.  Advective fluxes have units
!  of Tunits m3/s.  The new tracer has units of m Tunits.
# endif
!
# ifdef TS_MPDATA_NOT_YET
#  define I_RANGE IstrUm2,Iendp2i
#  define J_RANGE JstrVm2,Jendp2i
# else
#  define I_RANGE Istr,Iend
#  define J_RANGE Jstr,Jend
# endif
          DO j=J_RANGE
            DO i=I_RANGE
              cff=dt(ng)*pm(i,j)*pn(i,j)
!>            cff1=cff*(FX(i+1,j)-FX(i,j))
!>
              tl_cff1=cff*(tl_FX(i+1,j)-tl_FX(i,j))
!>            cff2=cff*(FE(i,j+1)-FE(i,j))
!>
              tl_cff2=cff*(tl_FE(i,j+1)-tl_FE(i,j))
!>            cff3=cff1+cff2
!>
              tl_cff3=tl_cff1+tl_cff2
# ifdef TS_MPDATA_NOT_YET
              Ta(i,j,k,itrc)=t(i,j,k,nnew,itrc)-cff3
              tl_Ta(i,j,k,itrc)=tl_t(i,j,k,nnew,itrc)-tl_cff3
# else
!>            t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-cff3
!>
              tl_t(i,j,k,nnew,itrc)=tl_t(i,j,k,nnew,itrc)-tl_cff3
# endif
# ifdef DIAGNOSTICS_TS
#  ifdef TS_MPDATA_NOT_YET
!!            Dhadv(i,j,iTxadv)=-cff1
!!            Dhadv(i,j,iTyadv)=-cff2
!!            Dhadv(i,j,iThadv)=-cff3
#  else
!!            DiaTwrk(i,j,k,itrc,iTxadv)=-cff1
!!            DiaTwrk(i,j,k,itrc,iTyadv)=-cff2
!!            DiaTwrk(i,j,k,itrc,iThadv)=-cff3
#  endif
# endif
            END DO
          END DO
# if defined DIAGNOSTICS_TS && defined TS_MPDATA_NOT_YET
!!        DO j=Jstr,Jend
!!          DO i=Istr,Iend
!!            DiaTwrk(i,j,k,itrc,iTxadv)=Dhadv(i,j,iTxadv)
!!            DiaTwrk(i,j,k,itrc,iTyadv)=Dhadv(i,j,iTyadv)
!!            DiaTwrk(i,j,k,itrc,iThadv)=Dhadv(i,j,iThadv)
!!          END DO
!!        END DO
# endif
        END DO K_LOOP
      END DO T_LOOP
!
!-----------------------------------------------------------------------
!  Time-step vertical advection term.
!-----------------------------------------------------------------------
!
      DO j=J_RANGE
        DO itrc=1,NT(ng)

# if defined TS_SVADVECTION_TL
!
!  Build conservative parabolic splines for the vertical derivatives
!  "FC" of the tracer.  Then, the interfacial "FC" values are
!  converted to vertical advective flux.
!
          DO i=Istr,Iend
#  ifdef NEUMANN
            FC(i,0)=1.5_r8*t(i,j,1,3,itrc)
            CF(i,1)=0.5_r8
#  else
            FC(i,0)=2.0_r8*t(i,j,1,3,itrc)
            CF(i,1)=1.0_r8
#  endif
          END DO
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(2.0_r8*Hz(i,j,k)+                             &
     &                    Hz(i,j,k+1)*(2.0_r8-CF(i,k)))
              CF(i,k+1)=cff*Hz(i,j,k)
              FC(i,k)=cff*(3.0_r8*(Hz(i,j,k  )*t(i,j,k+1,3,itrc)+       &
     &                             Hz(i,j,k+1)*t(i,j,k  ,3,itrc))-      &
     &                     Hz(i,j,k+1)*FC(i,k-1))
            END DO
          END DO
          DO i=Istr,Iend
#  ifdef NEUMANN
            FC(i,N(ng))=(3.0_r8*t(i,j,N(ng),3,itrc)-FC(i,N(ng)-1))/     &
     &                  (2.0_r8-CF(i,N(ng)))
#  else
            FC(i,N(ng))=(2.0_r8*t(i,j,N(ng),3,itrc)-FC(i,N(ng)-1))/     &
     &                  (1.0_r8-CF(i,N(ng)))
#  endif
          END DO
          DO k=N(ng)-1,0,-1
            DO i=Istr,Iend
              FC(i,k)=FC(i,k)-CF(i,k+1)*FC(i,k+1)
            END DO
          END DO
!
!  Now the tangent linear spline code.
!
          DO i=Istr,Iend
#  ifdef NEUMANN
!>          FC(i,0)=1.5_r8*t(i,j,1,3,itrc)
!>
            tl_FC(i,0)=1.5_r8*tl_t(i,j,1,3,itrc)
            CF(i,1)=0.5_r8
#  else
!>          FC(i,0)=2.0_r8*t(i,j,1,3,itrc)
!>
            tl_FC(i,0)=2.0_r8*tl_t(i,j,1,3,itrc)
            CF(i,1)=1.0_r8
#  endif
          END DO
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(2.0_r8*Hz(i,j,k)+                             &
     &                    Hz(i,j,k+1)*(2.0_r8-CF(i,k)))
              CF(i,k+1)=cff*Hz(i,j,k)
              tl_FC(i,k)=cff*                                           &
     &                   (3.0_r8*(Hz(i,j,k  )*tl_t(i,j,k+1,3,itrc)+     &
     &                            Hz(i,j,k+1)*tl_t(i,j,k  ,3,itrc)+     &
     &                            tl_Hz(i,j,k  )*t(i,j,k+1,3,itrc)+     &
     &                            tl_Hz(i,j,k+1)*t(i,j,k  ,3,itrc))-    &
     &                    (tl_Hz(i,j,k+1)*FC(i,k-1)+                    &
     &                     2.0_r8*(tl_Hz(i,j,k  )+                      &
     &                             tl_Hz(i,j,k+1))*FC(i,k)+             &
     &                     tl_Hz(i,j,k  )*FC(i,k+1))-                   &
     &                    Hz(i,j,k+1)*tl_FC(i,k-1))
            END DO
          END DO
          DO i=Istr,Iend
#  ifdef NEUMANN
!>          FC(i,N(ng))=(3.0_r8*t(i,j,N(ng),3,itrc)-FC(i,N(ng)-1))/     &
!>   &                  (2.0_r8-CF(i,N(ng)))
!>
            tl_FC(i,N(ng))=(3.0_r8*tl_t(i,j,N(ng),3,itrc)-              &
     &                      tl_FC(i,N(ng)-1))/                          &
     &                     (2.0_r8-CF(i,N(ng)))
#  else
!>          FC(i,N(ng))=(2.0_r8*t(i,j,N(ng),3,itrc)-FC(i,N(ng)-1))/     &
!>   &                  (1.0_r8-CF(i,N(ng)))
!>
            tl_FC(i,N(ng))=(2.0_r8*tl_t(i,j,N(ng),3,itrc)-              &
     &                      tl_FC(i,N(ng)-1))/                          &
     &                     (1.0_r8-CF(i,N(ng)))
#  endif
          END DO
          DO k=N(ng)-1,0,-1
            DO i=Istr,Iend
!>            FC(i,k)=FC(i,k)-CF(i,k+1)*FC(i,k+1)
!>
              tl_FC(i,k)=tl_FC(i,k)-CF(i,k+1)*tl_FC(i,k+1)
!>            FC(i,k+1)=W(i,j,k+1)*FC(i,k+1)
!>
              tl_FC(i,k+1)=tl_W(i,j,k+1)*FC(i,k+1)+                     &
     &                     W(i,j,k+1)*tl_FC(i,k+1)
            END DO
          END DO
          DO i=Istr,Iend
!>          FC(i,N(ng))=0.0_r8
!>
            tl_FC(i,N(ng))=0.0_r8
!>          FC(i,0)=0.0_r8
!>
            tl_FC(i,0)=0.0_r8
          END DO
!
!  Now complete the computation of the flux array FC.
!
          DO k=N(ng)-1,0,-1
            DO i=Istr,Iend
              FC(i,k+1)=W(i,j,k+1)*FC(i,k+1)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8
            FC(i,0)=0.0_r8
          END DO

# elif defined TS_A4VADVECTION_TL
!
!  Fourth-order, Akima vertical advective flux.
!
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=t(i,j,k+1,3,itrc)-                                &
     &                t(i,j,k  ,3,itrc)
              tl_FC(i,k)=tl_t(i,j,k+1,3,itrc)-                          &
     &                   tl_t(i,j,k  ,3,itrc)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=FC(i,1)
            tl_FC(i,0)=tl_FC(i,1)
            FC(i,N(ng))=FC(i,N(ng)-1)
            tl_FC(i,N(ng))=tl_FC(i,N(ng)-1)
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=2.0_r8*FC(i,k)*FC(i,k-1)
              tl_cff=2.0_r8*(tl_FC(i,k)*FC(i,k-1)+                      &
     &                       FC(i,k)*tl_FC(i,k-1))
              IF (cff.gt.eps) THEN
                CF(i,k)=cff/(FC(i,k)+FC(i,k-1))
                tl_CF(i,k)=((FC(i,k)+FC(i,k-1))*tl_cff-                 &
     &                      cff*(tl_FC(i,k)+tl_FC(i,k-1)))/             &
     &                     ((FC(i,k)+FC(i,k-1))*(FC(i,k)+FC(i,k-1)))
              ELSE
                CF(i,k)=0.0_r8
                tl_CF(i,k)=0.0_r8
              END IF
            END DO
          END DO
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=W(i,j,k)*                                         &
     &                0.5_r8*(t(i,j,k  ,3,itrc)+                        &
     &                        t(i,j,k+1,3,itrc)-                        &
     &                        cff1*(CF(i,k+1)-CF(i,k)))
              tl_FC(i,k)=0.5_r8*                                        &
     &                   (tl_W(i,j,k)*                                  &
     &                    (t(i,j,k  ,3,itrc)+                           &
     &                     t(i,j,k+1,3,itrc)-                           &
     &                     cff1*(CF(i,k+1)-CF(i,k)))+                   &
     &                    W(i,j,k)*                                     &
     &                    (tl_t(i,j,k  ,3,itrc)+                        &
     &                     tl_t(i,j,k+1,3,itrc)-                        &
     &                     cff1*(tl_CF(i,k+1)-tl_CF(i,k))))
            END DO
          END DO
          DO i=Istr,Iend
#  ifdef SED_MORPH
            FC(i,0)=W(i,j,0)*t(i,j,1,3,itrc)
            tl_FC(i,0)=tl_W(i,j,0)*t(i,j,1,3,itrc)+                     &
     &                 W(i,j,0)*tl_t(i,j,1,3,itrc)
#  else
            FC(i,0)=0.0_r8
            tl_FC(i,0)=0.0_r8
#  endif
            FC(i,N(ng))=0.0_r8
            tl_FC(i,N(ng))=0.0_r8
          END DO

# elif defined TS_C2VADVECTION_TL
!
!  Second-order, central differences vertical advective flux.
!
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=W(i,j,k)*                                         &
     &                0.5_r8*(t(i,j,k  ,3,itrc)+                        &
     &                        t(i,j,k+1,3,itrc))
              tl_FC(i,k)=0.5_r8*                                        &
     &                   (tl_W(i,j,k)*                                  &
     &                    (t(i,j,k  ,3,itrc)+                           &
     &                     t(i,j,k+1,3,itrc))+                          &
     &                    W(i,j,k)*                                     &
     &                    (tl_t(i,j,k  ,3,itrc)+                        &
     &                     tl_t(i,j,k+1,3,itrc)))
            END DO
          END DO
          DO i=Istr,Iend
#  ifdef SED_MORPH
            FC(i,0)=W(i,j,0)*t(i,j,1,3,itrc)
            tl_FC(i,0)=tl_W(i,j,0)*t(i,j,1,3,itrc)+                     &
     &                 W(i,j,0)*tl_t(i,j,1,3,itrc)
#  else
            FC(i,0)=0.0_r8
            tl_FC(i,0)=0.0_r8
#  endif
            FC(i,N(ng))=0.0_r8
            tl_FC(i,N(ng))=0.0_r8
          END DO

# elif defined TS_MPDATA_NOT_YET
!
!  First_order, upstream differences vertical advective flux.
!
          DO i=I_RANGE
            DO k=1,N(ng)-1
              cff1=MAX(W(i,j,k),0.0_r8)
              cff2=MIN(W(i,j,k),0.0_r8)
              tl_cff1=(0.5_r8+SIGN(0.5_r8, W(i,j,k)))*tl_W(i,j,k)
              tl_cff2=(0.5_r8+SIGN(0.5_r8,-W(i,j,k)))*tl_W(i,j,k)
              FC(i,k)=cff1*t(i,j,k  ,3,itrc)+                           &
     &                cff2*t(i,j,k+1,3,itrc)
              tl_FC(i,k)=tl_cff1*t(i,j,k  ,3,itrc)+                     &
     &                   cff1*tl_t(i,j,k  ,3,itrc)+                     &
     &                   tl_cff2*t(i,j,k+1,3,itrc)+                     &
     &                   cff2*tl_t(i,j,k+1,3,itrc)
            END DO
#  ifdef SED_MORPH
            FC(i,0)=W(i,j,0)*t(i,j,1,3,itrc)
            tl_FC(i,0)=tl_W(i,j,0)*t(i,j,1,3,itrc)+                     &
     &                 W(i,j,0)*tl_t(i,j,1,3,itrc)
#  else
            FC(i,0)=0.0_r8
            tl_FC(i,0)=0.0_r8
#  endif
            FC(i,N(ng))=0.0_r8
            tl_FC(i,N(ng))=0.0_r8
          END DO

# else
!
!  Fourth-order, central differences vertical advective flux.
!
          cff1=0.5_r8
          cff2=7.0_r8/12.0_r8
          cff3=1.0_r8/12.0_r8
          DO k=2,N(ng)-2
            DO i=Istr,Iend
              FC(i,k)=W(i,j,k)*                                         &
     &                (cff2*(t(i,j,k  ,3,itrc)+                         &
     &                       t(i,j,k+1,3,itrc))-                        &
     &                 cff3*(t(i,j,k-1,3,itrc)+                         &
     &                       t(i,j,k+2,3,itrc)))
              tl_FC(i,k)=tl_W(i,j,k)*                                   &
     &                   (cff2*(t(i,j,k  ,3,itrc)+                      &
     &                          t(i,j,k+1,3,itrc))-                     &
     &                    cff3*(t(i,j,k-1,3,itrc)+                      &
     &                          t(i,j,k+2,3,itrc)))+                    &
     &                   W(i,j,k)*                                      &
     &                   (cff2*(tl_t(i,j,k  ,3,itrc)+                   &
     &                          tl_t(i,j,k+1,3,itrc))-                  &
     &                    cff3*(tl_t(i,j,k-1,3,itrc)+                   &
     &                          tl_t(i,j,k+2,3,itrc)))
            END DO
          END DO
          DO i=Istr,Iend
#  ifdef SED_MORPH
            FC(i,0)=W(i,j,0)*2.0_r8*                                    &
     &              (cff2*t(i,j,1,3,itrc)-                              &
     &               cff3*t(i,j,2,3,itrc))
            tl_FC(i,0)=2.0_r8*                                          &
     &                 (tl_W(i,j,0)*                                    &
     &                  (cff2*t(i,j,1,3,itrc)-                          &
     &                   cff3*t(i,j,2,3,itrc))+                         &
     &                  W(i,j,0)*                                       &
     &                  (cff2*tl_t(i,j,1,3,itrc)-                       &
     &                   cff3*tl_t(i,j,2,3,itrc)))
#  else
            FC(i,0)=0.0_r8
            tl_FC(i,0)=0.0_r8
#  endif
            FC(i,1)=W(i,j,1)*                                           &
     &              (cff1*t(i,j,1,3,itrc)+                              &
     &               cff2*t(i,j,2,3,itrc)-                              &
     &               cff3*t(i,j,3,3,itrc))
            tl_FC(i,1)=tl_W(i,j,1)*                                     &
     &                 (cff1*t(i,j,1,3,itrc)+                           &
     &                  cff2*t(i,j,2,3,itrc)-                           &
     &                  cff3*t(i,j,3,3,itrc))+                          &
     &                 W(i,j,1)*                                        &
     &                 (cff1*tl_t(i,j,1,3,itrc)+                        &
     &                  cff2*tl_t(i,j,2,3,itrc)-                        &
     &                  cff3*tl_t(i,j,3,3,itrc))
            FC(i,N(ng)-1)=W(i,j,N(ng)-1)*                               &
     &                    (cff1*t(i,j,N(ng)  ,3,itrc)+                  &
     &                     cff2*t(i,j,N(ng)-1,3,itrc)-                  &
     &                     cff3*t(i,j,N(ng)-2,3,itrc))
            tl_FC(i,N(ng)-1)=tl_W(i,j,N(ng)-1)*                         &
     &                       (cff1*t(i,j,N(ng)  ,3,itrc)+               &
     &                        cff2*t(i,j,N(ng)-1,3,itrc)-               &
     &                        cff3*t(i,j,N(ng)-2,3,itrc))+              &
     &                       W(i,j,N(ng)-1)*                            &
     &                       (cff1*tl_t(i,j,N(ng)  ,3,itrc)+            &
     &                        cff2*tl_t(i,j,N(ng)-1,3,itrc)-            &
     &                        cff3*tl_t(i,j,N(ng)-2,3,itrc))
            FC(i,N(ng))=0.0_r8
            tl_FC(i,N(ng))=0.0_r8
          END DO
# endif
!
# ifdef TS_MPDATA_NOT_YET
!  Time-step vertical advection for intermediate diffusive tracer, Ta
!  (Tunits).
# else
#  ifdef SPLINES_VDIFF
!  Time-step vertical advection term (Tunits).
!  The BASIC STATE "t" used below must be in transport units, but "t"
!  retrived is in Tunits so we multiply by "Hz".
#  else
!  Time-step vertical advection term (m Tunits).
#  endif
# endif
!
          DO i=I_RANGE
            CF(i,0)=dt(ng)*pm(i,j)*pn(i,j)
          END DO
!
!  Apply mass point sources (volume vertical influx), if any.
!
          IF (LwSrc(ng).and.ANY(LtracerSrc(:,ng))) THEN
            DO is=1,Nsrc(ng)
              Isrc=SOURCES(ng)%Isrc(is)
              Jsrc=SOURCES(ng)%Jsrc(is)
              IF (LtracerSrc(itrc,ng).and.                              &
# ifdef TS_MPDATA_NOT_YET
     &            ((IstrUm2.le.Isrc).and.(Isrc.le.Iendp2i)).and.        &
# else
     &            ((IstrR.le.Isrc).and.(Isrc.le.IendR)).and.            &
# endif
     &            (j.eq.Jsrc)) THEN
                DO k=1,N(ng)-1
                  FC(Isrc,k)=FC(Isrc,k)+0.5_r8*                         &
     &                       (SOURCES(ng)%Qsrc(is,k  )*                 &
     &                        SOURCES(ng)%Tsrc(is,k  ,itrc)+            &
     &                        SOURCES(ng)%Qsrc(is,k+1)*                 &
     &                        SOURCES(ng)%Tsrc(is,k+1,itrc))
                  tl_FC(Isrc,k)=tl_FC(Isrc,k)+0.0_r8
                END DO
              END IF
            END DO
          END IF
!
          DO k=1,N(ng)
            DO i=I_RANGE
              cff1=CF(i,0)*(FC(i,k)-FC(i,k-1))
              tl_cff1=CF(i,0)*(tl_FC(i,k)-tl_FC(i,k-1))
# ifdef TS_MPDATA_NOT_YET
              Ta(i,j,k,itrc)=(Ta(i,j,k,itrc)-cff1)*oHz(i,j,k)
              tl_Ta(i,j,k,itrc)=(tl_Ta(i,j,k,itrc)-tl_cff1)*            &
     &                          oHz(i,j,k)+                             &
     &                          (Ta(i,j,k,itrc)-cff1)*                  &
     &                          tl_oHz(i,j,k)
#  ifdef DIAGNOSTICS_TS
!!            Dvadv(i,j,k,itrc)=-cff1
#  endif
# else
!>            t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-cff1
!>
              tl_t(i,j,k,nnew,itrc)=tl_t(i,j,k,nnew,itrc)-tl_cff1
#  ifdef SPLINES_VDIFF
!>            t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)*oHz(i,j,k)
!>
              tl_t(i,j,k,nnew,itrc)=tl_t(i,j,k,nnew,itrc)*              &
     &                              oHz(i,j,k)+                         &
     &                              (t(i,j,k,nnew,itrc)*Hz(i,j,k))*     &
     &                              tl_oHz(i,j,k)
#  endif
#  ifdef DIAGNOSTICS_TS
!!            DiaTwrk(i,j,k,itrc,iTvadv)=-cff1
!!            DO idiag=1,NDT
!!              DiaTwrk(i,j,k,itrc,idiag)=DiaTwrk(i,j,k,itrc,idiag)*    &
!!   &                                    oHz(i,j,k)
!!            END DO
#  endif
# endif
            END DO
          END DO
        END DO
# undef I_RANGE
# undef J_RANGE
# ifdef TS_MPDATA_NOT_YET
      END DO
!
!-----------------------------------------------------------------------
!  Compute anti-diffusive velocities to corrected advected tracers
!  using MPDATA recursive method.  Notice that pipelined J-loop ended.
!-----------------------------------------------------------------------
!
      DO itrc=1,NT(ng)
!>      CALL mpdata_adiff_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          IminS, ImaxS, JminS, JmaxS,             &
#  ifdef MASKING
!>   &                          rmask, umask, vmask,                    &
#  endif
#  ifdef WET_DRY
!>   &                          rmask_wet, umask_wet, vmask_wet,        &
#  endif
!>   &                          pm, pn, omn, om_u, on_v,                &
!>   &                          z_r, oHz,                               &
!>   &                          Huon, Hvom, W,                          &
!>   &                          t(:,:,:,3,itrc),                        &
!>   &                          Ta(:,:,:,itrc),  Ua, Va, Wa)
!>
        CALL tl_mpdata_adiff_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
#  ifdef MASKING
     &                             rmask, umask, vmask,                 &
#  endif
#  ifdef WET_DRY
     &                             rmask_wet, umask_wet, vmask_wet,     &
#  endif
     &                             pm, pn, omn, om_u, on_v,             &
     &                             z_r, tl_z_r,                         &
     &                             oHz, tl_oHz,                         &
     &                             Huon, tl_Huon,                       &
     &                             Hvom, tl_Hvom,                       &
     &                             W, tl_W,                             &
     &                             t(:,:,:,3,itrc), tl_t(:,:,:,3,itrc), &
     &                             Ta(:,:,:,itrc), tl_Ta(:,:,:,itrc),   &
     &                             Ua, tl_Ua,                           &
     &                             Va, tl_Va,                           &
     &                             Wa, tl_Wa)
!
!  Compute anti-diffusive corrected advection fluxes.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              cff1=MAX(Ua(i,j,k),0.0_r8)
              cff2=MIN(Ua(i,j,k),0.0_r8)
              tl_cff1=(0.5_r8+SIGN(0.5_r8, Ua(i,j,k)))*tl_Ua(i,j,k)
              tl_cff2=(0.5_r8+SIGN(0.5_r8,-Ua(i,j,k)))*tl_Ua(i,j,k)
!>            FX(i,j)=(cff1*Ta(i-1,j,k,itrc)+                           &
!>   &                 cff2*Ta(i  ,j,k,itrc))*                          &
!>   &                0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)
!>
              tl_FX(i,j)=0.5_r8*on_u(i,j)*                              &
     &                   ((tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*                &
     &                    (cff1*Ta(i-1,j,k,itrc)+                       &
     &                     cff2*Ta(i  ,j,k,itrc))+                      &
     &                    (Hz(i,j,k)+Hz(i-1,j,k))*                      &
     &                    (tl_cff1*Ta(i-1,j,k,itrc)+                    &
     &                     cff1*tl_Ta(i-1,j,k,itrc)+                    &
     &                     tl_cff2*Ta(i  ,j,k,itrc)+                    &
     &                     cff2*tl_Ta(i  ,j,k,itrc)))
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              cff1=MAX(Va(i,j,k),0.0_r8)
              cff2=MIN(Va(i,j,k),0.0_r8)
              tl_cff1=(0.5_r8+SIGN(0.5_r8, Va(i,j,k)))*tl_Va(i,j,k)
              tl_cff2=(0.5_r8+SIGN(0.5_r8,-Va(i,j,k)))*tl_Va(i,j,k)
!>            FE(i,j)=(cff1*Ta(i,j-1,k,itrc)+                           &
!>   &                 cff2*Ta(i,j  ,k,itrc))*                          &
!>   &                0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*om_v(i,j)
!>
              tl_FE(i,j)=0.5_r8*om_v(i,j)*                              &
     &                   ((tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*                &
     &                    (cff1*Ta(i,j-1,k,itrc)+                       &
     &                     cff2*Ta(i,j  ,k,itrc))+                      &
     &                    (Hz(i,j,k)+Hz(i,j-1,k))*                      &
     &                    (tl_cff1*Ta(i,j-1,k,itrc)+                    &
     &                     cff1*tl_Ta(i,j-1,k,itrc)+                    &
     &                     tl_cff2*Ta(i,j  ,k,itrc)+                    &
     &                     cff2*tl_Ta(i,j  ,k,itrc)))
            END DO
          END DO
!
!  Time-step corrected horizontal advection (Tunits m).
!
          DO j=Jstr,Jend
            DO i=Istr,Iend
              cff=dt(ng)*pm(i,j)*pn(i,j)
!>            cff1=cff*(FX(i+1,j)-FX(i,j))
!>
              tl_cff1=cff*(tl_FX(i+1,j)-tl_FX(i,j))
!>            cff2=cff*(FE(i,j+1)-FE(i,j))
!>
              tl_cff2=cff*(tl_FE(i,j+1)-tl_FE(i,j))
!>            cff3=cff1+cff2
!>
              tl_cff3=tl_cff1+tl_cff2
!>            t(i,j,k,nnew,itrc)=Ta(i,j,k,itrc)*Hz(i,j,k)-cff3
!>
              tl_t(i,j,k,nnew,itrc)=tl_Ta(i,j,k,itrc)*Hz(i,j,k)+        &
     &                              Ta(i,j,k,itrc)*tl_Hz(i,j,k)-tl_cff3
#  ifdef DIAGNOSTICS_TS
!!            DiaTwrk(i,j,k,itrc,iTxadv)=DiaTwrk(i,j,k,itrc,iTxadv)-    &
!!   &                                   cff1
!!            DiaTwrk(i,j,k,itrc,iTyadv)=DiaTwrk(i,j,k,itrc,iTyadv)-    &
!!   &                                   cff2
!!            DiaTwrk(i,j,k,itrc,iThadv)=DiaTwrk(i,j,k,itrc,iThadv)-    &
!!   &                                   cff3
#  endif
            END DO
          END DO
        END DO
!
!  Compute anti-diffusive corrected vertical advection flux.
!
        DO j=Jstr,Jend
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff1=MAX(Wa(i,j,k),0.0_r8)
              cff2=MIN(Wa(i,j,k),0.0_r8)
              tl_cff1=(0.5_r8+SIGN(0.5, Wa(i,j,k)))*tl_Wa(i,j,k)
              tl_cff2=(0.5_r8+SIGN(0.5,-Wa(i,j,k)))*tl_Wa(i,j,k)
              FC(i,k)=cff1*Ta(i,j,k  ,itrc)+                            &
     &                cff2*Ta(i,j,k+1,itrc)
              tl_FC(i,k)=tl_cff1*Ta(i,j,k  ,itrc)+                      &
     &                   cff1*tl_Ta(i,j,k  ,itrc)+                      &
     &                   tl_cff2*Ta(i,j,k+1,itrc)+                      &
     &                   cff2*tl_Ta(i,j,k+1,itrc)
            END DO
          END DO
          DO i=Istr,Iend
!>          FC(i,0)=0.0_r8
!>
            tl_FC(i,0)=0.0_r8
!>          FC(i,N(ng))=0.0_r8
!>
            tl_FC(i,N(ng))=0.0_r8
          END DO
!
!  Time-step corrected vertical advection (Tunits).
#  ifdef DIAGNOSTICS_TS
!  Convert units of tracer diagnostic terms to Tunits.
#  endif
!
          DO i=Istr,Iend
            CF(i,0)=dt(ng)*pm(i,j)*pn(i,j)
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
!>            cff1=CF(i,0)*(FC(i,k)-FC(i,k-1))
!>
              tl_cff1=CF(i,0)*(tl_FC(i,k)-tl_FC(i,k-1))
!>            t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-cff1
!>
              tl_t(i,j,k,nnew,itrc)=tl_t(i,j,k,nnew,itrc)-tl_cff1
#  ifdef DIAGNOSTICS_TS
!!            DiaTwrk(i,j,k,itrc,iTvadv)=Dvadv(i,j,k,itrc)-             &
!!   &                                   cff1
!!            DO idiag=1,NDT
!!              DiaTwrk(i,j,k,itrc,idiag)=DiaTwrk(i,j,k,itrc,idiag)*    &
!!   &                                    oHz(i,j,k)
!!            END DO
#  endif
            END DO
          END DO
        END DO
      END DO
!
!  Start pipelined J-loop.
!
      DO j=Jstr,Jend
# endif /* TS_MPDATA_NOT_YET */
!
!-----------------------------------------------------------------------
!  Time-step vertical diffusion term.
!-----------------------------------------------------------------------
!
        DO itrc=1,NT(ng)
          ltrc=MIN(NAT,itrc)

# if defined SPLINES_VDIFF && !defined TS_MPDATA_NOT_YET
!
!  Use conservative, parabolic spline reconstruction of BASIC STATE
!  vertical diffusion derivatives.  Solve BASIC STATE tridiagonal
!  system.
!
          cff1=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=cff1*Hz(i,j,k  )-                                 &
     &                dt(ng)*Akt(i,j,k-1,ltrc)*oHz(i,j,k  )
              CF(i,k)=cff1*Hz(i,j,k+1)-                                 &
     &                dt(ng)*Akt(i,j,k+1,ltrc)*oHz(i,j,k+1)
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,0)=0.0_r8
            DC(i,0)=0.0_r8
          END DO
!
!  LU decomposition and forward substitution.
!
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              BC(i,k)=cff1*(Hz(i,j,k)+Hz(i,j,k+1))+                     &
     &                dt(ng)*Akt(i,j,k,ltrc)*(oHz(i,j,k)+oHz(i,j,k+1))
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
              CF(i,k)=cff*CF(i,k)
              DC(i,k)=cff*(t(i,j,k+1,nnew,itrc)-t(i,j,k,nnew,itrc)-     &
     &                     FC(i,k)*DC(i,k-1))
            END DO
          END DO
!
!  Backward substitution. Save DC for the tangent linearization.
!  DC is scaled later by AKt.
!
          DO i=Istr,Iend
            DC(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
            END DO
          END DO
!
!  Use conservative, parabolic spline reconstruction of tangent linear
!  vertical diffusion derivatives.  Then, time step vertical diffusion
!  term implicitly.
!
!  Note that the BASIC STATE "t" must in Tunits when used in the
!  tangent spline routine below, which it does in the present code.
!
          cff1=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=cff1*Hz(i,j,k  )-                                 &
     &                dt(ng)*Akt(i,j,k-1,ltrc)*oHz(i,j,k  )
              tl_FC(i,k)=cff1*tl_Hz(i,j,k  )-                           &
     &                   dt(ng)*(tl_Akt(i,j,k-1,ltrc)*oHz(i,j,k  )+     &
     &                           Akt(i,j,k-1,ltrc)*tl_oHz(i,j,k  ))
              CF(i,k)=cff1*Hz(i,j,k+1)-                                 &
     &                dt(ng)*Akt(i,j,k+1,ltrc)*oHz(i,j,k+1)
              tl_CF(i,k)=cff1*tl_Hz(i,j,k+1)-                           &
     &                   dt(ng)*(tl_Akt(i,j,k+1,ltrc)*oHz(i,j,k+1)+     &
     &                           Akt(i,j,k+1,ltrc)*tl_oHz(i,j,k+1))
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,0)=0.0_r8
            tl_CF(i,0)=0.0_r8
            tl_DC(i,0)=0.0_r8
          END DO
!
!  Tangent linear LU decomposition and forward substitution.
!
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              BC(i,k)=cff1*(Hz(i,j,k)+Hz(i,j,k+1))+                     &
     &                dt(ng)*Akt(i,j,k,ltrc)*(oHz(i,j,k)+oHz(i,j,k+1))
              tl_BC(i,k)=cff1*(tl_Hz(i,j,k)+tl_Hz(i,j,k+1))+            &
     &                   dt(ng)*(tl_Akt(i,j,k,ltrc)*                    &
     &                           (oHz(i,j,k)+oHz(i,j,k+1))+             &
     &                           Akt(i,j,k,ltrc)*                       &
     &                           (tl_oHz(i,j,k)+tl_oHz(i,j,k+1)))
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
              CF(i,k)=cff*CF(i,k)
              tl_DC(i,k)=cff*(tl_t(i,j,k+1,nnew,itrc)-                  &
     &                        tl_t(i,j,k  ,nnew,itrc)-                  &
     &                        (tl_FC(i,k)*DC(i,k-1)+                    &
     &                         tl_BC(i,k)*DC(i,k  )+                    &
     &                         tl_CF(i,k)*DC(i,k+1))-                   &
     &                        FC(i,k)*tl_DC(i,k-1))
            END DO
          END DO
!
!  Tangent linear backward substitution.
!
          DO i=Istr,Iend
            tl_DC(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
            END DO
          END DO
!
!  Compute tl_DC before multiplying BASIC STATE spline gradients
!  DC by AKt.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              tl_DC(i,k)=tl_DC(i,k)*Akt(i,j,k,ltrc)+                    &
     &                   DC(i,k)*tl_Akt(i,j,k,ltrc)
              DC(i,k)=DC(i,k)*Akt(i,j,k,ltrc)
!>            cff1=dt(ng)*oHz(i,j,k)*(DC(i,k)-DC(i,k-1))
!>
              tl_cff1=dt(ng)*(tl_oHz(i,j,k)*(DC(i,k)-DC(i,k-1))+        &
     &                        oHz(i,j,k)*(tl_DC(i,k)-tl_DC(i,k-1)))
!>            t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+cff1
!>
              tl_t(i,j,k,nnew,itrc)=tl_t(i,j,k,nnew,itrc)+tl_cff1
#  ifdef DIAGNOSTICS_TS
!!            DiaTwrk(i,j,k,itrc,iTvdif)=DiaTwrk(i,j,k,itrc,iTvdif)+    &
!!   &                                   cff1
#  endif
            END DO
          END DO
# else
!
!  Compute off-diagonal coefficients FC [lambda*dt*Akt/Hz] for the
!  implicit vertical diffusion terms at future time step, located
!  at horizontal RHO-points and vertical W-points.
!  Also set FC at the top and bottom levels.
!
!  NOTE: The original code solves the tridiagonal system A*t=r where
!        A is a matrix and t and r are vectors. We need to solve the
!        tangent linear C of this system which is A*tl_t+tl_A*t=tl_r.
!        Here, tl_A*t and tl_r are known, so we must solve for tl_t
!        by inverting A*tl_t=tl_r-tl_A*t.
!
          cff=-dt(ng)*lambda
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff1=1.0_r8/(z_r(i,j,k+1)-z_r(i,j,k))
              tl_cff1=-cff1*cff1*(tl_z_r(i,j,k+1)-tl_z_r(i,j,k))
              FC(i,k)=cff*cff1*Akt(i,j,k,ltrc)
              tl_FC(i,k)=cff*(tl_cff1*Akt(i,j,k,ltrc)+                  &
     &                        cff1*tl_Akt(i,j,k,ltrc))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            tl_FC(i,0)=0.0_r8
            FC(i,N(ng))=0.0_r8
            tl_FC(i,N(ng))=0.0_r8
          END DO
!
!  Compute diagonal matrix coefficients BC.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              BC(i,k)=Hz(i,j,k)-FC(i,k)-FC(i,k-1)
              tl_BC(i,k)=tl_Hz(i,j,k)-tl_FC(i,k)-tl_FC(i,k-1)
            END DO
          END DO
!
!  Solve the tangent linear tridiagonal system.
!  (DC is a tangent linear variable here).
!
          DO k=2,N(ng)-1
            DO i=Istr,Iend
              DC(i,k)=tl_t(i,j,k,nnew,itrc)-                            &
     &                (tl_FC(i,k-1)*t(i,j,k-1,nnew,itrc)+               &
     &                 tl_BC(i,k  )*t(i,j,k  ,nnew,itrc)+               &
     &                 tl_FC(i,k  )*t(i,j,k+1,nnew,itrc))
            END DO
          END DO
          DO i=Istr,Iend
            DC(i,1)=tl_t(i,j,1,nnew,itrc)-                              &
     &              (tl_BC(i,1)*t(i,j,1,nnew,itrc)+                     &
     &               tl_FC(i,1)*t(i,j,2,nnew,itrc))
            DC(i,N(ng))=tl_t(i,j,N(ng),nnew,itrc)-                      &
     &                  (tl_FC(i,N(ng)-1)*t(i,j,N(ng)-1,nnew,itrc)+     &
     &                   tl_BC(i,N(ng)  )*t(i,j,N(ng)  ,nnew,itrc))
          END DO
!
          DO i=Istr,Iend
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FC(i,1)
            DC(i,1)=cff*DC(i,1)
          END DO
          DO k=2,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
              CF(i,k)=cff*FC(i,k)
              DC(i,k)=cff*(DC(i,k)-FC(i,k-1)*DC(i,k-1))
            END DO
          END DO
!
!  Compute new solution by back substitution.
!  (DC is a tangent linear variable here).
!
          DO i=Istr,Iend
#  ifdef DIAGNOSTICS_TS
!!          cff1=t(i,j,N(ng),nnew,itrc)*oHz(i,j,N(ng))
#  endif
            DC(i,N(ng))=(DC(i,N(ng))-FC(i,N(ng)-1)*DC(i,N(ng)-1))/      &
     &                  (BC(i,N(ng))-FC(i,N(ng)-1)*CF(i,N(ng)-1))
            tl_t(i,j,N(ng),nnew,itrc)=DC(i,N(ng))
#  ifdef DIAGNOSTICS_TS
!!          DiaTwrk(i,j,N(ng),itrc,iTvdif)=                             &
!!   &                             DiaTwrk(i,j,N(ng),itrc,iTvdif)+      &
!!   &                             t(i,j,N(ng),nnew,itrc)-cff1
#  endif
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
#  ifdef DIAGNOSTICS_TS
!!            cff1=t(i,j,k,nnew,itrc)*oHz(i,j,k)
#  endif
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
              tl_t(i,j,k,nnew,itrc)=DC(i,k)
#  ifdef DIAGNOSTICS_TS
!!            DiaTwrk(i,j,k,itrc,iTvdif)=DiaTwrk(i,j,k,itrc,iTvdif)+    &
!!   &                                   t(i,j,k,nnew,itrc)-cff1
#  endif
            END DO
          END DO
# endif
        END DO
      END DO

# if defined AGE_MEAN && defined T_PASSIVE
!
!-----------------------------------------------------------------------
!  If inert passive tracer and Mean Age, compute age concentration (even
!  inert index) forced by the right-hand-side term that is concentration
!  of an associated conservative passive tracer (odd inert index). Mean
!  Age is age concentration divided by conservative passive tracer
!  concentration. Code implements NPT/2 mean age tracer pairs.
!
!  Implemented and tested by W.G. Zhang and J. Wilkin. See following
!  reference for details.
!
!   Zhang et al. (2010): Simulation of water age and residence time in
!      the New York Bight, JPO, 40,965-982, doi:10.1175/2009JPO4249.1
!-----------------------------------------------------------------------
!
      DO itrc=1,NPT,2
        iage=inert(itrc+1)                     ! even inert tracer index
        DO k=1,N(ng)
          DO j=Jstr,Jend
            DO i=Istr,Iend
!>            t(i,j,k,nnew,iage)=t(i,j,k,nnew,iage)+                    &
!>   &                           dt(ng)*                                &
#  ifdef TS_MPDATA
!>   &                           t(i,j,k,nnew,inert(itrc))
#  else
!>   &                           t(i,j,k,3,inert(itrc))
#  endif
!>
              tl_t(i,j,k,nnew,iage)=tl_t(i,j,k,nnew,iage)+              &
     &                              dt(ng)*                             &
#  ifdef TS_MPDATA
     &                              tl_t(i,j,k,nnew,inert(itrc))
#  else
     &                              tl_t(i,j,k,3,inert(itrc))
#  endif
            END DO
          END DO
        END DO
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions and, if appropriate, nudge
!  to tracer data and apply Land/Sea mask.
!-----------------------------------------------------------------------
!
!  Initialize tracer counter index. The "tclm" array is only allocated
!  to the NTCLM fields that need to be processed. This is done to
!  reduce memory.
!
      ic=0
!
      DO itrc=1,NT(ng)
!
!  Set compact reduced memory tracer index for nudging coefficients and
!  climatology arrays.
!
        IF (LtracerCLM(itrc,ng).and.LnudgeTCLM(itrc,ng)) THEN
          ic=ic+1
        END IF
!
!  Set lateral boundary conditions.
!
!>      CALL t3dbc_tile (ng, tile, itrc, ic,                            &
!>   &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
!>   &                   IminS, ImaxS, JminS, JmaxS,                    &
!>   &                   nstp, nnew,                                    &
!>   &                   t)
!>
        CALL tl_t3dbc_tile (ng, tile, itrc, ic,                         &
     &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nstp, nnew,                                 &
     &                      tl_t)
!
!  Nudge towards tracer climatology.
!
        IF (LtracerCLM(itrc,ng).and.LnudgeTCLM(itrc,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
!>              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+                  &
!>   &                             dt(ng)*                              &
!>   &                             CLIMA(ng)%Tnudgcof(i,j,k,ic)*        &
!>   &                             (CLIMA(ng)%tclm(i,j,k,ic)-           &
!>   &                              t(i,j,k,nnew,itrc))
!>
                tl_t(i,j,k,nnew,itrc)=tl_t(i,j,k,nnew,itrc)-            &
     &                                dt(ng)*                           &
     &                                CLIMA(ng)%Tnudgcof(i,j,k,ic)*     &
     &                                tl_t(i,j,k,nnew,itrc)
              END DO
            END DO
          END DO
        END IF

# ifdef MASKING
!
!  Apply Land/Sea mask.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
!>            t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)*rmask(i,j)
!>
              tl_t(i,j,k,nnew,itrc)=tl_t(i,j,k,nnew,itrc)*rmask(i,j)
            END DO
          END DO
        END DO
# endif
# ifdef DIAGNOSTICS_TS
!!
!!  Compute time-rate-of-change diagnostic term.
!!
!!      DO k=1,N(ng)
!!        DO j=JstrR,JendR
!!          DO i=IstrR,IendR
!!            DiaTwrk(i,j,k,itrc,iTrate)=t(i,j,k,nnew,itrc)-            &
!!   &                                   DiaTwrk(i,j,k,itrc,iTrate)
!!            DiaTwrk(i,j,k,itrc,iTrate)=t(i,j,k,nnew,itrc)-            &
!!   &                                   t(i,j,k,nstp,itrc)
!!          END DO
!!        END DO
!!      END DO
# endif
!
!  Apply periodic boundary conditions.
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>        CALL exchange_r3d_tile (ng, tile,                             &
!>   &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
!>   &                            t(:,:,:,nnew,itrc))
!>
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            tl_t(:,:,:,nnew,itrc))
        END IF
      END DO
# ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
!>    CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    t(:,:,:,nnew,:))
!>
      CALL mp_exchange4d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_t(:,:,:,nnew,:))
# endif
# if defined FLOATS_NOT_YET && defined FLOAT_VWALK
!
!-----------------------------------------------------------------------
!  Compute vertical gradient in vertical T-diffusion coefficient for
!  floats random walk.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          DO k=1,N(ng)
            dAktdz(i,j,k)=(Akt(i,j,k,1)-Akt(i,j,k-1,1))/Hz(i,j,k)
          END DO
        END DO
      END DO
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          dAktdz)
      END IF

#  ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    dAktdz)
#  endif
# endif

      RETURN
      END SUBROUTINE tl_step3d_t_tile
#endif
      END MODULE tl_step3d_t_mod
