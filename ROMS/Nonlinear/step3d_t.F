#include "cppdefs.h"
      MODULE step3d_t_mod
#if !defined TS_FIXED && (defined NONLINEAR && defined SOLVE3D)
!
!svn $Id: step3d_t.F 1054 2021-03-06 19:47:12Z arango $
!=======================================================================
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine time-steps tracer equations.  Notice that advective    !
!  and diffusive terms are time-stepped differently. It applies the    !
!  corrector time-step for horizontal/vertical advection,  vertical    !
!  diffusion, nudging if necessary, and lateral boundary conditions.   !
!                                                                      !
!  A different horizontal/vertical advection scheme is allowed for     !
!  each tracer. If the MPDATA or HSIMT monotonic scheme, it is applied !
!  to both horizontal and vertical advective fluxes.                   !
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
      PUBLIC  :: step3d_t
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE step3d_t (ng, tile)
!***********************************************************************
!
      USE mod_param
# ifdef DIAGNOSTICS_TS
      USE mod_diags
# endif
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
# if defined TRC_PSOURCE
      USE mod_trc_sources
# endif
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 35, __LINE__, MyFile)
# endif
      CALL step3d_t_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nrhs(ng), nstp(ng), nnew(ng),                 &
# if defined TRC_PSOURCE
     &                    Msrcpt(ng), Nsrcpt(ng),                       &
     &                    TRC_SOURCES(ng) % Isrcpt,                     &
     &                    TRC_SOURCES(ng) % Jsrcpt,                     &
     &                    TRC_SOURCES(ng) % Tsrcpt,                     &
     &                    TRC_SOURCES(ng) % Lsrcpt,                     &
# endif
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
     &                    GRID(ng) % umask,                             &
     &                    GRID(ng) % vmask,                             &
# endif
# ifdef WET_DRY
     &                    GRID(ng) % rmask_wet,                         &
     &                    GRID(ng) % umask_wet,                         &
     &                    GRID(ng) % vmask_wet,                         &
# endif
     &                    GRID(ng) % omn,                               &
     &                    GRID(ng) % om_u,                              &
     &                    GRID(ng) % om_v,                              &
     &                    GRID(ng) % on_u,                              &
     &                    GRID(ng) % on_v,                              &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
     &                    GRID(ng) % Huon,                              &
     &                    GRID(ng) % Hvom,                              &
     &                    GRID(ng) % z_r,                               &
     &                    MIXING(ng) % Akt,                             &
#ifdef TS_VAR
     &                    MIXING(ng) % Akr,                             &
#endif
     &                    OCEAN(ng) % W,                                &
# ifdef WEC_VF
     &                    OCEAN(ng) % W_stokes,                         &
# endif
# if defined FLOATS && defined FLOAT_VWALK
     &                    MIXING(ng) % dAktdz,                          &
# endif
# ifdef DIAGNOSTICS_TS
     &                    DIAGS(ng) % DiaTwrk,                          &
# endif
     &                    OCEAN(ng) % t)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 35, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE step3d_t
!
!***********************************************************************
      SUBROUTINE step3d_t_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nrhs, nstp, nnew,                       &
# if defined TRC_PSOURCE
     &                          Msrcpt, Nsrcpt,                         &
     &                          Isrcpt, Jsrcpt, Tsrcpt, Lsrcpt,         &
# endif
# ifdef MASKING
     &                          rmask, umask, vmask,                    &
# endif
# ifdef WET_DRY
     &                          rmask_wet, umask_wet, vmask_wet,        &
# endif
     &                          omn, om_u, om_v, on_u, on_v,            &
     &                          pm, pn,                                 &
     &                          Hz, Huon, Hvom,                         &
     &                          z_r,                                    &
     &                          Akt,                                    &
#ifdef TS_VAR
     &                          Akr,                                    &
#endif
     &                          W,                                      &
# ifdef WEC_VF
     &                          W_stokes,                               &
# endif
# if defined FLOATS && defined FLOAT_VWALK
     &                          dAktdz,                                 &
# endif
# ifdef DIAGNOSTICS_TS
     &                          DiaTwrk,                                &
# endif
     &                          t)
!***********************************************************************
!
      USE mod_param
      USE mod_clima
      USE mod_ncparam
# if defined NESTING && !defined ONE_WAY
      USE mod_nesting
# endif
      USE mod_scalars
      USE mod_sources
!
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
      USE mp_exchange_mod, ONLY : mp_exchange4d
# endif
      USE mpdata_adiff_mod
# ifdef NESTING
      USE nesting_mod, ONLY : bry_fluxes
# endif
      USE t3dbc_mod, ONLY : t3dbc_tile
# ifdef T_PASSIVE
      USE pt3dbc_mod, ONLY : pt3dbc_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, nstp, nnew
# if defined TRC_PSOURCE
      integer, intent(in) :: Msrcpt, Nsrcpt
# endif

# ifdef ASSUMED_SHAPE
#  if defined TRC_PSOURCE
      logical, intent(in)  :: Lsrcpt(:,:)
      integer, intent(in)  :: Isrcpt(:)
      integer, intent(in)  :: Jsrcpt(:)
      real(r8), intent(in) :: Tsrcpt(:,:,:)
#  endif
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
#  ifdef SUN
      real(r8), intent(in) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
#  else
      real(r8), intent(in) :: Akt(LBi:,LBj:,0:,:)
#  endif
#  ifdef TS_VAR
      real(r8), intent(inout) :: Akr(LBi:,LBj:,0:)
#  endif
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
#  ifdef WEC_VF
      real(r8), intent(in) :: W_stokes(LBi:,LBj:,0:)
#  endif
#  ifdef DIAGNOSTICS_TS
      real(r8), intent(inout) :: DiaTwrk(LBi:,LBj:,:,:,:)
#  endif
#  ifdef SUN
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#  else
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#  endif
#  if defined FLOATS && defined FLOAT_VWALK
      real(r8), intent(out) :: dAktdz(LBi:,LBj:,:)
#  endif

# else

#  if defined TRC_PSOURCE
      logical, intent(in) :: Lsrcpt(Msrcpt,NT(ng))
      integer, intent(in) :: Isrcpt(Msrcpt)
      integer, intent(in) :: Jsrcpt(Msrcpt)
      real(r8), intent(in) :: Tsrcpt(Msrcpt,N(ng),NT(ng))
#  endif
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
#  ifdef TS_VAR
      real(r8), intent(inout) :: Akr(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
      real(r8), intent(in) :: W(LBi:UBi,LBj:UBj,0:N(ng))
#  ifdef WEC_VF
      real(r8), intent(in) :: W_stokes(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
#  ifdef DIAGNOSTICS_TS
      real(r8), intent(inout) :: DiaTwrk(LBi:UBi,LBj:UBj,N(ng),NT(ng),  &
     &                                   NDT)
#  endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))

#  if defined FLOATS && defined FLOAT_VWALK
      real(r8), intent(out) :: dAktdz(LBi:UBi,LBj:UBj,N(ng))
#  endif
# endif
!
!  Local variable declarations.
!
      logical :: LapplySrc, Lhsimt, Lmpdata
!
# ifdef NESTING
      integer :: ILB, IUB, JLB, JUB
      integer :: dg, cr, rg
# endif
      integer :: IminT, ImaxT, JminT, JmaxT
      integer :: Isrc, Jsrc
      integer :: i, ibt, ic, ii, is, itrc, j, jj, k, ltrc
# if defined AGE_MEAN && defined T_PASSIVE
      integer :: iage
# endif
# ifdef DIAGNOSTICS_TS
      integer :: idiag
# endif
      real(r8) :: eps = 1.0E-16_r8
      real(r8) :: eps1 = 1.0E-12_r8

      real(r8) :: cff, cff1, cff2, cff3
      real(r8) :: betaL, betaR, betaD, betaU
      real(r8) :: rL, rR, rD, rU, rkaL, rkaR, rkaD, rkaU
      real(r8) :: a1, b1, sw, sw_eta, sw_xi

      real(r8), dimension(IminS:ImaxS) :: gradX, KaX, oKaX
      real(r8), dimension(JminS:JmaxS) :: gradE, KaE, oKaE
      real(r8), dimension(0:N(ng))     :: gradZ, KaZ, oKaZ

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: BC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: curv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: oHz

# ifdef DIAGNOSTICS_TS
      real(r8), allocatable :: Dhadv(:,:,:)
      real(r8), allocatable :: Dvadv(:,:,:,:)
# endif
      real(r8) :: my_maxbio(15)
      real(r8), allocatable :: Ta(:,:,:,:)
      real(r8), allocatable :: Ua(:,:,:)
      real(r8), allocatable :: Va(:,:,:)
      real(r8), allocatable :: Wa(:,:,:)

# include "set_bounds.h"

# ifdef NESTING
!
!  Notice that the trace flux boundary arrays are dimensioned with the
!  global dimensions of grid to facilitate processing.
!
      ILB=BOUNDS(ng)%LBi(-1)
      IUB=BOUNDS(ng)%UBi(-1)
      JLB=BOUNDS(ng)%LBj(-1)
      JUB=BOUNDS(ng)%UBj(-1)
# endif
!
!-----------------------------------------------------------------------
!  Time-step horizontal advection term.
!-----------------------------------------------------------------------
!
      Lhsimt =ANY(Hadvection(:,ng)%HSIMT).and.                          &
     &        ANY(Vadvection(:,ng)%HSIMT)
      Lmpdata=ANY(Hadvection(:,ng)%MPDATA).and.                         &
     &        ANY(Vadvection(:,ng)%MPDATA)
!
!  Allocate local arrays for MPDATA.
!
      IF (Lmpdata) THEN
# ifdef DIAGNOSTICS_TS
        IF (.not.allocated(Dhadv)) THEN
          allocate ( Dhadv(IminS:ImaxS,JminS:JmaxS,3) )
          Dhadv=0.0_r8
        END IF
        IF (.not.allocated(Dvadv)) THEN
          allocate ( Dvadv(IminS:ImaxS,JminS:JmaxS,N(ng),NT(ng)) )
          Dvadv=0.0_r8
        END IF
# endif
        IF (.not.allocated(Ta)) THEN
          allocate ( Ta(IminS:ImaxS,JminS:JmaxS,N(ng),NT(ng)) )
          Ta=0.0_r8
        END IF
        IF (.not.allocated(Ua)) THEN
          allocate ( Ua(IminS:ImaxS,JminS:JmaxS,N(ng)) )
          Ua=0.0_r8
        END IF
        IF (.not.allocated(Va)) THEN
          allocate ( Va(IminS:ImaxS,JminS:JmaxS,N(ng)) )
          Va=0.0_r8
        END IF
        IF (.not.allocated(Wa)) THEN
          allocate ( Wa(IminS:ImaxS,JminS:JmaxS,0:N(ng)) )
          Wa=0.0_r8
        END IF
      END IF
!
!  Compute reciprocal thickness, 1/Hz.
!
      IF (Lmpdata.or.Lhsimt) THEN
        DO k=1,N(ng)
          DO j=Jstrm2,Jendp2
            DO i=Istrm2,Iendp2
              oHz(i,j,k)=1.0_r8/Hz(i,j,k)
            END DO
          END DO
        END DO
      ELSE
        DO k=1,N(ng)
          DO j=Jstr,Jend
            DO i=Istr,Iend
              oHz(i,j,k)=1.0_r8/Hz(i,j,k)
            END DO
          END DO
        END DO
      END IF
!
!  Horizontal tracer advection.  It is possible to have a different
!  advection schme for each tracer.
!
      T_LOOP1 : DO itrc=1,NT(ng)
!
!  The MPDATA and HSIMT algorithms requires a three-point footprint, so
!  exchange boundary data on t(:,:,:,nnew,:) so other processes computed
!  earlier (horizontal diffusion, biology, or sediment) are accounted.
!
        IF ((Hadvection(itrc,ng)%MPDATA).or.                            &
     &      (Hadvection(itrc,ng)%HSIMT)) THEN
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              t(:,:,:,nnew,itrc))
          END IF

# ifdef DISTRIBUTE
          CALL mp_exchange3d (ng, tile, iNLM, 1,                        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        t(:,:,:,nnew,itrc))
# endif
        END IF
!
!  Compute horizontal tracer advection fluxes.
!
        K_LOOP : DO k=1,N(ng)
!
          HADV_FLUX : IF (Hadvection(itrc,ng)%CENTERED2) THEN
!
!  Second-order, centered differences horizontal advective fluxes.
!
            DO j=Jstr,Jend
              DO i=Istr,Iend+1
                FX(i,j)=Huon(i,j,k)*                                    &
     &                  0.5_r8*(t(i-1,j,k,3,itrc)+                      &
     &                          t(i  ,j,k,3,itrc))
              END DO
            END DO
            DO j=Jstr,Jend+1
              DO i=Istr,Iend
                FE(i,j)=Hvom(i,j,k)*                                    &
     &                  0.5_r8*(t(i,j-1,k,3,itrc)+                      &
     &                          t(i,j  ,k,3,itrc))
              END DO
            END DO
!
          ELSE IF (Hadvection(itrc,ng)%MPDATA) THEN
!
!  First-order, upstream differences horizontal advective fluxes.
!
            DO j=JstrVm2,Jendp2i
              DO i=IstrUm2,Iendp3
                cff1=MAX(Huon(i,j,k),0.0_r8)
                cff2=MIN(Huon(i,j,k),0.0_r8)
                FX(i,j)=cff1*t(i-1,j,k,3,itrc)+                         &
     &                  cff2*t(i  ,j,k,3,itrc)
              END DO
            END DO
            DO j=JstrVm2,Jendp3
              DO i=IstrUm2,Iendp2i
                cff1=MAX(Hvom(i,j,k),0.0_r8)
                cff2=MIN(Hvom(i,j,k),0.0_r8)
                FE(i,j)=cff1*t(i,j-1,k,3,itrc)+                         &
     &                  cff2*t(i,j  ,k,3,itrc)
              END DO
            END DO
!
          ELSE IF (Hadvection(itrc,ng)%HSIMT) THEN
!
!  Third High-order Spatial Interpolation at the Middle Temporal level
!  (HSIMT; Wu and Zhu, 2010) with a Total Variation Diminishing (TVD)
!  limiter horizontal advection fluxes.
!
!  Hui Wu and Jianrong Zhu, 2010: Advection scheme with 3rd high-order
!    spatial interpolation at the middle temporal level and its
!    application to saltwater intrusion in the Changjiang Estuary,
!    Ocean Modelling 33, 33-51, doi:10.1016/j.ocemod.2009.12.001
!
            DO j=Jstr,Jend
              DO i=IstrU-1,Iendp2
                cff=0.125_r8*(pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))*   &
     &              dt(ng)
                cff1=cff*(oHz(i-1,j,k)+oHz(i,j,k))
                gradX(i)=t(i,j,k,3,itrc)-t(i-1,j,k,3,itrc)
                KaX(i)=1.0_r8-ABS(Huon(i,j,k)*cff1)
# ifdef MASKING
                gradX(i)=gradX(i)*umask(i,j)
                KaX(i)=KaX(i)*umask(i,j)
# endif
              END DO
              IF (.not.EWperiodic(ng)) THEN
                IF (DOMAIN(ng)%Western_Edge(tile)) THEN
                  IF (Huon(Istr,j,k).ge.0.0_r8) THEN
                    gradX(Istr-1)=0.0_r8
                    KaX(Istr-1)=0.0_r8
                  END IF
                END IF
                IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
                  IF (Huon(Iend+1,j,k).lt.0.0_r8) THEN
                    gradX(Iend+2)=0.0_r8
                    KaX(Iend+2)=0.0_r8
                  END IF
                END IF
              END IF
              DO i=Istr,Iend+1
                IF (KaX(i).le.eps1) THEN
                  oKaX(i)=0.0_r8
                ELSE
                  oKaX(i)=1.0_r8/MAX(KaX(i),eps1)
                END IF
                IF (Huon(i,j,k).ge.0.0_r8) THEN
                  IF (ABS(gradX(i)).le.eps1) THEN
                    rL=0.0_r8
                    rkaL=0.0_r8
                  ELSE
                    rL=gradX(i-1)/gradX(i)
                    rkaL=KaX(i-1)*oKaX(i)
                  END IF
                  a1= cc1*KaX(i)+cc2-cc3*oKaX(i)
                  b1=-cc1*KaX(i)+cc2+cc3*oKaX(i)
                  betaL=a1+b1*rL
                  cff=0.5_r8*MAX(0.0_r8,                                &
     &                           MIN(2.0_r8, 2.0_r8*rL*rkaL, betaL))*   &
     &                gradX(i)*KaX(i)
# ifdef MASKING
                  ii=MAX(i-2,0)
                  cff=cff*rmask(ii,j)
# endif
                  sw_xi=t(i-1,j,k,3,itrc)+cff
                ELSE
                  IF (ABS(gradX(i)).le.eps1) THEN
                    rR=0.0_r8
                    rkaR=0.0_r8
                  ELSE
                    rR=gradX(i+1)/gradX(i)
                    rkaR=KaX(i+1)*oKaX(i)
                  END IF
                  a1= cc1*KaX(i)+cc2-cc3*oKaX(i)
                  b1=-cc1*KaX(i)+cc2+cc3*oKaX(i)
                  betaR=a1+b1*rR
                  cff=0.5_r8*MAX(0.0_r8,                                &
     &                           MIN(2.0_r8, 2.0_r8*rR*rkaR, betaR))*   &
     &                gradX(i)*KaX(i)
# ifdef MASKING
                  ii=MIN(i+1,Lm(ng)+1)
                  cff=cff*rmask(ii,j)
# endif
                  sw_xi=t(i,j,k,3,itrc)-cff
                END IF
                FX(i,j)=sw_xi*Huon(i,j,k)
              END DO
            END DO
!
            DO i=Istr,Iend
              DO j=JstrV-1,Jendp2
                cff=0.125_r8*(pn(i,j)+pn(i,j-1))*(pm(i,j)+pm(i,j-1))*   &
     &              dt(ng)
                cff1=cff*(oHz(i,j,k)+oHz(i,j-1,k))
                gradE(j)=t(i,j,k,3,itrc)-t(i,j-1,k,3,itrc)
                KaE(j)=1.0_r8-ABS(Hvom(i,j,k)*cff1)
# ifdef MASKING
                gradE(j)=gradE(j)*vmask(i,j)
                KaE(j)=KaE(j)*vmask(i,j)
# endif
              END DO
              IF (.not.NSperiodic(ng)) THEN
                IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
                  IF (Hvom(i,Jstr,k).ge.0.0_r8) THEN
                    gradE(Jstr-1)=0.0_r8
                    KaE(Jstr-1)=0.0_r8
                  END IF
                END IF
                IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
                  IF (Hvom(i,Jend+1,k).lt.0.0_r8) THEN
                    gradE(Jend+2)=0.0_r8
                    KaE(Jend+2)=0.0_r8
                  END IF
                END IF
              END IF
              DO j=Jstr,Jend+1
                IF (KaE(j).le.eps1) THEN
                  oKaE(j)=0.0_r8
                ELSE
                  oKaE(j)=1.0_r8/MAX(KaE(j),eps1)
                END IF
                IF (Hvom(i,j,k).ge.0.0_r8) THEN
                  IF (ABS(gradE(j)).le.eps1) THEN
                    rD=0.0_r8
                    rkaD=0.0_r8
                  ELSE
                    rD=gradE(j-1)/gradE(j)
                    rkaD=KaE(j-1)*oKaE(j)
                  END IF
                  a1= cc1*KaE(j)+cc2-cc3*oKaE(j)
                  b1=-cc1*KaE(j)+cc2+cc3*oKaE(j)
                  betaD=a1+b1*rD
                  cff=0.5_r8*MAX(0.0_r8,                                &
     &                           MIN(2.0_r8, 2.0_r8*rD*rkaD, betaD))*   &
     &                gradE(j)*KaE(j)
# ifdef MASKING
                  jj=MAX(j-2,0)
                  cff=cff*rmask(i,jj)
# endif
                  sw_eta=t(i,j-1,k,3,itrc)+cff
                ELSE
                  IF (ABS(gradE(j)).le.eps1) THEN
                    rU=0.0_r8
                    rkaU=0.0_r8
                  ELSE
                    rU=gradE(j+1)/gradE(j)
                    rkaU=KaE(j+1)*oKaE(j)
                  END IF
                  a1= cc1*KaE(j)+cc2-cc3*oKaE(j)
                  b1=-cc1*KaE(j)+cc2+cc3*oKaE(j)
                  betaU=a1+b1*rU
                  cff=0.5*MAX(0.0_r8,                                   &
     &                        MIN(2.0_r8, 2.0_r8*rU*rkaU, betaU))*      &
     &                gradE(j)*KaE(j)
# ifdef MASKING
                  jj=MIN(j+1,Mm(ng)+1)
                  cff=cff*rmask(i,jj)
# endif
                  sw_eta=t(i,j,k,3,itrc)-cff
                END IF
                FE(i,j)=sw_eta*Hvom(i,j,k)
              END DO
            END DO
!
          ELSE IF ((Hadvection(itrc,ng)%AKIMA4).or.                     &
     &             (Hadvection(itrc,ng)%CENTERED4).or.                  &
     &             (Hadvection(itrc,ng)%SPLIT_U3).or.                   &
     &             (Hadvection(itrc,ng)%UPSTREAM3)) THEN
!
!  Fourth-order Akima, fourth-order centered differences, or third-order
!  upstream-biased horizontal advective fluxes.
!
            DO j=Jstr,Jend
              DO i=Istrm1,Iendp2
                FX(i,j)=t(i  ,j,k,3,itrc)-                              &
     &                  t(i-1,j,k,3,itrc)
# ifdef MASKING
                FX(i,j)=FX(i,j)*umask(i,j)
# endif
              END DO
            END DO
            IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
              IF (DOMAIN(ng)%Western_Edge(tile)) THEN
                DO j=Jstr,Jend
                  FX(Istr-1,j)=FX(Istr,j)
                END DO
              END IF
            END IF
            IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
              IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
                DO j=Jstr,Jend
                  FX(Iend+2,j)=FX(Iend+1,j)
                END DO
              END IF
            END IF
!
            DO j=Jstr,Jend
              DO i=Istr-1,Iend+1
                IF (Hadvection(itrc,ng)%UPSTREAM3) THEN
                  curv(i,j)=FX(i+1,j)-FX(i,j)
                ELSE IF (Hadvection(itrc,ng)%AKIMA4) THEN
                  cff=2.0_r8*FX(i+1,j)*FX(i,j)
                  IF (cff.gt.eps) THEN
                    grad(i,j)=cff/(FX(i+1,j)+FX(i,j))
                  ELSE
                    grad(i,j)=0.0_r8
                  END IF
                ELSE IF ((Hadvection(itrc,ng)%CENTERED4).or.            &
     &                   (Hadvection(itrc,ng)%SPLIT_U3)) THEN
                  grad(i,j)=0.5_r8*(FX(i+1,j)+FX(i,j))
                END IF
              END DO
            END DO
!
            cff1=1.0_r8/6.0_r8
            cff2=1.0_r8/3.0_r8
            DO j=Jstr,Jend
              DO i=Istr,Iend+1
                IF (Hadvection(itrc,ng)%UPSTREAM3) THEN
                  FX(i,j)=Huon(i,j,k)*0.5_r8*                           &
     &                    (t(i-1,j,k,3,itrc)+                           &
     &                     t(i  ,j,k,3,itrc))-                          &
     &                    cff1*(curv(i-1,j)*MAX(Huon(i,j,k),0.0_r8)+    &
     &                          curv(i  ,j)*MIN(Huon(i,j,k),0.0_r8))
                ELSE IF ((Hadvection(itrc,ng)%AKIMA4).or.               &
     &                   (Hadvection(itrc,ng)%CENTERED4).or.            &
     &                   (Hadvection(itrc,ng)%SPLIT_U3)) THEN
                  FX(i,j)=Huon(i,j,k)*0.5_r8*                           &
     &                    (t(i-1,j,k,3,itrc)+                           &
     &                     t(i  ,j,k,3,itrc)-                           &
     &                     cff2*(grad(i  ,j)-                           &
     &                           grad(i-1,j)))
                END IF
              END DO
            END DO
!
            DO j=Jstrm1,Jendp2
              DO i=Istr,Iend
                FE(i,j)=t(i,j  ,k,3,itrc)-                              &
     &                  t(i,j-1,k,3,itrc)
# ifdef MASKING
                FE(i,j)=FE(i,j)*vmask(i,j)
# endif
              END DO
            END DO
            IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
              IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
                DO i=Istr,Iend
                  FE(i,Jstr-1)=FE(i,Jstr)
                END DO
              END IF
            END IF
            IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
              IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
                DO i=Istr,Iend
                  FE(i,Jend+2)=FE(i,Jend+1)
                END DO
              END IF
            END IF
!
            DO j=Jstr-1,Jend+1
              DO i=Istr,Iend
                IF (Hadvection(itrc,ng)%UPSTREAM3) THEN
                  curv(i,j)=FE(i,j+1)-FE(i,j)
                ELSE IF (Hadvection(itrc,ng)%AKIMA4) THEN
                  cff=2.0_r8*FE(i,j+1)*FE(i,j)
                  IF (cff.gt.eps) THEN
                    grad(i,j)=cff/(FE(i,j+1)+FE(i,j))
                  ELSE
                    grad(i,j)=0.0_r8
                  END IF
                ELSE IF ((Hadvection(itrc,ng)%CENTERED4).or.            &
     &                   (Hadvection(itrc,ng)%SPLIT_U3)) THEN
                  grad(i,j)=0.5_r8*(FE(i,j+1)+FE(i,j))
                END IF
              END DO
            END DO
!
            cff1=1.0_r8/6.0_r8
            cff2=1.0_r8/3.0_r8
            DO j=Jstr,Jend+1
              DO i=Istr,Iend
                IF (Hadvection(itrc,ng)%UPSTREAM3) THEN
                  FE(i,j)=Hvom(i,j,k)*0.5_r8*                           &
     &                    (t(i,j-1,k,3,itrc)+                           &
     &                     t(i,j  ,k,3,itrc))-                          &
     &                    cff1*(curv(i,j-1)*MAX(Hvom(i,j,k),0.0_r8)+    &
     &                          curv(i,j  )*MIN(Hvom(i,j,k),0.0_r8))
                ELSE IF ((Hadvection(itrc,ng)%AKIMA4).or.               &
     &                   (Hadvection(itrc,ng)%CENTERED4).or.            &
     &                   (Hadvection(itrc,ng)%SPLIT_U3)) THEN
                  FE(i,j)=Hvom(i,j,k)*0.5_r8*                           &
     &                    (t(i,j-1,k,3,itrc)+                           &
     &                     t(i,j  ,k,3,itrc)-                           &
     &                     cff2*(grad(i,j  )-                           &
     &                           grad(i,j-1)))
                END IF
              END DO
            END DO
          END IF HADV_FLUX
!
!  Apply tracers point sources to the horizontal advection terms,
!  if any.
!
          IF (LuvSrc(ng)) THEN
            DO is=1,Nsrc(ng)
              Isrc=SOURCES(ng)%Isrc(is)
              Jsrc=SOURCES(ng)%Jsrc(is)
              IF (INT(SOURCES(ng)%Dsrc(is)).eq.0) THEN
                IF ((Hadvection(itrc,ng)%MPDATA).or.                    &
     &              (Hadvection(itrc,ng)%HSIMT)) THEN
                  LapplySrc=(IstrUm2.le.Isrc).and.                      &
     &                      (Isrc.le.Iendp3).and.                       &
     &                      (JstrVm2.le.Jsrc).and.                      &
     &                      (Jsrc.le.Jendp2i)
                ELSE
                  LapplySrc=(Istr.le.Isrc).and.                         &
     &                      (Isrc.le.Iend+1).and.                       &
     &                      (Jstr.le.Jsrc).and.                         &
     &                      (Jsrc.le.Jend)
                END IF
                IF (LapplySrc) THEN
                  IF (LtracerSrc(itrc,ng)) THEN
                    FX(Isrc,Jsrc)=Huon(Isrc,Jsrc,k)*                    &
# ifdef ONE_TRACER_SOURCE
     &                            SOURCES(ng)%Tsrc(itrc)
# elif defined TWO_D_TRACER_SOURCE
     &                            SOURCES(ng)%Tsrc(is,itrc)
# else
     &                            SOURCES(ng)%Tsrc(is,k,itrc)
# endif
# ifdef MASKING
                  ELSE
                    IF ((rmask(Isrc  ,Jsrc).eq.0.0_r8).and.             &
     &                  (rmask(Isrc-1,Jsrc).eq.1.0_r8)) THEN
                      FX(Isrc,Jsrc)=Huon(Isrc,Jsrc,k)*                  &
     &                              t(Isrc-1,Jsrc,k,3,itrc)
                    ELSE IF ((rmask(Isrc  ,Jsrc).eq.1.0_r8).and.        &
     &                       (rmask(Isrc-1,Jsrc).eq.0.0_r8)) THEN
                      FX(Isrc,Jsrc)=Huon(Isrc,Jsrc,k)*                  &
     &                              t(Isrc  ,Jsrc,k,3,itrc)
                    END IF
# endif
                  END IF
                END IF
              ELSE IF (INT(SOURCES(ng)%Dsrc(is)).eq.1) THEN
                IF ((Hadvection(itrc,ng)%MPDATA).or.                    &
     &              (Hadvection(itrc,ng)%HSIMT)) THEN
                  LapplySrc=(IstrUm2.le.Isrc).and.                      &
     &                      (Isrc.le.Iendp2i).and.                      &
     &                      (JstrVm2.le.Jsrc).and.                      &
     &                      (Jsrc.le.Jendp3)
                ELSE
                  LapplySrc=(Istr.le.Isrc).and.                         &
     &                      (Isrc.le.Iend).and.                         &
     &                      (Jstr.le.Jsrc).and.                         &
     &                      (Jsrc.le.Jend+1)
                END IF
                IF (LapplySrc) THEN
                  IF (LtracerSrc(itrc,ng)) THEN
                    FE(Isrc,Jsrc)=Hvom(Isrc,Jsrc,k)*                    &
# ifdef ONE_TRACER_SOURCE
     &                            SOURCES(ng)%Tsrc(itrc)
# elif defined TWO_D_TRACER_SOURCE
     &                            SOURCES(ng)%Tsrc(is,itrc)
# else
     &                            SOURCES(ng)%Tsrc(is,k,itrc)
# endif
# ifdef MASKING
                  ELSE
                    IF ((rmask(Isrc,Jsrc  ).eq.0.0_r8).and.             &
     &                  (rmask(Isrc,Jsrc-1).eq.1.0_r8)) THEN
                      FE(Isrc,Jsrc)=Hvom(Isrc,Jsrc,k)*                  &
     &                              t(Isrc,Jsrc-1,k,3,itrc)
                    ELSE IF ((rmask(Isrc,Jsrc  ).eq.1.0_r8).and.        &
     &                       (rmask(Isrc,Jsrc-1).eq.0.0_r8)) THEN
                      FE(Isrc,Jsrc)=Hvom(Isrc,Jsrc,k)*                  &
     &                              t(Isrc,Jsrc  ,k,3,itrc)
                    END IF
# endif
                  END IF
                END IF
              END IF
            END DO
          END IF

# if defined NESTING && !defined ONE_WAY
!
!  If refinement grids, extract tracer horizontal advection fluxes
!  (Hz*u*T/n, Hz*v*T/m) at the grid contact boundary (physical
!  domain perimeter) to be used in two-way nesting.
!
          IF (RefinedGrid(ng)) THEN
            DO cr=1,Ncontact
              dg=Rcontact(cr)%donor_grid
              rg=Rcontact(cr)%receiver_grid
              IF (ng.eq.rg) THEN
                CALL bry_fluxes (dg, rg, cr, iNLM, tile,                &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           ILB, IUB, JLB, JUB,                    &
     &                           dt(ng), FX, FE,                        &
     &                          BRY_CONTACT(iwest, cr)%Tflux(:,k,itrc), &
     &                          BRY_CONTACT(ieast, cr)%Tflux(:,k,itrc), &
     &                          BRY_CONTACT(isouth,cr)%Tflux(:,k,itrc), &
     &                          BRY_CONTACT(inorth,cr)%Tflux(:,k,itrc))
              END IF
            END DO
          END IF
# endif
!
!  If MPDATA, time-step horizontal advection for intermediate diffusive
!  tracer, Ta (m Tunits).
!
          HADV_STEPPING : IF (Hadvection(itrc,ng)%MPDATA) THEN
            DO j=JstrVm2,Jendp2i
              DO i=IstrUm2,Iendp2i
                cff=dt(ng)*pm(i,j)*pn(i,j)
                cff1=cff*(FX(i+1,j)-FX(i,j))
                cff2=cff*(FE(i,j+1)-FE(i,j))
                cff3=cff1+cff2
                Ta(i,j,k,itrc)=t(i,j,k,nnew,itrc)-cff3
# ifdef DIAGNOSTICS_TS
                Dhadv(i,j,iTxadv)=-cff1
                Dhadv(i,j,iTyadv)=-cff2
                Dhadv(i,j,iThadv)=-cff3
# endif
              END DO
            END DO
# ifdef DIAGNOSTICS_TS
!
            DO j=Jstr,Jend
              DO i=Istr,Iend
                DiaTwrk(i,j,k,itrc,iTxadv)=Dhadv(i,j,iTxadv)
                DiaTwrk(i,j,k,itrc,iTyadv)=Dhadv(i,j,iTyadv)
                DiaTwrk(i,j,k,itrc,iThadv)=Dhadv(i,j,iThadv)
              END DO
            END DO
# endif
!
!  OTHERWISE, time-step horizontal advection term.  Advective fluxes
!  have units of Tunits m3/s.  The new tracer has units of m Tunits.
!
          ELSE
            DO j=Jstr,Jend
              DO i=Istr,Iend
                cff=dt(ng)*pm(i,j)*pn(i,j)
                cff1=cff*(FX(i+1,j)-FX(i,j))
                cff2=cff*(FE(i,j+1)-FE(i,j))
                cff3=cff1+cff2
                t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-cff3
# ifdef DIAGNOSTICS_TS
                DiaTwrk(i,j,k,itrc,iTxadv)=-cff1
                DiaTwrk(i,j,k,itrc,iTyadv)=-cff2
                DiaTwrk(i,j,k,itrc,iThadv)=-cff3
# endif
              END DO
            END DO
          END IF HADV_STEPPING
        END DO K_LOOP
      END DO T_LOOP1
!
!-----------------------------------------------------------------------
!  Time-step vertical advection term.
!-----------------------------------------------------------------------
!
      T_LOOP2 : DO itrc=1,NT(ng)
        IF (Vadvection(itrc,ng)%MPDATA) THEN
          JminT=JstrVm2
          JmaxT=Jendp2i
        ELSE
          JminT=Jstr
          JmaxT=Jend
        END IF
!
        J_LOOP1 : DO j=JminT,JmaxT              ! start pipelined J-loop
!
          VADV_FLUX : IF (Vadvection(itrc,ng)%SPLINES) THEN
!
!  Build conservative parabolic splines for the vertical derivatives
!  "FC" of the tracer.  Then, the interfacial "FC" values are
!  converted to vertical advective flux (Tunits m3/s).
!
            DO i=Istr,Iend
# ifdef NEUMANN
              FC(i,0)=1.5_r8*t(i,j,1,3,itrc)
              CF(i,1)=0.5_r8
# else
              FC(i,0)=2.0_r8*t(i,j,1,3,itrc)
              CF(i,1)=1.0_r8
# endif
            END DO
            DO k=1,N(ng)-1
              DO i=Istr,Iend
                cff=1.0_r8/(2.0_r8*Hz(i,j,k)+                           &
     &                      Hz(i,j,k+1)*(2.0_r8-CF(i,k)))
                CF(i,k+1)=cff*Hz(i,j,k)
                FC(i,k)=cff*(3.0_r8*(Hz(i,j,k  )*t(i,j,k+1,3,itrc)+     &
     &                               Hz(i,j,k+1)*t(i,j,k  ,3,itrc))-    &
     &                       Hz(i,j,k+1)*FC(i,k-1))
              END DO
            END DO
            DO i=Istr,Iend
# ifdef NEUMANN
              FC(i,N(ng))=(3.0_r8*t(i,j,N(ng),3,itrc)-FC(i,N(ng)-1))/   &
     &                    (2.0_r8-CF(i,N(ng)))
# else
              FC(i,N(ng))=(2.0_r8*t(i,j,N(ng),3,itrc)-FC(i,N(ng)-1))/   &
     &                    (1.0_r8-CF(i,N(ng)))
# endif
            END DO
            DO k=N(ng)-1,0,-1
              DO i=Istr,Iend
                FC(i,k)=FC(i,k)-CF(i,k+1)*FC(i,k+1)
                FC(i,k+1)=W(i,j,k+1)*FC(i,k+1)
# if defined WEC_VF
                FC(i,k+1)=FC(i,k+1)+W_stokes(i,j,k+1)*FC(i,k+1)
# endif
              END DO
            END DO
            DO i=Istr,Iend
              FC(i,N(ng))=0.0_r8
              FC(i,0)=0.0_r8
            END DO
!
          ELSE IF (Vadvection(itrc,ng)%AKIMA4) THEN
!
!  Fourth-order, Akima vertical advective flux (Tunits m3/s).
!
            DO k=1,N(ng)-1
              DO i=Istr,Iend
                FC(i,k)=t(i,j,k+1,3,itrc)-                              &
     &                  t(i,j,k  ,3,itrc)
              END DO
            END DO
            DO i=Istr,Iend
              FC(i,0)=FC(i,1)
              FC(i,N(ng))=FC(i,N(ng)-1)
            END DO
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff=2.0_r8*FC(i,k)*FC(i,k-1)
                IF (cff.gt.eps) THEN
                  CF(i,k)=cff/(FC(i,k)+FC(i,k-1))
                ELSE
                  CF(i,k)=0.0_r8
                END IF
              END DO
            END DO
            cff1=1.0_r8/3.0_r8
            DO k=1,N(ng)-1
              DO i=Istr,Iend
# if defined WEC_VF
                FC(i,k)=(W(i,j,k)+W_stokes(i,j,k))*                     &
# else
                FC(i,k)=W(i,j,k)*                                       &
# endif
     &                  0.5_r8*(t(i,j,k  ,3,itrc)+                      &
     &                          t(i,j,k+1,3,itrc)-                      &
     &                          cff1*(CF(i,k+1)-CF(i,k)))
              END DO
            END DO
            DO i=Istr,Iend
# ifdef SED_MORPH
              FC(i,0)=W(i,j,0)*t(i,j,1,3,itrc)
# else
              FC(i,0)=0.0_r8
# endif
              FC(i,N(ng))=0.0_r8
            END DO
!
          ELSE IF (Vadvection(itrc,ng)%CENTERED2) THEN
!
!  Second-order, central differences vertical advective flux
!  (Tunits m3/s).
!
            DO k=1,N(ng)-1
              DO i=Istr,Iend
# if defined WEC_VF
                FC(i,k)=(W(i,j,k)+W_stokes(i,j,k))*                     &
# else
                FC(i,k)=W(i,j,k)*                                       &
# endif
     &                  0.5_r8*(t(i,j,k  ,3,itrc)+                      &
     &                          t(i,j,k+1,3,itrc))
              END DO
            END DO
            DO i=Istr,Iend
# ifdef SED_MORPH
              FC(i,0)=W(i,j,0)*t(i,j,1,3,itrc)
# else
              FC(i,0)=0.0_r8
# endif
              FC(i,N(ng))=0.0_r8
            END DO
!
          ELSE IF (Vadvection(itrc,ng)%MPDATA) THEN
!
!  First_order, upstream differences vertical advective flux
!  (Tunits m3/s).
!
            DO i=IstrUm2,Iendp2i
              DO k=1,N(ng)-1
# if defined WEC_VF
                cff1=MAX(W(i,j,k)+W_stokes(i,j,k),0.0_r8)
                cff2=MIN(W(i,j,k)+W_stokes(i,j,k),0.0_r8)
# else
                cff1=MAX(W(i,j,k),0.0_r8)
                cff2=MIN(W(i,j,k),0.0_r8)
# endif
                FC(i,k)=cff1*t(i,j,k  ,3,itrc)+                         &
     &                  cff2*t(i,j,k+1,3,itrc)
              END DO
# ifdef SED_MORPH
              FC(i,0)=W(i,j,0)*t(i,j,1,3,itrc)
# else
              FC(i,0)=0.0_r8
# endif
              FC(i,N(ng))=0.0_r8
            END DO
!
          ELSE IF (Vadvection(itrc,ng)%HSIMT) THEN
!
!  Third High-order Spatial Interpolation at the Middle Temporal level
!  (HSIMT; Wu and Zhu, 2010) with a Total Variation Diminishing (TVD)
!  limiter vertical advection flux (Tunits m3/s).
!
            DO i=Istr,Iend
              KaZ(0)=0.0_r8
              oKaZ(0)=0.0_r8
              gradZ(0)=0.0_r8
              DO k=1,N(ng)-1
                cff=pm(i,j)*pn(i,j)*dt(ng)
# if defined WEC_VF
                KaZ(k)=1.0_r8-ABS(cff*(W(i,j,k)+W_stokes(i,j,k))/       &
     &                            (z_r(i,j,k+1)-z_r(i,j,k)))
# else
                KaZ(k)=1.0_r8-ABS(cff*W(i,j,k)/                         &
     &                            (z_r(i,j,k+1)-z_r(i,j,k)))
# endif
                oKaZ(k)=1.0_r8/KaZ(k)
                gradZ(k)=t(i,j,k+1,3,itrc)-t(i,j,k,3,itrc)
              END DO
              KaZ(N(ng))=0.0_r8
              okaZ(N(ng))=0.0_r8
              gradZ(N(ng))=0.0_r8
!
              DO k=1,N(ng)-1
# if defined WEC_VF
                cff1=W(i,j,k)+W_stokes(i,j,k)
# else
                cff1=W(i,j,k)
# endif
                IF ((k.eq.1).and.(cff1.ge.0.0_r8)) THEN
                  FC(i,k)=cff1*t(i,j,k,3,itrc)
                ELSE IF ((k.eq.N(ng)-1).and.(cff1.lt.0.0_r8)) THEN
                  FC(i,k)=cff1*t(i,j,k+1,3,itrc)
                ELSE
                  IF (cff1.ge.0) THEN
                    IF (ABS(gradZ(k)).le.eps1) THEN
                      rD=0.0_r8
                      rkaD=0.0_r8
                    ELSE
                      rD=gradZ(k-1)/gradZ(k)
                      rkaD=KaZ(k-1)*oKaZ(k)
                    END IF
                    a1= cc1*KaZ(k)+cc2-cc3*oKaZ(k)
                    b1=-cc1*KaZ(k)+cc2+cc3*oKaZ(k)
                    betaD=a1+b1*rD
                    cff=0.5_r8*MAX(0.0_r8,                              &
     &                             MIN(2.0_r8, 2.0_r8*rD*rkaD, betaD))* &
     &                  gradZ(k)*KaZ(k)
                    sw=t(i,j,k,3,itrc)+cff
                  ELSE
                    IF (ABS(gradZ(k)).le.eps1) THEN
                      rU=0.0_r8
                      rkaU=0.0_r8
                    ELSE
                      rU=gradZ(k+1)/gradZ(k)
                      rkaU=KaZ(k+1)*oKaZ(k)
                    END IF
                    a1= cc1*KaZ(k)+cc2-cc3*oKaZ(k)
                    b1=-cc1*KaZ(k)+cc2+cc3*oKaZ(k)
                    betaU=a1+b1*rU
                    cff=0.5_r8*MAX(0.0_r8,                              &
     &                             MIN(2.0_r8, 2.0_r8*rU*rkaU, betaU))* &
     &                  gradZ(k)*KaZ(k)
                    sw=t(i,j,k+1,3,itrc)-cff
                  END IF
                  FC(i,k)=cff1*sw
                END IF
              END DO
# ifdef SED_MORPH
              FC(i,0)=W(i,j,0)*t(i,j,1,3,itrc)
# else
              FC(i,0)=0.0_r8
# endif
              FC(i,N(ng))=0.0_r8
            END DO
!
          ELSE IF ((Vadvection(itrc,ng)%CENTERED4).or.                  &
     &             (Vadvection(itrc,ng)%SPLIT_U3)) THEN
!
!  Fourth-order, central differences vertical advective flux
!  (Tunits m3/s).
!
            cff1=0.5_r8
            cff2=7.0_r8/12.0_r8
            cff3=1.0_r8/12.0_r8
            DO k=2,N(ng)-2
              DO i=Istr,Iend
# ifdef WEC_VF
                FC(i,k)=(W(i,j,k)+W_stokes(i,j,k))*                     &
# else
                FC(i,k)=W(i,j,k)*                                       &
# endif
     &                  (cff2*(t(i,j,k  ,3,itrc)+                       &
     &                         t(i,j,k+1,3,itrc))-                      &
     &                   cff3*(t(i,j,k-1,3,itrc)+                       &
     &                         t(i,j,k+2,3,itrc)))
              END DO
            END DO
            DO i=Istr,Iend
# ifdef SED_MORPH
              FC(i,0)=W(i,j,0)*2.0_r8*                                  &
     &                (cff2*t(i,j,1,3,itrc)-                            &
     &                 cff3*t(i,j,2,3,itrc))
# else
              FC(i,0)=0.0_r8
# endif
# ifdef WEC_VF
              FC(i,1)=(W(i,j,1)+W_stokes(i,j,1))*                       &
# else
              FC(i,1)=W(i,j,1)*                                         &
# endif
     &                (cff1*t(i,j,1,3,itrc)+                            &
     &                 cff2*t(i,j,2,3,itrc)-                            &
     &                 cff3*t(i,j,3,3,itrc))
# ifdef WEC_VF
              FC(i,N(ng)-1)=(W(i,j,N(ng)-1)+W_stokes(i,j,N(ng)-1))*     &
# else
              FC(i,N(ng)-1)=W(i,j,N(ng)-1)*                             &
# endif
     &                    (cff1*t(i,j,N(ng)  ,3,itrc)+                  &
     &                     cff2*t(i,j,N(ng)-1,3,itrc)-                  &
     &                     cff3*t(i,j,N(ng)-2,3,itrc))
              FC(i,N(ng))=0.0_r8
            END DO
          END IF VADV_FLUX
!
!  If MPDATA and cell-centered point source (LwSrc), augment the
!  intermediate tracer, Ta, with the tracer divergence. For other
!  advection schemes LwSrc is applied after the vertical advection
!  is completed (J. Wilkin).
!
          IF (LwSrc(ng)) THEN
            IF (Hadvection(itrc,ng)%MPDATA) THEN
              DO is=1,Nsrc(ng)
                Isrc=SOURCES(ng)%Isrc(is)
                Jsrc=SOURCES(ng)%Jsrc(is)
                IF (((Istr.le.Isrc).and.(Isrc.le.Iend+1)).and.          &
     &               ((Jstr.le.Jsrc).and.(Jsrc.le.Jend+1)).and.         &
     &               (j.eq.Jsrc)) THEN
                  DO k=1,N(ng)
                    cff=dt(ng)*pm(i,j)*pn(i,j)
                    IF (LtracerSrc(itrc,ng)) THEN
#  ifdef ONE_TRACER_SOURCE
                      cff3=SOURCES(ng)%Tsrc(itrc)
#  elif defined TWO_D_TRACER_SOURCE
                      cff3=SOURCES(ng)%Tsrc(is,itrc)
#  else
                      cff3=SOURCES(ng)%Tsrc(is,k,itrc)
#  endif
                    ELSE
                      cff3=t(Isrc,Jsrc,k,3,itrc)
                    END IF
                    Ta(Isrc,Jsrc,k,itrc)=Ta(Isrc,Jsrc,k,itrc)+          &
     &                                   cff*SOURCES(ng)%Qsrc(is,k)*    &
     &                                   cff3
                  END DO
                END IF
              END DO
            END IF
          END IF
!
!  If MPDATA, time-step vertical advection for intermediate diffusive
!  tracer, Ta (Tunits).
# ifdef DIAGNOSTICS_TS
!  Convert units of tracer diagnostic terms to Tunits.
# endif
!
          VADV_STEPPING : IF (Vadvection(itrc,ng)%MPDATA) THEN
            DO i=IstrUm2,Iendp2i
              CF(i,0)=dt(ng)*pm(i,j)*pn(i,j)
            END DO
            DO k=1,N(ng)
              DO i=IstrUm2,Iendp2i
                cff1=CF(i,0)*(FC(i,k)-FC(i,k-1))
                Ta(i,j,k,itrc)=(Ta(i,j,k,itrc)-cff1)*oHz(i,j,k)
# ifdef DIAGNOSTICS_TS
                Dvadv(i,j,k,itrc)=-cff1
# endif
              END DO
            END DO
!
!  OTHERWISE, time-step vertical advection term (m Tunits, or Tunits if
!  conservative, parabolic splines diffusion).
# ifdef DIAGNOSTICS_TS
!  Convert units of tracer diagnostic terms to Tunits.
# endif
!
          ELSE
            DO i=Istr,Iend
              CF(i,0)=dt(ng)*pm(i,j)*pn(i,j)
            END DO
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff1=CF(i,0)*(FC(i,k)-FC(i,k-1))
                t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-cff1
# ifdef SPLINES_VDIFF
                t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)*oHz(i,j,k)
# endif
# ifdef DIAGNOSTICS_TS
                DiaTwrk(i,j,k,itrc,iTvadv)=-cff1
                DO idiag=1,NDT
                  DiaTwrk(i,j,k,itrc,idiag)=DiaTwrk(i,j,k,itrc,idiag)*  &
     &                                      oHz(i,j,k)
                END DO
# endif
              END DO
            END DO
          END IF VADV_STEPPING
        END DO J_LOOP1
      END DO T_LOOP2
!
!-----------------------------------------------------------------------
!  Compute anti-diffusive velocities to corrected advected tracers
!  using MPDATA recursive method.  Notice that pipelined J-loop ended.
!-----------------------------------------------------------------------
!
      T_LOOP3 : DO itrc=1,NT(ng)
        MPDATA : IF ((Hadvection(itrc,ng)%MPDATA).and.                  &
     &               (Vadvection(itrc,ng)%MPDATA)) THEN
          CALL mpdata_adiff_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
# ifdef MASKING
     &                            rmask, umask, vmask,                  &
# endif
# ifdef WET_DRY
     &                            rmask_wet, umask_wet, vmask_wet,      &
# endif
     &                            pm, pn, omn, om_u, on_v,              &
     &                            z_r, oHz,                             &
     &                            Huon, Hvom, W,                        &
#  ifdef WEC_VF
     &                            W_stokes,                             &
#  endif
     &                            t(:,:,:,3,itrc),                      &
     &                            Ta(:,:,:,itrc),  Ua, Va, Wa)
!
!  Compute anti-diffusive corrected advection fluxes.
!
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend+1
                cff1=MAX(Ua(i,j,k),0.0_r8)
                cff2=MIN(Ua(i,j,k),0.0_r8)
                FX(i,j)=(cff1*Ta(i-1,j,k,itrc)+                         &
     &                   cff2*Ta(i  ,j,k,itrc))*                        &
     &                  0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)
              END DO
            END DO
            DO j=Jstr,Jend+1
              DO i=Istr,Iend
                cff1=MAX(Va(i,j,k),0.0_r8)
                cff2=MIN(Va(i,j,k),0.0_r8)
                FE(i,j)=(cff1*Ta(i,j-1,k,itrc)+                         &
     &                   cff2*Ta(i,j  ,k,itrc))*                        &
     &                  0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*om_v(i,j)
              END DO
            END DO
!
!  Time-step corrected horizontal advection (Tunits m).
!
            DO j=Jstr,Jend
              DO i=Istr,Iend
                cff=dt(ng)*pm(i,j)*pn(i,j)
                cff1=cff*(FX(i+1,j)-FX(i,j))
                cff2=cff*(FE(i,j+1)-FE(i,j))
                cff3=cff1+cff2
                t(i,j,k,nnew,itrc)=Ta(i,j,k,itrc)*Hz(i,j,k)-cff3
# ifdef DIAGNOSTICS_TS
                DiaTwrk(i,j,k,itrc,iTxadv)=DiaTwrk(i,j,k,itrc,iTxadv)-  &
     &                                     cff1
                DiaTwrk(i,j,k,itrc,iTyadv)=DiaTwrk(i,j,k,itrc,iTyadv)-  &
     &                                     cff2
                DiaTwrk(i,j,k,itrc,iThadv)=DiaTwrk(i,j,k,itrc,iThadv)-  &
     &                                     cff3
# endif
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
                FC(i,k)=cff1*Ta(i,j,k  ,itrc)+                          &
     &                  cff2*Ta(i,j,k+1,itrc)
              END DO
            END DO
            DO i=Istr,Iend
              FC(i,0)=0.0_r8
              FC(i,N(ng))=0.0_r8
            END DO
!
!  Time-step corrected vertical advection (Tunits).
# ifdef DIAGNOSTICS_TS
!  Convert units of tracer diagnostic terms to Tunits.
# endif
!
            DO i=Istr,Iend
              CF(i,0)=dt(ng)*pm(i,j)*pn(i,j)
            END DO
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff1=CF(i,0)*(FC(i,k)-FC(i,k-1))
                t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-cff1
# ifdef DIAGNOSTICS_TS
                DiaTwrk(i,j,k,itrc,iTvadv)=Dvadv(i,j,k,itrc)-           &
     &                                     cff1
                DO idiag=1,NDT
                  DiaTwrk(i,j,k,itrc,idiag)=DiaTwrk(i,j,k,itrc,idiag)*  &
     &                                      oHz(i,j,k)
                END DO
# endif
              END DO
            END DO
          END DO
        END IF MPDATA
      END DO T_LOOP3
!
!-----------------------------------------------------------------------
!  Add tracer divergence due to cell-centered (LwSrc) point sources.
!-----------------------------------------------------------------------
!
!  When LTracerSrc is .true. the inflowing concentration is Tsrc.
!  When LtracerSrc is .false. we add tracer mass to compensate for the
!  added volume to keep the receiving cell concentration unchanged.
!  J. Levin (Jupiter Intelligence Inc.) and J. Wilkin
!
      IF (LwSrc(ng)) THEN
        DO itrc=1,NT(ng)
          IF (.not.((Hadvection(itrc,ng)%MPDATA).and.                   &
     &              (Vadvection(itrc,ng)%MPDATA))) THEN
            DO is=1,Nsrc(ng)
              Isrc=SOURCES(ng)%Isrc(is)
              Jsrc=SOURCES(ng)%Jsrc(is)
              IF (((Istr.le.Isrc).and.(Isrc.le.Iend+1)).and.            &
     &            ((Jstr.le.Jsrc).and.(Jsrc.le.Jend+1))) THEN
                DO k=1,N(ng)
                  cff=dt(ng)*pm(i,j)*pn(i,j)
# ifdef SPLINES_VDIFF
                  cff=cff*oHz(Isrc,Jsrc,k)
# endif
                  IF (LtracerSrc(itrc,ng)) THEN
# ifdef ONE_TRACER_SOURCE
                    cff3=SOURCES(ng)%Tsrc(itrc)
# elif defined TWO_D_TRACER_SOURCE
                    cff3=SOURCES(ng)%Tsrc(is,itrc)
# else
                    cff3=SOURCES(ng)%Tsrc(is,k,itrc)
# endif
                  ELSE
                    cff3=t(Isrc,Jsrc,k,3,itrc)
                  END IF
                  t(Isrc,Jsrc,k,nnew,itrc)=t(Isrc,Jsrc,k,nnew,itrc)+    &
     &                                     cff*SOURCES(ng)%Qsrc(is,k)*  &
     &                                     cff3
                END DO
              END IF
            END DO
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Time-step vertical diffusion term.
!-----------------------------------------------------------------------
!
      J_LOOP2 : DO j=Jstr,Jend                  ! start pipelined J-loop
#  ifdef TS_VAR
        DO itrc=1,NAT
#  else
        DO itrc=1,NT(ng)
#  endif
          ltrc=MIN(NAT,itrc)

# ifdef SPLINES_VDIFF
          IF (.not.((Hadvection(itrc,ng)%MPDATA).and.                   &
     &              (Vadvection(itrc,ng)%MPDATA))) THEN
!
!  Use conservative, parabolic spline reconstruction of vertical
!  diffusion derivatives.  Then, time step vertical diffusion term
!  implicitly.
!
            cff1=1.0_r8/6.0_r8
            DO k=1,N(ng)-1
              DO i=Istr,Iend
                FC(i,k)=cff1*Hz(i,j,k  )-                               &
     &                  dt(ng)*Akt(i,j,k-1,ltrc)*oHz(i,j,k  )
                CF(i,k)=cff1*Hz(i,j,k+1)-                               &
     &                  dt(ng)*Akt(i,j,k+1,ltrc)*oHz(i,j,k+1)
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
                BC(i,k)=cff1*(Hz(i,j,k)+Hz(i,j,k+1))+                   &
     &                  dt(ng)*Akt(i,j,k,ltrc)*(oHz(i,j,k)+oHz(i,j,k+1))
                cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
                CF(i,k)=cff*CF(i,k)
                DC(i,k)=cff*(t(i,j,k+1,nnew,itrc)-t(i,j,k,nnew,itrc)-   &
     &                       FC(i,k)*DC(i,k-1))
              END DO
            END DO
!
!  Backward substitution.
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
            DO k=1,N(ng)
              DO i=Istr,Iend
                DC(i,k)=DC(i,k)*Akt(i,j,k,ltrc)
                cff1=dt(ng)*oHz(i,j,k)*(DC(i,k)-DC(i,k-1))
                t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+cff1
#  ifdef DIAGNOSTICS_TS
                DiaTwrk(i,j,k,itrc,iTvdif)=DiaTwrk(i,j,k,itrc,iTvdif)+  &
     &                                     cff1
#  endif
              END DO
            END DO
          ELSE
# endif
!
!  Compute off-diagonal coefficients FC [lambda*dt*Akt/Hz] for the
!  implicit vertical diffusion terms at future time step, located
!  at horizontal RHO-points and vertical W-points.
!  Also set FC at the top and bottom levels.
!
            cff=-dt(ng)*lambda
            DO k=1,N(ng)-1
              DO i=Istr,Iend
                cff1=1.0_r8/(z_r(i,j,k+1)-z_r(i,j,k))
                FC(i,k)=cff*cff1*Akt(i,j,k,ltrc)
              END DO
            END DO
            DO i=Istr,Iend
              FC(i,0)=0.0_r8
              FC(i,N(ng))=0.0_r8
            END DO
!
!  Compute diagonal matrix coefficients BC and load right-hand-side
!  terms for the tracer equation into DC.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                BC(i,k)=Hz(i,j,k)-FC(i,k)-FC(i,k-1)
                DC(i,k)=t(i,j,k,nnew,itrc)
              END DO
            END DO
!
!  Solve the tridiagonal system.
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
!
            DO i=Istr,Iend
# ifdef DIAGNOSTICS_TS
               cff1=t(i,j,N(ng),nnew,itrc)*oHz(i,j,N(ng))
# endif
               DC(i,N(ng))=(DC(i,N(ng))-FC(i,N(ng)-1)*DC(i,N(ng)-1))/   &
     &                     (BC(i,N(ng))-FC(i,N(ng)-1)*CF(i,N(ng)-1))
               t(i,j,N(ng),nnew,itrc)=DC(i,N(ng))
# ifdef DIAGNOSTICS_TS
               DiaTwrk(i,j,N(ng),itrc,iTvdif)=                          &
     &                              DiaTwrk(i,j,N(ng),itrc,iTvdif)+     &
     &                              t(i,j,N(ng),nnew,itrc)-cff1
# endif
            END DO
            DO k=N(ng)-1,1,-1
              DO i=Istr,Iend
# ifdef DIAGNOSTICS_TS
                cff1=t(i,j,k,nnew,itrc)*oHz(i,j,k)
# endif
                DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
                t(i,j,k,nnew,itrc)=DC(i,k)
# ifdef DIAGNOSTICS_TS
                DiaTwrk(i,j,k,itrc,iTvdif)=DiaTwrk(i,j,k,itrc,iTvdif)+  &
     &                                     t(i,j,k,nnew,itrc)-cff1
# endif
              END DO
            END DO
# ifdef SPLINES_VDIFF
          END IF
# endif
        END DO
      END DO J_LOOP2

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
              IF ((Hadvection(inert(itrc),ng)%MPDATA).and.              &
     &            (Vadvection(inert(itrc),ng)%MPDATA)) THEN
                t(i,j,k,nnew,iage)=t(i,j,k,nnew,iage)+                  &
     &                             dt(ng)*t(i,j,k,nnew,inert(itrc))
              ELSE
                t(i,j,k,nnew,iage)=t(i,j,k,nnew,iage)+                  &
     &                             dt(ng)*t(i,j,k,3,inert(itrc))
              END IF
            END DO
          END DO
        END DO
      END DO
# endif
# if defined TRC_PSOURCE
!
!  Apply tracers point sources to the horizontal advection terms.
!
        DO it=1,NPT
          itrc=inert(it)
          DO is=1,Nsrcpt
            i=Isrcpt(is)
            j=Jsrcpt(is)
            DO k=1,N(ng)
              IF ((Hadvection(inert(itrc),ng)%MPDATA).and.              &
     &            (Vadvection(inert(itrc),ng)%MPDATA)) THEN
                IF (((Istr-2.le.i).and.(i.le.Iend+2)).and.              &
     &              ((Jstr-2.le.j).and.(j.le.Jend+2))) THEN
              ELSE
                IF (((Istr.le.i).and.(i.le.Iend+1)).and.                &
     &              ((Jstr.le.j).and.(j.le.Jend))) THEN
              END IF
!
                IF (Lsrcpt(is,itrc)) THEN
                  t(i,j,k,3,itrc)=0.
                  t(i,j,k,nstp,itrc)=Tsrcpt(is,k,itrc)
                  t(i,j,k,nnew,itrc)=Tsrcpt(is,k,itrc)
                END IF
              END IF
            END DO
          END DO
        END DO
# endif
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
        CALL t3dbc_tile (ng, tile, itrc, ic,                            &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp, nnew,                                    &
     &                   t)
!
!  Nudge towards tracer climatology.
!
        IF (LtracerCLM(itrc,ng).and.LnudgeTCLM(itrc,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+                  &
     &                             dt(ng)*                              &
     &                             CLIMA(ng)%Tnudgcof(i,j,k,ic)*        &
     &                             (CLIMA(ng)%tclm(i,j,k,ic)-           &
     &                              t(i,j,k,nnew,itrc))
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
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)*rmask(i,j)
            END DO
          END DO
        END DO
# endif
# ifdef DIAGNOSTICS_TS
!
!  Compute time-rate-of-change diagnostic term.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              DiaTwrk(i,j,k,itrc,iTrate)=t(i,j,k,nnew,itrc)-            &
     &                                   DiaTwrk(i,j,k,itrc,iTrate)
!!            DiaTwrk(i,j,k,itrc,iTrate)=t(i,j,k,nnew,itrc)-            &
!!   &                                   t(i,j,k,nstp,itrc)
            END DO
          END DO
        END DO
# endif
!
!  Apply periodic boundary conditions.
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            t(:,:,:,nnew,itrc))
        END IF
      END DO
# ifdef TS_VAR
!
!  Compute tracer variance.
!
      itrc=inert(1)
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          DO k=1,N(ng)
#  if !defined SPLINES_VDIFF || defined TS_MPDATA
            t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)*oHz(i,j,k)
            t(i,j,k,nnew,itrc+1)=t(i,j,k,nnew,itrc+1)*oHz(i,j,k)
#  endif
            t(i,j,k,nnew,itrc+2)=(t(i,j,k,nnew,itrc+1)-                 &
     &                           t(i,j,k,nnew,itrc)*t(i,j,k,nnew,itrc))/&
     &                           dt(ng)
          END DO
          DO k=1,N(ng)-1
            cff=(t(i,j,k+1,nnew,isalt)-t(i,j,k,nnew,isalt))/            &
     &                         (z_r(i,j,k+1)-z_r(i,j,k))
            AKr(i,j,k)=2.0_r8*cff*cff*(Akt(i,j,k,itemp)-                &
     &                                 Akt_bak(itemp,ng))
          END DO
        END DO
      END DO
# endif
# ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    t(:,:,:,nnew,:))
# endif
# if defined FLOATS && defined FLOAT_VWALK
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
!  Apply periodic boundary conditions.
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
!
!-----------------------------------------------------------------------
!  If applicable, deallocate local arrays.
!-----------------------------------------------------------------------
!
      IF (Lmpdata) THEN
# ifdef DIAGNOSTICS_TS
        IF (allocated(Dhadv)) deallocate (Dhadv)
        IF (allocated(Dvadv)) deallocate (Dvadv)
# endif
        IF (allocated(Ta))    deallocate (Ta)
        IF (allocated(Ua))    deallocate (Ua)
        IF (allocated(Va))    deallocate (Va)
        IF (allocated(Wa))    deallocate (Wa)
      END IF
!
      RETURN
      END SUBROUTINE step3d_t_tile
#endif
      END MODULE step3d_t_mod
