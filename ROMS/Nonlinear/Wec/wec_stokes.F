#include "cppdefs.h"
      MODULE wec_stokes_mod
#if defined SOLVE3D && defined WEC
!
!svn $Id: wec_stokes.F 1428 2008-03-12 13:07:21Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!                                                   Nirnimesh Kumar    !
!================================================== John C. Warner ====!
!                                                                      !
!  This routine computes the stokes transport terms                    !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: wec_stokes
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wec_stokes (ng, tile)
!***********************************************************************
!
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_param
# if defined DIAGNOSTICS_UV
      USE mod_diags
# endif
!
      integer, intent(in) :: ng, tile
# include "tile.h"
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 21)
# endif
      CALL wec_stokes_tile (ng, tile, LBi, UBi, LBj, UBj, N(ng),        &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                            GRID(ng) % angler,                    &
     &                            GRID(ng) % h,                         &
# ifdef MASKING
     &                            GRID(ng) % rmask,                     &
     &                            GRID(ng) % umask,                     &
     &                            GRID(ng) % vmask,                     &
# endif
# ifdef WET_DRY
     &                            GRID(ng) % umask_wet,                 &
     &                            GRID(ng) % vmask_wet,                 &
# endif
     &                            GRID(ng) % on_u,                      &
     &                            GRID(ng) % on_v,                      &
     &                            GRID(ng) % om_u,                      &
     &                            GRID(ng) % om_v,                      &
     &                            GRID(ng) % Hz,                        &
     &                            GRID(ng) % z_r,                       &
     &                            GRID(ng) % z_w,                       &
     &                            FORCES(ng) % Hwave,                   &
     &                            FORCES(ng) % Dwave,                   &
     &                            FORCES(ng) % Lwave,                   &
# if defined WEC_ROLLER
     &                            FORCES(ng) % rollA,                   &
# endif
     &                            OCEAN(ng) % zeta,                     &
# ifdef WEC_MELLOR
     &                            OCEAN(ng) % rulag2d,                  &
     &                            OCEAN(ng) % rvlag2d,                  &
     &                            OCEAN(ng) % rulag3d,                  &
     &                            OCEAN(ng) % rvlag3d,                  &
# endif
     &                            OCEAN(ng) % ubar_stokes,              &
     &                            OCEAN(ng) % vbar_stokes,              &
     &                            OCEAN(ng) % u_stokes,                 &
     &                            OCEAN(ng) % v_stokes,                 &
     &                            OCEAN(ng) % W_stokes)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 21)
# endif
      RETURN
      END SUBROUTINE wec_stokes
!
!***********************************************************************
      SUBROUTINE wec_stokes_tile (ng, tile, LBi, UBi, LBj, UBj, UBk,    &
     &                                  IminS, ImaxS, JminS, JmaxS,     &
     &                                  angler, h,                      &
# ifdef MASKING
     &                                  rmask, umask, vmask,            &
# endif
# ifdef WET_DRY
     &                                  umask_wet, vmask_wet,           &
# endif
     &                                  on_u, on_v, om_u, om_v,         &
     &                                  Hz, z_r, z_w,                   &
     &                                  Hwave, Dwave, Lwave,            &
# if defined WEC_ROLLER
     &                                  rollA,                          &
# endif
     &                                  zeta,                           &
# ifdef WEC_MELLOR
     &                                  rulag2d, rvlag2d,               &
     &                                  rulag3d, rvlag3d,               &
# endif
     &                                  ubar_stokes, vbar_stokes,       &
     &                                  u_stokes, v_stokes, W_stokes)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE exchange_2d_mod
      USE exchange_3d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
# endif
      USE bc_2d_mod
      USE bc_3d_mod
      USE us2dbc_mod, ONLY : us2dbc_tile
      USE vs2dbc_mod, ONLY : vs2dbc_tile
# ifdef SOLVE3D
      USE us3dbc_mod, ONLY : us3dbc_tile
      USE vs3dbc_mod, ONLY : vs3dbc_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: angler(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: Hwave(LBi:,LBj:)
      real(r8), intent(in) :: Dwave(LBi:,LBj:)
      real(r8), intent(in) :: Lwave(LBi:,LBj:)
#  if defined WEC_ROLLER
      real(r8), intent(in) :: rollA(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#  ifdef WEC_MELLOR
      real(r8), intent(inout) :: rulag2d(LBi:,LBj:)
      real(r8), intent(inout) :: rvlag2d(LBi:,LBj:)
      real(r8), intent(inout) :: rulag3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: rvlag3d(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: ubar_stokes(LBi:,LBj:)
      real(r8), intent(inout) :: vbar_stokes(LBi:,LBj:)
      real(r8), intent(inout) :: u_stokes(LBi:,LBj:,:)
      real(r8), intent(inout) :: v_stokes(LBi:,LBj:,:)
      real(r8), intent(inout) :: W_stokes(LBi:,LBj:,0:)
# else
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: Hwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Dwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Lwave(LBi:UBi,LBj:UBj)
#  if defined WEC_ROLLER
      real(r8), intent(in) :: rollA(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#  ifdef WEC_MELLOR
      real(r8), intent(inout) :: rulag2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rvlag2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rulag3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: rvlag3d(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(inout) :: ubar_stokes(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vbar_stokes(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: u_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: v_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: W_stokes(LBi:UBi,LBj:UBj,0:N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8) :: cff, cff2, cff3, cff4, cff5, cff6, cff7
      real(r8) :: fac1, fac2, fac3, ofac3
      real(r8), parameter :: eps = 1.0E-14_r8
      real(r8), parameter :: kDmax = 200.0_r8
      real(r8), parameter :: Lwave_min = 1.0_r8
      real(r8), dimension(IminS:ImaxS) :: wrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Dstp
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: kD
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wavec
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waven
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: owaven
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wavenx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waveny
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waveE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: sigma
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Huons
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hvoms
# include "set_bounds.h"
      fac1=1.0_r8/dt(ng)
      fac2=1.0_r8/g
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
!
!  Compute total depth
!
          Dstp(i,j)=z_w(i,j,N(ng))-z_w(i,j,0)
!
!  Compute wave amplitude (0.5*Hrms), wave number, intrinsic frequency.
!
          waven(i,j)=2.0_r8*pi/MAX(Lwave(i,j),Lwave_min)
          owaven(i,j)=1.0_r8/waven(i,j)
          cff=1.5_r8*pi-Dwave(i,j)-angler(i,j)
          wavenx(i,j)=waven(i,j)*COS(cff)
          waveny(i,j)=waven(i,j)*SIN(cff)
          sigma(i,j)=SQRT(g*waven(i,j)*TANH(waven(i,j)*Dstp(i,j)))
          waveE(i,j)=0.0625_r8*g*Hwave(i,j)*Hwave(i,j)
!
!  Compute wave celerity and kD
!
          kD(i,j)=MIN(waven(i,j)*Dstp(i,j),kDmax)
          wavec(i,j)=SQRT(g*owaven(i,j)*TANH(kD(i,j)))
        END DO
      END DO
!
!---------------------------------------------------------------------------
! Stokes velocities.
!---------------------------------------------------------------------------
!
!  Compute u-stokes velocities, Eqn. 2.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
# ifdef WEC_MELLOR
          cff=fac1*om_u(i,j)*on_u(i,j)
# endif
          cff2=(waveE(i-1,j)+waveE(i,j))
# if defined WEC_ROLLER
          cff2=cff2+(rollA(i-1,j)*sigma(i-1,j)+                        &
     &               rollA(i,j)*sigma(i,j))
# endif
          cff3=(kD(i-1,j)+kD(i,j))
!         cff3=(waven(i-1,j)+waven(i,j))*(Dstp(i-1,j)+Dstp(i,j))
          cff4=wavenx(i-1,j)+wavenx(i,j)
          cff5=wavec(i-1,j)+wavec(i,j)
          fac3=Dstp(i-1,j)+Dstp(i,j)
          ofac3=1.0_r8/fac3
          DO k=1,N(ng)
            cff6=-1.0_r8+((z_w(i-1,j,k)+z_w(i,j,k))-                   &
     &                   (z_w(i-1,j,0)+z_w(i,j,0)))*ofac3
            cff7=-1.0_r8+((z_w(i-1,j,k-1)+z_w(i,j,k-1))-               &
     &                   (z_w(i-1,j,0)+z_w(i,j,0)))*ofac3
# ifdef WEC_MELLOR
!
!  Store old value to compute tendency term.
!
            rulag3d(i,j,k)=u_stokes(i,j,k)
# endif
            u_stokes(i,j,k)=0.25_r8*cff2*cff4*cff5*fac2/               &
     &                      (fac3*(cff6-cff7))/                        &
     &                      (SINH(0.5_r8*cff3)**2.0_r8)*               &
     &                      SINH(0.5_r8*cff3*(cff6-cff7))*             &
     &                      COSH(0.5_r8*cff3*(cff6+cff7+2.0_r8))
# ifdef MASKING
            u_stokes(i,j,k)=u_stokes(i,j,k)*umask(i,j)
# endif
# ifdef WET_DRY
            u_stokes(i,j,k)=u_stokes(i,j,k)*umask_wet(i,j)
# endif
# ifdef WEC_MELLOR
!
!  Finalize computation of stokes tendency term.
!
            rulag3d(i,j,k)=0.5_r8*cff*                                 &
     &                     (Hz(i,j,k)+Hz(i-1,j,k))*                    &
     &                     (u_stokes(i,j,k)-rulag3d(i,j,k))
# endif
          END DO
        END DO
      END DO
!
!  Compute v-stokes velocity, Eqn. 2.
!
      DO j=JstrV,Jend
        DO i=Istr,Iend
# ifdef WEC_MELLOR
          cff=fac1*om_v(i,j)*on_v(i,j)
# endif
          cff2=(waveE(i,j-1)+waveE(i,j))
# if defined WEC_ROLLER
          cff2=cff2+(rollA(i,j-1)*sigma(i,j-1)+                        &
     &               rollA(i,j)*sigma(i,j))
# endif
          cff3=(kD(i,j-1)+kD(i,j))
!         cff3=(waven(i,j-1)+waven(i,j))*(Dstp(i,j-1)+Dstp(i,j))
          cff4=waveny(i,j-1)+waveny(i,j)
          cff5=wavec(i,j-1)+wavec(i,j)
          fac3=Dstp(i,j-1)+Dstp(i,j)
          ofac3=1.0_r8/fac3
          DO k=1,N(ng)
            cff6=-1.0_r8+((z_w(i,j-1,k  )+z_w(i,j,k))-                  &
     &                    (z_w(i,j-1,0  )+z_w(i,j,0)))*ofac3
            cff7=-1.0_r8+((z_w(i,j-1,k-1)+z_w(i,j,k-1))-                &
     &                    (z_w(i,j-1,0  )+z_w(i,j,0)))*ofac3
# ifdef WEC_MELLOR
!
!  Store old value to compute tendency term.
!
            rvlag3d(i,j,k)=v_stokes(i,j,k)
# endif
            v_stokes(i,j,k)=0.25_r8*cff2*cff4*cff5*fac2/                &
     &                      (fac3*(cff6-cff7))/                         &
     &                      (SINH(0.5_r8*cff3)**2.0_r8)*                &
     &                      SINH(0.5_r8*cff3*(cff6-cff7))*              &
     &                      COSH(0.5_r8*cff3*(cff6+cff7+2.0_r8))
# ifdef MASKING
            v_stokes(i,j,k)=v_stokes(i,j,k)*vmask(i,j)
# endif
# ifdef WET_DRY
            v_stokes(i,j,k)=v_stokes(i,j,k)*vmask_wet(i,j)
# endif
# ifdef WEC_MELLOR
!
!  Finalize computation of stokes tendency term.
!
            rvlag3d(i,j,k)=0.5_r8*cff*                                  &
     &                     (Hz(i,j,k)+Hz(i,j-1,k))*                     &
     &                     (v_stokes(i,j,k)-rvlag3d(i,j,k))
# endif
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
! Set lateral boundary conditions.
!-----------------------------------------------------------------------
!
      CALL us3dbc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, N(ng),                      &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  u_stokes)
      CALL vs3dbc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, N(ng),                      &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  v_stokes)
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!
      CALL exchange_u3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        u_stokes(:,:,:))
      CALL exchange_v3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        v_stokes(:,:,:))
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    u_stokes, v_stokes)
# endif
!
!  Compute vertical stokes velocity, Eqn. 31.
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
             Huons(i,j)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*                 &
     &                  u_stokes(i,j,k)*on_u(i,j)
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
             Hvoms(i,j)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*                 &
     &                  v_stokes(i,j,k)*om_v(i,j)
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend
             W_stokes(i,j,k)=W_stokes(i,j,k-1)-                         &
     &                      (Huons(i+1,j)-Huons(i,j)+                   &
     &                       Hvoms(i,j+1)-Hvoms(i,j))
          END DO
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=Istr,Iend
          wrk(i)=W_stokes(i,j,N(ng))/(z_w(i,j,N(ng))-z_w(i,j,0))
        END DO
        DO k=N(ng)-1,1,-1
         DO i=Istr,Iend
            W_stokes(i,j,k)=W_stokes(i,j,k)-                            &
     &                      wrk(i)*(z_w(i,j,k)-z_w(i,j,0))
          END DO
        END DO
        DO i=Istr,Iend
          W_stokes(i,j,N(ng))=0.0_r8
        END DO
      END DO
!
!  For a 3D application, compute associated 2D fields by taking the
!  vertical integral of 3D fields.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff5=0.5_r8*(Hz(i-1,j,1)+Hz(i,j,1))
          ubar_stokes(i,j)=cff5*u_stokes(i,j,1)
# ifdef WEC_MELLOR
          rulag2d(i,j)=cff5*rulag3d(i,j,1)
# endif
          DO k=2,N(ng)
            cff5=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))
            ubar_stokes(i,j)=ubar_stokes(i,j)+cff5*u_stokes(i,j,k)
# ifdef WEC_MELLOR
            rulag2d(i,j)=rulag2d(i,j)+cff5*rulag3d(i,j,k)
# endif
          END DO
          cff4=2.0_r8/(Dstp(i-1,j)+Dstp(i,j))
          ubar_stokes(i,j)=ubar_stokes(i,j)*cff4
# ifdef WEC_MELLOR
          rulag2d(i,j)=rulag2d(i,j)*cff4
# endif
        END DO
      END DO
      DO i=Istr,Iend
        DO j=JstrV,Jend
          cff5=0.5_r8*(Hz(i,j-1,1)+Hz(i,j,1))
          vbar_stokes(i,j)=cff5*v_stokes(i,j,1)
# ifdef WEC_MELLOR
          rvlag2d(i,j)=cff5*rvlag3d(i,j,1)
# endif
          DO k=2,N(ng)
            cff5=0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))
            vbar_stokes(i,j)=vbar_stokes(i,j)+cff5*v_stokes(i,j,k)
# ifdef WEC_MELLOR
            rvlag2d(i,j)=rvlag2d(i,j)+cff5*rvlag3d(i,j,k)
# endif
          END DO
          cff4=2.0_r8/(Dstp(i,j-1)+Dstp(i,j))
          vbar_stokes(i,j)=vbar_stokes(i,j)*cff4
# ifdef WEC_MELLOR
          rvlag2d(i,j)=rvlag2d(i,j)*cff4
# endif
        END DO
      END DO
!
!  Apply boundary conditions.
!
        CALL us2dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    ubar_stokes)
        CALL vs2dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    vbar_stokes)
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ubar_stokes(:,:))
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        vbar_stokes(:,:))
      END IF
      CALL bc_w3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  W_stokes)
# ifdef WEC_MELLOR
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  rulag2d)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  rvlag2d)
      CALL bc_u3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  rulag3d)
      CALL bc_v3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  rvlag3d)
# endif
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    ubar_stokes, vbar_stokes)
#  ifdef WEC_MELLOR
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rulag2d, rvlag2d)
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rulag3d,  rvlag3d)
#  endif
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    W_stokes)
# endif
      RETURN
      END SUBROUTINE wec_stokes_tile
#endif
      END MODULE wec_stokes_mod
