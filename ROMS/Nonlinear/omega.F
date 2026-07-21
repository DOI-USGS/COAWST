#include "cppdefs.h"
      MODULE omega_mod
#ifdef SOLVE3D
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                            Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine computes S-coordinate vertical velocity (m^3/s),       !
!                                                                      !
!                  W=[Hz/(m*n)]*omega,                                 !
!                                                                      !
!  diagnostically at horizontal RHO-points and vertical W-points.      !
!                                                                      !
# ifdef OMEGA_IMPLICIT
!  Added implicit vertical advection following Shchepetkin (2015).     !
!  The vertical velocity is split into explicit (W) and implicit (Wi)  !
!  parts that adjust automatically to the local flow conditions based  !
!  on the Courant number for stability, allowing larger time steps.    !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Shchepetkin, A.F., 2015: An adaptive, Courant-number-dependent    !
!      implicit scheme for vertical advection in oceanic modeling,     !
!      Ocean Modelling, 91, 38-69, doi: 10.1016/j.ocemod.2015.03.006   !
!                                                                      !
# endif
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: omega, scale_omega
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE omega (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
# if defined SEDIMENT && defined SED_MORPH
      USE mod_sedbed
      USE mod_stepping
# endif
# ifdef EMINUSP_SHIMA
      USE mod_forces
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 13, __LINE__, MyFile)
# endif
      CALL omega_tile (ng, tile, model,                                 &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
# if defined SEDIMENT && defined SED_MORPH
     &                 nstp(ng), nnew(ng),                              &
     &                 GRID(ng) % omn,                                  &
     &                 SEDBED(ng) % bed_thick,                          &
# endif
# ifdef EMINUSP_SHIMA
     &                 FORCES(ng) % stflux,                             &
     &                 GRID(ng) % omn,                                  &
# endif
     &                 GRID(ng) % Huon,                                 &
     &                 GRID(ng) % Hvom,                                 &
# ifdef OMEGA_IMPLICIT
     &                 GRID(ng) % pm,                                   &
     &                 GRID(ng) % pn,                                   &
# endif
     &                 GRID(ng) % z_w,                                  &
# if defined WEC_VF
     &                 OCEAN(ng) % W_stokes,                            &
# endif
# ifdef OMEGA_IMPLICIT
     &                 OCEAN(ng) % Wi,                                  &
# endif
     &                 OCEAN(ng) % W)
# ifdef PROFILE
      CALL wclock_off (ng, model, 13, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE omega
!
!***********************************************************************
      SUBROUTINE omega_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
# if defined SEDIMENT && defined SED_MORPH
     &                       nstp, nnew,                                &
     &                       omn, bed_thick,                            &
# endif
# ifdef EMINUSP_SHIMA
     &                       stflux,                                    &
     &                       omn,                                       &
# endif
     &                       Huon, Hvom,                                &
# ifdef OMEGA_IMPLICIT
     &                       pm, pn,                                    &
# endif
     &                       z_w,                                       &
# if defined WEC_VF
     &                       W_stokes,                                  &
# endif
# ifdef OMEGA_IMPLICIT
     &                       Wi,                                        &
# endif
     &                       W)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sources
!
      USE bc_3d_mod, ONLY : bc_w3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
# if defined SEDIMENT && defined SED_MORPH
      integer, intent(in) :: nstp, nnew
# endif
!
# ifdef ASSUMED_SHAPE
#  if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in) :: bed_thick(LBi:,LBj:,:)
#  endif
#  ifdef EMINUSP_SHIMA
      real(r8), intent(in) :: stflux(LBi:,LBj:,:)
      real(r8), intent(in) :: omn(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
#  ifdef OMEGA_IMPLICIT
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
#  if defined WEC_VF
      real(r8), intent(in) :: W_stokes(LBi:,LBj:,0:)
#  endif
#  ifdef OMEGA_IMPLICIT
      real(r8), intent(out) :: Wi(LBi:,LBj:,0:)
#  endif
      real(r8), intent(out) :: W(LBi:,LBj:,0:)

# else

#  if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bed_thick(LBi:UBi,LBj:UBj,3)
#  endif
#  ifdef EMINUSP_SHIMA
      real(r8), intent(in) :: stflux(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
#  ifdef OMEGA_IMPLICIT
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#  if defined WEC_VF
      real(r8), intent(in) :: W_stokes(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
#  ifdef OMEGA_IMPLICIT
      real(r8), intent(out) :: Wi(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
      real(r8), intent(out) :: W(LBi:UBi,LBj:UBj,0:N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, ii, is, j, jj, k
# if defined SEDIMENT && defined SED_MORPH
      real(r8) :: cff1
# endif
# ifdef EMINUSP_SHIMA
      real(r8) :: fac2
# endif
      real(r8), dimension(IminS:ImaxS) :: wrk
# ifdef OMEGA_IMPLICIT
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: Cu_adv
      real(r8) :: cff
      real(r8) :: cw, c2d, dh, cw_max, cw_max2, cw_min
!
      real(r8), parameter :: amax = 0.75_r8    ! Maximum Courant number
      real(r8), parameter :: amin = 0.60_r8    ! Minimum Courant number
      real(r8), parameter :: cmnx_ratio = amin/amax
      real(r8), parameter :: cutoff = 2.0_r8-amin/amax
      real(r8), parameter :: r4cmx = 1.0_r8/(4.0_r8-4.0_r8*amin/amax)
# endif

# include "set_bounds.h"
!
!------------------------------------------------------------------------
!  Vertically integrate horizontal mass flux divergence:
!
!  Starting with zero vertical velocity at the bottom, integrate
!  from the bottom (k=0) to the free-surface (k=N).  The w(:,:,N(ng))
!  contains the vertical velocity at the free-surface, d(zeta)/d(t).
!  Notice that barotropic mass flux divergence is not used directly.
!------------------------------------------------------------------------
!
# if defined SEDIMENT && defined SED_MORPH
!  For sediment bed change, we need to include the mass change of
!  water volume due to change of the sea floor. This is similar to
!  the LwSrc point source approach.
!
      cff1=1.0_r8/(dt(ng)*N(ng))
# endif
# ifdef EMINUSP_SHIMA
      fac2=1.0_r8/N(ng)
# endif
      DO j=Jstr,Jend
        DO i=Istr,Iend
          W(i,j,0)=0.0_r8
# if defined SEDIMENT && defined SED_MORPH
          wrk(i)=cff1*(bed_thick(i,j,nstp)-                             &
     &                 bed_thick(i,j,nnew))*omn(i,j)
# endif
        END DO
        DO k=1,N(ng)
          DO i=Istr,Iend
            W(i,j,k)=W(i,j,k-1)-                                        &
# if defined SEDIMENT && defined SED_MORPH
     &               wrk(i)-                                            &
# endif
# ifdef EMINUSP_SHIMA
     &                stflux(i,j,isalt)*omn(i,j)*fac2-                  &
# endif
     &               (Huon(i+1,j,k)-Huon(i,j,k)+                        &
     &                Hvom(i,j+1,k)-Hvom(i,j,k))

# ifdef OMEGA_IMPLICIT
!
!  Compute the horizontal Courant number.
!
            Cu_adv(i,k)=MAX(Huon(i+1,j  ,k),0.0_r8)-                    &
     &                  MIN(Huon(i  ,j  ,k),0.0_r8)+                    &
     &                  MAX(Hvom(i  ,j+1,k),0.0_r8)-                    &
     &                  MIN(Hvom(i  ,j  ,k),0.0_r8)
# endif
          END DO
        END DO
!
!  Apply mass point sources (volume vertical influx), if any.
!
!  Overwrite W(Isrc,Jsrc,k) with the same divergence of Huon,Hvom as
!  above but add in point source Qsrc(k) and reaccumulate the vertical
!  sum to obtain the correct net Qbar given in user input - J. Levin
!  (Jupiter Intelligence Inc.) and J. Wilkin
!
!    Dsrc(is) = 2,  flow across grid cell w-face (positive or negative)
!
        IF (LwSrc(ng)) THEN
          DO is=1,Nsrc(ng)
            IF (INT(SOURCES(ng)%Dsrc(is)).eq.2) THEN
              ii=SOURCES(ng)%Isrc(is)
              jj=SOURCES(ng)%Jsrc(is)
              IF (((IstrR.le.ii).and.(ii.le.IendR)).and.                &
     &            ((JstrR.le.jj).and.(jj.le.JendR)).and.                &
     &            (j.eq.jj)) THEN
# if defined SEDIMENT && defined SED_MORPH
                wrk(ii)=cff1*(bed_thick(ii,jj,nstp)-                    &
     &                        bed_thick(ii,jj,nnew))*omn(ii,jj)
# endif
                DO k=1,N(ng)
                  W(ii,jj,k)=W(ii,jj,k-1)-                              &
# if defined SEDIMENT && defined SED_MORPH
     &                       wrk(ii)-                                   &
# endif
# ifdef EMINUSP_SHIMA
     &                       stflux(ii,jj,isalt)*omn(ii,jj)*fac2-       &
# endif
     &                       (Huon(ii+1,jj,k)-Huon(ii,jj,k)+            &
     &                        Hvom(ii,jj+1,k)-Hvom(ii,jj,k))+           &
     &                       SOURCES(ng)%Qsrc(is,k)
                END DO
              END IF
            END IF
          END DO
        END IF
!
        DO i=Istr,Iend
          wrk(i)=W(i,j,N(ng))/(z_w(i,j,N(ng))-z_w(i,j,0))
# ifdef OMEGA_IMPLICIT
          Cu_adv(i,0)=dt(ng)*pm(i,j)*pn(i,j)
# endif
        END DO
!
!  In order to insure zero vertical velocity at the free-surface,
!  subtract the vertical velocities of the moving S-coordinates
!  isosurfaces. These isosurfaces are proportional to d(zeta)/d(t).
!  The proportionally coefficients are a linear function of the
!  S-coordinate with zero value at the bottom (k=0) and unity at
!  the free-surface (k=N).
!
        DO k=N(ng)-1,1,-1
          DO i=Istr,Iend
            W(i,j,k)=W(i,j,k)-                                          &
# if defined WEC_VF
     &               W_stokes(i,j,k)-                                   &
# endif
     &               wrk(i)*(z_w(i,j,k)-z_w(i,j,0))

# ifdef OMEGA_IMPLICIT
!
!  Determine implicit part "Wi" of vertical advection: "W"  becomes the
!  explicit part "We".
!
            Wi(i,j,k)=W(i,j,k)
            IF (Wi(i,j,k).ge.0.0_r8) THEN     ! Three different variants
              c2d=Cu_adv(i,k)                 ! for computing 2D Courant
              dh=z_w(i,j,k)-z_w(i,j,k-1)      ! number at the interface:
            ELSE                              ! (1) use value from the
              c2d=Cu_adv(i,k+1)               !     grid box upstream in
              dh=z_w(i,j,k+1)-z_w(i,j,k)      !     vertical direction;
            END IF
!
!!          c2d=0.5_r8*(Cu_adv(i,k  )+                                  &
!!   &                  Cu_adv(i,k+1))        ! (2) average the two; or
!!          dh=0.5_r8*(z_w(i,j,k+1)-                                    &
!!   &                 z_w(i,j,k-1))
!!
!!          c2d=MAX(Cu_adv(i,k),                                        &
!!   &              Cu_adv(i,k+1))            ! (3) pick its maximum
!!          dh=MIN(z_w(i,j,k+1)-z_w(i,j,k),                             &
!!                 z_w(i,j,k)-z_w(i,j,k-1))
!!
            cw_max=amax*dh-c2d*Cu_adv(i,0)    ! compare the vertical
            IF (cw_max.ge.0.0_r8) THEN        ! displacement to dz*amax.
              cw_max2=cw_max*cw_max           ! Partition W into Wi and
              cw_min=cw_max*cmnx_ratio        ! We.
              cw=ABS(Wi(i,j,k))*Cu_adv(i,0)
              IF (cw.le.cw_min) THEN
                cff=cw_max2
              ELSE IF (cw.le.cutoff*cw_max) THEN
                cff=cw_max2+r4cmx*(cw-cw_min)**2
              ELSE
                cff=cw_max*cw
              END IF
!
              W(i,j,k)=cw_max2*Wi(i,j,k)/cff
              Wi(i,j,k)=Wi(i,j,k)-W(i,j,k)
            ELSE                              ! All the displacement is
              W(i,j,k)=0.0_r8                 ! greater than amax*dz, so
            END IF                            ! keep it all into Wi.
# endif
          END DO
        END DO
        DO i=Istr,Iend
          W(i,j,N(ng))=0.0_r8
        END DO
      END DO
!
!  Set lateral boundary conditions.
!
      CALL bc_w3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  W)
# ifdef OMEGA_IMPLICIT
      CALL bc_w3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  Wi)
# endif
# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    W)
#  ifdef OMEGA_IMPLICIT
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Wi)
#  endif
# endif
!
      RETURN
      END SUBROUTINE omega_tile
!
!***********************************************************************
      SUBROUTINE scale_omega (ng, tile, LBi, UBi, LBj, UBj, LBk, UBk,   &
     &                        pm, pn, W, Wscl)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_3d_mod, ONLY : exchange_w3d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: W(LBi:,LBj:,LBk:)
      real(r8), intent(out) :: Wscl(LBi:,LBj:,LBk:)
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: W(LBi:UBi,LBj:UBj,LBk:UBk)
      real(r8), intent(out) :: Wscl(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Scale omega vertical velocity to m/s.
!-----------------------------------------------------------------------
!
      DO k=LBk,UBk
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Wscl(i,j,k)=W(i,j,k)*pm(i,j)*pn(i,j)
          END DO
        END DO
      END DO
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          Wscl)
      END IF
!
      RETURN
      END SUBROUTINE scale_omega
#endif
      END MODULE omega_mod
