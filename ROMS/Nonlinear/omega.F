#include "cppdefs.h"
      MODULE omega_mod
#ifdef SOLVE3D
# if defined COMPOSED_GRID || defined REFINED_GRID
#  define I_RANGE IstrTU+1,MIN(IendT, Lm(ng)+2)
#  define J_RANGE JstrTV+1,MIN(JendT, Mm(ng)+2)
# else
#  define I_RANGE Istr,Iend
#  define J_RANGE Jstr,Jend
# endif
!
!svn $Id: omega.F 732 2008-09-07 01:55:51Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine computes S-coordinate vertical velocity (m^3/s),       !
!                                                                      !
!                  W=[Hz/(m*n)]*omega,                                 !
!                                                                      !
!  diagnostically at horizontal RHO-points and vertical W-points.      !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: omega, scale_omega

      CONTAINS
!
!***********************************************************************
      SUBROUTINE omega (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
# ifdef Q_PSOURCE
      USE mod_sources
# endif
# if defined Q_PSOURCE || (defined SEDIMENT && defined SED_MORPH)
      USE mod_stepping
# endif
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
      CALL wclock_on (ng, iNLM, 13)
# endif
      CALL omega_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
# if defined SEDIMENT && defined SED_MORPH
     &                 nstp(ng), nnew(ng),                              &
# endif
# ifdef Q_PSOURCE
     &                 Nsrc(ng),                                        &
     &                 SOURCES(ng) % Isrc,                              &
     &                 SOURCES(ng) % Jsrc,                              &
     &                 SOURCES(ng) % Dsrc,                              &
     &                 SOURCES(ng) % Qsrc,                              &
# endif
# if defined SEDIMENT && defined SED_MORPH
     &                 GRID(ng) % omn,                                  &
     &                 GRID(ng) % bed_thick,                            &
# endif
     &                 GRID(ng) % Huon,                                 &
     &                 GRID(ng) % Hvom,                                 &
     &                 GRID(ng) % z_w,                                  &
     &                 OCEAN(ng) % W)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 13)
# endif
      RETURN
      END SUBROUTINE omega
!
!***********************************************************************
      SUBROUTINE omega_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
# if defined SEDIMENT && defined SED_MORPH
     &                       nstp, nnew,                                &
# endif
# ifdef Q_PSOURCE
     &                       Nsrc, Isrc, Jsrc, Dsrc, Qsrc,              &
# endif
# if defined SEDIMENT && defined SED_MORPH
     &                       omn, bed_thick,                            &
# endif
     &                       Huon, Hvom,                                &
     &                       z_w,                                       &
     &                       W)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_3d_mod, ONLY : bc_w3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
# if defined SEDIMENT && defined SED_MORPH
      integer, intent(in) :: nstp, nnew
# endif
# ifdef Q_PSOURCE
      integer, intent(in) :: Nsrc
# endif
!
# ifdef ASSUMED_SHAPE
#  ifdef Q_PSOURCE
      integer, intent(in) :: Isrc(:)
      integer, intent(in) :: Jsrc(:)
      real(r8), intent(in) :: Dsrc(:)
      real(r8), intent(in) :: Qsrc(:,:)
#  endif
#  if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in):: bed_thick(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(out) :: W(LBi:,LBj:,0:)
# else
#  ifdef Q_PSOURCE
      integer, intent(in) :: Isrc(Nsrc)
      integer, intent(in) :: Jsrc(Nsrc)
      real(r8), intent(in) :: Dsrc(Nsrc)
      real(r8), intent(in) :: Qsrc(Nsrc,N(ng))
#  endif
#  if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
      real(r8), intent(in):: bed_thick(LBi:UBi,LBj:UBj,2)
#  endif
      real(r8), intent(in) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(out) :: W(LBi:UBi,LBj:UBj,0:N(ng))
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
# ifdef Q_PSOURCE
      integer :: ii, jj, is
# endif
# if defined SEDIMENT && defined SED_MORPH
      real(r8) :: cff1
# endif
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: wrk

# include "set_bounds.h"
!
!------------------------------------------------------------------------
!  Vertically integrate horizontal mass flux divergence.
!------------------------------------------------------------------------
!
!  Starting with zero vertical velocity at the bottom, integrate
!  from the bottom (k=0) to the free-surface (k=N).  The w(:,:,N(ng))
!  contains the vertical velocity at the free-surface, d(zeta)/d(t).
!  Notice that barotropic mass flux divergence is not used directly.
!
# if defined SEDIMENT && defined SED_MORPH
      cff1=1.0_r8/dt(ng)
# endif
      DO j=J_RANGE
        DO i=I_RANGE
# if defined SEDIMENT && defined SED_MORPH
          W(i,j,0)=-cff1*(bed_thick(i,j,nstp)-                          &
     &                    bed_thick(i,j,nnew))*omn(i,j)
# else
          W(i,j,0)=0.0_r8
# endif
        END DO
        DO k=1,N(ng)
          DO i=I_RANGE
            W(i,j,k)=W(i,j,k-1)-                                        &
     &               (Huon(i+1,j,k)-Huon(i,j,k)+                        &
     &                Hvom(i,j+1,k)-Hvom(i,j,k))
          END DO
        END DO
# ifdef Q_PSOURCE
!
!  Apply mass point sources - Volume influx.
!
        DO is=1,Nsrc
          ii=Isrc(is)
          jj=Jsrc(is)
          IF (((IstrR.le.ii).and.(ii.le.IendR)).and.                    &
     &        ((JstrR.le.jj).and.(jj.le.JendR)).and.                    &
     &        (j.eq.jj)) THEN
            DO k=1,N(ng)
              W(ii,jj,k)=W(ii,jj,k)+Qsrc(is,k)
            END DO
          END IF
        END DO
# endif
!
        DO i=I_RANGE
          wrk(i)=W(i,j,N(ng))/(z_w(i,j,N(ng))-z_w(i,j,0))
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
          DO i=I_RANGE
            W(i,j,k)=W(i,j,k)-                                          &
     &               wrk(i)*(z_w(i,j,k)-z_w(i,j,0))
          END DO
        END DO
        DO i=I_RANGE
          W(i,j,N(ng))=0.0_r8
        END DO
      END DO
!
!  Set lateral boundary conditions.
!
      CALL bc_w3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  W)
# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    W)
# endif
# undef I_RANGE
# undef J_RANGE

      RETURN
      END SUBROUTINE omega_tile
!
!***********************************************************************
      SUBROUTINE scale_omega (ng, tile, LBi, UBi, LBj, UBj, LBk, UBk,   &
     &                        pm, pn, W, Wscl)
!***********************************************************************
!
      USE mod_param
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
!      integer :: IstrR, IendR, JstrR, JendR
!      integer :: IstrT, IendT, JstrT, JendT
!      integer :: Istr, Iend, Jstr, Jend, IstrU, JstrV
      integer :: i, j, k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Scale omega vertical velocity to m/s.
!-----------------------------------------------------------------------
!
      DO k=LBk,UBk
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Wscl(i,j,k)=W(i,j,k)*pm(i,j)*pn(i,j)
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE scale_omega
#endif
      END MODULE omega_mod
