#include "cppdefs.h"
      MODULE tl_omega_mod
#if defined TANGENT && defined SOLVE3D
!
!svn $Id: tl_omega.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes S-coordinate vertical velocity (m^3/s),       !
!                                                                      !
!                  W=[Hz/(m*n)]*omega,                                 !
!                                                                      !
!  diagnostically at horizontal RHO-points and vertical W-points.      !
!                                                                      !
!  BASIC STATE variables needed: W, z_w.                               !
!                                                                      !
!  NOTE: We need to recompute  basic state W in this routine since     !
!  ----  intermediate values of W are needed by the tangent linear     !
!        and adjoint routines.                                         !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: tl_omega
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_omega (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
# if defined SEDIMENT && defined SED_MORPH
      USE mod_sedbed
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
      CALL wclock_on (ng, iTLM, 13)
# endif
      CALL tl_omega_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                    nstp(ng), nnew(ng),                           &
     &                    GRID(ng) % omn,                               &
     &                    SEDBED(ng) % bed_thick,                       &
     &                    SEDBED(ng) % tl_bed_thick,                    &
# endif
     &                    GRID(ng) % Huon,                              &
     &                    GRID(ng) % Hvom,                              &
     &                    GRID(ng) % z_w,                               &
     &                    GRID(ng) % tl_Huon,                           &
     &                    GRID(ng) % tl_Hvom,                           &
     &                    GRID(ng) % tl_z_w,                            &
     &                    OCEAN(ng) % W,                                &
     &                    OCEAN(ng) % tl_W)
# ifdef PROFILE
      CALL wclock_off (ng, iTLM, 13)
# endif

      RETURN
      END SUBROUTINE tl_omega
!
!***********************************************************************
      SUBROUTINE tl_omega_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                          nstp, nnew,                             &
     &                          omn,
     &                          bed_thick, tl_bed_thick,                &
# endif
     &                          Huon, Hvom, z_w,                        &
     &                          tl_Huon, tl_Hvom, tl_z_w,               &
     &                          W, tl_W)
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
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      integer, intent(in) :: nstp, nnew
# endif
!
# ifdef ASSUMED_SHAPE
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in):: bed_thick(LBi:,LBj:,:)
      real(r8), intent(in):: tl_bed_thick(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: tl_Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_z_w(LBi:,LBj:,0:)

      real(r8), intent(out) :: W(LBi:,LBj:,0:)
      real(r8), intent(out) :: tl_W(LBi:,LBj:,0:)

# else

#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
      real(r8), intent(in):: bed_thick(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in):: tl_bed_thick(LBi:UBi,LBj:UBj,2)
#  endif
      real(r8), intent(in) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: tl_Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_Hvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_z_w(LBi:UBi,LBj:UBj,0:N(ng))

      real(r8), intent(out) :: W(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(out) :: tl_W(LBi:UBi,LBj:UBj,0:N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, ii, is, j, jj, k
      real(r8) :: cff, tl_cff
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8) :: cff1
# endif
      real(r8), dimension(IminS:ImaxS) :: wrk
      real(r8), dimension(IminS:ImaxS) :: tl_wrk

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
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      cff1=1.0_r8/dt(ng)
# endif
      DO j=Jstr,Jend
        DO i=Istr,Iend
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
          W(i,j,0)=-cff1*(bed_thick(i,j,nstp)-                          &
     &                    bed_thick(i,j,nnew))*omn(i,j)
          tl_W(i,j,0)=-cff1*(tl_bed_thick(i,j,nstp)-                    &
     &                       tl_bed_thick(i,j,nnew))*omn(i,j)
# else
          W(i,j,0)=0.0_r8
          tl_W(i,j,0)=0.0_r8
# endif
        END DO
!
!  Code added to clear tl_W to be consistent with adjoint.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
           tl_W(i,j,k)=0.0_r8
          END DO
        END DO
        DO k=1,N(ng)
          DO i=Istr,Iend
            W(i,j,k)=W(i,j,k-1)-                                        &
     &               (Huon(i+1,j,k)-Huon(i,j,k)+                        &
     &                Hvom(i,j+1,k)-Hvom(i,j,k))
            tl_W(i,j,k)=tl_W(i,j,k-1)-                                  &
     &                  (tl_Huon(i+1,j,k)-tl_Huon(i,j,k)+               &
     &                   tl_Hvom(i,j+1,k)-tl_Hvom(i,j,k))
          END DO
        END DO
!
!  Apply mass point sources (volume vertical influx), if any.
!
        IF (LwSrc(ng)) THEN
          DO is=1,Nsrc(ng)
            ii=SOURCES(ng)%Isrc(is)
            jj=SOURCES(ng)%Jsrc(is)
            IF (((IstrR.le.ii).and.(ii.le.IendR)).and.                  &
     &          ((JstrR.le.jj).and.(jj.le.JendR)).and.                  &
     &          (j.eq.jj)) THEN
              DO k=1,N(ng)
                W(ii,jj,k)=W(ii,jj,k)+SOURCES(ng)%Qsrc(is,k)
!!              tl_W(ii,jj,k)=tl_W(ii,jj,k)+0.0_r8
              END DO
            END IF
          END DO
        END IF
!
        DO i=Istr,Iend
          cff=1.0_r8/(z_w(i,j,N(ng))-z_w(i,j,0))
          tl_cff=-cff*cff*(tl_z_w(i,j,N(ng))-tl_z_w(i,j,0))
          wrk(i)=cff*W(i,j,N(ng))
          tl_wrk(i)=tl_cff*W(i,j,N(ng))+cff*tl_W(i,j,N(ng))
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
     &               wrk(i)*(z_w(i,j,k)-z_w(i,j,0))
            tl_W(i,j,k)=tl_W(i,j,k)-                                    &
     &                  tl_wrk(i)*(z_w(i,j,k)-z_w(i,j,0))-              &
     &                  wrk(i)*(tl_z_w(i,j,k)-tl_z_w(i,j,0))
          END DO
        END DO
        DO i=Istr,Iend
          W(i,j,N(ng))=0.0_r8
          tl_W(i,j,N(ng))=0.0_r8
        END DO
      END DO
!
!  Set lateral boundary conditions.
!
      CALL bc_w3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  W)
      CALL bc_w3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  tl_W)

# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    W)
      CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_W)
# endif

      RETURN
      END SUBROUTINE tl_omega_tile
#endif
      END MODULE tl_omega_mod
