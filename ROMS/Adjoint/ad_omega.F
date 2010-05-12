#include "cppdefs.h"
      MODULE ad_omega_mod
#if defined ADJOINT && defined SOLVE3D
!
!svn $Id: ad_omega.F 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group       Andrew M. Moore   !
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
!  BASIC STATE variables needed: Huon, Hvom, z_w.                      !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ad_omega

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_omega (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
# if defined SEDIMENT && defined SED_MORPH
      USE mod_sedbed
# endif
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
      CALL wclock_on (ng, iADM, 13)
# endif
      CALL ad_omega_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                    nstp(ng), nnew(ng),                           &
# endif
# ifdef Q_PSOURCE
     &                    Msrc(ng), Nsrc(ng),                           &
     &                    SOURCES(ng) % Isrc,                           &
     &                    SOURCES(ng) % Jsrc,                           &
     &                    SOURCES(ng) % Dsrc,                           &
     &                    SOURCES(ng) % Qsrc,                           &
# endif
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                    GRID(ng) % omn,                               &
     &                    SEDBED(ng) % bed_thick,                       &
     &                    SEDBED(ng) % ad_bed_thick,                    &
# endif
     &                    GRID(ng) % Huon,                              &
     &                    GRID(ng) % Hvom,                              &
     &                    GRID(ng) % z_w,                               &
     &                    GRID(ng) % ad_Huon,                           &
     &                    GRID(ng) % ad_Hvom,                           &
     &                    GRID(ng) % ad_z_w,                            &
     &                    OCEAN(ng) % W,                                &
     &                    OCEAN(ng) % ad_W)
# ifdef PROFILE
      CALL wclock_off (ng, iADM, 13)
# endif
      RETURN
      END SUBROUTINE ad_omega
!
!***********************************************************************
      SUBROUTINE ad_omega_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                          nstp, nnew,                             &
# endif
# ifdef Q_PSOURCE
     &                          Msrc, Nsrc,                             &
     &                          Isrc, Jsrc, Dsrc, Qsrc,                 &
# endif
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                          omn,
     &                          bed_thick, ad_bed_thick,                &
# endif
     &                          Huon, Hvom, z_w,                        &
     &                          ad_Huon, ad_Hvom, ad_z_w,               &
     &                          W, ad_W)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE ad_bc_3d_mod, ONLY : ad_bc_w3d_tile
      USE bc_3d_mod, ONLY : bc_w3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange3d
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
# ifdef Q_PSOURCE
      integer, intent(in) :: Msrc, Nsrc
# endif
!
# ifdef ASSUMED_SHAPE
#  ifdef Q_PSOURCE
      integer, intent(in) :: Isrc(:)
      integer, intent(in) :: Jsrc(:)

      real(r8), intent(in) :: Dsrc(:)
      real(r8), intent(in) :: Qsrc(:,:)
#  endif
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in):: bed_thick(LBi:,LBj:,:)
      real(r8), intent(inout):: ad_bed_thick(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: ad_Huon(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_Hvom(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_w(LBi:,LBj:,0:)
      real(r8), intent(inout) :: ad_W(LBi:,LBj:,0:)

      real(r8), intent(out) :: W(LBi:,LBj:,0:)

# else

#  ifdef Q_PSOURCE
      integer, intent(in) :: Isrc(Msrc)
      integer, intent(in) :: Jsrc(Msrc)

      real(r8), intent(in) :: Dsrc(Msrc)
      real(r8), intent(in) :: Qsrc(Msrc,N(ng))
#  endif
      real(r8), intent(in) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
      real(r8), intent(in):: bed_thick(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout):: ad_bed_thick(LBi:UBi,LBj:UBj,2)
#  endif
      real(r8), intent(inout) :: ad_Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_Hvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: ad_W(LBi:UBi,LBj:UBj,0:N(ng))

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
      real(r8) :: cff
      real(r8) :: ad_cff, adfac
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8) :: cff1
# endif
      real(r8), dimension(IminS:ImaxS) :: wrk
      real(r8), dimension(IminS:ImaxS) :: ad_wrk

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_cff=0.0_r8
      DO i=IminS,ImaxS
        ad_wrk(i)=0.0_r8
      END DO
!
!-----------------------------------------------------------------------
!  Vertically integrage horizontal mass flux divergence.
!-----------------------------------------------------------------------
!
!  Set lateral boundary conditions.
!
# ifdef DISTRIBUTE
!>    CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    tl_W)
!>
      CALL ad_mp_exchange3d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj, 0, N(ng),              &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_W)
# endif
!>    CALL bc_w3d_tile (ng, tile,                                       &
!>   &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
!>   &                  tl_W)
!>
      CALL ad_bc_w3d_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj, 0, N(ng),                &
     &                     ad_W)
!
!  In order to insure zero vertical velocity at the free-surface,
!  subtract the vertical velocities of the moving S-coordinates
!  isosurfaces. These isosurfaces are proportional to d(zeta)/d(t).
!  The proportionaly coefficients are a linear function of the
!  S-coordinate with zero value at the bottom (k=0) and unity at
!  the free-surface (k=N).
!
!  Notice that here we need to recompute the intermediate value
!  of W which is needed for wrk.
!
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      cff1=1.0_r8/dt(ng)
# endif
      DO j=Jstr,Jend
        DO i=Istr,Iend
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
          W(i,j,0)=-cff1*(bed_thick(i,j,nstp)-                          &
     &                    bed_thick(i,j,nnew))*omn(i,j)
# else
          W(i,j,0)=0.0_r8
# endif

        END DO
        DO k=1,N(ng)
          DO i=Istr,Iend
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
        DO i=Istr,Iend
          wrk(i)=W(i,j,N(ng))/(z_w(i,j,N(ng))-z_w(i,j,0))
        END DO
!
!  Starting with zero vertical velocity at the bottom, integrate
!  from the bottom (k=0) to the free-surface (k=N).  The w(:,:,N)
!  contains the vertical velocity at the free-surface, d(zeta)/d(t).
!  Notice that barotropic mass flux divergence is not used directly.
!
        DO i=Istr,Iend
!>        tl_W(i,j,N(ng))=0.0_r8
!>
          ad_W(i,j,N(ng))=0.0_r8
        END DO
        DO k=N(ng)-1,1,-1
          DO i=Istr,Iend
!>          tl_W(i,j,k)=tl_W(i,j,k)-                                    &
!>   &                  tl_wrk(i)*(z_w(i,j,k)-z_w(i,j,0))-              &
!>   &                  wrk(i)*(tl_z_w(i,j,k)-tl_z_w(i,j,0))
!>
            adfac=wrk(i)*ad_W(i,j,k)
            ad_wrk(i)=ad_wrk(i)-                                        &
     &                ad_W(i,j,k)*(z_w(i,j,k)-z_w(i,j,0))
            ad_z_w(i,j,0)=ad_z_w(i,j,0)+adfac
            ad_z_w(i,j,k)=ad_z_w(i,j,k)-adfac
          END DO
        END DO
        DO i=Istr,Iend
          cff=1.0_r8/(z_w(i,j,N(ng))-z_w(i,j,0))
!>        tl_wrk(i)=tl_cff*W(i,j,N(ng))+cff*tl_W(i,j,N(ng))
!>
          ad_W(i,j,N(ng))=ad_W(i,j,N(ng))+cff*ad_wrk(i)
          ad_cff=ad_cff+W(i,j,N(ng))*ad_wrk(i)
          ad_wrk(i)=0.0_r8
!>        tl_cff=-cff*cff*(tl_z_w(i,j,N(ng))-tl_z_w(i,j,0))
!>
          adfac=-cff*cff*ad_cff
          ad_z_w(i,j,0    )=ad_z_w(i,j,0    )-adfac
          ad_z_w(i,j,N(ng))=ad_z_w(i,j,N(ng))+adfac
          ad_cff=0.0_r8
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
!>            tl_W(ii,jj,k)=tl_W(ii,jj,k)+0.0_r8
!>
!!            ad_W(ii,jj,k)=ad_W(ii,jj,k)+0.0_r8
            END DO
          END IF
        END DO
# endif
!>      DO k=1,N(ng)
!>
        DO k=N(ng),1,-1
          DO i=Istr,Iend
!>          tl_W(i,j,k)=tl_W(i,j,k-1)-                                  &
!>   &                  (tl_Huon(i+1,j,k)-tl_Huon(i,j,k)+               &
!>   &                   tl_Hvom(i,j+1,k)-tl_Hvom(i,j,k))
!>
            ad_W(i,j,k-1)=ad_W(i,j,k-1)+ad_W(i,j,k)
            ad_Huon(i  ,j,k)=ad_Huon(i  ,j,k)+ad_W(i,j,k)
            ad_Huon(i+1,j,k)=ad_Huon(i+1,j,k)-ad_W(i,j,k)
            ad_Hvom(i,j  ,k)=ad_Hvom(i,j  ,k)+ad_W(i,j,k)
            ad_Hvom(i,j+1,k)=ad_Hvom(i,j+1,k)-ad_W(i,j,k)
          END DO
        END DO
!
!  Clear ad_W here since it is a diagnostic variable.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
!>          tl_W(i,j,k)=0.0_r8
!>
            ad_W(i,j,k)=0.0_r8
          END DO
        END DO
        DO i=Istr,Iend
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
!>        tl_W(i,j,0)=-cff1*(tl_bed_thick(i,j,nstp)-                    &
!>   &                       tl_bed_thick(i,j,nnew))*omn(i,j)
!>
          adfac=cff1*omn(i,j)*ad_W(i,j,0)
          ad_bed_thick(i,j,nnew)=ad_bed_thick(i,j,nnew)+adfac
          ad_bed_thick(i,j,nstp)=ad_bed_thick(i,j,nstp)+adfac
          ad_W(i,j,0)=0.0_r8
# else
!>        tl_W(i,j,0)=0.0_r8
!>
          ad_W(i,j,0)=0.0_r8
# endif
        END DO
!
!  Complete the computation of BASIC STATE W here so that it is correct
!  for the remainder of the code.
!
        DO k=N(ng)-1,1,-1
          DO i=Istr,Iend
            W(i,j,k)=W(i,j,k)-wrk(i)*(z_w(i,j,k)-z_w(i,j,0))
          END DO
        END DO
        DO i=Istr,Iend
          W(i,j,N(ng))=0.0_r8
        END DO
      END DO
!
!  Set lateral boundary conditions for basic state.
!
      CALL bc_w3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  W)
# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iADM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    W)
# endif

      RETURN
      END SUBROUTINE ad_omega_tile
#endif
      END MODULE ad_omega_mod
