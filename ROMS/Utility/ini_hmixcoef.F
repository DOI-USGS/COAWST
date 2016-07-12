#include "cppdefs.h"
      MODULE ini_hmixcoef_mod
!
! svn $Id: ini_hmixcoef.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine initializes horizontal mixing coefficients arrays      !
!  according to the model flag.                                        !
!                                                                      !
!  WARNING:   All biharmonic coefficients are assumed to have the      !
!             square root taken and have  m^2 s^-1/2 units.  This      !
!             will allow multiplying the  biharmonic  coefficient      !
!             to harmonic operator.                                    !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ini_hmixcoef

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ini_hmixcoef (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
#ifdef SOLVE3D
      real(r8) :: diffusion2(MT), diffusion4(MT)
#endif
      real(r8) :: viscosity2, viscosity4
!
#include "tile.h"

      CALL ini_hmixcoef_tile (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        GRID(ng) % grdscl,                        &
#ifdef SOLVE3D
# ifdef TS_DIF2
     &                        MIXING(ng) % diff2,                       &
# endif
# ifdef TS_DIF4
     &                        MIXING(ng) % diff4,                       &
# endif
#endif
#ifdef UV_VIS2
     &                        MIXING(ng) % visc2_p,                     &
     &                        MIXING(ng) % visc2_r,                     &
#endif
#ifdef UV_VIS4
     &                        MIXING(ng) % visc4_p,                     &
     &                        MIXING(ng) % visc4_r,                     &
#endif
#ifdef SOLVE3D
     &                        diffusion2, diffusion4,                   &
#endif
     &                        viscosity2, viscosity4)

      RETURN
      END SUBROUTINE ini_hmixcoef
!
!***********************************************************************
      SUBROUTINE ini_hmixcoef_tile (ng, tile, model,                    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              grdscl,                             &
#ifdef SOLVE3D
# ifdef TS_DIF2
     &                              diff2,                              &
# endif
# ifdef TS_DIF4
     &                              diff4,                              &
# endif
#endif
#ifdef UV_VIS2
     &                              visc2_p,                            &
     &                              visc2_r,                            &
#endif
#ifdef UV_VIS4
     &                              visc4_p,                            &
     &                              visc4_r,                            &
#endif
#ifdef SOLVE3D
     &                              diffusion2, diffusion4,             &
#endif
     &                              viscosity2, viscosity4)
!***********************************************************************
!
      USE mod_param
      USE mod_mixing
      USE mod_scalars
!
      USE exchange_2d_mod
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#ifdef SOLVE3D
      real(r8), intent(out) :: diffusion2(MT), diffusion4(MT)
#endif
      real(r8), intent(out) :: viscosity2, viscosity4
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: grdscl(LBi:,LBj:)
# ifdef SOLVE3D
#  ifdef TS_DIF2
      real(r8), intent(inout) :: diff2(LBi:,LBj:,:)
#  endif
#  ifdef TS_DIF4
      real(r8), intent(inout) :: diff4(LBi:,LBj:,:)
#  endif
# endif
# if defined UV_VIS2
      real(r8), intent(inout) :: visc2_p(LBi:,LBj:)
      real(r8), intent(inout) :: visc2_r(LBi:,LBj:)
# endif
# ifdef UV_VIS4
      real(r8), intent(inout) :: visc4_p(LBi:,LBj:)
      real(r8), intent(inout) :: visc4_r(LBi:,LBj:)
# endif
#else
      real(r8), intent(in) :: grdscl(LBi:UBi,LBj:UBj)
# ifdef SOLVE3D
#  ifdef TS_DIF2
      real(r8), intent(inout) :: diff2(LBi:UBi,LBj:UBj,NT(ng))
#  endif
#  ifdef TS_DIF4
      real(r8), intent(inout) :: diff4(LBi:UBi,LBj:UBj,NT(ng))
#  endif
# endif
# if defined UV_VIS2
      real(r8), intent(inout) :: visc2_p(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: visc2_r(LBi:UBi,LBj:UBj)
# endif
# ifdef UV_VIS4
      real(r8), intent(inout) :: visc4_p(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: visc4_r(LBi:UBi,LBj:UBj)
# endif
#endif
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j
#ifdef SOLVE3D
      integer :: itrc
#endif
      real(r8) :: cff

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set horizontal, constant, mixing coefficient according to model flag.
!-----------------------------------------------------------------------
!
      IF (model.eq.iNLM) THEN
        viscosity2=nl_visc2(ng)
        viscosity4=nl_visc4(ng)
#ifdef SOLVE3D
        DO itrc=1,NT(ng)
          diffusion2(itrc)=nl_tnu2(itrc,ng)
          diffusion4(itrc)=nl_tnu4(itrc,ng)
        END DO
#endif
#if defined TANGENT || defined TL_IOMS
      ELSE IF ((model.eq.iTLM).or.(model.eq.iRPM)) THEN
        viscosity2=tl_visc2(ng)
        viscosity4=tl_visc4(ng)
# ifdef SOLVE3D
        DO itrc=1,NT(ng)
          diffusion2(itrc)=tl_tnu2(itrc,ng)
          diffusion4(itrc)=tl_tnu4(itrc,ng)
        END DO
# endif
#endif
#ifdef ADJOINT
      ELSE IF (model.eq.iADM) THEN
        viscosity2=ad_visc2(ng)
        viscosity4=ad_visc4(ng)
# ifdef SOLVE3D
        DO itrc=1,NT(ng)
          diffusion2(itrc)=ad_tnu2(itrc,ng)
          diffusion4(itrc)=ad_tnu4(itrc,ng)
        END DO
# endif
#endif
      END IF
!
!  Update generic values.
!
      IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
        visc2(ng)=viscosity2
        visc4(ng)=viscosity4
#ifdef SOLVE3D
        DO itrc=1,NT(ng)
          tnu2(itrc,ng)=diffusion2(itrc)
          tnu4(itrc,ng)=diffusion4(itrc)
        END DO
#endif
      END IF

#if defined UV_VIS2 || defined UV_VIS4  || \
  ((defined TS_DIF2 || defined TS_DIF4) && defined SOLVE3D)
!
!-----------------------------------------------------------------------
!  Initialize horizontal mixing arrays to constant mixing coefficient.
!-----------------------------------------------------------------------
!
# ifdef _OPENMP
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        Jmax=BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=Jend
      END IF
# else
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
# endif
!
# if defined UV_VIS2
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          visc2_p(i,j)=viscosity2
          visc2_r(i,j)=viscosity2
        END DO
      END DO
# endif
# ifdef UV_VIS4
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          visc4_p(i,j)=viscosity4
          visc4_r(i,j)=viscosity4
        END DO
      END DO
# endif
# ifdef SOLVE3D
#  ifdef TS_DIF2
      DO itrc=1,NT(ng)
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            diff2(i,j,itrc)=diffusion2(itrc)
          END DO
        END DO
      END DO
#  endif
#  ifdef TS_DIF4
      DO itrc=1,NT(ng)
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            diff4(i,j,itrc)=diffusion4(itrc)
          END DO
        END DO
      END DO
#  endif
# endif
#endif

#ifdef VISC_GRID
!
!-----------------------------------------------------------------------
!  Scale horizontal viscosity according to the grid size.  The values
!  used during initialization above are over-written.
!-----------------------------------------------------------------------
!
# ifdef UV_VIS2
      cff=viscosity2/grdmax(ng)
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          visc2_r(i,j)=cff*grdscl(i,j)
        END DO
      END DO
      cff=0.25_r8*cff
      DO j=JstrP,JendT
        DO i=IstrP,IendT
          visc2_p(i,j)=cff*(grdscl(i-1,j-1)+grdscl(i,j-1)+              &
     &                      grdscl(i-1,j  )+grdscl(i,j  ))
        END DO
      END DO
# endif
# ifdef UV_VIS4
      cff=viscosity4/(grdmax(ng)**3)
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          visc4_r(i,j)=cff*grdscl(i,j)**3
        END DO
      END DO
      cff=0.25_r8*cff
      DO j=JstrP,JendT
        DO i=IstrP,IendT
          visc4_p(i,j)=cff*(grdscl(i,j  )**3+grdscl(i-1,j  )**3+        &
     &                      grdscl(i,j-1)**3+grdscl(i-1,j-1)**3)
        END DO
      END DO
# endif
#endif

#if defined DIFF_GRID && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Scale horizontal diffusion according to the grid size.
!-----------------------------------------------------------------------
!
# ifdef TS_DIF2
      DO itrc=1,NT(ng)
        cff=diffusion2(itrc)/grdmax(ng)
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            diff2(i,j,itrc)=cff*grdscl(i,j)
          END DO
        END DO
      END DO
# endif
# ifdef TS_DIF4
      DO itrc=1,NT(ng)
        cff=diffusion4(itrc)/(grdmax(ng)**3)
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            diff4(i,j,itrc)=cff*grdscl(i,j)**3
          END DO
        END DO
      END DO
# endif
#endif

#if !defined ANA_SPONGE && \
    ( defined UV_VIS2   || defined UV_VIS4  || \
     ((defined TS_DIF2  || defined TS_DIF4) && defined SOLVE3D))
!
!-----------------------------------------------------------------------
!  Increase horizontal mixing coefficients in the sponge areas using
!  the nondimentional factors read from application Grid NetCDF file.
!-----------------------------------------------------------------------
!
      IF (LuvSponge(ng)) THEN
# ifdef UV_VIS2
        DO i=IstrT,IendT
          DO j=JstrT,JendT
            visc2_r(i,j)=ABS(MIXING(ng)%visc_factor(i,j))*              &
     &                   visc2_r(i,j)
          END DO
        END DO
        DO i=IstrP,IendT
          DO j=JstrP,JendT
            visc2_p(i,j)=0.25_r8*                                       &
     &                   ABS(MIXING(ng)%visc_factor(i-1,j-1)+           &
     &                       MIXING(ng)%visc_factor(i  ,j-1)+           &
     &                       MIXING(ng)%visc_factor(i-1,j  )+           &
     &                       MIXING(ng)%visc_factor(i  ,j  ))*          &
     &                   visc2_p(i,j)
          END DO
        END DO
# endif
# ifdef UV_VIS4
        DO i=IstrT,IendT
          DO j=JstrT,JendT
            visc4_r(i,j)=ABS(MIXING(ng)%visc_factor(i,j))*              &
     &                   visc4_r(i,j)
          END DO
        END DO
        DO i=IstrP,IendT
          DO j=JstrP,JendT
            visc4_p(i,j)=0.25_r8*                                       &
     &                   ABS(MIXING(ng)%visc_factor(i-1,j-1)+           &
     &                       MIXING(ng)%visc_factor(i  ,j-1)+           &
     &                       MIXING(ng)%visc_factor(i-1,j  )+           &
     &                       MIXING(ng)%visc_factor(i  ,j  ))*          &
     &                   visc4_p(i,j)
          END DO
        END DO
# endif
      END IF

# ifdef SOLVE3D
#  ifdef TS_DIF2
      DO itrc=1,NT(ng)
        IF (LtracerSponge(itrc,ng)) THEN
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              diff2(i,j,itrc)=ABS(MIXING(ng)%diff_factor(i,j))*         &
     &                        diff2(i,j,itrc)
            END DO
          END DO
        END IF
      END DO
#  endif
#  ifdef TS_DIF4
      DO itrc=1,NT(ng)
        IF (LtracerSponge(itrc,ng)) THEN
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              diff4(i,j,itrc)=ABS(MIXING(ng)%diff_factor(i,j))*         &
     &                        diff4(i,j,itrc)
            END DO
          END DO
        END IF
      END DO
#  endif
# endif
#endif
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
#ifdef UV_VIS2
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          visc2_r)
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          visc2_p)
      END IF
#endif
#ifdef UV_VIS4
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          visc4_r)
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          visc4_p)
      END IF
#endif

#ifdef SOLVE3D
# ifdef TS_DIF2
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        DO itrc=1,NT(ng)
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            diff2(:,:,itrc))
        END DO
      END IF
# endif
# ifdef TS_DIF4
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        DO itrc=1,NT(ng)
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            diff4(:,:,itrc))
        END DO
      END IF
# endif
#endif

#ifdef DISTRIBUTE
# ifdef UV_VIS2
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    visc2_r, visc2_p)
# endif
# ifdef UV_VIS4
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    visc4_r, visc4_p)
# endif
# ifdef SOLVE3D
#  ifdef TS_DIF2
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, NT(ng),                &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    diff2)
#  endif
#  ifdef TS_DIF4
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, NT(ng),                &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    diff4)
#  endif
# endif
#endif

      RETURN
      END SUBROUTINE ini_hmixcoef_tile

      END MODULE ini_hmixcoef_mod
