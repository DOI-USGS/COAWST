#include "cppdefs.h"
      MODULE ini_hmixcoef_mod
!
! svn $Id: ini_hmixcoef.F 442 2010-01-27 21:25:57Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
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
#ifdef SOLVE3D
# ifdef TS_DIF2
     &                        MIXING(ng) % diff2,                       &
# endif
# ifdef TS_DIF4
     &                        MIXING(ng) % diff4,                       &
# endif
#endif
#if defined UV_VIS2 || !defined SOLVE3D
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
#ifdef SOLVE3D
# ifdef TS_DIF2
     &                              diff2,                              &
# endif
# ifdef TS_DIF4
     &                              diff4,                              &
# endif
#endif
#if defined UV_VIS2 || !defined SOLVE3D
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
      USE mod_scalars
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
# ifdef SOLVE3D
#  ifdef TS_DIF2
      real(r8), intent(inout) :: diff2(LBi:,LBj:,:)
#  endif
#  ifdef TS_DIF4
      real(r8), intent(inout) :: diff4(LBi:,LBj:,:)
#  endif
# endif
# if defined UV_VIS2 || !defined SOLVE3D
      real(r8), intent(inout) :: visc2_p(LBi:,LBj:)
      real(r8), intent(inout) :: visc2_r(LBi:,LBj:)
# endif
# ifdef UV_VIS4
      real(r8), intent(inout) :: visc4_p(LBi:,LBj:)
      real(r8), intent(inout) :: visc4_r(LBi:,LBj:)
# endif
#else
# ifdef SOLVE3D
#  ifdef TS_DIF2
      real(r8), intent(inout) :: diff2(LBi:UBi,LBj:UBj,NT(ng))
#  endif
#  ifdef TS_DIF4
      real(r8), intent(inout) :: diff4(LBi:UBi,LBj:UBj,NT(ng))
#  endif
# endif
# if defined UV_VIS2 || !defined SOLVE3D
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
      IF (SOUTH_WEST_TEST) THEN
        visc2(ng)=viscosity2
        visc4(ng)=viscosity4
#ifdef SOLVE3D
        DO itrc=1,NT(ng)
          tnu2(itrc,ng)=diffusion2(itrc)
          tnu4(itrc,ng)=diffusion4(itrc)
        END DO
#endif
      END IF
!
!-----------------------------------------------------------------------
!  Initialize horizontal mixing arrays to constant mixing coefficient.
!-----------------------------------------------------------------------
!
#ifdef _OPENMP
      IF (WESTERN_EDGE) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (EASTERN_EDGE) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (SOUTHERN_EDGE) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (NORTHERN_EDGE) THEN
        Jmax=BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=Jend
      END IF
#else
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
#endif
!
#if defined UV_VIS2 || !defined SOLVE3D
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          visc2_p(i,j)=viscosity2
          visc2_r(i,j)=viscosity2
        END DO
      END DO
#endif
#ifdef UV_VIS4
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          visc4_p(i,j)=viscosity4
          visc4_r(i,j)=viscosity4
        END DO
      END DO
#endif
#ifdef SOLVE3D
# ifdef TS_DIF2
      DO itrc=1,NT(ng)
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            diff2(i,j,itrc)=diffusion2(itrc)
          END DO
        END DO
      END DO
# endif
# ifdef TS_DIF4
      DO itrc=1,NT(ng)
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            diff4(i,j,itrc)=diffusion4(itrc)
          END DO
        END DO
      END DO
# endif
#endif

      RETURN
      END SUBROUTINE ini_hmixcoef_tile

      END MODULE ini_hmixcoef_mod
