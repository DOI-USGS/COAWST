#include "cppdefs.h"
      MODULE ADfromTL_mod

#if defined SP4DVAR
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine initializes the adjoint state with the TL state        !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ADfromTL

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ADfromTL (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_clima
      USE mod_grid
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
      CALL ADfromTL_tile (ng, tile,                                     &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         Lold(ng), knew(ng),                      &
# ifdef SOLVE3D
     &                         nstp(ng),                                &
# endif
# ifdef MASKING
     &                         GRID(ng) % rmask,                        &
     &                         GRID(ng) % umask,                        &
     &                         GRID(ng) % vmask,                        &
# endif
# ifdef SOLVE3D
     &                         OCEAN(ng) % tl_u,                        &
     &                         OCEAN(ng) % tl_v,                        &
     &                         OCEAN(ng) % tl_t,                        &
# endif
# ifndef SOLVE3D
     &                         OCEAN(ng) % tl_ubar,                     &
     &                         OCEAN(ng) % tl_vbar,                     &
# endif
     &                         OCEAN(ng) % tl_zeta,                     &
# ifdef SOLVE3D
     &                         OCEAN(ng) % ad_u,                        &
     &                         OCEAN(ng) % ad_v,                        &
     &                         OCEAN(ng) % ad_t,                        &
# endif
# ifndef SOLVE3D
     &                         OCEAN(ng) % ad_ubar,                     &
     &                         OCEAN(ng) % ad_vbar,                     &
# endif
     &                         OCEAN(ng) % ad_zeta)

      RETURN
      END SUBROUTINE ADfromTL
!
!***********************************************************************
      SUBROUTINE ADfromTL_tile (ng, tile,                               &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               Linp, kout,                        &
# ifdef SOLVE3D
     &                               nout,                              &
# endif
# ifdef MASKING
     &                               rmask, umask, vmask,               &
# endif
# ifdef SOLVE3D
     &                               tl_u, tl_v, tl_t,                  &
# endif
# ifndef SOLVE3D
     &                               tl_ubar, tl_vbar,                  &
# endif
     &                               tl_zeta,                           &
# ifdef SOLVE3D
     &                               ad_u, ad_v, ad_t,                  &
# endif
# ifndef SOLVE3D
     &                               ad_ubar, ad_vbar,                  &
# endif
     &                               ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: kout, Linp
# ifdef SOLVE3D
      integer, intent(in) :: nout
# endif
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
#  endif
#  ifndef SOLVE3D
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
#  endif
#  ifndef SOLVE3D
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#  endif
#  ifndef SOLVE3D
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#  endif
#  ifndef SOLVE3D
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k

# include "set_bounds.h"
!
!------------------------------------------------------------------------
!  Initialize adjoint using the TL state in record 1.
!------------------------------------------------------------------------
!
!  Free-surface.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ad_zeta(i,j,kout)=tl_zeta(i,j,Linp)
#ifdef MASKING
          ad_zeta(i,j,kout)=ad_zeta(i,j,kout)*rmask(i,j)
#endif
        END DO
      END DO
#ifndef SOLVE3D
!
!  2D Momentum.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ad_ubar(i,j,kout)=tl_ubar(i,j,Linp)
#ifdef MASKING
          ad_ubar(i,j,kout)=ad_ubar(i,j,kout)*umask(i,j)
#endif
        END DO
      END DO
!
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          ad_vbar(i,j,kout)=tl_vbar(i,j,Linp)
#ifdef MASKING
          ad_vbar(i,j,kout)=ad_vbar(i,j,kout)*vmask(i,j)
#endif
        END DO
      END DO
#endif
# ifdef SOLVE3D
!
!  3D Momentum.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_u(i,j,k,nout)=tl_u(i,j,k,Linp)
#ifdef MASKING
            ad_u(i,j,k,nout)=ad_u(i,j,k,nout)*umask(i,j)
#endif
          END DO
        END DO
      END DO
!
      DO k=1,N(ng)
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_v(i,j,k,nout)=tl_v(i,j,k,Linp)
#ifdef MASKING
            ad_v(i,j,k,nout)=ad_v(i,j,k,nout)*vmask(i,j)
#endif
          END DO
        END DO
      END DO
!
!  Tracers.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              ad_t(i,j,k,nout,itrc)=tl_t(i,j,k,Linp,itrc)
#ifdef MASKING
              ad_t(i,j,k,nout,itrc)=ad_t(i,j,k,nout,itrc)*rmask(i,j)
#endif
            END DO
          END DO
        END DO
      END DO
# endif

      RETURN
      END SUBROUTINE ADfromTL_tile
#endif
      END MODULE ADfromTL_mod
