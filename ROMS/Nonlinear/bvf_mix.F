#include "cppdefs.h"
      MODULE bvf_mix_mod
#if defined NONLINEAR && defined BVF_MIXING
!
!svn $Id: bvf_mix.F 732 2008-09-07 01:55:51Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes tracer vertical mixing as a function of the   !
!  Brunt-Vaisala frequency.  The vertical mixing of momentum  is set   !
!  to its background value.  If static unstable regime, the vertical   !
!  mixing is set to "bvf_nu0c".                                        !
!                                                                      !
!======================================================================!
!
      implicit none
!
      PRIVATE
      PUBLIC  :: bvf_mix
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE bvf_mix (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_mixing
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL bvf_mix_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   MIXING(ng) % bvf,                              &
     &                   MIXING(ng) % AKt,                              &
     &                   MIXING(ng) % AKv)
      RETURN
      END SUBROUTINE bvf_mix
!
!***********************************************************************
      SUBROUTINE bvf_mix_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         bvf, Akt, Akv)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_3d_mod, ONLY : exchange_w3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: bvf(LBi:,LBj:,0:)
      real(r8), intent(out) :: Akt(LBi:,LBj:,0:,:)
      real(r8), intent(out) :: Akv(LBi:,LBj:,0:)
# else
      real(r8), intent(in) :: bvf(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(out) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
      real(r8), intent(out) :: Akv(LBi:UBi,LBj:UBj,0:N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k

      real(r8) :: cff

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set tracer diffusivity as function of the Brunt-vaisala frequency.
!  Set vertical viscosity to its background value.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)-1
        DO j=Jstr,Jend
          DO i=Istr,Iend
            Akv(i,j,k)=Akv_bak(ng)
            IF (bvf(i,j,k).lt.0.0_r8) THEN
              Akv(i,j,k)=bvf_nu0c
              Akt(i,j,k,itemp)=bvf_nu0c
# ifdef SALINITY
              Akt(i,j,k,isalt)=bvf_nu0c
# endif
            ELSE IF (bvf(i,j,k).eq.0.0_r8) THEN
              Akv(i,j,k)=Akv_bak(ng)
              Akt(i,j,k,itemp)=Akt_bak(itemp,ng)
# ifdef SALINITY
              Akt(i,j,k,isalt)=Akt_bak(isalt,ng)
# endif
            ELSE
              cff=bvf_nu0/SQRT(bvf(i,j,k))
              Akt(i,j,k,itemp)=MIN(bvf_numax,MAX(bvf_numin,cff))
              Akv(i,j,k)=Akt(i,j,k,itemp)
# ifdef SALINITY
              Akt(i,j,k,isalt)=Akt(i,j,k,itemp)
# endif
            END IF
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          Akv)
        DO itrc=1,NAT
          CALL exchange_w3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 0, N(ng),         &
     &                            Akt(:,:,:,itrc))
        END DO
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Akv)
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng), 1, NAT,         &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Akt)
# endif

      RETURN
      END SUBROUTINE bvf_mix_tile
#endif
      END MODULE bvf_mix_mod
