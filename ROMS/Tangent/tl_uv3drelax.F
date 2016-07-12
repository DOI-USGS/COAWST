#include "cppdefs.h"
      MODULE tl_uv3drelax_mod

#if defined TL_IOMS && defined RPM_RELAXATION && defined SOLVE3D && \
    defined W4DVAR_SENSITIVITY
!
!svn $Id: tl_uv3drelax.F 795 2016-05-11 01:42:43Z arango $
!=================================================== Andrew M. Moore ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group      Hernan G. Arango   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine relaxes current  representer  tangent linear 3D        !
!  momentum to previous Picard iteration solution (basic state)        !
!  to improve stability and convergence.                               !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE

      PUBLIC tl_uv3drelax

      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_uv3drelax (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
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
# ifdef PROFILE
      CALL wclock_on (ng, iTLM, 30)
# endif
      CALL tl_uv3drelax_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nrhs(ng), nnew(ng),                       &
# ifdef MASKING
     &                        GRID(ng) % pmask,                         &
# endif
     &                        GRID(ng) % Hz,                            &
     &                        GRID(ng) % pm,                            &
     &                        GRID(ng) % pmon_p,                        &
     &                        GRID(ng) % pmon_r,                        &
     &                        GRID(ng) % pn,                            &
     &                        GRID(ng) % pnom_p,                        &
     &                        GRID(ng) % pnom_r,                        &
     &                        OCEAN(ng) % tl_u,                         &
     &                        OCEAN(ng) % tl_v)
# ifdef PROFILE
      CALL wclock_off (ng, iTLM, 30)
# endif
      RETURN
      END SUBROUTINE tl_uv3drelax

!
!***********************************************************************
      SUBROUTINE tl_uv3drelax_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              nrhs, nnew,                         &
# ifdef MASKING
     &                              pmask,                              &
# endif
     &                              Hz,                                 &
     &                              pm, pmon_p, pmon_r,                 &
     &                              pn, pnom_p, pnom_r,                 &
     &                              tl_u, tl_v)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, nnew

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pmon_p(LBi:,LBj:)
      real(r8), intent(in) :: pmon_r(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: pnom_p(LBi:,LBj:)
      real(r8), intent(in) :: pnom_r(LBi:,LBj:)

      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_r(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8) :: cff

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFx

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute horizontal diffusion relaxation of 3D momentum between
!  current and previous representer tangent linear Picard iteration
!  trajectory (basic state).
!-----------------------------------------------------------------------
!
!  This is the tangent linear of the relaxation terms that appear in
!  the representer model. Since the representer model is linear, we
!  DO NOT need to include the tangent linear of the Hz terms below
!  because Hz is computed from the background trajectory (AMM).
!
      IF (tl_M3diff(ng).gt.0.0_r8) THEN
!
        K_LOOP : DO k=1,N(ng)
!
!  Compute flux-components of the diffusive relaxation (m4/s2) in XI-
!  and ETA-directions.
!
          DO j=Jstr,Jend
            DO i=IstrU-1,Iend
              UFx(i,j)=tl_M3diff(ng)*pmon_r(i,j)*                       &
     &                 Hz(i,j,k)*                                       &
     &                 (tl_u(i+1,j,k,nrhs)-                             &
     &                  tl_u(i  ,j,k,nrhs))
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=IstrU,Iend
              UFe(i,j)=tl_M3diff(ng)*pnom_p(i,j)*                       &
     &                 0.25_r8*(Hz(i,j  ,k)+Hz(i-1,j  ,k)+              &
     &                          Hz(i,j-1,k)+Hz(i-1,j-1,k))*             &
     &                 (tl_u(i,j  ,k,nrhs)-                             &
     &                  tl_u(i,j-1,k,nrhs))
# ifdef MASKING
              UFe(i,j)=UFe(i,j)*pmask(i,j)
# endif
            END DO
          END DO
          DO j=JstrV,Jend
            DO i=Istr,Iend+1
              VFx(i,j)=tl_M3diff(ng)*pmon_p(i,j)*                       &
     &                 0.25_r8*(Hz(i,j  ,k)+Hz(i-1,j  ,k)+              &
     &                          Hz(i,j-1,k)+Hz(i-1,j-1,k))*             &
     &                 (tl_v(i  ,j,k,nrhs)-                             &
     &                  tl_v(i-1,j,k,nrhs))
# ifdef MASKING
              VFx(i,j)=VFx(i,j)*pmask(i,j)
# endif
            END DO
          END DO
          DO j=JstrV-1,Jend
            DO i=Istr,Iend
              VFe(i,j)=tl_M3diff(ng)*pnom_r(i,j)*                       &
     &                 Hz(i,j,k)*                                       &
     &                 (tl_v(i,j+1,k,nrhs)-                             &
     &                  tl_v(i,j  ,k,nrhs))
            END DO
          END DO
!
! Time-step diffusive relaxation term.  Notice that momentum at this
! stage is HzU and HzV and has m2/s units.
!
          cff=dt(ng)*0.25_r8
          DO j=Jstr,Jend
            DO i=IstrU,Iend
              tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)+                        &
     &                         cff*(pm(i-1,j)+pm(i,j))*                 &
     &                             (pn(i-1,j)+pn(i,j))*                 &
     &                         (UFx(i,j)-UFx(i-1,j)+                    &
     &                          UFe(i,j+1)-UFe(i,j))

            END DO
          END DO
          DO j=JstrV,Jend
            DO i=Istr,Iend
              tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)+                        &
     &                         cff*(pm(i,j)+pm(i,j-1))*                 &
     &                             (pn(i,j)+pn(i,j-1))*                 &
     &                         (VFx(i+1,j)-VFx(i,j)+                    &
     &                          VFe(i,j)-VFe(i,j-1))
            END DO
          END DO
        END DO K_LOOP
      END IF
      RETURN
      END SUBROUTINE tl_uv3drelax_tile
#endif
      END MODULE tl_uv3drelax_mod
