#include "cppdefs.h"
      MODULE tl_t3drelax_mod

#if defined TL_IOMS && defined RPM_RELAXATION && defined SOLVE3D && \
    defined W4DVAR_SENSITIVITY
!
!svn $Id: tl_t3drelax.F 795 2016-05-11 01:42:43Z arango $
!=================================================== Andrew M. Moore ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group      Hernan G. Arango   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!
      implicit none

      PRIVATE

      PUBLIC :: tl_t3drelax

      CONTAINS

!
!***********************************************************************
      SUBROUTINE tl_t3drelax (ng, tile)
!***********************************************************************
!
      USE mod_param
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
      CALL wclock_on (ng, iTLM, 24)
# endif
      CALL tl_t3drelax_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nrhs(ng), nnew(ng),                        &
# ifdef MASKING
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
# endif
     &                       GRID(ng) % Hz,                             &
     &                       GRID(ng) % pmon_u,                         &
     &                       GRID(ng) % pnom_v,                         &
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pn,                             &
     &                       OCEAN(ng) % tl_t)
# ifdef PROFILE
      CALL wclock_off (ng, iTLM, 24)
# endif
      RETURN
      END SUBROUTINE tl_t3drelax
!
!***********************************************************************
      SUBROUTINE tl_t3drelax_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nrhs, nnew,                          &
# ifdef MASKING
     &                             umask, vmask,                        &
# endif
     &                             Hz,                                  &
     &                             pmon_u, pnom_v, pm, pn,              &
     &                             tl_t)
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
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: pmon_u(LBi:,LBj:)
      real(r8), intent(in) :: pnom_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)

      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX

# include "set_bounds.h"
!
!----------------------------------------------------------------------
!  Compute horizontal diffusion relaxation of tracers between current
!  and previous representer tangent linear Picard iteration trajectory
!  (basic state).
!----------------------------------------------------------------------
!
!  This is the tangent linear of the relaxation terms that appear in
!  the representer model. Since the representer model is linear, we
!  DO NOT need to include the tangent linear of the Hz terms below
!  because Hz is computed from the background trajectory (AMM).
!
      DO itrc=1,NT(ng)
        IF (tl_Tdiff(itrc,ng).gt.0.0_r8) THEN
          DO k=1,N(ng)
!
!  Compute XI- and ETA-components of diffusive tracer flux (T m3/s).
!
            DO j=Jstr,Jend
              DO i=Istr,Iend+1
                FX(i,j)=0.5_r8*tl_Tdiff(itrc,ng)*pmon_u(i,j)*           &
     &                  (Hz(i,j,k)+Hz(i-1,j,k))*                        &
     &                  (tl_t(i  ,j,k,nrhs,itrc)-                       &
     &                   tl_t(i-1,j,k,nrhs,itrc))
# ifdef MASKING
                FX(i,j)=FX(i,j)*umask(i,j)
# endif
              END DO
            END DO
            DO j=Jstr,Jend+1
              DO i=Istr,Iend
                FE(i,j)=0.5_r8*tl_Tdiff(itrc,ng)*pnom_v(i,j)*           &
     &                  (Hz(i,j,k)+Hz(i,j-1,k))*                        &
     &                  (tl_t(i,j  ,k,nrhs,itrc)-                       &
     &                   tl_t(i,j-1,k,nrhs,itrc))
# ifdef MASKING
                FE(i,j)=FE(i,j)*vmask(i,j)
# endif
              END DO
            END DO
!
! Time-step horizontal diffusion relaxation term (m Tunits).
!
            DO j=Jstr,Jend
              DO i=Istr,Iend
                tl_t(i,j,k,nnew,itrc)=tl_t(i,j,k,nnew,itrc)+            &
     &                                dt(ng)*pm(i,j)*pn(i,j)*           &
     &                                (FX(i+1,j)-FX(i,j)+               &
     &                                 FE(i,j+1)-FE(i,j))
              END DO
            END DO
          END DO
        END IF
      END DO
      RETURN
      END SUBROUTINE tl_t3drelax_tile
#endif
      END MODULE tl_t3drelax_mod
