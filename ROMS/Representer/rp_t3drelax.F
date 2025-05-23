#include "cppdefs.h"
      MODULE rp_t3drelax_mod

#if defined TL_IOMS && defined RPM_RELAXATION && defined SOLVE3D
!
!git $Id$
!svn $Id: rp_t3drelax.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine relaxes current representer tangent linear tracer      !
!  type variables to previous  Picard  iteration  solution (basic      !
!  state) to improve stability and convergence.                        !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
!
      PUBLIC :: rp_t3drelax
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE rp_t3drelax (ng, tile)
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
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iRPM, 24, __LINE__, MyFile)
# endif
      CALL rp_t3drelax_tile (ng, tile,                                  &
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
     &                       OCEAN(ng) % t,                             &
     &                       OCEAN(ng) % tl_t)
# ifdef PROFILE
      CALL wclock_off (ng, iRPM, 24, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE rp_t3drelax
!
!***********************************************************************
      SUBROUTINE rp_t3drelax_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nrhs, nnew,                          &
# ifdef MASKING
     &                             umask, vmask,                        &
# endif
     &                             Hz,                                  &
     &                             pmon_u, pnom_v, pm, pn,              &
     &                             t, tl_t)
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
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)

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
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))

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
     &                  (tl_t(i  ,j,k,nrhs,itrc)-t(i  ,j,k,nrhs,itrc)-  &
     &                   tl_t(i-1,j,k,nrhs,itrc)+t(i-1,j,k,nrhs,itrc))
# ifdef MASKING
                FX(i,j)=FX(i,j)*umask(i,j)
# endif
              END DO
            END DO
            DO j=Jstr,Jend+1
              DO i=Istr,Iend
                FE(i,j)=0.5_r8*tl_Tdiff(itrc,ng)*pnom_v(i,j)*           &
     &                  (Hz(i,j,k)+Hz(i,j-1,k))*                        &
     &                  (tl_t(i,j  ,k,nrhs,itrc)-t(i,j  ,k,nrhs,itrc)-  &
     &                   tl_t(i,j-1,k,nrhs,itrc)+t(i,j-1,k,nrhs,itrc))
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
!
      RETURN
      END SUBROUTINE rp_t3drelax_tile
#endif
      END MODULE rp_t3drelax_mod
