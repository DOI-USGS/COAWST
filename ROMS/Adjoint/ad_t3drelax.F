#include "cppdefs.h"
      MODULE ad_t3drelax_mod

#if defined TL_IOMS && defined RPM_RELAXATION && defined SOLVE3D
!
!svn $Id: ad_t3drelax.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine is the adjoint of relaxing current representer tangent !
!  linear tracer type variables to previous Picard iteration  solution !
!  (basic state) to improve stability and convergence.                 !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE

      PUBLIC :: ad_t3drelax

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_t3drelax (ng, tile)
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
      CALL wclock_on (ng, iADM, 24)
# endif
      CALL ad_t3drelax_tile (ng, tile,                                  &
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
     &                       OCEAN(ng) % ad_t)
# ifdef PROFILE
      CALL wclock_off (ng, iADM, 24)
# endif
      RETURN
      END SUBROUTINE ad_t3drelax
!
!***********************************************************************
      SUBROUTINE ad_t3drelax_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nrhs, nnew,                          &
# ifdef MASKING
     &                             umask, vmask,                        &
# endif
     &                             Hz,                                  &
     &                             pmon_u, pnom_v, pm, pn,              &
     &                             t, ad_t)
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

      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
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

      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k
      real(r8) :: adfac

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FX

# include "set_bounds.h"
!
!----------------------------------------------------------------------
!  Compute horizontal diffusion relaxation of tracers between current
!  and previous representer tangent linear Picard iteration trajectory
!  (basic state).
!----------------------------------------------------------------------
!
      ad_FX=0.0_r8
      ad_FE=0.0_r8
!
      DO itrc=1,NT(ng)
        IF (tl_Tdiff(itrc,ng).gt.0.0_r8) THEN
          DO k=1,N(ng)
!
! Time-step horizontal adjoint diffusion relaxation term (m Tunits).
!
            DO j=Jstr,Jend
              DO i=Istr,Iend
!>              tl_t(i,j,k,nnew,itrc)=tl_t(i,j,k,nnew,itrc)+            &
!>   &                                dt(ng)*pm(i,j)*pn(i,j)*           &
!>   &                                (FX(i+1,j)-FX(i,j)+               &
!>   &                                 FE(i,j+1)-FE(i,j))
!>
                adfac=dt(ng)*pm(i,j)*pn(i,j)*ad_t(i,j,k,nnew,itrc)
                ad_FX(i+1,j)=ad_FX(i+1,j)+adfac
                ad_FX(i,j)=ad_FX(i,j)-adfac
                ad_FE(i,j+1)=ad_FE(i,j+1)+adfac
                ad_FE(i,j)=ad_FE(i,j)+adfac
              END DO
            END DO
!
!  Compute XI- and ETA-components of adjoint diffusive tracer flux
!  (T m3/s).
!
            DO j=Jstr,Jend+1
              DO i=Istr,Iend
# ifdef MASKING
!>              FE(i,j)=FE(i,j)*vmask(i,j)
!>
                ad_FE(i,j)=ad_FE(i,j)*vmask(i,j)
# endif
!>              FE(i,j)=0.5_r8*tl_Tdiff(itrc,ng)*pnom_v(i,j)*           &
!>   &                  (Hz(i,j,k)+Hz(i,j-1,k))*                        &
!>   &                  (tl_t(i,j  ,k,nrhs,itrc)-t(i,j  ,k,nrhs,itrc)-  &
!>   &                   tl_t(i,j-1,k,nrhs,itrc)+t(i,j-1,k,nrhs,itrc))
!>
                adfac=0.5_r8*tl_Tdiff(itrc,ng)*pnom_v(i,j)*             &
     &                (Hz(i,j,k)+Hz(i,j-1,k))*ad_FE(i,j)
                ad_t(i,j  ,k,nrhs,itrc)=ad_t(i,j  ,k,nrhs,itrc)+adfac
                ad_t(i,j-1,k,nrhs,itrc)=ad_t(i,j-1,k,nrhs,itrc)-adfac
                ad_FE(i,j)=0.0_r8
              END DO
            END DO
!
            DO j=Jstr,Jend
              DO i=Istr,Iend+1
# ifdef MASKING
!>              FX(i,j)=FX(i,j)*umask(i,j)
!>
                ad_FX(i,j)=ad_FX(i,j)*umask(i,j)
# endif
!>              FX(i,j)=0.5_r8*tl_Tdiff(itrc,ng)*pmon_u(i,j)*           &
!>   &                  (Hz(i,j,k)+Hz(i-1,j,k))*                        &
!>   &                  (tl_t(i  ,j,k,nrhs,itrc)-t(i  ,j,k,nrhs,itrc)-  &
!>   &                   tl_t(i-1,j,k,nrhs,itrc)+t(i-1,j,k,nrhs,itrc))
!>
                adfac=0.5_r8*tl_Tdiff(itrc,ng)*pmon_u(i,j)*             &
     &                (Hz(i,j,k)+Hz(i-1,j,k))*ad_FX(i,j)
                ad_t(i  ,j,k,nrhs,itrc)=ad_t(i  ,j,k,nrhs,itrc)+adfac
                ad_t(i-1,j,k,nrhs,itrc)=ad_t(i-1,j,k,nrhs,itrc)-adfac
                ad_FX(i,j)=0.0_r8
              END DO
            END DO
          END DO
        END IF
      END DO

      RETURN
      END SUBROUTINE ad_t3drelax_tile
#endif
      END MODULE ad_t3drelax_mod
