#include "cppdefs.h"
      MODULE eikonal_inw_mod
#if defined INWAVE_MODEL
!
!svn $Id: eikonal_inw.F 732 2008-09-07 01:55:51Z jcwarner $

!======================================================================!
!                                                                      !
!  This routine computes the temporal change on the wave number        !
!  created by spatial changes of the absolute frequency.               !
!  @kxi/@t=- m*@wa/@xi                                                 !
!  @keta/@t=- m*@wa/@eta                                               !
!                                                                      !
!======================================================================!
!
      implicit none
      PRIVATE
      PUBLIC  :: eikonal_inw
      CONTAINS
!
!***********************************************************************
      SUBROUTINE eikonal_inw (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_grid
      USE mod_stepping
      USE mod_inwave_vars
      USE mod_inwave_params
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#  include "tile.h"
!
!#  ifdef PROFILE
!      CALL wclock_on (ng, iNLM, 35)
!#  endif

      CALL eikonal_inw_tile(ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
     &                      WAVEP(ng) % kwc,WAVEP(ng) % Ta,             &
     &                      WAVEG(ng) % wd,GRID(ng)%angler)

!#  ifdef PROFILE
!      CALL wclock_off (ng, iNLM, 35)
!#  endif

      RETURN
      END SUBROUTINE eikonal_inw
!
!***********************************************************************
      SUBROUTINE eikonal_inw_tile(ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            pm, pn, kwc, Ta,                      &
     &                            wd, angler)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_inwave_params
      USE bc_3d_mod
      USE exchange_3d_mod

#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      real(r8), intent(inout) :: wd(ND)

#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(inout) :: kwc(LBi:,LBj:,:)
      real(r8), intent(inout) :: Ta(LBi:,LBj:,:)
#  else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(inout) :: kwc(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(inout) :: Ta(LBi:UBi,LBj:UBj,ND)
#  endif
!
!  Local variable declarations.
!
      integer :: i, j, k, d

      real(r8) :: twopi, otwopi
      real(r8) :: alfa_wave, kx, ky

!
#  include "set_bounds.h"
!
      twopi=2.0_r8*pi
      otwopi=1.0_r8/twopi

!     DO j=Jstr+1,Jend-1
!       DO i=Istr+1,Iend-1
      DO j=Jstr,Jend
        DO i=Istr,Iend
          DO d=1,ND

!=======================================================================
! Compute xi and etai components of the wave number
!=======================================================================
!
            alfa_wave =(1.5_r8*pi-wd(d))-angler(i,j)
            kx= kwc(i,j,d)*cos(alfa_wave)
            ky= kwc(i,j,d)*sin(alfa_wave)
!
!=======================================================================
! Compute the time change of xi and etai components of the wave number
!=======================================================================
!
            kx=kx-pm(i,j)*dt(ng)*0.5_r8*                                &
     &         (twopi/Ta(i+1,j,d)-twopi/Ta(i-1,j,d))
            ky=ky-pn(i,j)*dt(ng)*0.5_r8*                                &
     &         (twopi/Ta(i,j+1,d)-twopi/Ta(i,j-1,d))
            kwc(i,j,d)=max(((kx*kx+ky*ky)**0.5_r8),0.000001_r8)
          END DO
        END DO
      END DO

!      DO j=Jstr+1,Jend-1
!        DO d=1,ND(ng)
!            kwc(Istr,j,d)=kwc(Istr+1,j,d)
!            kwc(Iend,j,d)=kwc(Iend-1,j,d)
!        END DO
!      END DO

!      DO i=Istr+1,Iend-1
!        DO d=1,ND(ng)
!          kwc(i,Jstr,d)=kwc(i,Jstr+1,d)
!          kwc(i,Jend,d)=kwc(i,Jend-1,d)
!        END DO
!      END DO

      CALL bc_r3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, ND,                      &
     &                  kwc)

!
# ifdef DISTRIBUTE
!
! Exchange boundary data.
!
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, ND,                    &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    kwc)
# endif

      RETURN
      END SUBROUTINE eikonal_inw_tile
#endif
      END MODULE eikonal_inw_mod
