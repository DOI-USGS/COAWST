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
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
     &                      GRID(ng) % om_u,                            &
     &                      GRID(ng) % on_v,                            &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask_wet,                       &
     &                      GRID(ng) % vmask_wet,                       &
# endif
     &                      WAVEP(ng) % cx,                             &
     &                      WAVEP(ng) % cy,                             &
     &                      WAVEP(ng) % kwc, WAVEP(ng) % Ta,            &
     &                      WAVEG(ng) % wd, GRID(ng)%angler)

!#  ifdef PROFILE
!      CALL wclock_off (ng, iNLM, 35)
!#  endif

      RETURN
      END SUBROUTINE eikonal_inw
!
!***********************************************************************
      SUBROUTINE eikonal_inw_tile(ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            pm, pn,                               &
     &                            om_u, on_v,                           &
# ifdef MASKING
     &                            rmask, umask, vmask,                  &
# endif
     &                            cx, cy, kwc, Ta,                      &
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
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      real(r8), intent(inout) :: wd(ND)

#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: cx(LBi:,LBj:,:)
      real(r8), intent(in) :: cy(LBi:,LBj:,:)
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(inout) :: kwc(LBi:,LBj:,:)
      real(r8), intent(in) :: Ta(LBi:,LBj:,:)
#  else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: cx(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(in) :: cy(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(inout) :: kwc(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(in) :: Ta(LBi:UBi,LBj:UBj,ND)
#  endif
!
!  Local variable declarations.
!
      integer :: i, j, k, d

      real(r8) :: twopi, otwopi, Ta_min
      real(r8) :: alfa_wave, kx, ky, cff, cff1, cff2, cff3

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: CF
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: mask_CF
!
#  include "set_bounds.h"
!
      twopi=2.0_r8*pi
      otwopi=1.0_r8/twopi
      Ta_min=1.0_r8
!
      DO d=1,ND
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
            IF ((Ta(i,j,d)-0.0001_r8).lt.(Ta_min)) THEN
              mask_CF(i,j)=0.0_r8
            ELSE
              mask_CF(i,j)=1.0_r8
            END IF
            cff=MAX(Ta_min,Ta(i,j,d))
            CF(i,j)=twopi/cff
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend
!
!=======================================================================
! Compute xi and etai components of the wave number
!=======================================================================
!
!           alfa_wave=(1.5_r8*pi-wd(d))-angler(i,j)
            IF (cx(i,j,d).eq.0.0_r8) THEN
              alfa_wave=0.5_r8*pi*SIGN(1.0_r8,cy(i,j,d))
            ELSE
              alfa_wave=ATAN2(cy(i,j,d),cx(i,j,d))
            ENDIF
            kx=kwc(i,j,d)*COS(alfa_wave)
            ky=kwc(i,j,d)*SIN(alfa_wave)
!
!=======================================================================
! Compute the time change of xi and etai components of the wave number
!=======================================================================
!
!           kx=kx-pm(i,j)*dt(ng)*0.5_r8*                                &
!    &         (twopi/Ta(i+1,j,d)-twopi/Ta(i-1,j,d))
!           ky=ky-pn(i,j)*dt(ng)*0.5_r8*                                &
!    &         (twopi/Ta(i,j+1,d)-twopi/Ta(i,j-1,d))
!
            cff=0.5_r8*pm(i,j)*dt(ng)*                                  &
     &          (CF(i+1,j)-CF(i-1,j))*mask_CF(i+1,j)*mask_CF(i-1,j)
            kx=kx-cff*mask_CF(i,j)
!
            cff=0.5_r8*pn(i,j)*dt(ng)*                                  &
     &         (CF(i,j+1)-CF(i,j-1))*mask_CF(i,j+1)*mask_CF(i,j-1)
            ky=ky-cff*mask_CF(i,j)
!
            kwc(i,j,d)=MAX(SQRT(kx*kx+ky*ky),0.000001_r8)
          END DO
        END DO
      END DO
!
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
