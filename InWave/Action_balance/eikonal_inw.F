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
     &                      GRID(ng)%angler,                            &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
# endif
     &                      WAVEG(ng) % wd,                             &
     &                      WAVEP(ng) % cx,                             &
     &                      WAVEP(ng) % cy,                             &
     &                      WAVEP(ng) % h_tot,                          &
     &                      WAVEP(ng) % u_rho,                          &
     &                      WAVEP(ng) % v_rho,                          &
     &                      WAVEP(ng) % kwc)

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
     &                            angler, pm, pn,                       &
     &                            on_u, om_v,                           &
# ifdef MASKING
     &                            rmask,                                &
# endif
# ifdef WET_DRY
     &                            rmask_wet,                            &
# endif
     &                            wd, cx, cy, h_tot,                    &
     &                            u_rho, v_rho, kwc)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_inwave_params
      USE mod_inwave_swan
      USE bc_3d_mod
      USE exchange_3d_mod
      USE kwc3dbc_mod
#  ifdef DOPPLER_SHEAR
      USE mod_grid
      USE mod_ocean
#  endif

#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS

#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
#   ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: wd(:)
      real(r8), intent(in) :: cx(LBi:,LBj:,:)
      real(r8), intent(in) :: cy(LBi:,LBj:,:)
      real(r8), intent(in) :: h_tot(LBi:,LBj:)
      real(r8), intent(in) :: u_rho(LBi:,LBj:)
      real(r8), intent(in) :: v_rho(LBi:,LBj:)
      real(r8), intent(inout) :: kwc(LBi:,LBj:,:)
#  else
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: wd(ND)
      real(r8), intent(in) :: cx(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(in) :: cy(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(in) :: h_tot(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: u_rho(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: v_rho(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: kwc(LBi:UBi,LBj:UBj,ND)
#  endif
!
!  Local variable declarations.
!
      integer :: i, j, k, d

      real(r8) :: alfa_wave, kx, ky, cff, cff1, cff2
      real(r8) :: crk, crx, cry, kh, S
!     real(r8), parameter :: kwc_max = 10.0_r8
      real(r8), parameter :: kwc_max = 1.0_r8
      real(r8), parameter :: kwc_min = 0.015_r8
      real(r8), parameter :: eps = 1.0e-10_r8

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
!
#  include "set_bounds.h"
!
!=======================================================================
! Compute xi and etai components of the wave number
! Compute the time change of xi and etai components of the wave number
!=======================================================================
      DO d=1,ND
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            alfa_wave=(1.5_r8*pi-wd(d))-angler(i,j)
            cff1=kwc(i-1,j,d)*h_tot(i-1,j)
            cff2=kwc(i  ,j,d)*h_tot(i  ,j)
            crk=0.5_r8*(sqrt(MAX(g/kwc(i-1,j,d)*tanh(cff1),eps))+       &
     &                  sqrt(MAX(g/kwc(i  ,j,d)*tanh(cff2),eps)))
            cff1=0.5_r8*(u_rho(i-1,j)+u_rho(i,j))
#   ifdef DOPPLER_SHEAR
            cff=GRID(ng)%Hz(i,j,1)
            DO k=2,N(ng)
              cff=cff+GRID(ng)%Hz(i,j,k)
            END DO
            S=OCEAN(ng)%u(i,j,N(ng),1)-OCEAN(ng)%u(i,j,1,1)
!           crk=crk*sqrt(1.0_r8+S/cff/g*(cff1-cx(i,j,d)))
            crk=crk*sqrt(1.0_r8+S/cff/g*cx(i,j,d))
#   endif
            cff=(crk*cos(alfa_wave)+cff1)*on_u(i,j)
            cff1=MAX(cff,0.0_r8)
            cff2=MIN(cff,0.0_r8)
            FX(i,j)=cff1*kwc(i-1,j,d)+                                  &
     &              cff2*kwc(i  ,j,d)
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            alfa_wave=(1.5_r8*pi-wd(d))-angler(i,j)
            cff1=kwc(i,j-1,d)*h_tot(i,j-1)
            cff2=kwc(i,j  ,d)*h_tot(i,j  )
            crk=0.5_r8*(sqrt(MAX(g/kwc(i,j-1,d)*tanh(cff1),eps))+       &
     &                  sqrt(MAX(g/kwc(i,j  ,d)*tanh(cff2),eps)))
            cff1=0.5_r8*(v_rho(i,j-1)+v_rho(i,j))
            cff=(crk*sin(alfa_wave)+cff1)*om_v(i,j)
            cff1=MAX(cff,0.0_r8)
            cff2=MIN(cff,0.0_r8)
            FE(i,j)=cff1*kwc(i,j-1,d)+                                  &
     &              cff2*kwc(i,j  ,d)
          END DO
        END DO
!
!  Time-step horizontal advection.
!
        cff=dt(ng)*WAVES(ng)%dtifac
        DO j=Jstr,Jend
          DO i=Istr,Iend
            kwc(i,j,d)=kwc(i,j,d)-                                      &
     &                   cff*pm(i,j)*pn(i,j)*                           &
     &                  (FX(i+1,j)-FX(i,j)+                             &
     &                   FE(i,j+1)-FE(i,j))
            kwc(i,j,d)=MIN(kwc(i,j,d),kwc_max)
            kwc(i,j,d)=MAX(kwc_min,kwc(i,j,d))

!#  ifdef MASKING
!            kwc(i,j,d)=kwc(i,j,d)*rmask(i,j)
!#  endif
!#  ifdef WET_DRY
!            kwc(i,j,d)=kwc(i,j,d)*rmask_wet(i,j)
!#  endif
         END DO
        END DO
      END DO
!
!  Apply nonperiodic boundary conditions in xi and eta space.
!
      CALL kwc3dbc_tile (ng, tile,                                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     kwc)
!
!  Apply periodic boundary conditions.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, ND,              &
     &                          kwc)
      END IF
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
