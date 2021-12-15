#include "cppdefs.h"
      MODULE dissip_inw_mod

#if defined INWAVE_MODEL
# if defined WDISS_ROELVINK || defined WDISS_GAMMA
!
!=======================================================================
!                                                                      !
!  This routine computes the energy dissipation                        !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: dissip_inw_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE dissip_inw (ng, tile, nout)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_inwave_params
      USE mod_inwave_vars
      USE mod_ocean
      USE mod_coupling
      USE mod_stepping
      USE mod_forces
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, nout
!
!  Local variable declarations.
!
# include "tile.h"
!
!#  ifdef PROFILE
!      CALL wclock_on (ng, iNLM, 35)
!#  endif

      CALL dissip_inw_tile(ng, tile,                                    &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 nstp(ng), nout,                                  &
# ifdef MASKING
     &                 GRID(ng) % rmask,                                &
# endif
# ifdef WET_DRY
     &                 GRID(ng) % rmask_wet,                            &
     &                 GRID(ng) % umask_wet,                            &
     &                 GRID(ng) % vmask_wet,                            &
# endif
     &                 FORCES(ng) % Dissip_break,                       &
     &                 FORCES(ng) % Dissip_wcap,                        &
     &                 FORCES(ng) % Lwave,                              &
# ifdef WAVES_BOT_PERIOD
     &                 FORCES(ng) % Pwave_bot,                          &
# endif
     &                 WAVEP(ng) % h_tot,                               &
     &                 WAVEP(ng) % AC,                                  &
     &                 WAVEP(ng) % Tr,                                  &
     &                 WAVEP(ng) % kwc)
!#  ifdef PROFILE
!      CALL wclock_off (ng, iNLM, 35)
!#  endif
      RETURN
      END SUBROUTINE dissip_inw
!
!***********************************************************************
      SUBROUTINE dissip_inw_tile(ng, tile,                              &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nstp, nout,                                &
# ifdef MASKING
     &                       rmask,                                     &
# endif
# ifdef WET_DRY
     &                       rmask_wet, umask_wet, vmask_wet,           &
# endif
     &                       Dissip_break,                              &
     &                       Dissip_wcap,                               &
     &                       Lwave,                                     &
# ifdef WAVES_BOT_PERIOD
     &                       Pwave_bot,                                 &
# endif
     &                       h_tot,                                     &
     &                       AC, Tr, kwc)
!***********************************************************************
!
      USE mod_param
      USE mod_inwave_params
      USE mod_inwave_swan
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
      USE exchange_2d_mod

#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  endif
      USE bc_2d_mod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nout
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
#   ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
#   endif
      real(r8), intent(inout) :: Dissip_break(LBi:,LBj:)
      real(r8), intent(inout) :: Dissip_wcap(LBi:,LBj:)
      real(r8), intent(in) :: Lwave(LBi:,LBj:)
#   ifdef WAVES_BOT_PERIOD
      real(r8), intent(in) :: Pwave_bot(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: h_tot(LBi:,LBj:)
      real(r8), intent(inout) :: AC(LBi:,LBj:,:,:)
      real(r8), intent(in) :: Tr(LBi:,LBj:,:)
      real(r8), intent(in) :: kwc(LBi:,LBj:,:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(inout) :: Dissip_break(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Dissip_wcap(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Lwave(LBi:UBi,LBj:UBj)
#   ifdef WAVES_BOT_PERIOD
      real(r8), intent(in) :: Pwave_bot(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: h_tot(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: AC(LBi:UBi,LBj:UBj,ND,3)
      real(r8), intent(in) :: Tr(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(in) :: kwc(LBi:UBi,LBj:UBj,ND)
#  endif
!
!  Local variable declarations.
!
      integer :: i, j, d
      real(r8) :: EW, oEW, TRM, H, Qb, Hmax_r, diff, Emax_r
      real(r8) :: twopi, otwopi, ogrho0, cff, orho0, fac
#  ifdef WDISS_FRIC
      real(r8) :: gamma, kd, dissfric
      real(r8), parameter :: eps=1.0e-08_r8
#  endif
      real(r8), parameter :: Trmin=1.0_r8
      real(r8), parameter :: alfa=1.0_r8
      real(r8), parameter :: breakr=0.45_r8
      real(r8), parameter :: breakg=0.45_r8
      real(r8), parameter :: n_r=15.0_r8
      real(r8), parameter :: EWlim=0.00001_r8
      real(r8):: EN(ND)

# include "set_bounds.h"

      twopi=2.0_r8*pi
      otwopi=1.0_r8/twopi
      ogrho0=1.0_r8/(g*rho0)
      orho0=1.0_r8/rho0
      fac=1.0_r8/(dt(ng)*WAVES(ng)%dtifac)
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          EW=0.0_r8
#  ifdef WDISS_ROELVINK
          TRM=0.0_r8
#  endif
          DO d=1,ND
!=======================================================================
!  Compute the energy from action balance and wave heigth
!=======================================================================
            EN(d)=AC(i,j,d,nout)*twopi/MAX(Trmin,Tr(i,j,d))
!=======================================================================
!  Compute the total energy
!=======================================================================
            EW=EW+EN(d)
!=======================================================================
!  Compute the mean wave number and intrinsic periods
!  What we do is give more importance to those wave 
!  numbers with more energy
!=======================================================================
#  ifdef WDISS_ROELVINK
            TRM=TRM+Tr(i,j,d)*EN(d)
#  endif
          ENDDO
#  ifdef WDISS_ROELVINK
          cff=1.0_r8/(max(EW,EWlim))
          TRM=TRM*cff
#  endif
!         EW=MAX(EW,EWlim)
          EW=MAX(EW,0.0_r8)  !this was needed
#  ifdef WDISS_ROELVINK
!=======================================================================
!  Compute the wave height. This is based on Hrms.
!=======================================================================
          H=(8.0_r8*EW*ogrho0)**0.5_r8
#  endif
!=======================================================================
!  Compute the energy dissipation
!=======================================================================
          IF (h_tot(i,j).ge.Dcrit(ng)) THEN
#  ifdef WDISS_FRIC
            gamma=0.067_r8/g**2
            kD=twopi/(Lwave(i,j)+eps)*h_tot(i,j)
            dissfric=dt(ng)*gamma*EW/((sinh(kD))**2)*twopi*twopi/       &
     &         MAX(Trmin,Pwave_bot(i,j)*Pwave_bot(i,j))*WAVES(ng)%dtifac
#  endif
#  ifdef WDISS_ROELVINK
            Hmax_r=breakr*(MAX(h_tot(i,j),0.0_r8))
            Qb=MIN(1.0_r8,1.0_r8-EXP(-(H/Hmax_r)**n_r))
            IF (TRM.gt.0.0001_r8) THEN
              Dissip_break(i,j)=2.0_r8*alfa/TRM*EW*Qb*dt(ng)*           &
     &                          WAVES(ng)%dtifac
            ELSE
              Dissip_break(i,j)=0.0_r8
            END IF
#  elif defined WDISS_GAMMA
            Hmax_r=breakg*(MAX(h_tot(i,j),0.0_r8))
            Emax_r=0.125_r8*g*rho0*Hmax_r**2.0_r8
            diff=EW-Emax_r
            Dissip_break(i,j)=MAX(0.0_r8,diff)
#  endif
          ELSE
            Dissip_break(i,j)=0.0_r8
#  ifdef WDISS_FRIC
            dissfric=0.0_r8
#  endif
          END IF
#  ifdef MASKING
          Dissip_break(i,j)=Dissip_break(i,j)*rmask(i,j)
#   ifdef WDISS_FRIC
          dissfric=dissfric*rmask(i,j)
#   endif
#  endif
#  ifdef WET_DRY
          Dissip_break(i,j)=Dissip_break(i,j)*rmask_wet(i,j)
#   ifdef WDISS_FRIC
          dissfric=dissfric*rmask_wet(i,j)
#   endif
#  endif
!=======================================================================
!  Distribute dissipation over directions and recompute Ac
!=======================================================================
          oEW=1.0_r8/MAX(EW,EWlim)
          DO d=1,ND
            IF ((h_tot(i,j).ge.Dcrit(ng)).and.(EW.gt.EWlim)) THEN
#   ifdef WDISS_FRIC
              EN(d)=MAX(0.0_r8,EN(d)-(Dissip_break(i,j)+dissfric)*      &
     &                                EN(d)*oEW)
#   else
              EN(d)=MAX(0.0_r8,EN(d)-Dissip_break(i,j)*EN(d)*oEW)
#   endif
              AC(i,j,d,nout)=EN(d)*Tr(i,j,d)*otwopi
            ELSE
              AC(i,j,d,nout)=0.0_r8
            ENDIF
#  ifdef MASKING
            AC(i,j,d,nout)=AC(i,j,d,nout)*rmask(i,j)
#  endif
#  ifdef WET_DRY
            AC(i,j,d,nout)=AC(i,j,d,nout)*rmask_wet(i,j)
#  endif
          ENDDO
          Dissip_wcap(i,j)=0.0_r8
          Dissip_break(i,j)=Dissip_break(i,j)*orho0*fac
        ENDDO
      ENDDO
!
!  Apply boundary conditions.
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Dissip_break)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Dissip_wcap)
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Dissip_break, Dissip_wcap)
#  endif
      RETURN
      END SUBROUTINE dissip_inw_tile
# endif
#endif
      END MODULE dissip_inw_mod
