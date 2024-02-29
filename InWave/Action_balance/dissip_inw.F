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
#  ifdef DIAGNOSTICS_INWAVE
      USE mod_diags
#  endif
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
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
#  ifdef WET_DRY
     &                 GRID(ng) % rmask_wet,                            &
     &                 GRID(ng) % umask_wet,                            &
     &                 GRID(ng) % vmask_wet,                            &
#  endif
#  ifdef DIAGNOSTICS_INWAVE
     &                 DIAGS(ng) % DiaInwrk,                            &
#  endif
#  if defined WAVE_PARTITION
     &                 FORCES(ng) % Dissip_break_dir,                   &
#  endif
     &                 FORCES(ng) % Dissip_break,                       &
     &                 FORCES(ng) % Dissip_wcap,                        &
     &                 FORCES(ng) % Lwave,                              &
#  ifdef WAVES_BOT_PERIOD
     &                 FORCES(ng) % Pwave_bot,                          &
#  endif
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
#  ifdef MASKING
     &                       rmask,                                     &
#  endif
#  ifdef WET_DRY
     &                       rmask_wet, umask_wet, vmask_wet,           &
#  endif
#  ifdef DIAGNOSTICS_INWAVE
     &                       DiaInwrk,                                  &
#  endif
#  if defined WAVE_PARTITION
     &                       Dissip_break_dir,                          &
#  endif
     &                       Dissip_break,                              &
     &                       Dissip_wcap,                               &
     &                       Lwave,                                     &
#  ifdef WAVES_BOT_PERIOD
     &                       Pwave_bot,                                 &
#  endif
     &                       h_tot,                                     &
     &                       AC, Tr, kwc)
!***********************************************************************
!
      USE mod_param
      USE mod_inwave_params
      USE mod_inwave_vars
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
      USE exchange_2d_mod
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  endif
      USE bc_2d_mod
#  if defined WAVE_PARTITION
      USE exchange_3d_mod
#   ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
#   endif
      USE bc_3d_mod
#  endif
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
#   ifdef DIAGNOSTICS_INWAVE
      real(r8), intent(inout) :: DiaInwrk(LBi:,LBj:,:,:)
#   endif
#   if defined WAVE_PARTITION
      real(r8), intent(inout) :: Dissip_break_dir(LBi:,LBj:,:)
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
#   ifdef DIAGNOSTICS_INWAVE
      real(r8), intent(inout) :: DiaInwrk(LBi:UBi,LBj:UBj,1:ND,NDInw)
#   endif
#   if defined WAVE_PARTITION
      real(r8), intent(inout) :: Dissip_break_dir(LBi:UBi,LBj:UBj,ND)
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
      real(r8) :: twopi, otwopi, ogrho0, cff, cff1, orho0, fac
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
      real(r8) :: EN(ND)
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: disstot

# include "set_bounds.h"

      twopi=2.0_r8*pi
      otwopi=1.0_r8/twopi
      ogrho0=1.0_r8/(g*rho0)
      orho0=1.0_r8/rho0
      fac=1.0_r8/(dt(ng)*WAVEG(ng)%dtifac)
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
     &         MAX(Trmin,Pwave_bot(i,j)*Pwave_bot(i,j))*WAVEG(ng)%dtifac
#  endif
#  ifdef WDISS_ROELVINK
            Hmax_r=breakr*(MAX(h_tot(i,j),0.0_r8))
            Qb=MIN(1.0_r8,1.0_r8-EXP(-(H/Hmax_r)**n_r))
            IF (TRM.gt.0.0001_r8) THEN
              disstot(i,j)=2.0_r8*alfa/TRM*EW*Qb*dt(ng)*                &
     &                          WAVEG(ng)%dtifac
            ELSE
              disstot(i,j)=0.0_r8
            END IF
#  elif defined WDISS_GAMMA
            Hmax_r=breakg*(MAX(h_tot(i,j),0.0_r8))
            Emax_r=0.125_r8*g*rho0*Hmax_r**2.0_r8
            diff=EW-Emax_r
            disstot(i,j)=MAX(0.0_r8,diff)
#  else
            disstot(i,j)=0.0_r8
#  endif
          ELSE
            disstot(i,j)=0.0_r8
#  ifdef WDISS_FRIC
            dissfric=0.0_r8
#  endif
          END IF
#  ifdef MASKING
          disstot(i,j)=disstot(i,j)*rmask(i,j)
#   ifdef WDISS_FRIC
          dissfric=dissfric*rmask(i,j)
#   endif
#  endif
#  ifdef WET_DRY
          disstot(i,j)=disstot(i,j)*rmask_wet(i,j)
#   ifdef WDISS_FRIC
          dissfric=dissfric*rmask_wet(i,j)
#   endif
#  endif
          Dissip_break(i,j)=disstot(i,j)*orho0*fac
!=======================================================================
!  Distribute dissipation over directions and recompute AC
!=======================================================================
          oEW=1.0_r8/MAX(EW,EWlim)
          DO d=1,ND
!           IF ((h_tot(i,j).ge.Dcrit(ng)).and.(EW.gt.EWlim)) THEN
            IF ((h_tot(i,j).ge.1.5_r8*Dcrit(ng)).and.(EW.gt.EWlim)) THEN
              cff1=EN(d)*oEW
              cff=disstot(i,j)*cff1
#  if defined WAVE_PARTITION
              Dissip_break_dir(i,j,d)=MAX(0.0_r8,cff)
#  endf
              EN(d)=MAX(0.0_r8,EN(d)-cff)
#  ifdef DIAGNOSTICS_INWAVE
              DiaInwrk(i,j,d,iIbrak)=cff*Tr(i,j,d)*otwopi
#  endif
#  ifdef WDISS_FRIC
              cff=dissfric*cff1
              EN(d)=MAX(0.0_r8,EN(d)-cff)
#   ifdef DIAGNOSTICS_INWAVE
              DiaInwrk(i,j,d,iIfric)=cff*Tr(i,j,d)*otwopi
#   endif
#  endif
              AC(i,j,d,nout)=EN(d)*Tr(i,j,d)*otwopi
            ELSE
#  if defined WAVE_PARTITION
              Dissip_break_dir(i,j,d)=0.0_r8
#  endif
              AC(i,j,d,nout)=0.0_r8
            ENDIF
#  ifdef MASKING
            AC(i,j,d,nout)=AC(i,j,d,nout)*rmask(i,j)
#   if defined WAVE_PARTITION
            Dissip_break_dir(i,j,d)=Dissip_break_dir(i,j,d)*rmask(i,j)
#   endif
#  endif
#  ifdef WET_DRY
            AC(i,j,d,nout)=AC(i,j,d,nout)*rmask_wet(i,j)
#   if defined WAVE_PARTITION
            Dissip_break_dir(i,j,d)=Dissip_break_dir(i,j,d)*            &
     &                              rmask_wet(i,j)
#   endif
#  endif
#  if defined DIAGNOSTICS_INWAVE
            DiaInwrk(i,j,d,iIrate)=DiaInwrk(i,j,d,iIrate)-              &
     &                             AC(i,j,d,nout)
#  endif
#  if defined WAVE_PARTITION
            Dissip_break_dir(i,j,d)=Dissip_break_dir(i,j,d)*orho0*fac
#  endif
          ENDDO
          Dissip_wcap(i,j)=0.0_r8
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
#  if defined WAVE_PARTITION
      CALL bc_r3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, ND,                      &
     &                  Dissip_break_dir)
#  endif
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Dissip_wcap)
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Dissip_break)
#   if defined WAVE_PARTITION
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, ND,                    &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Dissip_break_dir)
#   endif
#  endif
      RETURN
      END SUBROUTINE dissip_inw_tile
# endif
#endif
      END MODULE dissip_inw_mod
