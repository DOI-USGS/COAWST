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
      PUBLIC  :: dissip_inw

      CONTAINS
!
!***********************************************************************
      SUBROUTINE dissip_inw (ng, tile)
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
      integer, intent(in) :: ng, tile
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
     &                 nstp(ng), nnew(ng),                              &
#  ifdef DISSIP_BREAK_DIR
     &                 GRID(ng) % angler,                               &
#  endif
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
     &                 GRID(ng) % umask,                                &
     &                 GRID(ng) % vmask,                                &
#  endif
#  ifdef WET_DRY
     &                 GRID(ng) % rmask_wet,                            &
     &                 GRID(ng) % umask_wet,                            &
     &                 GRID(ng) % vmask_wet,                            &
#  endif
#  ifdef DIAGNOSTICS_INWAVE
     &                 DIAGS(ng) % DiaInwrk,                            &
#  endif
# ifdef DISSIP_BREAK_DIR
     &                 FORCES(ng) % Dissip_breakx,                      &
     &                 FORCES(ng) % Dissip_breaky,                      &
     &                 FORCES(ng) % Dissip_wcapx,                       &
     &                 FORCES(ng) % Dissip_wcapy,                       &
     &                 FORCES(ng) % Dwave,                              &
# else
     &                 FORCES(ng) % Dissip_break,                       &
     &                 FORCES(ng) % Dissip_wcap,                        &
# endif
     &                 FORCES(ng) % Lwave,                              &
#  ifdef WAVES_BOT_PERIOD
     &                 FORCES(ng) % Pwave_bot,                          &
#  endif
     &                 FORCES(ng) % Wave_break,                         &
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
#  ifdef DISSIP_BREAK_DIR
     &                       angler,                                    &
#  endif
#  ifdef MASKING
     &                       rmask, umask, vmask,                       &
#  endif
#  ifdef WET_DRY
     &                       rmask_wet, umask_wet, vmask_wet,           &
#  endif
#  ifdef DIAGNOSTICS_INWAVE
     &                       DiaInwrk,                                  &
#  endif
# ifdef DISSIP_BREAK_DIR
     &                       Dissip_breakx, Dissip_breaky,              &
     &                       Dissip_wcapx, Dissip_wcapy,                &
     &                       Dwave,                                     &
# else
     &                       Dissip_break, Dissip_wcap,                 &
# endif
     &                       Lwave,                                     &
#  ifdef WAVES_BOT_PERIOD
     &                       Pwave_bot,                                 &
#  endif
     &                       Wave_break,                                &
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
      USE exchange_3d_mod
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
      USE mp_exchange_mod, ONLY : mp_exchange3d
#  endif
      USE bc_2d_mod
      USE AC3dbc_mod, ONLY : AC3dbc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nout
!
#  ifdef ASSUMED_SHAPE
#   ifdef DISSIP_BREAK_DIR
      real(r8), intent(in) :: angler(LBi:,LBj:)
#   endif
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
#   ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
#   endif
#   ifdef DIAGNOSTICS_INWAVE
      real(r8), intent(inout) :: DiaInwrk(LBi:,LBj:,:,:)
#   endif
#   ifdef DISSIP_BREAK_DIR
      real(r8), intent(inout) :: Dissip_breakx(LBi:,LBj:)
      real(r8), intent(inout) :: Dissip_breaky(LBi:,LBj:)
      real(r8), intent(inout) :: Dissip_wcapx(LBi:,LBj:)
      real(r8), intent(inout) :: Dissip_wcapy(LBi:,LBj:)
      real(r8), intent(in)    :: Dwave(LBi:,LBj:)
#   else
      real(r8), intent(inout) :: Dissip_break(LBi:,LBj:)
      real(r8), intent(inout) :: Dissip_wcap(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: Lwave(LBi:,LBj:)
#   ifdef WAVES_BOT_PERIOD
      real(r8), intent(in) :: Pwave_bot(LBi:,LBj:)
#   endif
      real(r8), intent(inout) :: Wave_break(LBi:,LBj:)
      real(r8), intent(in) :: h_tot(LBi:,LBj:)
      real(r8), intent(inout) :: AC(LBi:,LBj:,:,:)
      real(r8), intent(in) :: Tr(LBi:,LBj:,:)
      real(r8), intent(in) :: kwc(LBi:,LBj:,:)
#  else
#   ifdef DISSIP_BREAK_DIR
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
#   endif
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
#   endif
#   ifdef DIAGNOSTICS_INWAVE
      real(r8), intent(inout) :: DiaInwrk(LBi:UBi,LBj:UBj,1:ND,NDInw)
#   endif
#   ifdef DISSIP_BREAK_DIR
      real(r8), intent(inout) :: Dissip_breakx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Dissip_breaky(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Dissip_wcapx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Dissip_wcapy(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: Dwave(LBi:UBi,LBj:UBj)
#   else
      real(r8), intent(inout) :: Dissip_break(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Dissip_wcap(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: Lwave(LBi:UBi,LBj:UBj)
#   ifdef WAVES_BOT_PERIOD
      real(r8), intent(in) :: Pwave_bot(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(inout) :: Wave_break(LBi:UBi,LBj:UBj)
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
      real(r8) :: ETOT, B, B2, QO, Z
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
      real(r8) :: disstot
      real(r8) :: EN(ND)

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
              disstot=2.0_r8*alfa/TRM*EW*Qb*dt(ng)*                     &
     &                          WAVEG(ng)%dtifac
            ELSE
              disstot=0.0_r8
            END IF
#  elif defined WDISS_GAMMA
            Hmax_r=breakg*(MAX(h_tot(i,j),0.0_r8))
            Emax_r=0.125_r8*g*rho0*Hmax_r**2.0_r8
            diff=EW-Emax_r
            disstot=MAX(0.0_r8,diff)
#  else
            disstot=0.0_r8
!  set an Hmax_r for the break percent.
            Hmax_r=breakr*(MAX(h_tot(i,j),0.0_r8))
#  endif
!
! Compute percent wave breaking.
! Copied from swancom2.ftn
!
            ETOT=EW*ogrho0
            IF ( (Hmax_r.gt.0.0_r8).and.(ETOT.ge.0.0_r8) ) THEN
              B=SQRT(8.0_r8 * ETOT / (Hmax_r*Hmax_r) )
            ELSE
              B=0.0_r8
            END IF
            IF ( B.le.0.5_r8 ) THEN
              QO=0.0_r8
            ELSE IF (B.le.1.0_r8) THEN
              QO=(2.0_r8*B-1.0_r8)**2
            END IF
!
            IF (B.le.0.2_r8) THEN
              Wave_break(i,j)=0.0_r8
            ELSE IF (B.lt.1.0_r8) THEN
!
!       *** second iteration to find Qb ***
!
              B2 = B*B
              Z  = EXP((QO-1.0_r8)/B2)
              Wave_break(i,j)=QO - B2 * (QO-Z)/(B2-Z)
            ELSE
              Wave_break(i,j)=1.0_r8
            END IF
          ELSE
            disstot=0.0_r8
#  ifdef WDISS_FRIC
            dissfric=0.0_r8
#  endif
            Wave_break(i,j)=0.0_r8
          END IF
#  ifdef MASKING
          disstot=disstot*rmask(i,j)
#   ifdef WDISS_FRIC
          dissfric=dissfric*rmask(i,j)
#   endif
          Wave_break(i,j)=Wave_break(i,j)*rmask(i,j)
#  endif
#  ifdef WET_DRY
          disstot=disstot*rmask_wet(i,j)
#   ifdef WDISS_FRIC
          dissfric=dissfric*rmask_wet(i,j)
#   endif
          Wave_break(i,j)=Wave_break(i,j)*rmask_wet(i,j)
#  endif
#  ifdef DISSIP_BREAK_DIR
          cff=1.5_r8*pi-Dwave(i,j)-angler(i,j)
          Dissip_breakx(i,j)=disstot*orho0*fac*COS(cff)
          Dissip_breaky(i,j)=disstot*orho0*fac*SIN(cff)
#  else
          Dissip_break(i,j)=disstot*orho0*fac
#  endif
!=======================================================================
!  Distribute dissipation over directions and recompute AC
!=======================================================================
          oEW=1.0_r8/MAX(EW,EWlim)
          DO d=1,ND
!           IF ((h_tot(i,j).ge.Dcrit(ng)).and.(EW.gt.EWlim)) THEN
            IF ((h_tot(i,j).ge.1.5_r8*Dcrit(ng)).and.(EW.gt.EWlim)) THEN
              cff1=EN(d)*oEW
              cff=disstot*cff1
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
              AC(i,j,d,nout)=0.0_r8
            ENDIF
#  ifdef MASKING
            AC(i,j,d,nout)=AC(i,j,d,nout)*rmask(i,j)
#  endif
#  ifdef WET_DRY
            AC(i,j,d,nout)=AC(i,j,d,nout)*rmask_wet(i,j)
#  endif
#  if defined DIAGNOSTICS_INWAVE
            DiaInwrk(i,j,d,iIrate)=DiaInwrk(i,j,d,iIrate)-              &
     &                             AC(i,j,d,nout)
#  endif
          ENDDO
#  ifdef DISSIP_BREAK_DIR
          Dissip_wcapx(i,j)=0.0_r8
          Dissip_wcapy(i,j)=0.0_r8
#  else
          Dissip_wcap(i,j)=0.0_r8
#  endif
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions.
!-----------------------------------------------------------------------
!
      CALL AC3dbc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  nstp, nout,                                     &
     &                  AC)
#  ifdef DISSIP_BREAK_DIR
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Dissip_breakx)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Dissip_breaky)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Dissip_wcapx)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Dissip_wcapy)
#  else
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Dissip_break)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Dissip_wcap)
#  endif
!  Apply periodic boundary conditions.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_AC3d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, 1, ND,             &
     &                           AC(:,:,:,nout))
#  ifdef DISSIP_BREAK_DIR
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Dissip_breakx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Dissip_breaky)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Dissip_wcapx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Dissip_wcapx)
#  else
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Dissip_break)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Dissip_wcap)
#  endif
      END IF
#  ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, ND,                    &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    AC(:,:,:,nout))
#   ifdef DISSIP_BREAK_DIR
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Dissip_breakx, Dissip_breaky)
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Dissip_wcapx, Dissip_wcapy)
#   else
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Dissip_break)
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Dissip_wcap)
#   endif
#  endif
      RETURN
      END SUBROUTINE dissip_inw_tile
# endif
#endif
      END MODULE dissip_inw_mod
