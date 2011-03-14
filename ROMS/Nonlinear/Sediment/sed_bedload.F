#include "cppdefs.h"
#define SLOPE_NEMETH
#undef  SLOPE_LESSER
#define BSTRESS_UPWIND

      MODULE sed_bedload_mod

#if defined NONLINEAR && defined SEDIMENT && defined BEDLOAD
!
!svn $Id: sed_bedload.F 429 2009-12-20 17:30:26Z arango $
!==================================================== John C. Warner ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group      Hernan G. Arango   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes sediment bedload transport using the Meyer-   !
!  Peter and Muller (1948) formulation  for unidirectional flow and    !
!  Souksby and Damgaard (2005) algorithm that accounts for combined    !
!  effect of currents and waves.                                       !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Meyer-Peter, E. and R. Muller, 1948: Formulas for bedload transport !
!    In: Report on the 2nd Meeting International Association Hydraulic !
!    Research, Stockholm, Sweden, pp 39-64.                            !
!                                                                      !
!  Soulsby, R.L. and J.S. Damgaard, 2005: Bedload sediment transport   !
!    in coastal waters, Coastal Engineering, 52 (8), 673-689.          !
!                                                                      !
!  Warner, J.C., C.R. Sherwood, R.P. Signell, C.K. Harris, and H.G.    !
!    Arango, 2008:  Development of a three-dimensional,  regional,     !
!    coupled wave, current, and sediment-transport model, Computers    !
!    & Geosciences, 34, 1284-1306.                                     !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: sed_bedload

      CONTAINS
!
!***********************************************************************
      SUBROUTINE sed_bedload (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_sedbed
      USE mod_stepping
# ifdef BBL_MODEL
      USE mod_bbl
# endif
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
      CALL wclock_on (ng, iNLM, 16)
# endif
      CALL sed_bedload_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nstp(ng), nnew(ng),                        &
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pn,                             &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
# endif
# ifdef WET_DRY
     &                       GRID(ng) % rmask_wet,                      &
# endif
     &                       GRID(ng) % z_w,                            &
# ifdef BBL_MODEL
     &                       BBL(ng) % bustrc,                          &
     &                       BBL(ng) % bvstrc,                          &
     &                       BBL(ng) % bustrw,                          &
     &                       BBL(ng) % bvstrw,                          &
     &                       BBL(ng) % bustrcwmax,                      &
     &                       BBL(ng) % bvstrcwmax,                      &
     &                       FORCES(ng) % Dwave,                        &
     &                       FORCES(ng) % Pwave_bot,                    &
# endif
     &                       FORCES(ng) % bustr,                        &
     &                       FORCES(ng) % bvstr,                        &
# if defined BEDLOAD_SOULSBY
     &                   FORCES(ng) % Hwave,                            &
     &                   FORCES(ng) % Lwave,                            &
     &                   GRID(ng) % angler,                             &
# endif
# if defined SED_MORPH
     &                   SEDBED(ng) % bed_thick,                        &
# endif
# if defined BEDLOAD_MPM || defined BEDLOAD_SOULSBY
     &                   GRID(ng) % h,                                  &
     &                   GRID(ng) % om_r,                               &
     &                   GRID(ng) % om_u,                               &
     &                   GRID(ng) % on_r,                               &
     &                   GRID(ng) % on_v,                               &
     &                   SEDBED(ng) % bedldu,                           &
     &                   SEDBED(ng) % bedldv,                           &
# endif
     &                   SEDBED(ng) % bed,                              &
     &                   SEDBED(ng) % bed_frac,                         &
     &                   SEDBED(ng) % bed_mass,                         &
     &                   SEDBED(ng) % bottom)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 16)
# endif
      RETURN
      END SUBROUTINE sed_bedload
!
!***********************************************************************
      SUBROUTINE sed_bedload_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nstp, nnew,                          &
     &                             pm, pn,                              &
# ifdef MASKING
     &                             rmask, umask, vmask,                 &
# endif
# ifdef WET_DRY
     &                             rmask_wet,                           &
# endif
     &                             z_w,                                 &
# ifdef BBL_MODEL
     &                             bustrc, bvstrc,                      &
     &                             bustrw, bvstrw,                      &
     &                             bustrcwmax, bvstrcwmax,              &
     &                             Dwave, Pwave_bot,                    &
# endif
     &                             bustr, bvstr,                        &
# if defined BEDLOAD_SOULSBY
     &                             Hwave, Lwave,                        &
     &                             angler,                              &
# endif
# if defined SED_MORPH
     &                             bed_thick,                           &
# endif
# if defined BEDLOAD_MPM || defined BEDLOAD_SOULSBY
     &                             h, om_r, om_u, on_r, on_v,           &
     &                             bedldu, bedldv,                      &
# endif
     &                             bed, bed_frac, bed_mass,             &
     &                             bottom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sediment
!
      USE bc_3d_mod, ONLY : bc_r3d_tile
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# endif
# ifdef BEDLOAD
#  if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod, ONLY : exchange_u2d_tile, exchange_v2d_tile
#  endif
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
#  ifdef BBL_MODEL
      real(r8), intent(in) :: bustrc(LBi:,LBj:)
      real(r8), intent(in) :: bvstrc(LBi:,LBj:)
      real(r8), intent(in) :: bustrw(LBi:,LBj:)
      real(r8), intent(in) :: bvstrw(LBi:,LBj:)
      real(r8), intent(in) :: bustrcwmax(LBi:,LBj:)
      real(r8), intent(in) :: bvstrcwmax(LBi:,LBj:)
      real(r8), intent(in) :: Dwave(LBi:,LBj:)
      real(r8), intent(in) :: Pwave_bot(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
#  if defined BEDLOAD_SOULSBY
      real(r8), intent(in) :: Hwave(LBi:,LBj:)
      real(r8), intent(in) :: Lwave(LBi:,LBj:)
      real(r8), intent(in) :: angler(LBi:,LBj:)
#  endif
#  if defined SED_MORPH
      real(r8), intent(inout):: bed_thick(LBi:,LBj:,:)
#  endif
#  if defined BEDLOAD_MPM || defined BEDLOAD_SOULSBY
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: om_r(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: on_r(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(inout) :: bedldu(LBi:,LBj:,:)
      real(r8), intent(inout) :: bedldv(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: bed(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_frac(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_mass(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: bottom(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#  ifdef BBL_MODEL
      real(r8), intent(in) :: bustrc(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstrc(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bustrw(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstrw(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bustrcwmax(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstrcwmax(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Dwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Pwave_bot(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstr(LBi:UBi,LBj:UBj)
#  if defined BEDLOAD_SOULSBY
      real(r8), intent(in) :: Hwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Lwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
#  endif
#  if defined SED_MORPH
      real(r8), intent(inout):: bed_thick(LBi:UBi,LBj:UBj,2)
#  endif
#  if defined BEDLOAD_MPM || defined BEDLOAD_SOULSBY
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: bedldu(LBi:UBi,LBj:UBj,NST)
      real(r8), intent(inout) :: bedldv(LBi:UBi,LBj:UBj,NST)
#  endif
      real(r8), intent(inout) :: bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(inout) :: bed_frac(LBi:UBi,LBj:UBj,Nbed,NST)
      real(r8), intent(inout) :: bed_mass(LBi:UBi,LBj:UBj,Nbed,1:2,NST)
      real(r8), intent(inout) :: bottom(LBi:UBi,LBj:UBj,MBOTP)
# endif
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      integer :: i, ised, j, k

      real(r8), parameter :: eps = 1.0E-14_r8

      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5

# if defined BEDLOAD_MPM
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tau_w
# endif
# ifdef BSTRESS_UPWIND
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tau_wX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tau_wE
# endif
# ifdef BEDLOAD
      real(r8) :: a_slopex, a_slopey, sed_angle
      real(r8) :: bedld, bedld_mass, dzdx, dzdy
      real(r8) :: smgd, smgdr, osmgd, Umag
      real(r8) :: rhs_bed, Ua, Ra, phi, Clim

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX_r
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE_r
# endif
# if defined BEDLOAD_MPM
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: angleu
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: anglev
# endif
# if defined BEDLOAD_SOULSBY
      real(r8) :: theta_mean, theta_wav, w_asym
      real(r8) :: theta_max, theta_max1, theta_max2
      real(r8) :: phi_x1, phi_x2, phi_x, phi_y, Dstp
      real(r8) :: bedld_x, bedld_y, tau_cur, waven, wavec

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: phic
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: phicw
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tau_wav
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tau_mean
      real(r8), parameter :: kdmax = 100.0_r8
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
! Compute maximum bottom stress for MPM bedload or suspended load.
!-----------------------------------------------------------------------
!
# if defined BEDLOAD_MPM
#  ifdef BBL_MODEL
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          tau_w(i,j)=SQRT(bustrcwmax(i,j)*bustrcwmax(i,j)+              &
     &                    bvstrcwmax(i,j)*bvstrcwmax(i,j))
#   ifdef WET_DRY
          tau_w(i,j)=tau_w(i,j)*rmask_wet(i,j)
#   endif
        END DO
      END DO
#  else
#   ifdef EW_PERIODIC
#    define I_RANGE Istr-1,Iend+1
#   else
#    define I_RANGE MAX(Istr-1,1),MIN(Iend+1,Lm(ng))
#   endif
#   ifdef NS_PERIODIC
#    define J_RANGE Jstr-1,Jend+1
#   else
#    define J_RANGE MAX(Jstr-1,1),MIN(Jend+1,Mm(ng))
#   endif
      DO i=I_RANGE
        DO j=J_RANGE
#   ifdef BSTRESS_UPWIND
          cff1=0.5_r8*(1.0_r8+SIGN(1.0_r8,bustr(i+1,j)))
          cff2=0.5_r8*(1.0_r8-SIGN(1.0_r8,bustr(i+1,j)))
          cff3=0.5_r8*(1.0_r8+SIGN(1.0_r8,bustr(i  ,j)))
          cff4=0.5_r8*(1.0_r8-SIGN(1.0_r8,bustr(i  ,j)))
          tau_wX(i,j)=cff3*(cff1*bustr(i,j)+                            &
     &                cff2*0.5_r8*(bustr(i,j)+bustr(i+1,j)))+           &
     &                cff4*(cff2*bustr(i+1,j)+                          &
     &                cff1*0.5_r8*(bustr(i,j)+bustr(i+1,j)))
          cff1=0.5_r8*(1.0_r8+SIGN(1.0_r8,bvstr(i,j+1)))
          cff2=0.5_r8*(1.0_r8-SIGN(1.0_r8,bvstr(i,j+1)))
          cff3=0.5_r8*(1.0_r8+SIGN(1.0_r8,bvstr(i,j)))
          cff4=0.5_r8*(1.0_r8-SIGN(1.0_r8,bvstr(i,j)))
          tau_wE(i,j)=cff3*(cff1*bvstr(i,j)+                            &
     &                cff2*0.5_r8*(bvstr(i,j)+bvstr(i,j+1)))+           &
     &                cff4*(cff2*bvstr(i,j+1)+                          &
     &                cff1*0.5_r8*(bvstr(i,j)+bvstr(i,j+1)))
#   endif
          tau_w(i,j)=0.5_r8*SQRT((bustr(i,j)+bustr(i+1,j))*             &
     &                           (bustr(i,j)+bustr(i+1,j))+             &
     &                           (bvstr(i,j)+bvstr(i,j+1))*             &
     &                           (bvstr(i,j)+bvstr(i,j+1)))
#   ifdef WET_DRY
          tau_w(i,j)=tau_w(i,j)*rmask_wet(i,j)
#   endif
        END DO
      END DO
#   undef I_RANGE
#   undef J_RANGE
#  endif
# endif

# ifdef BEDLOAD
!
!-----------------------------------------------------------------------
!  Compute bedload sediment transport.
!-----------------------------------------------------------------------
!
#  ifdef EW_PERIODIC
#   define I_RANGE Istr-1,Iend+1
#  else
#   define I_RANGE MAX(Istr-1,1),MIN(Iend+1,Lm(ng))
#  endif
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr-1,Jend+1
#  else
#   define J_RANGE MAX(Jstr-1,1),MIN(Jend+1,Mm(ng))
#  endif
!
! Compute some constant bed slope parameters.
!
      sed_angle=DTAN(33.0_r8*pi/180.0_r8)
!
!  Compute angle between currents and waves (radians).
!
      DO j=J_RANGE
        DO i=I_RANGE
#  if defined BEDLOAD_SOULSBY
!
! Compute angle between currents and waves, measure CCW from current
! direction toward wave vector.
!
          IF (bustrc(i,j).eq.0.0_r8) THEN
            phic(i,j)=0.5_r8*pi*SIGN(1.0_r8,bvstrc(i,j))
          ELSE
            phic(i,j)=ATAN2(bvstrc(i,j),bustrc(i,j))
          ENDIF
          phicw(i,j)=1.5_r8*pi-Dwave(i,j)-phic(i,j)-angler(i,j)
!
! Compute stress components at rho points.
!
          tau_cur=SQRT(bustrc(i,j)*bustrc(i,j)+                         &
     &                 bvstrc(i,j)*bvstrc(i,j))
          tau_wav(i,j)=SQRT(bustrw(i,j)*bustrw(i,j)+                    &
     &                      bvstrw(i,j)*bvstrw(i,j))
          tau_mean(i,j)=tau_cur*(1.0_r8+1.2_r8*((tau_wav(i,j)/          &
     &                  (tau_cur+tau_wav(i,j)+eps))**3.2_r8))
!
#  elif defined BEDLOAD_MPM
          cff1=0.5_r8*(bustr(i,j)+bustr(i+1,j))
          cff2=0.5_r8*(bvstr(i,j)+bvstr(i,j+1))
          Umag=SQRT(cff1*cff1+cff2*cff2)+eps
          angleu(i,j)=cff1/Umag
          anglev(i,j)=cff2/Umag
#  endif
        END DO
      END DO
!
      DO ised=NCS+1,NST
        smgd=(Srho(ised,ng)/rho0-1.0_r8)*g*Sd50(ised,ng)
        osmgd=1.0_r8/smgd
        smgdr=SQRT(smgd)*Sd50(ised,ng)*Srho(ised,ng)
!
        DO j=J_RANGE
          DO i=I_RANGE
#  ifdef BEDLOAD_SOULSBY
!
! Compute wave asymmetry factor, based on Fredosoe and Deigaard.
!
            Dstp=z_w(i,j,N(ng))+h(i,j)
            waven=2.0_r8*pi/(Lwave(i,j)+eps)
            wavec=SQRT(g/waven*tanh(waven*Dstp))
            cff4=MIN(waven*Dstp,kdmax)
            cff1=-0.1875_r8*wavec*(waven*Dstp)**2/(SINH(cff4))**4
            cff2=0.125_r8*g*Hwave(i,j)**2/(wavec*Dstp+eps)
            cff3=pi*Hwave(i,j)/(Pwave_bot(i,j)*SINH(cff4)+eps)
            w_asym=MAX(MIN((cff1-cff2)/cff3,0.2_r8),0.0_r8)
            w_asym=0.0_r8
!
! Compute nondimensional stresses.
!
            theta_wav=tau_wav(i,j)*osmgd+eps
            theta_mean=tau_mean(i,j)*osmgd
!
            cff1=theta_wav*(1.0_r8+w_asym)
            cff2=theta_wav*(1.0_r8-w_asym)
            theta_max1=SQRT((theta_mean+                                &
     &                       cff1*COS(phicw(i,j)))**2+                  &
     &                      (cff1*SIN(phicw(i,j)))**2)
            theta_max2=SQRT((theta_mean+                                &
     &                       cff2*COS(phicw(i,j)+pi))**2+               &
     &                      (cff2*SIN(phicw(i,j)+pi))**2)
            theta_max=MAX(theta_max1,theta_max2)
!
! Motion initiation factor.
!
            cff3=0.5_r8*(1.0_r8+SIGN(1.0_r8,                            &
     &                               theta_max/tau_ce(ised,ng)-1.0_r8))
!
! Calculate bed loads in direction of current and perpendicular
! direction.
!
            phi_x1=12.0_r8*SQRT(theta_mean)*                            &
     &             MAX((theta_mean-tau_ce(ised,ng)),0.0_r8)
            phi_x2=12.0_r8*(0.9534_r8+0.1907*COS(2.0_r8*phicw(i,j)))*   &
     &             SQRT(theta_wav)*theta_mean+                          &
     &             12.0_r8*(0.229_r8*w_asym*theta_wav**1.5_r8*          &
     &                      COS(phicw(i,j)))
!           phi_x=MAX(phi_x1,phi_x2) !  <- original
            IF (ABS(phi_x2).gt.phi_x1) THEN
              phi_x=phi_x2
            ELSE
              phi_x=phi_x1
            END IF
            bedld_x=phi_x*smgdr*cff3
!
            cff5=theta_wav**1.5_r8+1.5_r8*(theta_mean**1.5_r8)
            phi_y=12.0_r8*0.1907_r8*theta_wav*theta_wav*                &
     &            (theta_mean*SIN(2.0_r8*phicw(i,j))+1.2_r8*w_asym*     &
     &            theta_wav*SIN(phicw(i,j)))/cff5*cff3
            bedld_y=phi_y*smgdr
!
! Partition bedld into xi and eta directions, still at rho points.
! (FX_r and FE_r have dimensions of kg).
!
            FX_r(i,j)=(bedld_x*COS(phic(i,j))-bedld_y*SIN(phic(i,j)))*  &
     &                on_r(i,j)*dt(ng)
            FE_r(i,j)=(bedld_x*SIN(phic(i,j))+bedld_y*COS(phic(i,j)))*  &
     &                om_r(i,j)*dt(ng)
#  elif defined BEDLOAD_MPM
#   ifdef BSTRESS_UPWIND
!
! Magnitude of bed load at rho points. Meyer-Peter Muller formulation.
! bedld has dimensions of kg m-1 s-1. Use partitions of stress
! from upwind direction, still at rho points.
! (FX_r and FE_r have dimensions of kg).
!
            bedld=8.0_r8*(MAX((ABS(tau_wX(i,j))*osmgd-0.047_r8),        &
     &                        0.0_r8)**1.5_r8)*smgdr*                   &
     &                        SIGN(1.0_r8,tau_wX(i,j))
            FX_r(i,j)=bedld*on_r(i,j)*dt(ng)
            bedld=8.0_r8*(MAX((ABS(tau_wE(i,j))*osmgd-0.047_r8),        &
     &                        0.0_r8)**1.5_r8)*smgdr*                   &
     &                        SIGN(1.0_r8,tau_wE(i,j))
            FE_r(i,j)=bedld*om_r(i,j)*dt(ng)
#   else
!
! Magnitude of bed load at rho points. Meyer-Peter Muller formulation.
! (BEDLD has dimensions of kg m-1 s-1).
!
            bedld=8.0_r8*(MAX((tau_w(i,j)*osmgd-0.047_r8),              &
     &                        0.0_r8)**1.5_r8)*smgdr
!
! Partition bedld into xi and eta directions, still at rho points.
! (FX_r and FE_r have dimensions of kg).
!
            FX_r(i,j)=angleu(i,j)*bedld*on_r(i,j)*dt(ng)
            FE_r(i,j)=anglev(i,j)*bedld*om_r(i,j)*dt(ng)
#   endif
#  endif
!
! Correct for along-direction slope. Limit slope to 0.9*sed angle.
!
            cff1=0.5_r8*(1.0_r8+SIGN(1.0_r8,FX_r(i,j)))
            cff2=0.5_r8*(1.0_r8-SIGN(1.0_r8,FX_r(i,j)))
#  if defined SLOPE_NEMETH
            dzdx=(h(i+1,j)-h(i  ,j))/om_u(i+1,j)*cff1+                  &
     &           (h(i-1,j)-h(i  ,j))/om_u(i  ,j)*cff2
            dzdy=(h(i,j+1)-h(i,j  ))/on_v(i,j+1)*cff1+                  &
     &           (h(i,j-1)-h(i,j  ))/on_v(i  ,j)*cff2
#   ifdef BEDLOAD_MPM
            cff=ABS(tau_w(i,j))
#   else
            cff=ABS(tau_mean(i,j))
#   endif
            a_slopex=0.3_r8*cff**0.5_r8*0.002_r8*dzdx+                  &
     &               0.3_r8*cff**1.5_r8*3.330_r8*dzdx
            a_slopey=0.3_r8*cff**0.5_r8*0.002_r8*dzdy+                  &
     &               0.3_r8*cff**1.5_r8*3.330_r8*dzdy
!
! Add contriubiton of bed slope to bed load transport fluxes.
!
            FX_r(i,j)=FX_r(i,j)+a_slopex
            FE_r(i,j)=FE_r(i,j)+a_slopey
#  elif defined SLOPE_LESSER
            dzdx=MIN(((h(i+1,j)-h(i  ,j))/om_u(i+1,j)*cff1+             &
     &                (h(i  ,j)-h(i-1,j))/om_u(i  ,j)*cff2),0.52_r8)*   &
     &                SIGN(1.0_r8,FX_r(i,j))
            dzdy=MIN(((h(i,j+1)-h(i,j  ))/on_v(i,j+1)*cff1+             &
     &                (h(i,j  )-h(i,j-1))/on_v(i  ,j)*cff2),0.52_r8)*   &
     &                SIGN(1.0_r8,FE_r(i,j))
            cff=DATAN(dzdx)
            a_slopex=sed_angle/(COS(cff)*(sed_angle-dzdx))
            cff=DATAN(dzdy)
            a_slopey=sed_angle/(COS(cff)*(sed_angle-dzdy))
!
! Add contriubiton of bed slope to bed load transport fluxes.
!
            FX_r(i,j)=FX_r(i,j)*a_slopex
            FE_r(i,j)=FE_r(i,j)*a_slopey
#  endif
!
!
#  ifdef SED_MORPH
!
! Apply morphology factor.
!
            FX_r(i,j)=FX_r(i,j)*morph_fac(ised,ng)
            FE_r(i,j)=FE_r(i,j)*morph_fac(ised,ng)

#  endif
!
! Apply bedload transport rate coefficient. Also limit
! bedload to the fraction of each sediment class.
!
            FX_r(i,j)=FX_r(i,j)*bedload_coeff(ng)*bed_frac(i,j,1,ised)
            FE_r(i,j)=FE_r(i,j)*bedload_coeff(ng)*bed_frac(i,j,1,ised)
!
! Limit bed load to available bed mass.
!
            bedld_mass=ABS(FX_r(i,j))+ABS(FE_r(i,j))+eps
            FX_r(i,j)=MIN(ABS(FX_r(i,j)),                               &
     &                    bed_mass(i,j,1,nstp,ised)*                    &
     &                    om_r(i,j)*on_r(i,j)*ABS(FX_r(i,j))/           &
     &                    bedld_mass)*                                  &
     &                    SIGN(1.0_r8,FX_r(i,j))
            FE_r(i,j)=MIN(ABS(FE_r(i,j)),                               &
     &                    bed_mass(i,j,1,nstp,ised)*                    &
     &                    om_r(i,j)*on_r(i,j)*ABS(FE_r(i,j))/           &
     &                    bedld_mass)*                                  &
     &                    SIGN(1.0_r8,FE_r(i,j))
          END DO
        END DO
#  ifndef EW_PERIODIC
        IF (WESTERN_EDGE) THEN
          DO j=J_RANGE
            FX_r(Istr-1,j)=FX_r(Istr,j)
            FE_r(Istr-1,j)=FE_r(Istr,j)
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO j=J_RANGE
            FX_r(Iend+1,j)=FX_r(Iend,j)
            FE_r(Iend+1,j)=FE_r(Iend,j)
          END DO
        END IF
#  endif
#  ifndef NS_PERIODIC
        IF (SOUTHERN_EDGE) THEN
          DO i=I_RANGE
            FX_r(i,Jstr-1)=FX_r(i,Jstr)
            FE_r(i,Jstr-1)=FE_r(i,Jstr)
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
          DO i=I_RANGE
            FX_r(i,Jend+1)=FX_r(i,Jend)
            FE_r(i,Jend+1)=FE_r(i,Jend)
          END DO
        END IF
#  endif
#  undef I_RANGE
#  undef J_RANGE
#  if !defined EW_PERIODIC && !defined NS_PERIODIC
        IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
          FX_r(Istr-1,Jstr-1)=0.5_r8*(FX_r(Istr  ,Jstr-1)+              &
     &                                FX_r(Istr-1,Jstr  ))
          FE_r(Istr-1,Jstr-1)=0.5_r8*(FE_r(Istr  ,Jstr-1)+              &
     &                                FE_r(Istr-1,Jstr  ))
        END IF
        IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
          FX_r(Iend+1,Jstr-1)=0.5_r8*(FX_r(Iend  ,Jstr-1)+              &
     &                                FX_r(Iend+1,Jstr  ))
          FE_r(Iend+1,Jstr-1)=0.5_r8*(FE_r(Iend  ,Jstr-1)+              &
     &                                FE_r(Iend+1,Jstr  ))
        END IF
        IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
          FX_r(Istr-1,Jend+1)=0.5_r8*(FX_r(Istr-1,Jend  )+              &
     &                                FX_r(Istr  ,Jend+1))
          FE_r(Istr-1,Jend+1)=0.5_r8*(FE_r(Istr-1,Jend  )+              &
     &                                FE_r(Istr  ,Jend+1))
        END IF
        IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
          FX_r(Iend+1,Jend+1)=0.5_r8*(FX_r(Iend+1,Jend  )+              &
     &                                FX_r(Iend  ,Jend+1))
          FE_r(Iend+1,Jend+1)=0.5_r8*(FE_r(Iend+1,Jend  )+              &
     &                                FE_r(Iend  ,Jend+1))
        END IF
#  endif
!
! Upwind shift FX_r and FE_r to u and v points.
!
        DO j=Jstr-1,Jend+1
          DO i=Istr,Iend+1
            cff1=0.5_r8*(1.0_r8+SIGN(1.0_r8,FX_r(i,j)))
            cff2=0.5_r8*(1.0_r8-SIGN(1.0_r8,FX_r(i,j)))
            FX(i,j)=0.5_r8*(1.0_r8+SIGN(1.0_r8,FX_r(i-1,j)))*           &
     &              (cff1*FX_r(i-1,j)+                                  &
     &               cff2*0.5_r8*(FX_r(i-1,j)+FX_r(i,j)))+              &
     &              0.5_r8*(1.0_r8-SIGN(1.0_r8,FX_r(i-1,j)))*           &
     &              (cff2*FX_r(i  ,j)+                                  &
     &               cff1*0.5_r8*(FX_r(i-1,j)+FX_r(i,j)))
#  ifdef MASKING
            FX(i,j)=FX(i,j)*umask(i,j)
#  endif
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr-1,Iend+1
            cff1=0.5_r8*(1.0_r8+SIGN(1.0_r8,FE_r(i,j)))
            cff2=0.5_r8*(1.0_r8-SIGN(1.0_r8,FE_r(i,j)))
            FE(i,j)=0.5_r8*(1.0_r8+SIGN(1.0_r8,FE_r(i,j-1)))*           &
     &              (cff1*FE_r(i,j-1)+                                  &
     &               cff2*0.5_r8*(FE_r(i,j-1)+FE_r(i,j)))+              &
     &              0.5_r8*(1.0_r8-SIGN(1.0_r8,FE_r(i,j-1)))*           &
     &              (cff2*FE_r(i  ,j)+                                  &
     &               cff1*0.5_r8*(FE_r(i,j-1)+FE_r(i,j)))
#  ifdef MASKING
            FE(i,j)=FE(i,j)*vmask(i,j)
#  endif
          END DO
        END DO
!
! Limit fluxes to prevent bottom from breaking thru water surface.
!
!        DO j=Jstr,Jend
!          DO i=Istr,Iend
!            cff1=1.0_r8/(Srho(ised,ng)*(1.0_r8-bed(i,j,1,iporo)))
!            rhs_bed=(FX(i+1,j)-FX(i,j)+                                 &
!     &               FE(i,j+1)-FE(i,j))*pm(i,j)*pn(i,j)
!            cff2=MAX(rhs_bed*cff1+h(i,j)-Dcrit(ng),0.0_r8)
!            cff=cff2/ABS(cff2+eps)
!            FX(i  ,j  )=MAX(FX(i  ,j  ),0.0_r8)*cff+                    &
!     &                  MIN(FX(i  ,j  ),0.0_r8)
!            FX(i+1,j  )=MAX(FX(i+1,j  ),0.0_r8)+                        &
!     &                  MIN(FX(i+1,j  ),0.0_r8)*cff
!            FE(i  ,j  )=MAX(FE(i  ,j  ),0.0_r8)*cff+                    &
!     &                  MIN(FE(i  ,j  ),0.0_r8)
!            FE(i  ,j+1)=MAX(FE(i  ,j+1),0.0_r8)+                        &
!     &                  MIN(FE(i  ,j+1),0.0_r8)*cff
!          END DO
!        END DO
!
#  ifndef EW_PERIODIC
#   ifdef WESTERN_WALL
        IF (WESTERN_EDGE) THEN
          DO j=Jstr-1,Jend+1
            FX(Istr,j)=0.0_r8
          END DO
        END IF
#   endif
#   ifdef EASTERN_WALL
        IF (EASTERN_EDGE) THEN
          DO j=Jstr-1,Jend+1
            FX(Iend+1,j)=0.0_r8
          END DO
        END IF
#   endif
#  endif
#  ifndef NS_PERIODIC
#   ifdef SOUTHERN_WALL
        IF (SOUTHERN_EDGE) THEN
          DO i=Istr-1,Iend+1
            FE(i,Jstr)=0.0_r8
          END DO
        END IF
#   endif
#   ifdef NORTHERN_WALL
        IF (NORTHERN_EDGE) THEN
          DO i=Istr-1,Iend+1
            FE(i,Jend+1)=0.0_r8
          END DO
        END IF
#   endif
#  endif
!
!  Determine flux divergence and evaluate change in bed properties.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff=(FX(i+1,j)-FX(i,j)+                                     &
     &           FE(i,j+1)-FE(i,j))*pm(i,j)*pn(i,j)
            bed_mass(i,j,1,nnew,ised)=MAX(bed_mass(i,j,1,nstp,ised)-    &
     &                                    cff,0.0_r8)
#  if !defined SUSPLOAD
            DO k=2,Nbed
              bed_mass(i,j,k,nnew,ised)=bed_mass(i,j,k,nstp,ised)
            END DO
#  endif
            bed(i,j,1,ithck)=MAX(bed(i,j,1,ithck)-                      &
     &                           cff/(Srho(ised,ng)*                    &
     &                                (1.0_r8-bed(i,j,1,iporo))),       &
     &                           0.0_r8)
#  ifdef MASKING
            bed(i,j,1,ithck)=bed(i,j,1,ithck)*rmask(i,j)
#  endif
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Output bedload fluxes.
!-----------------------------------------------------------------------
!
        cff=0.5_r8/dt(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            bedldu(i,j,ised)=FX(i,j)*(pn(i-1,j)+pn(i,j))*cff
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            bedldv(i,j,ised)=FE(i,j)*(pm(i,j-1)+pm(i,j))*cff
          END DO
        END DO
      END DO
!
!  Need to update bed mass for the non-cohesive sediment types, becasue 
!  they did not partake in the bedload transport.
!
      DO ised=1,NCS
        DO j=Jstr,Jend
          DO i=Istr,Iend
            bed_mass(i,j,1,nnew,ised)=bed_mass(i,j,1,nstp,ised)
#  if !defined SUSPLOAD
            DO k=2,Nbed
              bed_mass(i,j,k,nnew,ised)=bed_mass(i,j,k,nstp,ised)
            END DO
#  endif
          END DO
        END DO
      END DO
!
!  Update mean surface properties.
!  Sd50 must be positive definite, due to BBL routines.
!  Srho must be >1000, due to (s-1) in BBL routines.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          cff3=0.0_r8
          DO ised=1,NST
            cff3=cff3+bed_mass(i,j,1,nnew,ised)
          END DO
          IF (cff3.eq.0.0_r8) THEN
            cff3=eps
          END IF
          DO ised=1,NST
            bed_frac(i,j,1,ised)=bed_mass(i,j,1,nnew,ised)/cff3
          END DO
!
          cff1=1.0_r8
          cff2=1.0_r8
          cff3=1.0_r8
          cff4=1.0_r8
          DO ised=1,NST
            cff1=cff1*tau_ce(ised,ng)**bed_frac(i,j,1,ised)
            cff2=cff2*Sd50(ised,ng)**bed_frac(i,j,1,ised)
            cff3=cff3*(wsed(ised,ng)+eps)**bed_frac(i,j,1,ised)
            cff4=cff4*Srho(ised,ng)**bed_frac(i,j,1,ised)
          END DO
          bottom(i,j,itauc)=cff1
          bottom(i,j,isd50)=MIN(cff2,Zob(ng))
          bottom(i,j,iwsed)=cff3
          bottom(i,j,idens)=MAX(cff4,1050.0_r8)
        END DO
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Apply periodic or gradient boundary conditions to property arrays.
!-----------------------------------------------------------------------
!
      DO ised=1,NST
        CALL bc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    bed_frac(:,:,:,ised))
        CALL bc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    bed_mass(:,:,:,nnew,ised))
# if defined BEDLOAD && (defined EW_PERIODIC || defined NS_PERIODIC)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          bedldu(:,:,ised))
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          bedldv(:,:,ised))
# endif
      END DO
# ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, Nbed, 1, NST,          &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    bed_frac,                                     &
     &                    bed_mass(:,:,:,nnew,:))
#  if defined BEDLOAD && (defined EW_PERIODIC || defined NS_PERIODIC)
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, NST,                   &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    bedldu, bedldv)
#  endif
# endif

      DO i=1,MBEDP
        CALL bc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    bed(:,:,:,i))
      END DO
# ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, Nbed, 1, MBEDP,        &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    bed)
# endif

      CALL bc_r3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, MBOTP,                   &
     &                  bottom)
# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, MBOTP,                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    bottom)
# endif

      RETURN
      END SUBROUTINE sed_bedload_tile
#endif
      END MODULE sed_bedload_mod
