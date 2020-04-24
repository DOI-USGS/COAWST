#include "cppdefs.h"
 
       MODULE marsh_wave_thrust_mod
#if defined MARSH_DYNAMICS && defined MARSH_WAVE_THRUST 
!
!svn $Id: marsh_wave_thrust.F 429 2019-04-20 17:30:26Z arango $
!======================================================================!
!  Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!==============================================Tarandeep S. Kalra======!
!================================================Neil K. Ganju  =======!
!================================================John C. Warner========!
!==============================================Julia M. Moriarty=======!
!                                                                      ! 
!  This routine computes the wave thrust on marsh edge from wave       !
!  climate from the adjacent water cells.                              !
!  Thrust is modified based on water level of the adjacent water cell  !
!  providing wave climate and bathymetry of marsh cell.                !
!                                                                      !
!  References:                                                         !   
!                                                                      !
!======================================================================!
!  Dean, R.G. and Dalrymple, R.A., 1991: Water Wave Mechanics for      !
!  Engineers and Scientists, World Scientific Publications.            !
!                                                                      !
!  Tonelli, M., Fagherazzi, Sergio., and Petti., M., 2010:             !
!  Modeling wave impact on salt marsh boundaries, Journal of           !
!  Geophysical Research, 115, 0148-0227.                               !
!======================================================================!
!                                                                      ! 
      implicit none

      PRIVATE
      PUBLIC  :: marsh_wave_thrust

      CONTAINS
!
!***********************************************************************
      SUBROUTINE marsh_wave_thrust (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ocean 
      USE mod_stepping
      USE mod_vegarr
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL marsh_wave_thrust_tile  (ng, tile,                           &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nstp(ng),                                 &
     &                        GRID(ng) % h,                             &
     &                        GRID(ng) % angler,                        &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
# endif
# ifdef WET_DRY
     &                       GRID(ng) % rmask_wet,                      &
# endif
# if defined WAVES_HEIGHT 
     &                        FORCES(ng) % Hwave,                       &
# endif 
# if defined WAVES_LENGTH 
     &                        FORCES(ng) % Lwave,                       &
# endif 
# if defined WAVES_DIR 
     &                        FORCES(ng) % Dwave,                       &
# endif 
     &                        VEG(ng) % marsh_mask,                     &
     &                        VEG(ng) % umask_marsh,                    &
     &                        VEG(ng) % vmask_marsh,                    &
     &                        VEG(ng) % Thrust_xi,                      &
     &                        VEG(ng) % Thrust_eta,                     &
     &                        VEG(ng) % Thrust_total,                   &
     &                        OCEAN(ng)  % zeta)
              
      RETURN
      END SUBROUTINE marsh_wave_thrust
!
!***********************************************************************
!
      SUBROUTINE marsh_wave_thrust_tile  (ng, tile,                     &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              nstp,                               &
     &                              h,angler,                           &
# ifdef MASKING
     &                             rmask, umask, vmask,                 &
# endif
# ifdef WET_DRY
     &                             rmask_wet,                           &
# endif
# if defined WAVES_HEIGHT  
     &                              Hwave,                              &
# endif 
# if defined WAVES_LENGTH  
     &                              Lwave,                              &
# endif 
# if defined WAVES_DIR  
     &                              Dwave,                              &
# endif 
     &                              marsh_mask,                         &
     &                              umask_marsh,                        &
     &                              vmask_marsh,                        &
     &                              Thrust_xi,                          &
     &                              Thrust_eta,                         &
     &                              Thrust_total,                       &
     &                              zeta)           
!
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_forces
      USE mod_ocean 
      USE mod_scalars
      USE mod_vegetation
      USE mod_vegarr
      USE bc_2d_mod
      USE exchange_2d_mod
!
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp 
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)  :: h(LBi:,LBj:)
      real(r8), intent(in)  :: angler(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
#  endif
#  ifdef WAVES_HEIGHT  
      real(r8), intent(in)  :: Hwave(LBi:,LBj:)
#  endif 
#  ifdef WAVES_LENGTH
      real(r8), intent(in)  :: Lwave(LBi:,LBj:)
#  endif 
#  ifdef WAVES_DIR
      real(r8), intent(in)  :: Dwave(LBi:,LBj:)
#  endif 
      real(r8), intent(in)       :: marsh_mask(LBi:,LBj:)
      real(r8), intent(inout)    :: umask_marsh(LBi:,LBj:)
      real(r8), intent(inout)    :: vmask_marsh(LBi:,LBj:)
!
      real(r8), intent(inout)    :: Thrust_xi(LBi:,LBj:)
      real(r8), intent(inout)    :: Thrust_eta(LBi:,LBj:)
      real(r8), intent(inout)    :: Thrust_total(LBi:,LBj:)
      real(r8), intent(in)       :: zeta(LBi:,LBj:,:)
# else
      real(r8), intent(in)  :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: angler(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
#  endif 
#  ifdef WAVES_HEIGHT  
      real(r8), intent(in)  :: Hwave(LBi:UBi,LBj:UBj)
#  endif 
#  ifdef WAVES_HEIGHT  
      real(r8), intent(in)  :: Lwave(LBi:UBi,LBj:UBj)
#  endif 
#  ifdef WAVES_DIR
      real(r8), intent(in)  :: Dwave(LBi:UBi,LBj:UBj)
#  endif 
      real(r8), intent(in)     :: marsh_mask(LBi:UBi,LBj:UBj)
      real(r8), intent(inout)    :: umask_marsh(LBi:UBi,LBj:UBj)
      real(r8), intent(inout)    :: vmask_marsh(LBi:UBi,LBj:UBj)
     
      real(r8), intent(inout) :: Thrust_xi(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Thrust_eta(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Thrust_total(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i,j

      real(r8), parameter :: Inival=0.0_r8
      real(r8), parameter :: half=0.5_r8
      real(r8), parameter :: N_kN=0.001_r8
      real(r8), parameter :: eps=1.0e-12_r8
!      real(r8), parameter :: depth_const=0.6_r8

      real(r8) :: Kw, Integral_Kp 
      real(r8) :: depth_all

      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5, cff6, cff7, cff8

      real(r8) :: F_asl, F_bsl, F_msl

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FY

# include "set_bounds.h"
!
!----------------------------------------------------------------------
!  Executing the code
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!  Marsh mask at U-points and V-points.
!----------------------------------------------------------------------
!
      DO j=JstrR,JendR
         DO i=Istr,IendR
          umask_marsh(i,j)=marsh_mask(i-1,j)+marsh_mask(i,j)
          IF (umask_marsh(i,j).eq.1.0_r8) THEN
            umask_marsh(i,j)=marsh_mask(i-1,j)-marsh_mask(i,j)
          ELSE
            umask_marsh(i,j)=0.0_r8 
          END IF
        END DO
      END DO
!
!
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          vmask_marsh(i,j)=marsh_mask(i,j-1)+marsh_mask(i,j)
          IF (vmask_marsh(i,j).eq.1.0_r8) THEN
            vmask_marsh(i,j)=marsh_mask(i,j-1)-marsh_mask(i,j)
          ELSE
            vmask_marsh(i,j)=0.0_r8
          END IF
        END DO
      END DO
!
! Exchange boundary data 
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          umask_marsh(:,:))
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vmask_marsh(:,:))
      END IF
!
# ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    umask_marsh, vmask_marsh)
# endif
!
!----------------------------------------------------------------------
!  Compute total thrust on all cell centers
!  where wave data is available at RHO points.
!----------------------------------------------------------------------
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          cff=1.5_r8*pi-Dwave(i,j)-angler(i,j)
!
!----------------------------------------------------------------------
!  Calculate Thrust below and above sea level.
!----------------------------------------------------------------------
!
          kw=2.0_r8*pi/Lwave(i,j)
!          cff1=MAX(0.0_r8,h(i,j)+zeta(i,j,1))
!          cff1=h(i,j)
!          Integral_kp=ABS(sinh(kw*cff1)/(kw*cosh(kw*cff1)))
          cff1=h(i,j)+zeta(i,j,1)
          Integral_kp=(sinh(kw*cff1)/(kw*cosh(kw*h(i,j))))
!
          cff2=rho0*g*Hwave(i,j)*N_kN
          F_asl=half*cff2*Hwave(i,j)
          F_bsl=cff2*Integral_kp
!
!---------------------------------------------------------------------
!  Total wave thrust at mean sea level.
!---------------------------------------------------------------------
!
          F_msl=F_asl+F_bsl
          FX(i,j)=F_msl*cos(cff)
          FY(i,j)=F_msl*sin(cff)
! 
        END DO
      END DO
!
!---------------------------------------------------------------------
! Get the wave thrust on cell faces
!---------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
!
          cff1=SIGN(1.0_r8,umask_marsh(i,j))
          cff2=half*(1.0_r8+cff1)
          cff3=half*(1.0_r8-cff1)
!
          cff4=SIGN(1.0_r8,FX(i,j))
          cff5=half*(1.0_r8-cff4)
!
          cff6=SIGN(1.0_r8,FX(i-1,j))
          cff7=half*(1.0_r8+cff6)
!
          Thrust_xi(i,j)=ABS(umask_marsh(i,j))*                         &
     &                      (cff2*cff5*FX(i,j)+cff3*cff7*FX(i-1,j))
!
! Account for marsh cell depth by exponentially decreasing 
! thrust if marsh is fully submerged. Otherwise, constant thrust.
!
          depth_all=cff2*(h(i-1,j  )+zeta(i,j,1)  )+                    &
     &              cff3*(h(i,j    )+zeta(i-1,j,1))
          IF(depth_all.lt.0.0_r8) THEN
            cff8=1.0_r8
          ELSEIF(depth_all.ge.0.0_r8) THEN
            cff8=exp(-depth_all*3.0_r8)
          ENDIF
!
!---------------------------------------------------------------------
! Modify xi cell face thrust based on water depth.
!---------------------------------------------------------------------
!
          Thrust_xi(i,j)=cff8*Thrust_xi(i,j)
!
        END DO
      END DO
!
      DO j=Jstr,JendR
        DO i=IstrR,IendR
!
          cff1=SIGN(1.0_r8,vmask_marsh(i,j))
          cff2=half*(1.0_r8+cff1)
          cff3=half*(1.0_r8-cff1)
!
          cff4=SIGN(1.0_r8,FY(i,j))
          cff5=half*(1.0_r8-cff4)
!
          cff6=SIGN(1.0_r8,FY(i,j-1))
          cff7=half*(1.0_r8+cff6)
!
          Thrust_eta(i,j)=ABS(vmask_marsh(i,j))*                        &
     &                       (cff3*cff7*FY(i,j-1)+cff2*cff5*FY(i,j))
!
! Account for marsh cell depth by exponentially decreasing 
! thrust if marsh is fully submerged. Otherwise, constant thrust.
!
          depth_all=cff2*(h(i,j-1)+zeta(i,j,1  ))+                      &
     &              cff3*(h(i,j  )+zeta(i,j-1,1))
!
          IF(depth_all.lt.0.0_r8) THEN
            cff8=1.0_r8
          ELSEIF(depth_all.ge.0.0_r8) THEN
            cff8=exp(-depth_all*3.0_r8)
          ENDIF
!
!---------------------------------------------------------------------
! Modify eta cell face thrust based on water depth.
!---------------------------------------------------------------------
!
          Thrust_eta(i,j)=cff8*Thrust_eta(i,j)
        END DO
      END DO
!
!
! Exchange boundary data 
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Thrust_xi(:,:))
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Thrust_eta(:,:))
      END IF
!
# ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Thrust_xi, Thrust_eta)
# endif

!
!---------------------------------------------------------------------
! Sum the thrust from all cell faces to get thrust at cell center
!---------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          cff=ABS(Thrust_xi(i,j ))+ABS(Thrust_xi(i+1,j ))+              &
     &        ABS(Thrust_eta(i,j))+ABS(Thrust_eta(i,j+1))
          Thrust_total(i,j)=cff*marsh_mask(i,j)
        END DO
      END DO
!
!---------------------------------------------------------------------
!  Apply periodic or gradient boundary conditions for output
!  purposes only.
!---------------------------------------------------------------------
!
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Thrust_xi)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Thrust_eta)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Thrust_total)
! 
# ifdef DISTRIBUTE
!
! Exchange boundary data 
!
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Thrust_total)
!
# endif
!
      END SUBROUTINE marsh_wave_thrust_tile
# endif 
      END MODULE marsh_wave_thrust_mod
