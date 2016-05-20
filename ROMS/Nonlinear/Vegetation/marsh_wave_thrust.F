#include "cppdefs.h"
 
       MODULE marsh_wave_thrust_mod
#if defined VEGETATION && defined MARSH_WAVE_THRUST
!
!svn $Id: marsh_wave_thrust.F 429 2015-04-20 17:30:26Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!================================================John C. Warner=========
!================================================Neil K. Ganju  ========
!==============================================Tarandeep S. Kalra=======
!                                                                      ! 
!  This routine computes the wave thrust on marshes. Marsh thrust      !
!  values are computed with correction from the wave angle. For each   !
!  cell if one side is sheltered from other cells, that side is not    !
!  exposed to waves. Each cell has four cell normals directed towards  !
!  the center of the cell. The angle of the normals is with respect to !
!  the North and clockwise direction. For a submerged marsh,           !
!  "Tonelli mask" is used to reduce the value of the wave thrust.      !
!                                                                      !
!  References:                                                         !   
!                                                                      !
!=======================================================================
!                                                                      !
!  Tonelli, M., Fagherazzi, Sergio., and Petti., M., 2010: Modeling    !
!  wave impact on salt marsh boundaries, Journal of Geophysical        !
!  Research, 115, 0148-0227.                                           !   
!                                                                      !
!  Dean, R.G. and Dalrymple, R.A., 1991: Water Wave Mechanics for      !
!  Engineers and Scientists, World Scientific Publications             !
!                                                                      !
!=======================================================================

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

# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 16)
# endif
      CALL marsh_wave_thrust_tile  (ng, tile,                           &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nstp(ng),                                 &
     &                        GRID(ng) % h,                             &
     &                        GRID(ng) % angler,                        &
     &                        GRID(ng) % om_r,                          &
     &                        GRID(ng) % on_r,                          &
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
     &                        VEG(ng) % mask_thrust,                    &
     &                        VEG(ng) % Thrust_max,                     &
     &                        VEG(ng) % Thrust_tonelli,                 &
     &                        OCEAN(ng)  % zeta)
              
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 16)
# endif
      RETURN
      END SUBROUTINE marsh_wave_thrust

!***********************************************************************
      SUBROUTINE marsh_wave_thrust_tile  (ng, tile,                     &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              nstp,                               &
     &                              h,angler,                           &
     &                              om_r, on_r ,                        &
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
     &                              mask_thrust,                        &
     &                              Thrust_max, Thrust_tonelli,         &
     &                              zeta)           
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_forces
      USE mod_ocean 
      USE mod_scalars
      USE bc_2d_mod
      USE mod_vegetation
      USE mod_vegarr
!#ifdef DISTRIBUTE
!      USE mp_exchange_mod, ONLY : mp_exchange2d
!#endif

!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp 
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)  :: om_r(LBi:,LBj:)
      real(r8), intent(in)  :: on_r(LBi:,LBj:)
      real(r8), intent(in)  :: h(LBi:,LBj:)
      real(r8), intent(in)  :: angler(LBi:,LBj:)
#  ifdef WAVES_HEIGHT  
      real(r8), intent(in)  :: Hwave(LBi:,LBj:)
#  endif 
#  ifdef WAVES_LENGTH
      real(r8), intent(in)  :: Lwave(LBi:,LBj:)
#  endif 
#  ifdef WAVES_DIR
      real(r8), intent(in)  :: Dwave(LBi:,LBj:)
#  endif 
      real(r8), intent(in)    :: marsh_mask(LBi:,LBj:)
      real(r8), intent(inout) :: mask_thrust(LBi:,LBj:)
      real(r8), intent(inout) :: Thrust_max(LBi:,LBj:)
      real(r8), intent(inout) :: Thrust_tonelli(LBi:,LBj:) 
      real(r8), intent(in)    :: zeta(LBi:,LBj:,:)

# else
      real(r8), intent(in)  :: om_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: on_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: angler(LBi:UBi,LBj:UBj)
#  ifdef WAVES_HEIGHT  
      real(r8), intent(in)  :: Hwave(LBi:UBi,LBj:UBj)
#  endif 
#  ifdef WAVES_HEIGHT  
      real(r8), intent(in)  :: Lwave(LBi:UBi,LBj:UBj)
#  endif 
#  ifdef WAVES_DIR
      real(r8), intent(in)  :: Dwave(LBi:UBi,LBj:UBj)
#  endif 
      real(r8), intent(in)    :: marsh_mask(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: mask_thrust(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Thrust_max(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Thrust_tonelli(LBi:UBi,LBj:UBj) 
      real(r8), intent(in)    :: zeta(LBi:UBi,LBj:UBj,3)
# endif

!  Local variable declarations.
!
      integer :: i,j
      integer, parameter :: zero = 0

      integer :: bound_wd1,bound_wd2,cffi
      integer :: sum_l,sum_r,sum_u,sum_d
      integer :: Isolated_cell_l,Isolated_cell_r
      integer :: Isolated_cell_u,Isolated_cell_d
      integer :: Isolated_cell

      real(r8) :: Kw,Integral_Kp,Fw1,Fw2,Fw,cff
      real(r8) :: depth_all,mask_local_tonelli

      real(r8) :: angler_deg,eft_angle
      real(r8) :: dir_l,dir_r,dir_u,dir_d
      real(r8) :: exp_l,exp_r,exp_u,exp_d
      real(r8) :: thrust_l,thrust_r,thrust_u,thrust_d
      real(r8) :: thrust_w

      integer, dimension(IminS:ImaxS,JminS:JmaxS) :: mask_wd
      integer, dimension(IminS:ImaxS,JminS:JmaxS) :: bound_wdv
      integer, dimension(IminS:ImaxS,JminS:JmaxS) :: bound_wdh
      integer, dimension(IminS:ImaxS,JminS:JmaxS) :: bound_wd
      integer, dimension(IminS:ImaxS,JminS:JmaxS) :: bound_wd_nan

      integer, dimension(IminS:ImaxS,JminS:JmaxS) :: bound_wd_l
      integer, dimension(IminS:ImaxS,JminS:JmaxS) :: bound_wd_r
      integer, dimension(IminS:ImaxS,JminS:JmaxS) :: bound_wd_u
      integer, dimension(IminS:ImaxS,JminS:JmaxS) :: bound_wd_d

# include "set_bounds.h"
!----------------------------------------------------------------------
!----------Executing the code---------------------------------------
!----------------------------------------------------------------------
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          mask_wd(i,j)     = marsh_mask(i,j) 
          bound_wdv(i,j)   = zero          
          bound_wdh(i,j)   = zero
          bound_wd(i,j)    = zero   
          bound_wd_nan(i,j)= zero  
        END DO
      END DO 
!----------------------------------------------------------------------
!             FIND THE LAST WET AND DRY BOUNDARY POINT
!----------------------------------------------------------------------
!----Point=1 is last wet point; point=-1 are last dry point------------
!----------------------------------------------------------------------
! 
!----------------------------------------------------------------------
!----Do the sweep for vertical direction first ------------------------
!----------------------------------------------------------------------
      DO j=Jstr,Jend 
        DO i=Istr,Iend
          cffi=mask_wd(i,j)-mask_wd(i,j+1)
          bound_wd1=MAX(0,cffi)
          bound_wd2=MAX(0,-cffi)
          bound_wdv(i,j)=bound_wd1+bound_wd2
        END DO
      END DO

!----------------------------------------------------------------------
! Do the sweep for horizontal direction first 
!---------------------------------------------------------------------
      DO j=Jstr,Jend
        DO i=Istr,Iend  
          cffi=mask_wd(i,j)-mask_wd(i+1,j)
          bound_wd1=MAX(0,cffi)
          bound_wd2=MAX(0,-cffi)
          bound_wdh(i,j) = bound_wd1 + bound_wd2
        END DO
      END DO

!----------------------------------------------------------------------
!-------These are the points used to compute thrust -------------------
!----------------------------------------------------------------------
      DO j=Jstr,Jend
        DO i=Istr,Iend
          bound_wd(i,j) = bound_wdv(i,j) + bound_wdh(i,j) 
          IF(bound_wd(i,j).ne.0) THEN 
            bound_wd(i,j)=1 
          ENDIF 
          bound_wd_nan(i,j) = bound_wd(i,j) 

!----------------------------------------------------------------------
!------"10" is corresponding to a flag to not include these pts.-------
!------ it can be changed to any number greater than 4 ----------------
!----------------------------------------------------------------------
          IF(bound_wd(i,j).eq.0) THEN     
            bound_wd_nan(i,j)=10          
          ENDIF
        END DO
      END DO

!----------------------------------------------------------------------
!                              FIND NEIGHBOURS
!----------------------------------------------------------------------
      DO j=Jstr,Jend 
        DO i=Istr,Iend 
          bound_wd_l(i,j)=mask_wd(i-1,j)
          bound_wd_r(i,j)=mask_wd(i+1,j)
        END DO 
      END DO
      DO j=Jstr,Jend 
        DO i=Istr,Iend
          bound_wd_u(i,j)=mask_wd(i,j-1)
          bound_wd_d(i,j)=mask_wd(i,j+1)
        END DO 
      END DO
      DO j=Jstr,Jend 
        DO i=Istr,Iend
          Isolated_cell_l = zero
          Isolated_cell_r = zero
          Isolated_cell_u = zero
          Isolated_cell_d = zero
!
          sum_l = mask_wd(i,j)+bound_wd_l(i,j)
          sum_r = mask_wd(i,j)+bound_wd_r(i,j)
          sum_u = mask_wd(i,j)+bound_wd_u(i,j)
          sum_d = mask_wd(i,j)+bound_wd_d(i,j)
!----------------------------------------------------------------------
!---------------------left dir ---------------------------------------
!          if it is not equal to 2 then it is an isolated cell
!----------------------------------------------------------------------
          IF (sum_l.ne.2) THEN
            Isolated_cell_l = 1
            sum_l = 10            
          ELSE 
!----------------------------------------------------------------------
! if the cell has a wet neighbor make sum_l = 1 
!----------------------------------------------------------------------
            sum_l = 1
          ENDIF
!----------------------------------------------------------------------
!-----------------------right dir-------------------------------------
!----------------------------------------------------------------------
          IF (sum_r.ne.2) THEN
            Isolated_cell_r = 1
            sum_r = 10         
          ELSE 
            sum_r = 1
          ENDIF
!----------------------------------------------------------------------
!------------------- up dir-------------------------------------------
!----------------------------------------------------------------------
          IF (sum_u.ne.2) THEN
            Isolated_cell_u = 1
            sum_u = 10         
          ELSE 
            sum_u = 1
          ENDIF
!----------------------------------------------------------------------
!-------------------bottom dir----------------------------------------
!----------------------------------------------------------------------
          IF (sum_d.ne.2) THEN
            Isolated_cell_d = 1
            sum_d = 10            
          ELSE 
            sum_d = 1
          ENDIF
          Isolated_cell=Isolated_cell_l+Isolated_cell_r+                &
     &                  Isolated_cell_u+Isolated_cell_d
!----------------------------------------------------------------------
!---- Have different conditions for isolated cell masking 
!----------------------------------------------------------------------
          IF (Isolated_cell.ne.4) THEN 
            Isolated_cell = 10 
          ELSE
            Isolated_cell = 1 
          ENDIF 

          Isolated_cell = Isolated_cell + mask_wd(i,j)
!----------------------------------------------------------------------
!------if isolated cell=2 this means cell is isolated + wet 
!----------------------------------------------------------------------
          IF (Isolated_cell.ne.2) THEN 
            Isolated_cell = 0
          ELSE 
            Isolated_cell = 1 
          ENDIF 

          IF (Isolated_cell.eq.1) THEN 
            bound_wd(i,j) = 1 
            sum_l  = 1
            sum_r  = 1
            sum_u  = 1
            sum_d  = 1
          ENDIF 
!---------------------------------------------------------------------
! angler is in radians between xi and actual East ! convert to degrees
!---------------------------------------------------------------------
          angler_deg = angler(i,j)*180.0_r8/pi
          dir_l = 270.0_r8 - angler_deg 
          dir_r = 90.0_r8 - angler_deg 
          dir_u = angler_deg
          dir_d = 180.0_r8 - angler_deg

!---------------------------------------------------------------------
!                     CALCULATE THE MARSH THRUST
!---------------------------------------------------------------------
          kw = 2.0_r8*pi/Lwave(i,j)              
          Integral_kp = sinh(kw*h(i,j))/(kw*cosh(h(i,j)*kw))   
          Fw1 = rho0*g*Hwave(i,j)*Integral_kp*0.001_r8   
          Fw2 = (rho0*g*Hwave(i,j))*Hwave(i,j)*0.5_r8*0.001_r8  
!---------------------------------------------------------------------
!         Total wave thrust at mean sea level 
!---------------------------------------------------------------------
          Fw  = Fw1 + Fw2 

          IF (sum_l.ne.10) THEN             
            exp_l=dir_l*sum_l*bound_wd(i,j)
            eft_angle = exp_l - Dwave(i,j)*rad2deg
            thrust_l=abs(Fw*cos(eft_angle*deg2rad))
!---------------------------------------------------------------------
!       If 90>angle>270 then waves arrive from opposite direction with 
!       respect to the cell side so impact other side of cell 
!---------------------------------------------------------------------
            IF (abs(eft_angle).ge.90.0_r8.and.                          &
     &          abs(eft_angle).le.270.0_r8) THEN
              thrust_l = 0.0_r8 
            ENDIF 
          ELSE  
            thrust_l = 0.0_r8
          ENDIF 

          IF (sum_r.ne.10) THEN           
            exp_r=dir_r*sum_r*bound_wd(i,j)
            eft_angle = exp_r - Dwave(i,j)*rad2deg 
            thrust_r=abs(Fw*cos(eft_angle*deg2rad))
            IF (abs(eft_angle).ge.90.0_r8.and.                          &
     &          abs(eft_angle).le.270.0_r8) THEN           
              thrust_r = 0.0_r8 
            ENDIF 
          ELSE  
            thrust_r = 0.0_r8
          ENDIF     
 
          IF (sum_u.ne.10) THEN           
            exp_u=dir_u*sum_u*bound_wd(i,j)
            eft_angle = exp_u - Dwave(i,j)*rad2deg 
            thrust_u=abs(Fw*cos(eft_angle*deg2rad))
            IF (abs(eft_angle).ge.90.0_r8.and.                          & 
     &          abs(eft_angle).le.270.0_r8) THEN           
              thrust_u = 0.0_r8 
            ENDIF 
          ELSE  
              thrust_u = 0.0_r8
          ENDIF 

          IF (sum_d.ne.10) THEN           
            exp_d=dir_d*sum_d*bound_wd(i,j)
            eft_angle = exp_d - Dwave(i,j)*rad2deg 
            thrust_d=abs(Fw*cos(eft_angle*deg2rad))
            IF (abs(eft_angle).ge.90.0_r8.and.                          &
     &          abs(eft_angle).le.270.0_r8) THEN           
              thrust_d = 0.0_r8 
            ENDIF 
          ELSE  
            thrust_d = 0.0_r8
          ENDIF 
          thrust_w = thrust_l+thrust_r+thrust_u+thrust_d
!---------------------------------------------------------------------
!   if marsh is submerged in water depending on depth, wave thrust to 
!   be reduced 
!---------------------------------------------------------------------
          depth_all=h(i,j)+zeta(i,j,3)
          IF (depth_all.lt.0.2_r8) THEN
            cff=1.0_r8-0.45_r8*depth_all*5.0_r8
          ELSEIF (0.2_r8.lt.depth_all.and.depth_all.lt.0.4_r8) THEN
            cff=0.55_r8*(1.0_r8-2.5_r8*(depth_all-0.2_r8))
          ELSE 
            cff=0.275_r8 
          ENDIF
          mask_thrust(i,j)=cff
!
          IF (bound_wd_nan(i,j).ne.10) THEN   
            Thrust_tonelli(i,j) = bound_wd_nan(i,j)*thrust_w*           &
     &                            mask_thrust(i,j)  

            Thrust_max(i,j) = bound_wd_nan(i,j)*thrust_w 
          ELSE
            Thrust_tonelli(i,j) = 0.0_r8  
            Thrust_max(i,j)     = 0.0_r8 
          ENDIF
        END DO
      END DO
!---------------------------------------------------------------------
!  Apply periodic or gradient boundary conditions for output
!  purposes only.
!---------------------------------------------------------------------
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  mask_thrust)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Thrust_max)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Thrust_tonelli)
!# ifdef DISTRIBUTE
!      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
!     &                    LBi, UBi, LBj, UBj,                           &
!     &                    NghostPoints,                                 &
!     &                    EWperiodic(ng), NSperiodic(ng),               &
!     &                    mask_thrust, Thrust_max,                      &
!     &                    Thrust_tonelli)
!# endif
      END SUBROUTINE marsh_wave_thrust_tile
#endif
      END MODULE marsh_wave_thrust_mod
