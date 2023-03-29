#include "cppdefs.h"
 
       MODULE marsh_tidal_range_calc_mod
#if defined MARSH_TIDAL_RANGE_CALC
!
!svn $Id: marsh_tidal_range.F 429 2015-04-20 17:30:26Z arango $
!======================================================================!
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!==============================================Tarandeep S. Kalra======!
!=============================================Alfredo Aretxabaleta ====!
!================================================Neil K. Ganju  =======!
!================================================John C. Warner========!
!                                                                      ! 
!  Compute the tidal range and mean high water that are later used     !
!  for biomass production and vertical growth of marsh.                !
!======================================================================!
!                                                                      !

      implicit none

      PRIVATE
      PUBLIC  :: marsh_tidal_range_calc

      CONTAINS
!
!***********************************************************************
      SUBROUTINE marsh_tidal_range_calc (ng, tile)
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
      CALL marsh_tidal_range_calc_tile (ng, tile,                       &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nstp(ng),                                 &
     &                        OCEAN(ng)  % zeta,                        &
     &                        VEG(ng)    % zeta_max1,                   &
     &                        VEG(ng)    % zeta_min1,                   &
     &                        VEG(ng)    % zeta_max_rec,                &
     &                        VEG(ng)    % zeta_min_rec,                &
     &                        VEG(ng)    % counter_loc_rl,              &
     &                        VEG(ng)    % marsh_mask,                  &
     &                        VEG(ng)    % marsh_high_water,            &
     &                        VEG(ng)    % marsh_tidal_range) 
              
      RETURN
      END SUBROUTINE marsh_tidal_range_calc
!
!***********************************************************************
!
      SUBROUTINE marsh_tidal_range_calc_tile  (ng, tile,                &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              nstp,                               &
     &                              zeta,                               &
     &                              zeta_max1,                          &
     &                              zeta_min1,                          &
     &                              zeta_max_rec, zeta_min_rec,         &
     &                              counter_loc_rl,                     &
     &                              marsh_mask,                         &
     &                              marsh_high_water,                   &
     &                              marsh_tidal_range) 
!
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_forces
      USE mod_ocean 
      USE mod_scalars
      USE mod_vegetation
      USE bc_2d_mod, ONLY : bc_r2d_tile
      USE bc_3d_mod, ONLY : bc_r3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp 
!
      real(r8), intent(inout):: counter_loc_rl
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)   :: zeta(LBi:,LBj:,:)
      real(r8), intent(inout):: zeta_max1(LBi:,LBj:)
      real(r8), intent(inout):: zeta_min1(LBi:,LBj:)
      real(r8), intent(inout):: zeta_max_rec(LBi:,LBj:,:)
      real(r8), intent(inout):: zeta_min_rec(LBi:,LBj:,:)
      real(r8), intent(inout):: marsh_mask(LBi:,LBj:)
      real(r8), intent(inout):: marsh_high_water(LBi:,LBj:)
      real(r8), intent(inout):: marsh_tidal_range(LBi:,LBj:)
# else
      real(r8), intent(in)   :: zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout):: zeta_max1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout):: zeta_min1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout)::                                         &
     &                   zeta_max_rec(LBi:UBi,LBj:UBj,NTIMES_MARSH)
      real(r8), intent(inout)::                                         &
     &                    zeta_min_rec(LBi:UBi,LBj:UBj,NTIMES_MARSH)
      real(r8), intent(inout):: marsh_mask(LBi:UBj,LBj:UBj)
      real(r8), intent(inout):: marsh_high_water(LBi:UBj,LBj:UBj)
      real(r8), intent(inout):: marsh_tidal_range(LBi:UBj,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer  :: i, j, it, iic_loc
      integer  :: counter_new_time, counter_loc
      integer  :: freq 
      real(r8) :: ntimes_marsh_rl
      real(r8) :: zeta_loc, zeta_max2, zeta_min, tdays_marsh
      real(r8) :: tot_zeta_max, tot_zeta_min
      real(r8) :: mean_zeta_max, mean_zeta_min
      real(r8) :: time_loc
      real(r8) :: modm_freq, day2sec_loc, day_inp_loc, tdays_loc
!
# include "set_bounds.h"
!
!----------------------------------------------------------------------
!  Executing the code
!----------------------------------------------------------------------
! 
       day2sec_loc=86400.0_r8 ! 86400.0_r8 
       day_inp_loc=1.0_r8     ! record zeta max, min per day 
!       
       time_loc=dt(ng)*iic(ng)                            ! ocean time 
       tdays_loc=MAX((FLOOR(time_loc/(day_inp_loc*day2sec_loc))),1)! time in days 
       freq=CEILING(day_inp_loc*day2sec_loc/dt(ng))        ! time step freq. to record zeta max,min

       modm_freq=MOD(iic(ng),freq)
!
       DO j=Jstr,Jend
        DO i=Istr,Iend
         
          zeta_loc=zeta(i,j,3)     ! store zeta locally first
!          
          ntimes_marsh_rl=REAL(NTIMES_MARSH)
          zeta_max1(i,j)=MAX(zeta_max1(i,j),zeta_loc)
          zeta_min1(i,j)=MIN(zeta_min1(i,j),zeta_loc)
!
! Retain the value of counter_loc_rl in the next time step. 
!
          IF(modm_freq==0.and.counter_loc_rl<=ntimes_marsh_rl) THEN 
            counter_loc_rl=MOD(tdays_loc,ntimes_marsh_rl)
            IF(counter_loc_rl.eq.0.0) THEN 
              counter_loc_rl=ntimes_marsh_rl
            END IF
!
! this has correct value of counter. 
!
            zeta_max_rec(i,j,INT(counter_loc_rl))=(zeta_max1(i,j))
            zeta_min_rec(i,j,INT(counter_loc_rl))=(zeta_min1(i,j))
!
! Once get zeta_max_rec for a day, then reinitialize zeta_max1 
! for the next cycle of saving values.
!
            zeta_max1(i,j)=-10.0_r8
            zeta_min1(i,j)= 10.0_r8
         END IF
          
!          
! Find the mean for the tidal range.
!  
          tot_zeta_max=0.0_r8
          tot_zeta_min=0.0_r8 
          DO it=1,NTIMES_MARSH
            tot_zeta_max=zeta_max_rec(i,j,it)+tot_zeta_max
            tot_zeta_min=zeta_min_rec(i,j,it)+tot_zeta_min
          END DO 
          mean_zeta_max=tot_zeta_max/NTIMES_MARSH
          mean_zeta_min=tot_zeta_min/NTIMES_MARSH
!
          marsh_high_water(i,j)=mean_zeta_max 
          marsh_tidal_range(i,j)=mean_zeta_max-mean_zeta_min 
          marsh_tidal_range(i,j)=marsh_tidal_range(i,j)
!
        END DO
      END DO
!
        CALL bc_r2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    marsh_tidal_range(:,:))
        CALL bc_r2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    marsh_high_water(:,:))

      DO it=1,NTIMES_MARSH
        CALL bc_r2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    zeta_max_rec(:,:,it))
        CALL bc_r2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    zeta_min_rec(:,:,it))
      END DO
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        marsh_tidal_range, marsh_high_water)
!
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, NTIMES_MARSH,          &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    zeta_max_rec, zeta_min_rec)
# endif 
      END SUBROUTINE marsh_tidal_range_calc_tile
#endif 
!        
      END MODULE marsh_tidal_range_calc_mod
