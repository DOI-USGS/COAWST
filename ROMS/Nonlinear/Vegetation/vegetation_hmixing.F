#include "cppdefs.h"
 
       MODULE vegetation_hmixing_mod
#if defined VEGETATION && defined VEG_HMIXING
!
!svn $Id: vegetation_hmixing.F 429 2015-04-20 17:30:26Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!================================================John C. Warner=========
!================================================Neil K. Ganju  ========
!=============================================Alexis  Beudin============
!==============================================Tarandeep S. Kalra=======
!                                                                      !
!  Calculate viscosity change at vegetation iterface and add in        !
!  hmixing.F                                                           ! 
!                                                                      !
!  References:                                                         !   
!                                                                      !
!=======================================================================
!                                                                      !
!                                                                      !
!=======================================================================

      implicit none

      PRIVATE
      PUBLIC  :: vegetation_hmixing_cal

      CONTAINS
!
!***********************************************************************
      SUBROUTINE vegetation_hmixing_cal (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_stepping
      USE mod_grid 
      USE mod_ocean
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
      CALL vegetation_hmixing_tile  (ng, tile,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        VEG(ng) % Lveg,                           &
     &                        VEG(ng) % plant,                          &
     &                        VEG(ng) % visc2d_r_veg,                   &
     &                        VEG(ng) % visc3d_r_veg)
             
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 16)
# endif
      RETURN
      END SUBROUTINE vegetation_hmixing_cal

!***********************************************************************
      SUBROUTINE vegetation_hmixing_tile (ng, tile,                     &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              Lveg,                               &
     &                              plant,                              &
     &                              visc2d_r_veg,                       &
     &                              visc3d_r_veg)          
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_vegetation
      USE mod_vegarr
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Lveg(LBi:,LBj:,:)
      real(r8), intent(in) :: plant(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: visc2d_r_veg(LBi:,LBj:)
      real(r8), intent(inout) :: visc3d_r_veg(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: Lveg(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: plant(LBi:UBi,LBj:UBj,NVEG,NVEGP)
      real(r8), intent(inout) :: visc2d_r_veg(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: visc3d_r_veg(LBi:UBi,LBj:UBj,N(ng))
# endif

!  Local variable declarations.
!
      integer :: i, j, k, iveg
      real(r8) :: visc3d_r_veg_loc, cff 
      real(r8) :: maxgradvegdens, gradu, gradv, wrku, wrkv
      real(r8), parameter :: Inival=0.0_r8
!
# include "set_bounds.h"
!----------------------------------------------------------------------
!----------Executing the code------------------------------------------
!----------------------------------------------------------------------
! 
!   This is currently assuming that vegetation types do not overlap
! 
      visc2d_r_veg=Inival
      visc3d_r_veg=Inival
      visc3d_r_veg_loc=Inival

      DO iveg=1,NVEG
        DO k=1,N(ng)  
          DO j=Jstr,Jend
            DO i=Istr,Iend
              maxgradvegdens=plant(i,j,iveg,pdens)
              maxgradvegdens=MAX(maxgradvegdens,plant(i,j,iveg,pdens))
            END DO
          END DO
        ENDDO
      ENDDO
!
      DO iveg=1,NVEG
        DO k=1,N(ng)  
          DO j=Jstr,Jend
            DO i=Istr,Iend
              cff=VegHMixCoef(iveg,ng)*Lveg(i,j,k)/maxgradvegdens
              gradu=ABS(plant(i,j,iveg,pdens)-plant(i-1,j,iveg,pdens))
              wrku=cff*gradu
              gradv=ABS(plant(i,j,iveg,pdens)-plant(i,j-1,iveg,pdens))
              wrkv=cff*gradv
!
              visc3d_r_veg_loc=SQRT(wrku**2+wrkv**2)
!
!   Adding from multiple veg types
! 
              visc3d_r_veg(i,j,k)=visc3d_r_veg_loc+visc3d_r_veg(i,j,k)
!
!   For overlapping veg types, maximum limit 
!
              visc3d_r_veg(i,j,k)=MIN(visc3d_r_veg(i,j,k),              &
                                      cff*maxgradvegdens)
            END DO
          END DO
        END DO
      END DO 
      END SUBROUTINE vegetation_hmixing_tile
#endif
      END MODULE vegetation_hmixing_mod
