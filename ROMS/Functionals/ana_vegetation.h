      SUBROUTINE ana_vegetation (ng, tile, model)
!                                                                      ! 
!! svn $Id: ana_vegetation.h 429 2015-18-05 17:00:25 Z arango $        !
!!=====================================================================!
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!================================================== John C. Warner ====!
!==================================================== Neil K. Ganju  ==!
!==================================================== Alexis Beudin  ==!
!==================================================Tarandeep S. Kalra==!
!                                                                      !
!  Vegetation Model Kernel Variables:                                  !
!  NVEG          Number of vegetation types                            !
!  NVEGP         Varying vegetation properties                         !
!  plant         Vegetation variable properties:                       !
!                   plant(:,:,:,pdiam) => diameter                     !
!                   plant(:,:,:,phght) => height                       !
!                   plant(:,:,:,pdens) => density                      !
!                   plant(:,:,:,pthck) => thickness                    !
!                   plant(:,:,:,pupbm) => above ground biomass         !
!                   plant(:,:,:,pdwbm) => below ground biomass         !
!  marsh_mask    Initialize the mask to get wave thrust on marsh       ! 
!                                                                      !
!  This routine sets initial conditions for vegetation fields          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ncparam
      USE mod_vegarr
      
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_vegetation_tile (ng, tile, model,                        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
#ifdef VEG_DRAG
     &                        VEG(ng) % plant,                          &
#endif
#ifdef MARSH_WAVE_THRUST
     &                        VEG(ng) % marsh_mask                     )
#endif
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(48)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_vegetation
!
!***********************************************************************
      SUBROUTINE ana_vegetation_tile (ng, tile, model,                  &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
#ifdef VEG_DRAG
     &                              plant,                              &
#endif
#ifdef WAVE_THRUST_MARSH
     &                              marsh_mask                          )
#endif
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_vegetation
      USE mod_vegarr
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
#ifdef VEG_DRAG
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: plant(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: plant(LBi:UBi,LBj:UBj,NVEG,NVEGP)
# endif
#endif 
!
#ifdef WAVE_THRUST_MARSH
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: marsh_mask(LBi:,LBj:)
# else
      real(r8), intent(inout) :: marsh_mask(LBi:UBi,LBj:UBj)
# endif 
#endif
!
!  Local variable declarations.
!
#ifdef DISTRIBUTE
      integer :: Tstr, Tend
#endif
      integer :: i, j, k, iveg

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set initial properties for each plant 
!  To have variable properties in array->plant(x,y,iveg,iprop)
!----------------------------------------------------------------------- 
!
#ifdef VEG_DRAG
      DO iveg=1,NVEG
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            plant(i,j,iveg,pdiam)=0.01_r8        !Diameter
            plant(i,j,iveg,phght)=2.0_r8         !Height
            plant(i,j,iveg,pdens)=800.0_r8       !Density
            plant(i,j,iveg,pthck)=0.0005_r8      !Thickness
#ifdef VEGETATION_BIOMASS
            plant(i,j,iveg,pagbm)=0.0_r8         !Above ground Biomass
            plant(i,j,iveg,pbgbm)=0.0_r8         !Below ground Biomass
#endif  
          END DO
        END DO
      END DO
#endif
!      
#ifdef WAVE_THRUST_MARSH
      DO j=Jstr,JendT
        DO i=IstrT,IendT
          marsh_mask(i,j)=1.0_r8 
        END DO 
      END DO  
#endif            
!                                        
      RETURN

      END SUBROUTINE ana_vegetation_tile
