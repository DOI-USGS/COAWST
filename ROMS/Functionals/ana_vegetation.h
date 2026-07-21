!
      SUBROUTINE ana_vegetation (ng, tile, model)
!
!! git $Id$
!!================================================== John C. Warner ====
!! Copyright (c) 2002-2026 The ROMS Group             Neil K. Ganju    !
!!   Licensed under a MIT/X style license             Alexis Beudin    !
!!   See License_ROMS.md                         Tarandeep S. Kalra    !
!=======================================================================
!                                                                      !
!  It sets analytical initial condition for the submerge aquatic       !
!  vegetation model.                                                   !
!                                                                      !
!  Basic model parameters and variables:                               !
!                                                                      !
!  NVEG          Number of kernel aquatic vegetation types             !
!                                                                      !
!  NVEGP         Number of vegetation properties, usually NVEGP=4      !
!                                                                      !
!  plant         Vegetation variable properties:                       !
!                   plant(:,:,:,isDens) => density, individuals per m2 !
!                   plant(:,:,:,isDiam) => mean diameter (m)           !
!                   plant(:,:,:,isHght) => mean height (m)             !
!                   plant(:,:,:,isThck) => mean thickness (m)          !
!                                                                      !
!  marsh_mask    Marsh mask: [1] marsh cell, [0] non-marsh cell        ! 
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ncparam
      USE mod_vegetation
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_vegetation_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
#ifdef MARSH_WAVE_THRUST
     &                          VEG(ng) % marsh_mask,                   &
#endif
     &                          VEG(ng) % plant)
!
!  Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(48)=__FILE__
      END IF
!
      RETURN
      END SUBROUTINE ana_vegetation
!
!***********************************************************************
      SUBROUTINE ana_vegetation_tile (ng, tile, model,                  &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
#ifdef WAVE_THRUST_MARSH
     &                                marsh_mask,                       &
#endif
     &                                plant)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_vegetation
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: plant(LBi:,LBj:,:,:)
# ifdef WAVE_THRUST_MARSH
      real(r8), intent(inout) :: marsh_mask(LBi:,LBj:)
# endif
#else
      real(r8), intent(inout) :: plant(LBi:UBi,LBj:UBj,NVEG,NVEGP)
# ifdef WAVE_THRUST_MARSH
      real(r8), intent(inout) :: marsh_mask(LBi:UBi,LBj:UBj)
# endif
#endif 
!
!  Local variable declarations.
!
      integer :: i, j, k, iveg
#ifdef DISTRIBUTE
      integer :: Tstr, Tend
#endif
!
#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set initial dominant plant properties for each vegetation type.
!----------------------------------------------------------------------- 
!
#ifdef VEG_DRAG
      DO iveg=1,NVEG
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            plant(i,j,iveg,isDens)=800.0_r8       ! Density
            plant(i,j,iveg,isDiam)=0.01_r8        ! Diameter
            plant(i,j,iveg,isHght)=2.0_r8         ! Height
            plant(i,j,iveg,isThck)=0.0005_r8      ! Thickness
          END DO
        END DO
      END DO
#endif

#ifdef WAVE_THRUST_MARSH
!
!-----------------------------------------------------------------------
!  Set initial marsh mask.
!----------------------------------------------------------------------- 
!
      DO j=Jstr,JendT
        DO i=IstrT,IendT
          marsh_mask(i,j)=1.0_r8 
        END DO 
      END DO  
#endif            
!                                        
      RETURN
      END SUBROUTINE ana_vegetation_tile
