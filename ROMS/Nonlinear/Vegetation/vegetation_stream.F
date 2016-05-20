#include "cppdefs.h"
 
       MODULE vegetation_stream_mod
#if defined VEGETATION && defined VEG_STREAMING 
!
!svn $Id: vegetation_stream.F 429 2015-04-20 17:30:26Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!================================================John C. Warner=========
!================================================Neil K. Ganju  ========
!================================================Alexis Beudin  ========
!==============================================Tarandeep S. Kalra=======
!                                                                      ! 
!  References:                                                         !   
!                                                                      !
!=======================================================================
!                                                                      !
!=======================================================================

      implicit none

      PRIVATE
      PUBLIC  :: vegetation_stream_cal

      CONTAINS
!
!***********************************************************************
      SUBROUTINE vegetation_stream_cal (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
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
      CALL vegetation_stream_tile  (ng, tile,                           &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                         GRID(ng) % angler,                       &
# ifdef SOLVE3D
     &                         GRID(ng) % z_w,                          &
# endif
     &                         FORCES(ng) % Dwave,                      &
     &                         FORCES(ng) % Lwave,                      &
     &                         VEG(ng) % dissip_veg,                    &
     &                         VEG(ng) % Lveg,                          &
     &                         VEG(ng) % BWDXL_veg,                     &
     &                         VEG(ng) % BWDYL_veg)
              
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 16)
# endif
      RETURN
      END SUBROUTINE vegetation_stream_cal 

!***********************************************************************
      SUBROUTINE vegetation_stream_tile  (ng, tile,                     &
     &                              LBi, UBi, LBj, UBj,                 &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                              angler,                             &
# ifdef SOLVE3D
     &                              z_w,                                &
# endif
     &                              Dwave,                              &
     &                              Lwave,                              &
     &                              dissip_veg, Lveg,                   &
     &                              BWDXL_veg, BWDYL_veg)  
!***********************************************************************
!
      USE mod_param
      USE mod_grid
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
      real(r8), intent(in) :: angler(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
#  endif
      real(r8), intent(in) :: Lwave(LBi:,LBj:)
      real(r8), intent(in) :: Dwave(LBi:,LBj:)
      real(r8), intent(in) :: dissip_veg(LBi:,LBj:)
      real(r8), intent(in) :: Lveg(LBi:,LBj:,:)
      real(r8), intent(inout) :: BWDXL_veg(LBi:,LBj:,:)
      real(r8), intent(inout) :: BWDYL_veg(LBi:,LBj:,:)
# else
      real(r8), intent(in)  :: angler(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,UBk)
#  endif
      real(r8), intent(in)  :: Lwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: Dwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: dissip_veg(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Lveg(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: BWDXL_veg(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: BWDYL_veg(LBi:UBi,LBj:UBj,N(ng))
# endif

!  Local variable declarations.
!
      integer :: i, j, k, iveg
      real(r8) :: cff1, cff2
      real(r8) :: EWD_veg
      real(r8), parameter :: Lwave_min = 1.0_r8

      real(r8) :: Dstp
      real(r8) :: waven, wavenx, waveny
      real(r8) :: sigma, osigma

# include "set_bounds.h"
!----------------------------------------------------------------------
!----------Executing the code------------------------------------------
!----------------------------------------------------------------------
!     
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
            Dstp=z_w(i,j,N(ng))-z_w(i,j,0)
!----------------------------------------------------------------------
!  Compute wave amplitude (0.5*Hrms), wave number, intrinsic frequency.
!----------------------------------------------------------------------
            waven=2.0_r8*pi/MAX(Lwave(i,j),Lwave_min)
            cff1=1.5_r8*pi-Dwave(i,j)-angler(i,j)
            wavenx=waven*COS(cff1)
            waveny=waven*SIN(cff1)
            sigma=MIN(SQRT(g*waven*TANH(waven*Dstp)),2.0_r8)
            osigma=1.0_r8/sigma
! 
!----------------------------------------------------------------------
!   Note: Alexis - check if we need a local dissip_veg here 
!   Also Lveg is for 1 veg type only 
!----------------------------------------------------------------------
!
            EWD_veg=dissip_veg(i,j)
            cff2=EWD_veg*osigma*Lveg(i,j,k)
            BWDXL_veg(i,j,k)=cff2*wavenx
            BWDYL_veg(i,j,k)=cff2*waveny
!           
           END DO
        END DO
      END DO    
!
      END SUBROUTINE vegetation_stream_tile
#endif
      END MODULE vegetation_stream_mod
