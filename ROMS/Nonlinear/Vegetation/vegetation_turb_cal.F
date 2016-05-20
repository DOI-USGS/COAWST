#include "cppdefs.h"

      MODULE vegetation_turb_mod

#if defined NONLINEAR && defined VEGETATION && defined VEG_TURB
!
!svn $Id: vegetation_turb_cal.F 429 2015-06-10 12:30:26Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                                              !
!==================================================== John C. Warner ===
!==================================================== Neil K. Ganju  ===
!==================================================== Alexis Beudin  ===
!==================================================Tarandeep S. Kalra===
!                                                                      !
!  This routine computes the turbulent kinetic energy and length scale !
!  modifications due to vegetation for gls_corstep.F                   !
!                                                                      !
!  References:                                                         !
!                                                                      !
!   Uittenbogaard R. (2003): Modelling turbulence in vegetated aquatic !
!   flows. International workshop on RIParian FORest vegetated         !
!   channels: hydraulic, morphological and ecological aspects,         !
!   20-22 February 2003, Trento, Italy.                                !
!                                                                      !
!   Warner J.C., C.R. Sherwood, H.G. Arango, and R.P. Signell (2005):  !
!   Performance of four turbulence closure models implemented using a  !
!   generic length scale method, Ocean Modelling 8: 81-113.            !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: vegetation_turb_cal

      CONTAINS
!
!***********************************************************************
      SUBROUTINE vegetation_turb_cal (ng, tile)
!***********************************************************************
!
      USE mod_stepping 
      USE mod_grid
      USE mod_ocean
      USE mod_mixing 
      USE mod_vegarr
      USE vegetation_drag_mod
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
      CALL vegetation_turb_tile  ( ng, tile,                            &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nstp(ng), nnew(ng),                       &
     &                        OCEAN(ng) % u,                            &
     &                        OCEAN(ng) % v,                            &
     &                        VEG(ng) % ru_loc_veg,                     &
     &                        VEG(ng) % rv_loc_veg,                     &
     &                        VEG(ng) % plant,                          &
# ifdef VEG_FLEX 
     &                        VEG(ng) % bend,                           &
# endif 
     &                        MIXING(ng) % gls,                         &
     &                        MIXING(ng) % tke,                         &
     &                        VEG(ng) % gls_veg,                        &
     &                        VEG(ng) % tke_veg )
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 16)
# endif
      RETURN
      END SUBROUTINE vegetation_turb_cal 
!
!***********************************************************************
      SUBROUTINE vegetation_turb_tile ( ng, tile,                       &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              nstp, nnew,                         &
     &                              u, v,                               &
     &                              ru_loc_veg, rv_loc_veg,             &
     &                              plant,                              &
# ifdef VEG_FLEX 
     &                              bend,                               &
# endif 
     &                              gls, tke,                           &
     &                              gls_veg, tke_veg )
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_vegetation 
      USE mod_vegarr
      USE vegetation_drag_mod
!
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew 
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ru_loc_veg(LBi:,LBj:,:,:)
      real(r8), intent(in) :: rv_loc_veg(LBi:,LBj:,:,:)
      real(r8), intent(in) :: plant(LBi:,LBj:,:,:)
# ifdef VEG_FLEX 
      real(r8), intent(in) :: bend(LBi:,LBj:,:)
# endif 
      real(r8), intent(in) :: gls(LBi:,LBj:,0:,:)
      real(r8), intent(in) :: tke(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: gls_veg(LBi:,LBj:,0:)
      real(r8), intent(inout) :: tke_veg(LBi:,LBj:,0:)
# else
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),nstp)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),nstp)
      real(r8), intent(in) :: ru_loc_veg(LBi:UBi,LBj:UBj,N(ng),NVEG)
      real(r8), intent(in) :: rv_loc_veg(LBi:UBi,LBj:UBj,N(ng),NVEG)
      real(r8), intent(in) :: plant(LBi:UBi,LBj:UBj,NVEG,NVEGP)
# ifdef VEG_FLEX 
      real(r8), intent(in) :: bend(LBi:UBi,LBj:UBj,NVEG)
# endif 
      real(r8), intent(in) :: gls(LBi:UBi,LBj:UBj,0:N(ng),nnew)
      real(r8), intent(in) :: tke(LBi:UBi,LBj:UBj,0:N(ng),nnew)
      real(r8), intent(inout) :: gls_veg(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: tke_veg(LBi:UBi,LBj:UBj,0:N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k, iveg
!
      real(r8), parameter :: one_half=1.0_r8/2.0_r8 
      real(r8), parameter :: one_third=1.0_r8/3.0_r8 
      real(r8), parameter :: Inival=0.0_r8
      real(r8), parameter :: cl_veg=1.0_r8, ck=0.09_r8
      real(r8), parameter :: max_L=10.0e10_r8 
      real(r8), parameter :: eps=1.0e-12_r8 
      real(r8) :: wrku1, wrku2, wrku3, wrku4, wrku
      real(r8) :: wrkv1, wrkv2, wrkv3, wrkv4, wrkv
      real(r8) :: wrk, cff1, cff2, cff3, dissip, inverse_dissip
      real(r8) :: solid, L, eqvegT
      real(r8) :: taufree, tauveg, taueff
      real(r8) :: tke_loc_veg, gls_loc_veg
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vegu
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vegv 

# include "set_bounds.h"
!
      DO k=1,N(ng)-1
        DO j=Jstr,Jend
          DO i=Istr,Iend
            tke_veg(i,j,k)=Inival
            gls_veg(i,j,k)=Inival
          END DO 
        END DO 
      END DO 

      cff1=3.0_r8+gls_p(ng)/gls_n(ng)
      cff2=1.5_r8+gls_m(ng)/gls_n(ng)
      cff3=-1.0_r8/gls_n(ng)

      VEG_LOOP: DO iveg=1,NVEG
        DO k=1,N(ng)-1
          DO j=Jstr,Jend
            DO i=Istr,Iend
!
!-----------------------------------------------------------------------
! Additional turbulence generated by the vegetation = 
! work spent by the fluid against the plants (in m3/s3)
!-----------------------------------------------------------------------
!
              wrku1=ru_loc_veg(i,j,k,iveg)*u(i,j,k,nstp)
              wrku2=ru_loc_veg(i,j,k+1,iveg)*u(i,j,k+1,nstp)
              wrku3=ru_loc_veg(i+1,j,k,iveg)*u(i+1,j,k,nstp)
              wrku4=ru_loc_veg(i+1,j,k+1,iveg)*u(i+1,j,k+1,nstp)
              wrku=0.25_r8*(wrku1+wrku2+wrku3+wrku4)
              wrkv1=rv_loc_veg(i,j,k,iveg)*v(i,j,k,nstp)
              wrkv2=rv_loc_veg(i,j,k+1,iveg)*v(i,j,k+1,nstp)
              wrkv3=rv_loc_veg(i,j+1,k,iveg)*v(i,j+1,k,nstp)
              wrkv4=rv_loc_veg(i,j+1,k+1,iveg)*v(i,j+1,k+1,nstp)
              wrkv=0.25_r8*(wrkv1+wrkv2+wrkv3+wrkv4)
              tke_loc_veg=sqrt(wrku*wrku+wrkv*wrkv)
!
!-----------------------------------------------------------------------
! Dissipation due to vegetation
!-----------------------------------------------------------------------
! Dissipation in GLS (Eq. 12 in Warner et al., 2005)
!
              wrk=MAX(tke(i,j,k,nstp),gls_Kmin(ng))
              dissip=(gls_cmu0(ng)**cff1)*(wrk**cff2)*                  &
     &                 (gls(i,j,k,nstp)**cff3)
              inverse_dissip=1.0_r8/MAX(dissip,eps)
!
! Dissipation time-scale for free turbulence
!
              taufree=wrk*inverse_dissip
!
!# ifdef VEG_FLEX
!
! Equivalent thickness: horizontal projection of the bending plant
! 
!              eqvegT=plant(i,j,iveg,pthck)+sin(bend(i,j,iveg))*         &
!     &                                       plant(i,j,iveg,phght)
!# else
              eqvegT=plant(i,j,iveg,pthck)
!# endif                 
!
!
! Solidity:cross-sectional area of a plant the number of plants per m2
!
!
              solid=plant(i,j,iveg,pdiam)*eqvegT*plant(i,j,iveg,pdens)
!
! Eddies typical size constrained by distance in between the plants
!
              L=cl_veg*((1.0_r8-MIN(solid,1.0_r8))/                     &
     &                 plant(i,j,iveg,pdens))**one_half
              L=MIN(L,max_L)
!
! Dissipation time-scale of eddies in between the plants
!
              tauveg=(L**2.0_r8/(ck**2.0_r8*tke_loc_veg))**one_third
!
! Effective dissipation time-scale
! 
              taueff=MIN(taufree,tauveg)
              gls_loc_veg=gls_c2(ng)*tke_loc_veg/taueff
!
!-----------------------------------------------------------------------
! Add the tke and gls changes from all vegetation types
!-----------------------------------------------------------------------
! 
              tke_veg(i,j,k)=tke_loc_veg + tke_veg(i,j,k)
              gls_veg(i,j,k)=gls_loc_veg + gls_veg(i,j,k)
                                            
            END DO 
          END DO
        END DO 
      END DO VEG_LOOP
!
      RETURN
      END SUBROUTINE vegetation_turb_tile
#endif
      END MODULE vegetation_turb_mod
