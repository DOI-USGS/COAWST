#include "cppdefs.h"
      MODULE ad_wvelocity_mod
#if defined ADJOINT && defined SOLVE3D
!
!=======================================================================
!  Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                            Hernan G. Arango  !
!=================================================== Andrew M. Moore ===
!                                                                      !
!  This subroutines computes vertical velocity (m/s) at W-points       !
!  from the vertical mass flux (omega*hz/m*n).  This computation       !
!  is done solely for output purposes.                                 !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ad_wvelocity

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_wvelocity (ng, tile, Ninp)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Ninp
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL ad_wvelocity_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        Ninp,                                     &
     &                        GRID(ng) % pm,                            &
     &                        GRID(ng) % pn,                            &
     &                        GRID(ng) % z_r,                           &
     &                        GRID(ng) % z_w,                           &
     &                        GRID(ng) % ad_z_r,                        &
     &                        GRID(ng) % ad_z_w,                        &
     &                        COUPLING(ng) % DU_avg1,                   &
     &                        COUPLING(ng) % DV_avg1,                   &
     &                        COUPLING(ng) % ad_DU_avg1,                &
     &                        COUPLING(ng) % ad_DV_avg1,                &
     &                        OCEAN(ng) % u,                            &
     &                        OCEAN(ng) % v,                            &
     &                        OCEAN(ng) % W,                            &
     &                        OCEAN(ng) % ad_u,                         &
     &                        OCEAN(ng) % ad_v,                         &
     &                        OCEAN(ng) % ad_W,                         &
     &                        OCEAN(ng) % ad_wvel)
      RETURN
      END SUBROUTINE ad_wvelocity
!
!***********************************************************************
      SUBROUTINE ad_wvelocity_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              Ninp,                               &
     &                              pm, pn, z_r, z_w,                   &
     &                              ad_z_r, ad_z_w,                     &
     &                              DU_avg1, DV_avg1,                   &
     &                              ad_DU_avg1, ad_DV_avg1,             &
     &                              u, v, W,                            &
     &                              ad_u, ad_v, ad_W,                   &
     &                              ad_wvel)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_3d_mod,    ONLY : bc_w3d_tile
      USE ad_bc_3d_mod, ONLY : ad_bc_w3d_tile
      USE ad_exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d, ad_mp_exchange3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Ninp
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: DU_avg1(LBi:,LBj:)
      real(r8), intent(in) :: DV_avg1(LBi:,LBj:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
      real(r8), intent(inout) :: ad_z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_w(LBi:,LBj:,0:)
      real(r8), intent(inout) :: ad_DU_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: ad_DV_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_W(LBi:,LBj:,0:)
      real(r8), intent(inout) :: ad_wvel(LBi:,LBj:,0:)
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: DU_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: DV_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: W(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: ad_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: ad_DU_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_DV_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_W(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: ad_wvel(LBi:UBi,LBj:UBj,0:N(ng))
# endif

!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8) :: cff1, cff2, cff3, cff4, cff5, slope , ad_slope
      real(r8) :: adfac, adfac1, adfac2, adfac3

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: vert
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: ad_vert

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_wrk

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute adjoint "true" vertical velocity (m/s).
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
      CALL ad_mp_exchange3d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj, 0, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng), &
     &                       ad_wvel)
# endif
!
!  Set lateral boundary conditions.
!
      CALL ad_bc_w3d_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj, 0, N(ng),                &
     &                     ad_wvel)
!
!  Initialize local adjoint variables and arrays.
!
      ad_slope=0.0_r8

      DO j=JminS,JmaxS
        DO i=IminS,ImaxS
          ad_wrk(i,j)=0.0_r8
        END DO
      END DO
      DO k=1,N(ng)
        DO j=JminS,JmaxS
          DO i=IminS,ImaxS
            ad_vert(i,j,k)=0.0_r8
          END DO
        END DO
      END DO
!
!  Compute vert.
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            wrk(i,j)=u(i,j,k,Ninp)*(z_r(i,j,k)-z_r(i-1,j,k))*           &
     &                             (pm(i-1,j)+pm(i,j))
          END DO
          DO i=Istr,Iend
            vert(i,j,k)=0.25_r8*(wrk(i,j)+wrk(i+1,j))
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            wrk(i,j)=v(i,j,k,Ninp)*(z_r(i,j,k)-z_r(i,j-1,k))*           &
     &                             (pn(i,j-1)+pn(i,j))
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend
            vert(i,j,k)=vert(i,j,k)+0.25_r8*(wrk(i,j)+wrk(i,j+1))
          END DO
        END DO
      END DO
!
      cff1=3.0_r8/8.0_r8
      cff2=3.0_r8/4.0_r8
      cff3=1.0_r8/8.0_r8
      cff4=9.0_r8/16.0_r8
      cff5=1.0_r8/16.0_r8

      J_LOOP : DO j=Jstr,Jend
        DO i=Istr,Iend
          wrk(i,j)=(DU_avg1(i,j)-DU_avg1(i+1,j)+                        &
     &              DV_avg1(i,j)-DV_avg1(i,j+1))/                       &
     &             (z_w(i,j,N(ng))-z_w(i,j,0))
        END DO
        DO i=Istr,Iend
          slope=(z_w(i,j,N(ng))-z_r(i,j,N(ng)  ))/                      &
     &          (z_r(i,j,N(ng))-z_r(i,j,N(ng)-1))  ! extrapolation slope
!>        tl_wvel(i,j,N(ng)-1)=pm(i,j)*pn(i,j)*                         &
!>   &                         (tl_W(i,j,N(ng)-1)+                      &
!>   &                          tl_wrk(i,j)*                            &
!>   &                          (z_w(i,j,N(ng)-1)-z_w(i,j,0))+          &
!>   &                          wrk(i,j)*                               &
!>   &                          (tl_z_w(i,j,N(ng)-1)-tl_z_w(i,j,0)))+   &
!>   &                         cff1*tl_vert(i,j,N(ng)  )+               &
!>   &                         cff2*tl_vert(i,j,N(ng)-1)-               &
!>   &                         cff3*tl_vert(i,j,N(ng)-2)
!>
          adfac=pm(i,j)*pn(i,j)*ad_wvel(i,j,N(ng)-1)
          adfac1=wrk(i,j)*adfac
          ad_W(i,j,N(ng)-1)=ad_W(i,j,N(ng)-1)+adfac
          ad_wrk(i,j)=ad_wrk(i,j)+(z_w(i,j,N(ng)-1)-z_w(i,j,0))*adfac
          ad_z_w(i,j,N(ng)-1)=ad_z_w(i,j,N(ng)-1)+adfac1
          ad_z_w(i,j,0)=ad_z_w(i,j,0)-adfac1
          ad_vert(i,j,N(ng)  )=ad_vert(i,j,N(ng))+                      &
     &                         cff1*ad_wvel(i,j,N(ng)-1)
          ad_vert(i,j,N(ng)-1)=ad_vert(i,j,N(ng)-1)+                    &
     &                         cff2*ad_wvel(i,j,N(ng)-1)
          ad_vert(i,j,N(ng)-2)=ad_vert(i,j,N(ng)-2)-                    &
     &                         cff3*ad_wvel(i,j,N(ng)-1)
          ad_wvel(i,j,N(ng)-1)=0.0_r8
!>        tl_wvel(i,j,N(ng))=pm(i,j)*pn(i,j)*                           &
!>   &                       (tl_wrk(i,j)*                              &
!>   &                        (z_w(i,j,N(ng))-z_w(i,j,0))+              &
!>   &                        wrk(i,j)*                                 &
!>   &                        (tl_z_w(i,j,N(ng))-tl_z_w(i,j,0)))+       &
!>   &                       cff1*(tl_vert(i,j,N(ng))+                  &
!>   &                             tl_slope*(vert(i,j,N(ng)  )-         &
!>   &                                       vert(i,j,N(ng)-1))+        &
!>   &                             slope*(tl_vert(i,j,N(ng)  )-         &
!>   &                                    tl_vert(i,j,N(ng)-1)))+       &
!>   &                       cff2*tl_vert(i,j,N(ng)  )-                 &
!>   &                       cff3*tl_vert(i,j,N(ng)-1)
!>
          adfac=pm(i,j)*pn(i,j)*ad_wvel(i,j,N(ng))
          adfac1=wrk(i,j)*adfac
          adfac2=cff1*ad_wvel(i,j,N(ng))
          adfac3=slope*adfac2
          ad_wrk(i,j)=ad_wrk(i,j)+(z_w(i,j,N(ng))-z_w(i,j,0))*adfac
          ad_z_w(i,j,N(ng))=ad_z_w(i,j,N(ng))+adfac1
          ad_z_w(i,j,0    )=ad_z_w(i,j,0    )-adfac1
          ad_vert(i,j,N(ng))=ad_vert(i,j,N(ng))+adfac2
          ad_slope=ad_slope+                                            &
     &             (vert(i,j,N(ng)  )-vert(i,j,N(ng)-1))*adfac2
          ad_vert(i,j,N(ng)  )=ad_vert(i,j,N(ng)  )+                    &
     &                         adfac3+cff2*ad_wvel(i,j,N(ng))
          ad_vert(i,j,N(ng)-1)=ad_vert(i,j,N(ng)-1)-                    &
     &                         adfac3-cff3*ad_wvel(i,j,N(ng))
          ad_wvel(i,j,N(ng))=0.0_r8
!>        tl_slope=(tl_z_w(i,j,N(ng))-tl_z_r(i,j,N(ng)  ))/             &
!>   &             (z_r(i,j,N(ng))-z_r(i,j,N(ng)-1))-                   &
!>   &             (tl_z_r(i,j,N(ng))-tl_z_r(i,j,N(ng)-1))*slope/       &
!>   &             (z_r(i,j,N(ng))-z_r(i,j,N(ng)-1))
!>
          adfac1=ad_slope/(z_r(i,j,N(ng))-z_r(i,j,N(ng)-1))
          adfac2=slope*adfac1
          ad_z_w(i,j,N(ng))=ad_z_w(i,j,N(ng))+adfac1
          ad_z_r(i,j,N(ng)  )=ad_z_r(i,j,N(ng)  )-adfac1-adfac2
          ad_z_r(i,j,N(ng)-1)=ad_z_r(i,j,N(ng)-1)+adfac2
          ad_slope=0.0_r8
        END DO
        DO k=2,N(ng)-2
          DO i=Istr,Iend
!>          tl_wvel(i,j,k)=pm(i,j)*pn(i,j)*                             &
!>   &                     (tl_W(i,j,k)+                                &
!>   &                      tl_wrk(i,j)*(z_w(i,j,k)-z_w(i,j,0))+        &
!>   &                      wrk(i,j)*(tl_z_w(i,j,k)-tl_z_w(i,j,0)))+    &
!>   &                     cff4*(tl_vert(i,j,k  )+tl_vert(i,j,k+1))-    &
!>   &                     cff5*(tl_vert(i,j,k-1)+tl_vert(i,j,k+2))
!>
            adfac=pm(i,j)*pn(i,j)*ad_wvel(i,j,k)
            adfac1=wrk(i,j)*adfac
            adfac2=cff4*ad_wvel(i,j,k)
            adfac3=cff5*ad_wvel(i,j,k)
            ad_W(i,j,k)=ad_W(i,j,k)+adfac
            ad_wrk(i,j)=ad_wrk(i,j)+(z_w(i,j,k)-z_w(i,j,0))*adfac
            ad_z_w(i,j,k)=ad_z_w(i,j,k)+adfac1
            ad_z_w(i,j,0)=ad_z_w(i,j,0)-adfac1
            ad_vert(i,j,k  )=ad_vert(i,j,k  )+adfac2
            ad_vert(i,j,k+1)=ad_vert(i,j,k+1)+adfac2
            ad_vert(i,j,k-1)=ad_vert(i,j,k-1)-adfac3
            ad_vert(i,j,k+2)=ad_vert(i,j,k+2)-adfac3
            ad_wvel(i,j,k)=0.0_r8
          END DO
        END DO
        DO i=Istr,Iend
          slope=(z_r(i,j,1)-z_w(i,j,0))/                                &
     &          (z_r(i,j,2)-z_r(i,j,1))            ! extrapolation slope
!>        tl_wvel(i,j,1)=pm(i,j)*pn(i,j)*                               &
!>   &                   (tl_W(i,j,1)+                                  &
!>   &                    tl_wrk(i,j)*(z_w(i,j,1)-z_w(i,j,0))+          &
!>   &                    wrk(i,j)*(tl_z_w(i,j,1)-tl_z_w(i,j,0)))+      &
!>   &                   cff1*tl_vert(i,j,1)+                           &
!>   &                   cff2*tl_vert(i,j,2)-                           &
!>   &                   cff3*tl_vert(i,j,3)
!>
          adfac=pm(i,j)*pn(i,j)*ad_wvel(i,j,1)
          adfac1=wrk(i,j)*adfac
          ad_W(i,j,1)=ad_W(i,j,1)+adfac
          ad_wrk(i,j)=ad_wrk(i,j)+(z_w(i,j,1)-z_w(i,j,0))*adfac
          ad_z_w(i,j,1)=ad_z_w(i,j,1)+adfac1
          ad_z_w(i,j,0)=ad_z_w(i,j,0)-adfac1
          ad_vert(i,j,1)=ad_vert(i,j,1)+cff1*ad_wvel(i,j,1)
          ad_vert(i,j,2)=ad_vert(i,j,2)+cff2*ad_wvel(i,j,1)
          ad_vert(i,j,3)=ad_vert(i,j,3)-cff3*ad_wvel(i,j,1)
          ad_wvel(i,j,1)=0.0_r8
!>        tl_wvel(i,j,0)=cff1*(tl_vert(i,j,1)-                          &
!>   &                         tl_slope*(vert(i,j,2)-                   &
!>   &                                   vert(i,j,1))-                  &
!>   &                         slope*(tl_vert(i,j,2)-                   &
!>   &                                tl_vert(i,j,1)))+                 &
!>   &                   cff2*tl_vert(i,j,1)-                           &
!>   &                   cff3*tl_vert(i,j,2)
!>
          adfac=cff1*ad_wvel(i,j,0)
          adfac1=slope*adfac
          ad_vert(i,j,1)=ad_vert(i,j,1)+adfac
          ad_slope=ad_slope-(vert(i,j,2)-vert(i,j,1))*adfac
          ad_vert(i,j,2)=ad_vert(i,j,2)-adfac1
          ad_vert(i,j,1)=ad_vert(i,j,1)+adfac1
          ad_vert(i,j,1)=ad_vert(i,j,1)+cff2*ad_wvel(i,j,0)
          ad_vert(i,j,2)=ad_vert(i,j,2)-cff3*ad_wvel(i,j,0)
          ad_wvel(i,j,0)=0.0_r8
!>        tl_slope=(tl_z_r(i,j,1)-tl_z_w(i,j,0))/                       &
!>   &             (z_r(i,j,2)-z_r(i,j,1))-                             &
!>   &             (tl_z_r(i,j,2)-tl_z_r(i,j,1))*slope/                 &
!>   &             (z_r(i,j,2)-z_r(i,j,1))
!>
          adfac=ad_slope/(z_r(i,j,2)-z_r(i,j,1))
          adfac1=slope*adfac
          ad_z_r(i,j,1)=ad_z_r(i,j,1)+adfac
          ad_z_w(i,j,0)=ad_z_w(i,j,0)-adfac
          ad_z_r(i,j,2)=ad_z_r(i,j,2)-adfac1
          ad_z_r(i,j,1)=ad_z_r(i,j,1)+adfac1
          ad_slope=0.0_r8
        END DO
        DO i=Istr,Iend
!>        tl_wrk(i,j)=(tl_DU_avg1(i,j)-tl_DU_avg1(i+1,j)+               &
!>   &                 tl_DV_avg1(i,j)-tl_DV_avg1(i,j+1))/              &
!>   &                (z_w(i,j,N(ng))-z_w(i,j,0))-                      &
!>   &                (tl_z_w(i,j,N(ng))-tl_z_w(i,j,0))*wrk(i,j)/       &
!>   &                (z_w(i,j,N(ng))-z_w(i,j,0))
!>
          adfac=ad_wrk(i,j)/(z_w(i,j,N(ng))-z_w(i,j,0))
          adfac1=wrk(i,j)*adfac
          ad_DU_avg1(i  ,j)=ad_DU_avg1(i  ,j)+adfac
          ad_DU_avg1(i+1,j)=ad_DU_avg1(i+1,j)-adfac
          ad_DV_avg1(i,j  )=ad_DV_avg1(i,j  )+adfac
          ad_DV_avg1(i,j+1)=ad_DV_avg1(i,j+1)-adfac
          ad_z_w(i,j,N(ng))=ad_z_w(i,j,N(ng))-adfac1
          ad_z_w(i,j,0)=ad_z_w(i,j,0)+adfac1
          ad_wrk(i,j)=0.0_r8
        END DO

      END DO J_LOOP
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
!>          tl_vert(i,j,k)=tl_vert(i,j,k)+                              &
!>   &                      0.25_r8*(tl_wrk(i,j)+tl_wrk(i,j+1))
!>
            adfac=0.25_r8*ad_vert(i,j,k)
            ad_wrk(i,j  )=ad_wrk(i,j  )+adfac
            ad_wrk(i,j+1)=ad_wrk(i,j+1)+adfac
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
!>          tl_wrk(i,j)=(pn(i,j-1)+pn(i,j))*                            &
!>   &                  (tl_v(i,j,k,Ninp)*(z_r(i,j,k)-z_r(i,j-1,k))+    &
!>   &                   v(i,j,k,Ninp)*(tl_z_r(i,j,k)-tl_z_r(i,j-1,k)))
!>
            adfac=(pn(i,j-1)+pn(i,j))*ad_wrk(i,j)
            adfac1=v(i,j,k,Ninp)*adfac
            ad_v(i,j,k,Ninp)=ad_v(i,j,k,Ninp)+                          &
     &                       (z_r(i,j,k)-z_r(i,j-1,k))*adfac
            ad_z_r(i,j  ,k)=ad_z_r(i,j  ,k)+adfac1
            ad_z_r(i,j-1,k)=ad_z_r(i,j-1,k)-adfac1
            ad_wrk(i,j)=0.0_r8
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend
!>          tl_vert(i,j,k)=0.25_r8*(tl_wrk(i,j)+tl_wrk(i+1,j))
!>
            adfac=0.25_r8*ad_vert(i,j,k)
            ad_wrk(i  ,j)=ad_wrk(i  ,j)+adfac
            ad_wrk(i+1,j)=ad_wrk(i+1,j)+adfac
            ad_vert(i,j,k)=0.0_r8
          END DO
          DO i=Istr,Iend+1
!>          tl_wrk(i,j)=(pm(i-1,j)+pm(i,j))*                            &
!>   &                  (tl_u(i,j,k,Ninp)*(z_r(i,j,k)-z_r(i-1,j,k))+    &
!>   &                   u(i,j,k,Ninp)*(tl_z_r(i,j,k)-tl_z_r(i-1,j,k)))
!>
            adfac=(pm(i-1,j)+pm(i,j))*ad_wrk(i,j)
            adfac1=u(i,j,k,Ninp)*adfac
            ad_u(i,j,k,Ninp)=ad_u(i,j,k,Ninp)+                          &
     &                       (z_r(i,j,k)-z_r(i-1,j,k))*adfac
            ad_z_r(i  ,j,k)=ad_z_r(i  ,j,k)+adfac1
            ad_z_r(i-1,j,k)=ad_z_r(i-1,j,k)-adfac1
            ad_wrk(i,j)=0.0_r8
          END DO
        END DO
      END DO
!
!  Exchange time-averaged fields.
!
#  ifdef DISTRIBUTE
      CALL ad_mp_exchange2d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_DU_avg1, ad_DV_avg1)
#  endif

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL ad_exchange_u2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_DU_avg1)
        CALL ad_exchange_v2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_DV_avg1)
      END IF

      RETURN
      END SUBROUTINE ad_wvelocity_tile
#endif
      END MODULE ad_wvelocity_mod
