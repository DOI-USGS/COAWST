#include "cppdefs.h"
      MODULE step3d_uv_mod
#if defined NONLINEAR && defined SOLVE3D
!
!svn $Id: step3d_uv.F 732 2008-09-07 01:55:51Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!==================================================== John C. Warner ===
!                                                                      !
!  This subroutine time-steps the nonlinear  horizontal  momentum      !
!  equations. The vertical viscosity terms are time-stepped using      !
!  an implicit algorithm.                                              !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: step3d_uv
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE step3d_uv (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
# ifdef DIAGNOSTICS_UV
      USE mod_diags
# endif
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
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
      CALL wclock_on (ng, iNLM, 34)
# endif
      CALL step3d_uv_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nrhs(ng), nstp(ng), nnew(ng),                &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
# endif
# ifdef WET_DRY
     &                     GRID(ng) % umask_wet,                        &
     &                     GRID(ng) % vmask_wet,                        &
# endif
     &                     GRID(ng) % om_v,                             &
     &                     GRID(ng) % on_u,                             &
     &                     GRID(ng) % pm,                               &
     &                     GRID(ng) % pn,                               &
     &                     GRID(ng) % Hz,                               &
     &                     GRID(ng) % z_r,                              &
     &                     GRID(ng) % z_w,                              &
# ifdef UV_WAVEDRAG
     &                     GRID(ng) % wavedrag,                         &
# endif
     &                     MIXING(ng) % Akv,                            &
     &                     COUPLING(ng) % DU_avg1,                      &
     &                     COUPLING(ng) % DV_avg1,                      &
     &                     COUPLING(ng) % DU_avg2,                      &
     &                     COUPLING(ng) % DV_avg2,                      &
     &                     OCEAN(ng) % ru,                              &
     &                     OCEAN(ng) % rv,                              &
# ifdef DIAGNOSTICS_UV
     &                     DIAGS(ng) % DiaU2wrk,                        &
     &                     DIAGS(ng) % DiaV2wrk,                        &
     &                     DIAGS(ng) % DiaU2int,                        &
     &                     DIAGS(ng) % DiaV2int,                        &
     &                     DIAGS(ng) % DiaU3wrk,                        &
     &                     DIAGS(ng) % DiaV3wrk,                        &
     &                     DIAGS(ng) % DiaRU,                           &
     &                     DIAGS(ng) % DiaRV,                           &
# endif
     &                     OCEAN(ng) % u,                               &
     &                     OCEAN(ng) % v,                               &
     &                     OCEAN(ng) % ubar,                            &
     &                     OCEAN(ng) % vbar,                            &
# ifdef WEC
     &                     OCEAN(ng) % ubar_stokes,                     &
     &                     OCEAN(ng) % vbar_stokes,                     &
     &                     OCEAN(ng) % u_stokes,                        &
     &                     OCEAN(ng) % v_stokes,                        &
# endif
     &                     GRID(ng) % Huon,                             &
     &                     GRID(ng) % Hvom)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 34)
# endif

      RETURN
      END SUBROUTINE step3d_uv
!
!***********************************************************************
      SUBROUTINE step3d_uv_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nrhs, nstp, nnew,                      &
# ifdef MASKING
     &                           umask, vmask,                          &
# endif
# ifdef WET_DRY
     &                           umask_wet, vmask_wet,                  &
# endif
     &                           om_v, on_u, pm, pn,                    &
     &                           Hz, z_r, z_w,                          &
# ifdef UV_WAVEDRAG
     &                           wavedrag,                              &
# endif
     &                           Akv,                                   &
     &                           DU_avg1, DV_avg1,                      &
     &                           DU_avg2, DV_avg2,                      &
     &                           ru, rv,                                &
# ifdef DIAGNOSTICS_UV
     &                           DiaU2wrk, DiaV2wrk,                    &
     &                           DiaU2int, DiaV2int,                    &
     &                           DiaU3wrk, DiaV3wrk,                    &
     &                           DiaRU,    DiaRV,                       &
# endif
     &                           u, v,                                  &
     &                           ubar, vbar,                            &
# ifdef WEC
     &                           ubar_stokes, vbar_stokes,              &
     &                           u_stokes, v_stokes,                    &
# endif
     &                           Huon, Hvom)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_sources
!
      USE exchange_2d_mod
      USE exchange_3d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
# endif
      USE u3dbc_mod, ONLY : u3dbc_tile
      USE v3dbc_mod, ONLY : v3dbc_tile
# ifdef UV_WAVEDRAG
      USE mod_tides, only : TIDES, drag_scale
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, nstp, nnew
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
#  ifdef UV_WAVEDRAG
      real(r8), intent(in) :: wavedrag(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Akv(LBi:,LBj:,0:)
      real(r8), intent(in) :: DU_avg1(LBi:,LBj:)
      real(r8), intent(in) :: DV_avg1(LBi:,LBj:)
      real(r8), intent(in) :: DU_avg2(LBi:,LBj:)
      real(r8), intent(in) :: DV_avg2(LBi:,LBj:)
      real(r8), intent(inout) :: ru(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: rv(LBi:,LBj:,0:,:)
#  ifdef WEC
      real(r8), intent(in) :: ubar_stokes(LBi:,LBj:)
      real(r8), intent(in) :: vbar_stokes(LBi:,LBj:)
#  endif
#  ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaU2wrk(LBi:,LBj:,:)
      real(r8), intent(inout) :: DiaV2wrk(LBi:,LBj:,:)
      real(r8), intent(inout) :: DiaU2int(LBi:,LBj:,:)
      real(r8), intent(inout) :: DiaV2int(LBi:,LBj:,:)
      real(r8), intent(inout) :: DiaU3wrk(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: DiaV3wrk(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: DiaRU(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: DiaRV(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
#  ifdef WEC
      real(r8), intent(inout) :: u_stokes(LBi:,LBj:,:)
      real(r8), intent(inout) :: v_stokes(LBi:,LBj:,:)
#  endif
      real(r8), intent(out) :: ubar(LBi:,LBj:,:)
      real(r8), intent(out) :: vbar(LBi:,LBj:,:)
      real(r8), intent(out) :: Huon(LBi:,LBj:,:)
      real(r8), intent(out) :: Hvom(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#  ifdef UV_WAVEDRAG
      real(r8), intent(in) :: wavedrag(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Akv(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: DU_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: DV_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: DU_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: DV_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ru(LBi:UBi,LBj:UBj,0:N(ng),2)
      real(r8), intent(inout) :: rv(LBi:UBi,LBj:UBj,0:N(ng),2)
#  ifdef WEC
      real(r8), intent(in) :: ubar_stokes(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vbar_stokes(LBi:UBi,LBj:UBj)
#  endif
#  ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaU2wrk(LBi:UBi,LBj:UBj,NDM2d)
      real(r8), intent(inout) :: DiaV2wrk(LBi:UBi,LBj:UBj,NDM2d)
      real(r8), intent(inout) :: DiaU2int(LBi:UBi,LBj:UBj,NDM2d)
      real(r8), intent(inout) :: DiaV2int(LBi:UBi,LBj:UBj,NDM2d)
      real(r8), intent(inout) :: DiaU3wrk(LBi:UBi,LBj:UBj,N(ng),NDM3d)
      real(r8), intent(inout) :: DiaV3wrk(LBi:UBi,LBj:UBj,N(ng),NDM3d)
      real(r8), intent(inout) :: DiaRU(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
      real(r8), intent(inout) :: DiaRV(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
#  endif
      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#  ifdef WEC
      real(r8), intent(inout) :: u_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: v_stokes(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(out) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(out) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(out) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, idiag, is, j, k

      real(r8) :: cff, cff1, cff2

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: AK
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: BC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
# ifdef WEC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CFs
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DCs
# endif
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hzk
      real(r8), dimension(IminS:ImaxS,N(ng)) :: oHz
# ifdef DIAGNOSTICS_UV
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: wrk
      real(r8), dimension(IminS:ImaxS,1:NDM2d-1) :: Dwrk
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Time step momentum equation in the XI-direction.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          AK(i,0)=0.5_r8*(Akv(i-1,j,0)+                                 &
     &                    Akv(i  ,j,0))
          DO k=1,N(ng)
            AK(i,k)=0.5_r8*(Akv(i-1,j,k)+                               &
     &                      Akv(i  ,j,k))
            Hzk(i,k)=0.5_r8*(Hz(i-1,j,k)+                               &
     &                       Hz(i  ,j,k))
# if defined SPLINES_VVISC || defined DIAGNOSTICS_UV
            oHz(i,k)=1.0_r8/Hzk(i,k)
# endif
          END DO
        END DO
!
!  Time step right-hand-side terms.
!
        IF (iic(ng).eq.ntfirst(ng)) THEN
          cff=0.25_r8*dt(ng)
        ELSE IF (iic(ng).eq.(ntfirst(ng)+1)) THEN
          cff=0.25_r8*dt(ng)*3.0_r8/2.0_r8
        ELSE
          cff=0.25_r8*dt(ng)*23.0_r8/12.0_r8
        END IF
        DO i=IstrU,Iend
          DC(i,0)=cff*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
        END DO
        DO k=1,N(ng)
          DO i=IstrU,Iend
            u(i,j,k,nnew)=u(i,j,k,nnew)+                                &
     &                    DC(i,0)*ru(i,j,k,nrhs)
# ifdef SPLINES_VVISC
            u(i,j,k,nnew)=u(i,j,k,nnew)*oHz(i,k)
# endif
# ifdef DIAGNOSTICS_UV
            DO idiag=1,M3pgrd
              DiaU3wrk(i,j,k,idiag)=(DiaU3wrk(i,j,k,idiag)+             &
     &                               DC(i,0)*DiaRU(i,j,k,nrhs,idiag))*  &
     &                              oHz(i,k)
            END DO
#  if defined UV_VIS2 || defined UV_VIS4
            DiaU3wrk(i,j,k,M3xvis)=DiaU3wrk(i,j,k,M3xvis)*oHz(i,k)
            DiaU3wrk(i,j,k,M3yvis)=DiaU3wrk(i,j,k,M3yvis)*oHz(i,k)
            DiaU3wrk(i,j,k,M3hvis)=DiaU3wrk(i,j,k,M3hvis)*oHz(i,k)
#  endif
            DiaU3wrk(i,j,k,M3vvis)=DiaU3wrk(i,j,k,M3vvis)*oHz(i,k)
#  ifdef BODYFORCE
            DiaU3wrk(i,j,k,M3vvis)=DiaU3wrk(i,j,k,M3vvis)+              &
     &                             DC(i,0)*DiaRU(i,j,k,nrhs,M3vvis)*    &
     &                             oHz(i,k)
#  endif
            DiaU3wrk(i,j,k,M3rate)=DiaU3wrk(i,j,k,M3rate)*oHz(i,k)
# endif
          END DO
        END DO

# ifdef SPLINES_VVISC
!
!  Use conservative, parabolic spline reconstruction of vertical
!  viscosity derivatives.  Then, time step vertical viscosity term
!  implicitly by solving a tridiagonal system.
!
        cff1=1.0_r8/6.0_r8
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            FC(i,k)=cff1*Hzk(i,k  )-dt(ng)*AK(i,k-1)*oHz(i,k  )
            CF(i,k)=cff1*Hzk(i,k+1)-dt(ng)*AK(i,k+1)*oHz(i,k+1)
          END DO
        END DO
        DO i=IstrU,Iend
          CF(i,0)=0.0_r8
          DC(i,0)=0.0_r8
        END DO
!
!  LU decomposition and forward substitution.
!
        cff1=1.0_r8/3.0_r8
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            BC(i,k)=cff1*(Hzk(i,k)+Hzk(i,k+1))+                         &
     &              dt(ng)*AK(i,k)*(oHz(i,k)+oHz(i,k+1))
            cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
            CF(i,k)=cff*CF(i,k)
            DC(i,k)=cff*(u(i,j,k+1,nnew)-u(i,j,k,nnew)-                 &
     &                   FC(i,k)*DC(i,k-1))
          END DO
        END DO
!
!  Backward substitution.
!
        DO i=IstrU,Iend
          DC(i,N(ng))=0.0_r8
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrU,Iend
            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
          END DO
        END DO
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            DC(i,k)=DC(i,k)*AK(i,k)
            cff=dt(ng)*oHz(i,k)*(DC(i,k)-DC(i,k-1))
            u(i,j,k,nnew)=u(i,j,k,nnew)+cff
#  ifdef DIAGNOSTICS_UV
            DiaU3wrk(i,j,k,M3vvis)=DiaU3wrk(i,j,k,M3vvis)+cff
#  endif
          END DO
        END DO
# else
!
!  Compute off-diagonal coefficients [lambda*dt*Akv/Hz] for the
!  implicit vertical viscosity term at horizontal U-points and
!  vertical W-points.
!
        cff=-lambda*dt(ng)/0.5_r8
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            cff1=1.0_r8/(z_r(i,j,k+1)+z_r(i-1,j,k+1)-                   &
     &                   z_r(i,j,k  )-z_r(i-1,j,k  ))
            FC(i,k)=cff*cff1*AK(i,k)
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,0)=0.0_r8
          FC(i,N(ng))=0.0_r8
        END DO
!
!  Solve the tridiagonal system.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            DC(i,k)=u(i,j,k,nnew)
            BC(i,k)=Hzk(i,k)-FC(i,k)-FC(i,k-1)
          END DO
        END DO
        DO i=IstrU,Iend
          cff=1.0_r8/BC(i,1)
          CF(i,1)=cff*FC(i,1)
          DC(i,1)=cff*DC(i,1)
        END DO
        DO k=2,N(ng)-1
          DO i=IstrU,Iend
            cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
            CF(i,k)=cff*FC(i,k)
            DC(i,k)=cff*(DC(i,k)-FC(i,k-1)*DC(i,k-1))
          END DO
        END DO
!
!  Compute new solution by back substitution.
!
        DO i=IstrU,Iend
#  ifdef DIAGNOSTICS_UV
          wrk(i,N(ng))=u(i,j,N(ng),nnew)*oHz(i,N(ng))
#  endif
          DC(i,N(ng))=(DC(i,N(ng))-FC(i,N(ng)-1)*DC(i,N(ng)-1))/        &
     &                (BC(i,N(ng))-FC(i,N(ng)-1)*CF(i,N(ng)-1))
          u(i,j,N(ng),nnew)=DC(i,N(ng))
#  ifdef DIAGNOSTICS_UV
          DiaU3wrk(i,j,N(ng),M3vvis)=DiaU3wrk(i,j,N(ng),M3vvis)+        &
     &                               u(i,j,N(ng),nnew)-wrk(i,N(ng))
#  endif
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrU,Iend
#  ifdef DIAGNOSTICS_UV
            wrk(i,k)=u(i,j,k,nnew)*oHz(i,k)
#  endif
            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
            u(i,j,k,nnew)=DC(i,k)
#  ifdef DIAGNOSTICS_UV
            DiaU3wrk(i,j,k,M3vvis)=DiaU3wrk(i,j,k,M3vvis)+              &
     &                             u(i,j,k,nnew)-wrk(i,k)
#  endif
          END DO
        END DO
# endif
!
!  Replace INTERIOR POINTS incorrect vertical mean with more accurate
!  barotropic component, ubar=DU_avg1/(D*on_u). Recall that, D=CF(:,0).
!
        DO i=IstrU,Iend
          CF(i,0)=Hzk(i,1)
          DC(i,0)=u(i,j,1,nnew)*Hzk(i,1)
# ifdef WEC
          DCs(i,0)=u_stokes(i,j,1)*Hzk(i,1)
# endif
# ifdef DIAGNOSTICS_UV
          Dwrk(i,M2pgrd)=DiaU3wrk(i,j,1,M3pgrd)*Hzk(i,1)
          Dwrk(i,M2bstr)=DiaU3wrk(i,j,1,M3vvis)*Hzk(i,1)
#  ifdef UV_COR
          Dwrk(i,M2fcor)=DiaU3wrk(i,j,1,M3fcor)*Hzk(i,1)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
          Dwrk(i,M2xvis)=DiaU3wrk(i,j,1,M3xvis)*Hzk(i,1)
          Dwrk(i,M2yvis)=DiaU3wrk(i,j,1,M3yvis)*Hzk(i,1)
          Dwrk(i,M2hvis)=DiaU3wrk(i,j,1,M3hvis)*Hzk(i,1)
#  endif
#  ifdef UV_ADV
          Dwrk(i,M2xadv)=DiaU3wrk(i,j,1,M3xadv)*Hzk(i,1)
          Dwrk(i,M2yadv)=DiaU3wrk(i,j,1,M3yadv)*Hzk(i,1)
          Dwrk(i,M2hadv)=DiaU3wrk(i,j,1,M3hadv)*Hzk(i,1)
#  endif
#  ifdef WEC_MELLOR
          Dwrk(i,M2hrad)=DiaU3wrk(i,j,1,M3hrad)*Hzk(i,1)
#  endif
#  ifdef WEC_VF
          Dwrk(i,M2hjvf)=DiaU3wrk(i,j,1,M3hjvf)*Hzk(i,1)
          Dwrk(i,M2kvrf)=DiaU3wrk(i,j,1,M3kvrf)*Hzk(i,1)
#   ifdef UV_COR
          Dwrk(i,M2fsco)=DiaU3wrk(i,j,1,M3fsco)*Hzk(i,1)
#   endif
#   ifdef BOTTOM_STREAMING
          Dwrk(i,M2bstm)=DiaU3wrk(i,j,1,M3bstm)*Hzk(i,1)
#   endif
#   ifdef SURFACE_STREAMING
          Dwrk(i,M2sstm)=DiaU3wrk(i,j,1,M3sstm)*Hzk(i,1)
#   endif
          Dwrk(i,M2wrol)=DiaU3wrk(i,j,1,M3wrol)*Hzk(i,1)
          Dwrk(i,M2wbrk)=DiaU3wrk(i,j,1,M3wbrk)*Hzk(i,1)
#  endif
# endif
        END DO
        DO k=2,N(ng)
          DO i=IstrU,Iend
            CF(i,0)=CF(i,0)+Hzk(i,k)
            DC(i,0)=DC(i,0)+u(i,j,k,nnew)*Hzk(i,k)
# ifdef WEC
            DCs(i,0)=DCs(i,0)+u_stokes(i,j,k)*Hzk(i,k)
# endif
# ifdef DIAGNOSTICS_UV
            Dwrk(i,M2pgrd)=Dwrk(i,M2pgrd)+                              &
     &                     DiaU3wrk(i,j,k,M3pgrd)*Hzk(i,k)
            Dwrk(i,M2bstr)=Dwrk(i,M2bstr)+                              &
     &                     DiaU3wrk(i,j,k,M3vvis)*Hzk(i,k)
#  ifdef UV_COR
            Dwrk(i,M2fcor)=Dwrk(i,M2fcor)+                              &
     &                     DiaU3wrk(i,j,k,M3fcor)*Hzk(i,k)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
            Dwrk(i,M2xvis)=Dwrk(i,M2xvis)+                              &
     &                     DiaU3wrk(i,j,k,M3xvis)*Hzk(i,k)
            Dwrk(i,M2yvis)=Dwrk(i,M2yvis)+                              &
     &                     DiaU3wrk(i,j,k,M3yvis)*Hzk(i,k)
            Dwrk(i,M2hvis)=Dwrk(i,M2hvis)+                              &
     &                     DiaU3wrk(i,j,k,M3hvis)*Hzk(i,k)
#  endif
#  ifdef UV_ADV
            Dwrk(i,M2xadv)=Dwrk(i,M2xadv)+                              &
     &                     DiaU3wrk(i,j,k,M3xadv)*Hzk(i,k)
            Dwrk(i,M2yadv)=Dwrk(i,M2yadv)+                              &
     &                     DiaU3wrk(i,j,k,M3yadv)*Hzk(i,k)
            Dwrk(i,M2hadv)=Dwrk(i,M2hadv)+                              &
     &                     DiaU3wrk(i,j,k,M3hadv)*Hzk(i,k)
#  endif
#  ifdef WEC_MELLOR
            Dwrk(i,M2hrad)=Dwrk(i,M2hrad)+                              &
     &                     DiaU3wrk(i,j,k,M3hrad)*Hzk(i,k)
#  endif
#  ifdef WEC_VF
            Dwrk(i,M2hjvf)=Dwrk(i,M2hjvf)+                              &
     &                     DiaU3wrk(i,j,k,M3hjvf)*Hzk(i,k)
            Dwrk(i,M2kvrf)=Dwrk(i,M2kvrf)+                              &
     &                     DiaU3wrk(i,j,k,M3kvrf)*Hzk(i,k)
#   ifdef UV_COR
            Dwrk(i,M2fsco)=Dwrk(i,M2fsco)+                              &
     &                     DiaU3wrk(i,j,k,M3fsco)*Hzk(i,k)
#   endif
#   ifdef BOTTOM_STREAMING
            Dwrk(i,M2bstm)=Dwrk(i,M2bstm)+                              &
     &                     DiaU3wrk(i,j,k,M3bstm)*Hzk(i,k)
#   endif
#   ifdef SURFACE_STREAMING
            Dwrk(i,M2sstm)=Dwrk(i,M2sstm)+                              &
     &                     DiaU3wrk(i,j,k,M3sstm)*Hzk(i,k)
#   endif
            Dwrk(i,M2wrol)=Dwrk(i,M2wrol)+                              &
     &                     DiaU3wrk(i,j,k,M3wrol)*Hzk(i,k)
            Dwrk(i,M2wbrk)=Dwrk(i,M2wbrk)+                              &
     &                     DiaU3wrk(i,j,k,M3wbrk)*Hzk(i,k)
#  endif
# endif
          END DO
        END DO
        DO i=IstrU,Iend
          cff1=1.0_r8/(CF(i,0)*on_u(i,j))
          DC(i,0)=(DC(i,0)*on_u(i,j)-DU_avg1(i,j))*cff1      ! recursive
# ifdef WEC
          cff2=1.0_r8/CF(i,0)
          DCs(i,0)=DCs(i,0)*cff2-ubar_stokes(i,j)            ! recursive
# endif
# ifdef DIAGNOSTICS_UV
          DO idiag=1,M2pgrd
            Dwrk(i,idiag)=(Dwrk(i,idiag)*on_u(i,j)-                     &
     &                     DiaU2wrk(i,j,idiag))*cff1
          END DO
          Dwrk(i,M2bstr)=(Dwrk(i,M2bstr)*on_u(i,j)-                     &
     &                    DiaU2wrk(i,j,M2bstr)-DiaU2wrk(i,j,M2sstr))*   &
     &                   cff1
# endif
        END DO
!
!  Couple and update new solution.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            u(i,j,k,nnew)=u(i,j,k,nnew)-DC(i,0)
# ifdef MASKING
            u(i,j,k,nnew)=u(i,j,k,nnew)*umask(i,j)
# endif
# ifdef WET_DRY
            u(i,j,k,nnew)=u(i,j,k,nnew)*umask_wet(i,j)
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)*umask_wet(i,j)
# endif
# ifdef WEC
            u_stokes(i,j,k)=u_stokes(i,j,k)-DCs(i,0)
#  ifdef MASKING
            u_stokes(i,j,k)=u_stokes(i,j,k)*umask(i,j)
#  endif
#  ifdef WET_DRY
            u_stokes(i,j,k)=u_stokes(i,j,k)*umask_wet(i,j)
#  endif
# endif
# ifdef DIAGNOSTICS_UV
            DiaU3wrk(i,j,k,M3pgrd)=DiaU3wrk(i,j,k,M3pgrd)-              &
     &                             Dwrk(i,M2pgrd)
            DiaU3wrk(i,j,k,M3vvis)=DiaU3wrk(i,j,k,M3vvis)-              &
     &                             Dwrk(i,M2bstr)
#  ifdef UV_COR
            DiaU3wrk(i,j,k,M3fcor)=DiaU3wrk(i,j,k,M3fcor)-              &
     &                             Dwrk(i,M2fcor)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
            DiaU3wrk(i,j,k,M3xvis)=DiaU3wrk(i,j,k,M3xvis)-              &
     &                             Dwrk(i,M2xvis)
            DiaU3wrk(i,j,k,M3yvis)=DiaU3wrk(i,j,k,M3yvis)-              &
     &                             Dwrk(i,M2yvis)
            DiaU3wrk(i,j,k,M3hvis)=DiaU3wrk(i,j,k,M3hvis)-              &
     &                             Dwrk(i,M2hvis)
#  endif
#  ifdef UV_ADV
            DiaU3wrk(i,j,k,M3xadv)=DiaU3wrk(i,j,k,M3xadv)-              &
     &                             Dwrk(i,M2xadv)
            DiaU3wrk(i,j,k,M3yadv)=DiaU3wrk(i,j,k,M3yadv)-              &
     &                             Dwrk(i,M2yadv)
            DiaU3wrk(i,j,k,M3hadv)=DiaU3wrk(i,j,k,M3hadv)-              &
     &                             Dwrk(i,M2hadv)
#  endif
#  ifdef WEC_MELLOR
            DiaU3wrk(i,j,k,M3hrad)=DiaU3wrk(i,j,k,M3hrad)-              &
     &                             Dwrk(i,M2hrad)
#  endif
#  ifdef WEC_VF
            DiaU3wrk(i,j,k,M3hjvf)=DiaU3wrk(i,j,k,M3hjvf)-              &
     &                             Dwrk(i,M2hjvf)
            DiaU3wrk(i,j,k,M3kvrf)=DiaU3wrk(i,j,k,M3kvrf)-              &
     &                             Dwrk(i,M2kvrf)
#   ifdef UV_COR
            DiaU3wrk(i,j,k,M3fsco)=DiaU3wrk(i,j,k,M3fsco)-              &
     &                             Dwrk(i,M2fsco)
#   endif
#   ifdef BOTTOM_STREAMING
            DiaU3wrk(i,j,k,M3bstm)=DiaU3wrk(i,j,k,M3bstm)-              &
     &                             Dwrk(i,M2bstm)
#   endif
#   ifdef SURFACE_STREAMING
            DiaU3wrk(i,j,k,M3sstm)=DiaU3wrk(i,j,k,M3sstm)-              &
     &                             Dwrk(i,M2sstm)
#   endif
            DiaU3wrk(i,j,k,M3wrol)=DiaU3wrk(i,j,k,M3wrol)-              &
     &                             Dwrk(i,M2wrol)
            DiaU3wrk(i,j,k,M3wbrk)=DiaU3wrk(i,j,k,M3wbrk)-              &
     &                             Dwrk(i,M2wbrk)
#  endif
# endif
          END DO
        END DO

# if defined DIAGNOSTICS_UV && defined MASKING
        DO k=1,N(ng)
          DO i=IstrU,Iend
            DO idiag=1,NDM3d
              DiaU3wrk(i,j,k,idiag)=DiaU3wrk(i,j,k,idiag)*umask(i,j)
            END DO
          END DO
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Time step momentum equation in the ETA-direction.
!-----------------------------------------------------------------------
!
        IF (j.ge.JstrV) THEN
          DO i=Istr,Iend
            AK(i,0)=0.5_r8*(Akv(i,j-1,0)+                               &
     &                      Akv(i,j  ,0))
            DO k=1,N(ng)
              AK(i,k)=0.5_r8*(Akv(i,j-1,k)+                             &
     &                        Akv(i,j  ,k))
              Hzk(i,k)=0.5_r8*(Hz(i,j-1,k)+                             &
     &                         Hz(i,j  ,k))
# if defined SPLINES_VVISC || defined DIAGNOSTICS_UV
              oHz(i,k)=1.0_r8/Hzk(i,k)
# endif
            END DO
          END DO
!
!  Time step right-hand-side terms.
!
          IF (iic(ng).eq.ntfirst(ng)) THEN
            cff=0.25_r8*dt(ng)
          ELSE IF (iic(ng).eq.(ntfirst(ng)+1)) THEN
            cff=0.25_r8*dt(ng)*3.0_r8/2.0_r8
          ELSE
            cff=0.25_r8*dt(ng)*23.0_r8/12.0_r8
          END IF
          DO i=Istr,Iend
            DC(i,0)=cff*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              v(i,j,k,nnew)=v(i,j,k,nnew)+DC(i,0)*rv(i,j,k,nrhs)
# ifdef SPLINES_VVISC
              v(i,j,k,nnew)=v(i,j,k,nnew)*oHz(i,k)
# endif
# ifdef DIAGNOSTICS_UV
              DO idiag=1,M3pgrd
                DiaV3wrk(i,j,k,idiag)=(DiaV3wrk(i,j,k,idiag)+           &
     &                                 DC(i,0)*                         &
     &                                 DiaRV(i,j,k,nrhs,idiag))*        &
     &                                oHz(i,k)
              END DO
#  if defined UV_VIS2 || defined UV_VIS4
              DiaV3wrk(i,j,k,M3xvis)=DiaV3wrk(i,j,k,M3xvis)*oHz(i,k)
              DiaV3wrk(i,j,k,M3yvis)=DiaV3wrk(i,j,k,M3yvis)*oHz(i,k)
              DiaV3wrk(i,j,k,M3hvis)=DiaV3wrk(i,j,k,M3hvis)*oHz(i,k)
#  endif
              DiaV3wrk(i,j,k,M3vvis)=DiaV3wrk(i,j,k,M3vvis)*oHz(i,k)
#  ifdef BODYFORCE
              DiaV3wrk(i,j,k,M3vvis)=DiaV3wrk(i,j,k,M3vvis)+            &
     &                               DC(i,0)*DiaRV(i,j,k,nrhs,M3vvis)*  &
     &                               oHz(i,k)
#  endif
              DiaV3wrk(i,j,k,M3rate)=DiaV3wrk(i,j,k,M3rate)*oHz(i,k)
# endif
            END DO
          END DO

# ifdef SPLINES_VVISC
!
!  Use conservative, parabolic spline reconstruction of vertical
!  viscosity derivatives.  Then, time step vertical viscosity term
!  implicitly by solving a tridiagonal system.
!
          cff1=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=cff1*Hzk(i,k  )-dt(ng)*AK(i,k-1)*oHz(i,k  )
              CF(i,k)=cff1*Hzk(i,k+1)-dt(ng)*AK(i,k+1)*oHz(i,k+1)
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,0)=0.0_r8
            DC(i,0)=0.0_r8
          END DO
!
!  LU decomposition and forward substitution.
!
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              BC(i,k)=cff1*(Hzk(i,k)+Hzk(i,k+1))+                       &
     &                dt(ng)*AK(i,k)*(oHz(i,k)+oHz(i,k+1))
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
              CF(i,k)=cff*CF(i,k)
              DC(i,k)=cff*(v(i,j,k+1,nnew)-v(i,j,k,nnew)-               &
     &                     FC(i,k)*DC(i,k-1))
            END DO
          END DO
!
!  Backward substitution.
!
          DO i=Istr,Iend
            DC(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
            END DO
          END DO
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)*AK(i,k)
              cff=dt(ng)*oHz(i,k)*(DC(i,k)-DC(i,k-1))
              v(i,j,k,nnew)=v(i,j,k,nnew)+cff
#  ifdef DIAGNOSTICS_UV
              DiaV3wrk(i,j,k,M3vvis)=DiaV3wrk(i,j,k,M3vvis)+cff
#  endif
            END DO
          END DO
# else
!
!  Compute off-diagonal coefficients [lambda*dt*Akv/Hz] for the
!  implicit vertical viscosity term at horizontal V-points and
!  vertical W-points.
!
          cff=-lambda*dt(ng)/0.5_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff1=1.0_r8/(z_r(i,j,k+1)+z_r(i,j-1,k+1)-                 &
     &                     z_r(i,j,k  )-z_r(i,j-1,k  ))
              FC(i,k)=cff*cff1*AK(i,k)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            FC(i,N(ng))=0.0_r8
          END DO
!
!  Solve the tridiagonal system.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              DC(i,k)=v(i,j,k,nnew)
              BC(i,k)=Hzk(i,k)-FC(i,k)-FC(i,k-1)
            END DO
          END DO
          DO i=Istr,Iend
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FC(i,1)
            DC(i,1)=cff*DC(i,1)
          END DO
          DO k=2,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
              CF(i,k)=cff*FC(i,k)
              DC(i,k)=cff*(DC(i,k)-FC(i,k-1)*DC(i,k-1))
            END DO
          END DO
!
!  Compute new solution by back substitution.
!
          DO i=Istr,Iend
#  ifdef DIAGNOSTICS_UV
            wrk(i,N(ng))=v(i,j,N(ng),nnew)*oHz(i,N(ng))
#  endif
            DC(i,N(ng))=(DC(i,N(ng))-FC(i,N(ng)-1)*DC(i,N(ng)-1))/      &
     &                  (BC(i,N(ng))-FC(i,N(ng)-1)*CF(i,N(ng)-1))
            v(i,j,N(ng),nnew)=DC(i,N(ng))
#  ifdef DIAGNOSTICS_UV
            DiaV3wrk(i,j,N(ng),M3vvis)=DiaV3wrk(i,j,N(ng),M3vvis)+      &
     &                                 v(i,j,N(ng),nnew)-wrk(i,N(ng))
#  endif
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
#  ifdef DIAGNOSTICS_UV
              wrk(i,k)=v(i,j,k,nnew)*oHz(i,k)
#  endif
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
              v(i,j,k,nnew)=DC(i,k)
#  ifdef DIAGNOSTICS_UV
              DiaV3wrk(i,j,k,M3vvis)=DiaV3wrk(i,j,k,M3vvis)+            &
     &                               v(i,j,k,nnew)-wrk(i,k)
#  endif
            END DO
          END DO
# endif
!
!  Replace INTERIOR POINTS incorrect vertical mean with more accurate
!  barotropic component, vbar=DV_avg1/(D*om_v). Recall that, D=CF(:,0).
!
          DO i=Istr,Iend
            CF(i,0)=Hzk(i,1)
            DC(i,0)=v(i,j,1,nnew)*Hzk(i,1)
# ifdef WEC
            DCs(i,0)=v_stokes(i,j,1)*Hzk(i,1)
# endif
# ifdef DIAGNOSTICS_UV
            Dwrk(i,M2pgrd)=DiaV3wrk(i,j,1,M3pgrd)*Hzk(i,1)
            Dwrk(i,M2bstr)=DiaV3wrk(i,j,1,M3vvis)*Hzk(i,1)
#  ifdef UV_COR
            Dwrk(i,M2fcor)=DiaV3wrk(i,j,1,M3fcor)*Hzk(i,1)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
            Dwrk(i,M2xvis)=DiaV3wrk(i,j,1,M3xvis)*Hzk(i,1)
            Dwrk(i,M2yvis)=DiaV3wrk(i,j,1,M3yvis)*Hzk(i,1)
            Dwrk(i,M2hvis)=DiaV3wrk(i,j,1,M3hvis)*Hzk(i,1)
#  endif
#  ifdef UV_ADV
            Dwrk(i,M2xadv)=DiaV3wrk(i,j,1,M3xadv)*Hzk(i,1)
            Dwrk(i,M2yadv)=DiaV3wrk(i,j,1,M3yadv)*Hzk(i,1)
            Dwrk(i,M2hadv)=DiaV3wrk(i,j,1,M3hadv)*Hzk(i,1)
#  endif
#  ifdef WEC_MELLOR
            Dwrk(i,M2hrad)=DiaV3wrk(i,j,1,M3hrad)*Hzk(i,1)
#  endif
#  ifdef WEC_VF
            Dwrk(i,M2hjvf)=DiaV3wrk(i,j,1,M3hjvf)*Hzk(i,1)
            Dwrk(i,M2kvrf)=DiaV3wrk(i,j,1,M3kvrf)*Hzk(i,1)
#   ifdef UV_COR
            Dwrk(i,M2fsco)=DiaV3wrk(i,j,1,M3fsco)*Hzk(i,1)
#   endif
#   ifdef BOTTOM_STREAMING
            Dwrk(i,M2bstm)=DiaV3wrk(i,j,1,M3bstm)*Hzk(i,1)
#   endif
#   ifdef SURFACE_STREAMING
            Dwrk(i,M2sstm)=DiaV3wrk(i,j,1,M3sstm)*Hzk(i,1)
#   endif
            Dwrk(i,M2wrol)=DiaV3wrk(i,j,1,M3wrol)*Hzk(i,1)
            Dwrk(i,M2wbrk)=DiaV3wrk(i,j,1,M3wbrk)*Hzk(i,1)
#  endif
# endif
          END DO
          DO k=2,N(ng)
            DO i=Istr,Iend
              CF(i,0)=CF(i,0)+Hzk(i,k)
              DC(i,0)=DC(i,0)+v(i,j,k,nnew)*Hzk(i,k)
# ifdef WEC
              DCs(i,0)=DCs(i,0)+v_stokes(i,j,k)*Hzk(i,k)
# endif
# ifdef DIAGNOSTICS_UV
              Dwrk(i,M2pgrd)=Dwrk(i,M2pgrd)+                            &
     &                       DiaV3wrk(i,j,k,M3pgrd)*Hzk(i,k)
              Dwrk(i,M2bstr)=Dwrk(i,M2bstr)+                            &
     &                       DiaV3wrk(i,j,k,M3vvis)*Hzk(i,k)
#  ifdef UV_COR
              Dwrk(i,M2fcor)=Dwrk(i,M2fcor)+                            &
     &                       DiaV3wrk(i,j,k,M3fcor)*Hzk(i,k)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
              Dwrk(i,M2xvis)=Dwrk(i,M2xvis)+                            &
     &                       DiaV3wrk(i,j,k,M3xvis)*Hzk(i,k)
              Dwrk(i,M2yvis)=Dwrk(i,M2yvis)+                            &
     &                       DiaV3wrk(i,j,k,M3yvis)*Hzk(i,k)
              Dwrk(i,M2hvis)=Dwrk(i,M2hvis)+                            &
     &                       DiaV3wrk(i,j,k,M3hvis)*Hzk(i,k)
#  endif
#  ifdef UV_ADV
              Dwrk(i,M2xadv)=Dwrk(i,M2xadv)+                            &
     &                       DiaV3wrk(i,j,k,M3xadv)*Hzk(i,k)
              Dwrk(i,M2yadv)=Dwrk(i,M2yadv)+                            &
     &                       DiaV3wrk(i,j,k,M3yadv)*Hzk(i,k)
              Dwrk(i,M2hadv)=Dwrk(i,M2hadv)+                            &
     &                       DiaV3wrk(i,j,k,M3hadv)*Hzk(i,k)
#  endif
#  ifdef WEC_MELLOR
              Dwrk(i,M2hrad)=Dwrk(i,M2hrad)+                            &
     &                       DiaV3wrk(i,j,k,M3hrad)*Hzk(i,k)
#  endif
#  ifdef WEC_VF
              Dwrk(i,M2hjvf)=Dwrk(i,M2hjvf)+                            &
     &                       DiaV3wrk(i,j,k,M3hjvf)*Hzk(i,k)
              Dwrk(i,M2kvrf)=Dwrk(i,M2kvrf)+                            &
     &                       DiaV3wrk(i,j,k,M3kvrf)*Hzk(i,k)
#   ifdef UV_COR
              Dwrk(i,M2fsco)=Dwrk(i,M2fsco)+                            &
     &                       DiaV3wrk(i,j,k,M3fsco)*Hzk(i,k)
#   endif
#   ifdef BOTTOM_STREAMING
              Dwrk(i,M2bstm)=Dwrk(i,M2bstm)+                            &
     &                       DiaV3wrk(i,j,k,M3bstm)*Hzk(i,k)
#   endif
#   ifdef SURFACE_STREAMING
              Dwrk(i,M2sstm)=Dwrk(i,M2sstm)+                            &
     &                       DiaV3wrk(i,j,k,M3sstm)*Hzk(i,k)
#   endif
              Dwrk(i,M2wrol)=Dwrk(i,M2wrol)+                            &
     &                       DiaV3wrk(i,j,k,M3wrol)*Hzk(i,k)
              Dwrk(i,M2wbrk)=Dwrk(i,M2wbrk)+                            &
     &                       DiaV3wrk(i,j,k,M3wbrk)*Hzk(i,k)
#  endif
# endif
            END DO
          END DO
          DO i=Istr,Iend
            cff1=1.0_r8/(CF(i,0)*om_v(i,j))
            DC(i,0)=(DC(i,0)*om_v(i,j)-DV_avg1(i,j))*cff1    ! recursive
# ifdef WEC
            cff2=1.0_r8/CF(i,0)
            DCs(i,0)=DCs(i,0)*cff2-vbar_stokes(i,j)          ! recursive
# endif
# ifdef DIAGNOSTICS_UV
            DO idiag=1,M2pgrd
              Dwrk(i,idiag)=(Dwrk(i,idiag)*om_v(i,j)-                   &
     &                       DiaV2wrk(i,j,idiag))*cff1
            END DO
            Dwrk(i,M2bstr)=(Dwrk(i,M2bstr)*om_v(i,j)-                   &
     &                      DiaV2wrk(i,j,M2bstr)-DiaV2wrk(i,j,M2sstr))* &
     &                     cff1
# endif
          END DO
!
!  Couple and update new solution.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              v(i,j,k,nnew)=v(i,j,k,nnew)-DC(i,0)
# ifdef MASKING
              v(i,j,k,nnew)=v(i,j,k,nnew)*vmask(i,j)
# endif
# ifdef WET_DRY
              v(i,j,k,nnew)=v(i,j,k,nnew)*vmask_wet(i,j)
              rv(i,j,k,nrhs)=rv(i,j,k,nrhs)*vmask_wet(i,j)
# endif
# ifdef WEC
              v_stokes(i,j,k)=v_stokes(i,j,k)-DCs(i,0)
#  ifdef MASKING
              v_stokes(i,j,k)=v_stokes(i,j,k)*vmask(i,j)
#  endif
#  ifdef WET_DRY
              v_stokes(i,j,k)=v_stokes(i,j,k)*vmask_wet(i,j)
#  endif
# endif
# ifdef DIAGNOSTICS_UV
              DiaV3wrk(i,j,k,M3pgrd)=DiaV3wrk(i,j,k,M3pgrd)-            &
     &                               Dwrk(i,M2pgrd)
              DiaV3wrk(i,j,k,M3vvis)=DiaV3wrk(i,j,k,M3vvis)-            &
     &                               Dwrk(i,M2bstr)
#  ifdef UV_COR
              DiaV3wrk(i,j,k,M3fcor)=DiaV3wrk(i,j,k,M3fcor)-            &
     &                               Dwrk(i,M2fcor)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
              DiaV3wrk(i,j,k,M3xvis)=DiaV3wrk(i,j,k,M3xvis)-            &
     &                               Dwrk(i,M2xvis)
              DiaV3wrk(i,j,k,M3yvis)=DiaV3wrk(i,j,k,M3yvis)-            &
     &                               Dwrk(i,M2yvis)
              DiaV3wrk(i,j,k,M3hvis)=DiaV3wrk(i,j,k,M3hvis)-            &
     &                               Dwrk(i,M2hvis)
#  endif
#  ifdef UV_ADV
              DiaV3wrk(i,j,k,M3xadv)=DiaV3wrk(i,j,k,M3xadv)-            &
     &                               Dwrk(i,M2xadv)
              DiaV3wrk(i,j,k,M3yadv)=DiaV3wrk(i,j,k,M3yadv)-            &
     &                               Dwrk(i,M2yadv)
              DiaV3wrk(i,j,k,M3hadv)=DiaV3wrk(i,j,k,M3hadv)-            &
     &                               Dwrk(i,M2hadv)
#  endif
#  ifdef WEC_MELLOR
              DiaV3wrk(i,j,k,M3hrad)=DiaV3wrk(i,j,k,M3hrad)-            &
     &                               Dwrk(i,M2hrad)
#  endif
#  ifdef WEC_VF
              DiaV3wrk(i,j,k,M3hjvf)=DiaV3wrk(i,j,k,M3hjvf)-            &
     &                               Dwrk(i,M2hjvf)
              DiaV3wrk(i,j,k,M3kvrf)=DiaV3wrk(i,j,k,M3kvrf)-            &
     &                               Dwrk(i,M2kvrf)
#   ifdef UV_COR
              DiaV3wrk(i,j,k,M3fsco)=DiaV3wrk(i,j,k,M3fsco)-            &
     &                               Dwrk(i,M2fsco)
#   endif
#   ifdef BOTTOM_STREAMING
              DiaV3wrk(i,j,k,M3bstm)=DiaV3wrk(i,j,k,M3bstm)-            &
     &                               Dwrk(i,M2bstm)
#   endif
#   ifdef SURFACE_STREAMING
              DiaV3wrk(i,j,k,M3sstm)=DiaV3wrk(i,j,k,M3sstm)-            &
     &                               Dwrk(i,M2sstm)
#   endif
              DiaV3wrk(i,j,k,M3wrol)=DiaV3wrk(i,j,k,M3wrol)-            &
     &                               Dwrk(i,M2wrol)
              DiaV3wrk(i,j,k,M3wbrk)=DiaV3wrk(i,j,k,M3wbrk)-            &
     &                               Dwrk(i,M2wbrk)
#  endif
# endif
            END DO
          END DO
# if defined DIAGNOSTICS_UV && defined MASKING
          DO k=1,N(ng)
            DO idiag=1,NDM3d
              DO i=Istr,Iend
                DiaV3wrk(i,j,k,idiag)=DiaV3wrk(i,j,k,idiag)*vmask(i,j)
              END DO
            END DO
          END DO
# endif
        END IF
      END DO
# ifdef UV_WAVEDRAG
! Apply wave drag on tidal velocities only. filt_ubot is an estimate of
! the detided velocity. Apply only in the bottom 500 m. See Arbic et
! al., 2010.
      IF (iic(ng)-ntstart(ng) >= 86400) THEN
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff = 0.5*drag_scale*(wavedrag(i,j)+wavedrag(i-1,j))
            IF (cff > 0.0) THEN
              DO k=1,N(ng)
                IF (0.5*(z_r(i,j,k)+z_r(i-1,j,k)) <=                    &
     &             (0.5*(z_w(i,j,0)+z_w(i-1,j,0))+500.)) THEN
                  u(i,j,k,nnew)=u(i,j,k,nstp)-dt(ng)*cff*               &
     &                (u(i,j,k,nnew)-TIDES(ng)%filt_ubot(i,j))
                END IF
              END DO
            END IF
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff = 0.5*drag_scale*(wavedrag(i,j)+wavedrag(i,j-1))
            IF (cff > 0.0) THEN
              DO k=1,N(ng)
                IF (0.5*(z_r(i,j,k)+z_r(i,j-1,k)) <=                    &
     &             (0.5*(z_w(i,j,0)+z_w(i,j-1,0))+500.)) THEN
                  v(i,j,k,nnew)=v(i,j,k,nstp)-dt(ng)*cff*               &
     &                  (v(i,j,k,nnew)-TIDES(ng)%filt_vbot(i,j))
                END IF
              END DO
            END IF
          END DO
        END DO
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Set lateral boundary conditions.
!-----------------------------------------------------------------------
!
      CALL u3dbc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 nstp, nnew,                                      &
     &                 u)
      CALL v3dbc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 nstp, nnew,                                      &
     &                 v)
!
!-----------------------------------------------------------------------
!  Apply momentum transport point sources (like river runoff), if any.
!-----------------------------------------------------------------------
!
      IF (LuvSrc(ng)) THEN
        DO is=1,Nsrc(ng)
          i=SOURCES(ng)%Isrc(is)
          j=SOURCES(ng)%Jsrc(is)
          IF (((IstrR.le.i).and.(i.le.IendR)).and.                      &
     &        ((JstrR.le.j).and.(j.le.JendR))) THEN
            IF (INT(SOURCES(ng)%Dsrc(is)).eq.0) THEN
              DO k=1,N(ng)
                cff1=1.0_r8/(on_u(i,j)*                                 &
     &                       0.5_r8*(z_w(i-1,j,k)-z_w(i-1,j,k-1)+       &
     &                               z_w(i  ,j,k)-z_w(i  ,j,k-1)))
                u(i,j,k,nnew)=SOURCES(ng)%Qsrc(is,k)*cff1
              END DO
            ELSE
              DO k=1,N(ng)
                cff1=1.0_r8/(om_v(i,j)*                                 &
     &                       0.5_r8*(z_w(i,j-1,k)-z_w(i,j-1,k-1)+       &
     &                               z_w(i,j  ,k)-z_w(i,j  ,k-1)))
                v(i,j,k,nnew)=SOURCES(ng)%Qsrc(is,k)*cff1
              END DO
            END IF
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Couple 2D and 3D momentum equations.
!-----------------------------------------------------------------------
!
!  Couple velocity component in the XI-direction.
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          DC(i,0)=0.0_r8
          CF(i,0)=0.0_r8
# ifdef WEC
          CFs(i,0)=0.0_r8
# endif
          FC(i,0)=0.0_r8
        END DO
!
!  Compute thicknesses of U-boxes DC(i,1:N), total depth of the water
!  column DC(i,0), and incorrect vertical mean CF(i,0).  Notice that
!  barotropic component is replaced with its fast-time averaged
!  values.
!
        DO k=1,N(ng)
          DO i=IstrP,IendT
            cff=0.5_r8*on_u(i,j)
            DC(i,k)=cff*(Hz(i,j,k)+Hz(i-1,j,k))
            DC(i,0)=DC(i,0)+DC(i,k)
            CF(i,0)=CF(i,0)+                                            &
     &              DC(i,k)*u(i,j,k,nnew)
# ifdef WEC
            CFs(i,0)=CFs(i,0)+                                          &
     &               DC(i,k)*u_stokes(i,j,k)
# endif
          END DO
        END DO
        DO i=IstrP,IendT
          cff1=DC(i,0)                                  ! intermediate
# ifdef WEC
!          cff2=on_u(i,j)*DC(i,0)*ubar_stokes(i,j)
          cff2=DC(i,0)*ubar_stokes(i,j)
# endif
          DC(i,0)=1.0_r8/DC(i,0)                        ! recursive
          CF(i,0)=DC(i,0)*(CF(i,0)-DU_avg1(i,j))        ! recursive
# ifdef WEC
          CFs(i,0)=DC(i,0)*(CFs(i,0)-cff2)              ! recursive
# endif
          ubar(i,j,1)=DC(i,0)*DU_avg1(i,j)
# ifdef WET_DRY
          ubar(i,j,1)=ubar(i,j,1)*umask_wet(i,j)
# endif
          ubar(i,j,2)=ubar(i,j,1)
# ifdef DIAGNOSTICS_UV
          DiaU2wrk(i,j,M2rate)=ubar(i,j,1)-DiaU2int(i,j,M2rate)*DC(i,0)
          DiaU2int(i,j,M2rate)=ubar(i,j,1)*cff1
# endif
        END DO
# ifdef DIAGNOSTICS_UV
!
!  Convert the units of the fast-time integrated change in mass flux
!  diagnostic values to change in velocity (m s-1).
!
        DO idiag=1,NDM2d-1
          DO i=IstrP,IendT
            DiaU2wrk(i,j,idiag)=DC(i,0)*DiaU2wrk(i,j,idiag)
#  ifdef MASKING
            DiaU2wrk(i,j,idiag)=DiaU2wrk(i,j,idiag)*umask(i,j)
#  endif
          END DO
        END DO
# endif
!
!  Replace only BOUNDARY POINTS incorrect vertical mean with more
!  accurate barotropic component, ubar=DU_avg1/(D*on_u). Recall that,
!  D=CF(:,0).
!
!  NOTE:  Only the BOUNDARY POINTS need to be replaced. Avoid redundant
!         update in the interior again for computational purposes which
!         will not affect the nonlinear code.  However, the adjoint
!         code is wrong because the interior solution is corrected
!         twice. The replacement is avoided if the boundary edge is
!         periodic. The J-loop is pipelined, so we need to do a special
!         test for the southern and northern domain edges.
!
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO k=1,N(ng)
              u(Istr,j,k,nnew)=u(Istr,j,k,nnew)-CF(Istr,0)
# ifdef MASKING
              u(Istr,j,k,nnew)=u(Istr,j,k,nnew)*                        &
     &                         umask(Istr,j)
# endif
# ifdef WET_DRY
              u(Istr,j,k,nnew)=u(Istr,j,k,nnew)*                        &
     &                         umask_wet(Istr,j)
# endif
# ifdef WEC
              u_stokes(Istr,j,k)=u_stokes(Istr,j,k)-CFs(Istr,0)
#  ifdef MASKING
              u_stokes(Istr,j,k)=u_stokes(Istr,j,k)*                    &
     &                           umask(Istr,j)
#  endif
#  ifdef WET_DRY
              u_stokes(Istr,j,k)=u_stokes(Istr,j,k)*                    &
     &                           umask_wet(Istr,j)
#  endif
# endif
            END DO
          END IF
        END IF
!
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO k=1,N(ng)
              u(Iend+1,j,k,nnew)=u(Iend+1,j,k,nnew)-CF(Iend+1,0)
# ifdef MASKING
              u(Iend+1,j,k,nnew)=u(Iend+1,j,k,nnew)*                    &
     &                           umask(Iend+1,j)
# endif
# ifdef WET_DRY
              u(Iend+1,j,k,nnew)=u(Iend+1,j,k,nnew)*                    &
     &                           umask_wet(Iend+1,j)
# endif
# ifdef WEC
              u_stokes(Iend+1,j,k)=u_stokes(Iend+1,j,k)-CFs(Iend+1,0)
#  ifdef MASKING
              u_stokes(Iend+1,j,k)=u_stokes(Iend+1,j,k)*                &
     &                             umask(Iend+1,j)
#  endif
#  ifdef WET_DRY
              u_stokes(Iend+1,j,k)=u_stokes(Iend+1,j,k)*                &
     &                             umask_wet(Iend+1,j)
#  endif
# endif
            END DO
          END IF
        END IF
!
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (j.eq.0) THEN                           ! southern boundary
            DO k=1,N(ng)                             ! J-loop pipelined
              DO i=IstrU,Iend
                u(i,j,k,nnew)=u(i,j,k,nnew)-CF(i,0)
# ifdef MASKING
                u(i,j,k,nnew)=u(i,j,k,nnew)*                            &
     &                        umask(i,j)
# endif
# ifdef WET_DRY
                u(i,j,k,nnew)=u(i,j,k,nnew)*                            &
     &                        umask_wet(i,j)
# endif
# ifdef WEC
                u_stokes(i,j,k)=u_stokes(i,j,k)-CFs(i,0)
#  ifdef MASKING
                u_stokes(i,j,k)=u_stokes(i,j,k)*                        &
     &                          umask(i,j)
#  endif
#  ifdef WET_DRY
                u_stokes(i,j,k)=u_stokes(i,j,k)*                        &
     &                          umask_wet(i,j)
#  endif
# endif
              END DO
            END DO
          END IF
        END IF
!
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (j.eq.Mm(ng)+1) THEN                    ! northern boundary
            DO k=1,N(ng)                             ! J-loop pipelined
              DO i=IstrU,Iend
                u(i,j,k,nnew)=u(i,j,k,nnew)-CF(i,0)
# ifdef MASKING
                u(i,j,k,nnew)=u(i,j,k,nnew)*                            &
     &                        umask(i,j)
# endif
# ifdef WET_DRY
                u(i,j,k,nnew)=u(i,j,k,nnew)*                            &
     &                        umask_wet(i,j)
# endif
# ifdef WEC
                u_stokes(i,j,k)=u_stokes(i,j,k)-CFs(i,0)
#  ifdef MASKING
                u_stokes(i,j,k)=u_stokes(i,j,k)*                        &
     &                          umask(i,j)
#  endif
#  ifdef WET_DRY
                u_stokes(i,j,k)=u_stokes(i,j,k)*                        &
     &                          umask_wet(i,j)
#  endif
# endif
              END DO
            END DO
          END IF
        END IF
!
!  Compute correct mass flux, Hz*u/n.
!
        DO k=N(ng),1,-1
          DO i=IstrP,IendT
            Huon(i,j,k)=0.5_r8*(Huon(i,j,k)+u(i,j,k,nnew)*DC(i,k))
# ifdef WEC
            Huon(i,j,k)=Huon(i,j,k)+0.5_r8*u_stokes(i,j,k)*DC(i,k)
# endif
            FC(i,0)=FC(i,0)+Huon(i,j,k)
# ifdef DIAGNOSTICS_UV
            DiaU3wrk(i,j,k,M3rate)=u(i,j,k,nnew)-DiaU3wrk(i,j,k,M3rate)
# endif
          END DO
        END DO
        DO i=IstrP,IendT
          FC(i,0)=DC(i,0)*(FC(i,0)-DU_avg2(i,j))        ! recursive
        END DO
        DO k=1,N(ng)
          DO i=IstrP,IendT
            Huon(i,j,k)=Huon(i,j,k)-DC(i,k)*FC(i,0)
          END DO
        END DO
!
!  Couple velocity component in the ETA-direction.
!
        IF (j.ge.Jstr) THEN
          DO i=IstrT,IendT
            DC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
# ifdef WEC
            CFs(i,0)=0.0_r8
# endif
            FC(i,0)=0.0_r8
          END DO
!
!  Compute thicknesses of V-boxes DC(i,1:N), total depth of the water
!  column DC(i,0), and incorrect vertical mean CF(i,0).  Notice that
!  barotropic component is replaced with its fast-time averaged
!  values.
!
          DO k=1,N(ng)
            DO i=IstrT,IendT
              cff=0.5_r8*om_v(i,j)
              DC(i,k)=cff*(Hz(i,j,k)+Hz(i,j-1,k))
              DC(i,0)=DC(i,0)+DC(i,k)
              CF(i,0)=CF(i,0)+                                          &
     &                DC(i,k)*v(i,j,k,nnew)
# ifdef WEC
              CFs(i,0)=CFs(i,0)+                                        &
     &                 DC(i,k)*v_stokes(i,j,k)
# endif
            END DO
          END DO
          DO i=IstrT,IendT
            cff1=DC(i,0)                                 ! Intermediate
# ifdef WEC
!            cff2=om_v(i,j)*DC(i,0)*vbar_stokes(i,j)
            cff2=DC(i,0)*vbar_stokes(i,j)
# endif
            DC(i,0)=1.0_r8/DC(i,0)                       ! recursive
            CF(i,0)=DC(i,0)*(CF(i,0)-DV_avg1(i,j))       ! recursive
# ifdef WEC
            CFs(i,0)=DC(i,0)*(CFs(i,0)-cff2)             ! recursive
# endif
            vbar(i,j,1)=DC(i,0)*DV_avg1(i,j)
# ifdef WET_DRY
            vbar(i,j,1)=vbar(i,j,1)*vmask_wet(i,j)
# endif
            vbar(i,j,2)=vbar(i,j,1)
# ifdef DIAGNOSTICS_UV
            DiaV2wrk(i,j,M2rate)=vbar(i,j,1)-                           &
     &                           DiaV2int(i,j,M2rate)*DC(i,0)
            DiaV2int(i,j,M2rate)=vbar(i,j,1)*cff1
!!          DiaV2wrk(i,j,M2rate)=vbar(i,j,1)-DiaV2int(i,j,M2rate)
!!          DiaV2int(i,j,M2rate)=vbar(i,j,1)
# endif
          END DO
# ifdef DIAGNOSTICS_UV
!
!  Convert the units of the fast-time integrated change in mass flux
!  diagnostic values to change in velocity (m s-1).
!
          DO idiag=1,NDM2d-1
            DO i=IstrT,IendT
              DiaV2wrk(i,j,idiag)=DC(i,0)*DiaV2wrk(i,j,idiag)
#  ifdef MASKING
              DiaV2wrk(i,j,idiag)=DiaV2wrk(i,j,idiag)*vmask(i,j)
#  endif
            END DO
          END DO
# endif
!
!  Replace only BOUNDARY POINTS incorrect vertical mean with more
!  accurate barotropic component, vbar=DV_avg1/(D*om_v).  Recall that,
!  D=CF(:,0).
!
!  NOTE:  Only the BOUNDARY POINTS need to be replaced. Avoid redundant
!         update in the interior again for computational purposes which
!         will not affect the nonlinear code.  However, the adjoint
!         code is wrong because the interior solution is corrected
!         twice. The replacement is avoided if the boundary edge is
!         periodic. The J-loop is pipelined, so we need to do a special
!         test for the southern and northern domain edges.
!
          IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Western_Edge(tile)) THEN
              DO k=1,N(ng)
                v(Istr-1,j,k,nnew)=v(Istr-1,j,k,nnew)-CF(Istr-1,0)
# ifdef MASKING
                v(Istr-1,j,k,nnew)=v(Istr-1,j,k,nnew)*                  &
     &                             vmask(Istr-1,j)
# endif
# ifdef WET_DRY
                v(Istr-1,j,k,nnew)=v(Istr-1,j,k,nnew)*                  &
     &                             vmask_wet(Istr-1,j)
# endif
# ifdef WEC
                v_stokes(Istr-1,j,k)=v_stokes(Istr-1,j,k)-              &
     &                               CFs(Istr-1,0)
#  ifdef MASKING
                v_stokes(Istr-1,j,k)=v_stokes(Istr-1,j,k)*              &
     &                               vmask(Istr-1,j)
#  endif
#  ifdef WET_DRY
                v_stokes(Istr-1,j,k)=v_stokes(Istr-1,j,k)*              &
     &                               vmask_wet(Istr-1,j)
#  endif
# endif
              END DO
            END IF
          END IF
!
          IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
              DO k=1,N(ng)
                v(Iend+1,j,k,nnew)=v(Iend+1,j,k,nnew)-CF(Iend+1,0)
# ifdef MASKING
                v(Iend+1,j,k,nnew)=v(Iend+1,j,k,nnew)*                  &
     &                             vmask(Iend+1,j)
# endif
# ifdef WET_DRY
                v(Iend+1,j,k,nnew)=v(Iend+1,j,k,nnew)*                  &
     &                             vmask_wet(Iend+1,j)
# endif
# ifdef WEC
                v_stokes(Iend+1,j,k)=v_stokes(Iend+1,j,k)-              &
     &                               CFs(Iend+1,0)
#  ifdef MASKING
                v_stokes(Iend+1,j,k)=v_stokes(Iend+1,j,k)*              &
     &                               vmask(Iend+1,j)
#  endif
#  ifdef WET_DRY
                v_stokes(Iend+1,j,k)=v_stokes(Iend+1,j,k)*              &
     &                               vmask_wet(Iend+1,j)
#  endif
# endif
              END DO
            END IF
          END IF
!
          IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
            IF (j.eq.1) THEN                         ! southern boundary
              DO k=1,N(ng)                           ! J-loop pipelined
                DO i=Istr,Iend
                  v(i,j,k,nnew)=v(i,j,k,nnew)-CF(i,0)
# ifdef MASKING
                  v(i,j,k,nnew)=v(i,j,k,nnew)*                          &
     &                          vmask(i,j)
# endif
# ifdef WET_DRY
                  v(i,j,k,nnew)=v(i,j,k,nnew)*                          &
     &                          vmask_wet(i,j)
# endif
# ifdef WEC
                  v_stokes(i,j,k)=v_stokes(i,j,k)-CFs(i,0)
#  ifdef MASKING
                  v_stokes(i,j,k)=v_stokes(i,j,k)*                      &
     &                            vmask(i,j)
#  endif
#  ifdef WET_DRY
                  v_stokes(i,j,k)=v_stokes(i,j,k)*                      &
     &                            vmask_wet(i,j)
#  endif
# endif
                END DO
              END DO
            END IF
          END IF
!
          IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
            IF (j.eq.Mm(ng)+1) THEN                  ! northern boundary
              DO k=1,N(ng)                           ! J-loop pipelined
                DO i=Istr,Iend
                  v(i,j,k,nnew)=v(i,j,k,nnew)-CF(i,0)
# ifdef MASKING
                  v(i,j,k,nnew)=v(i,j,k,nnew)*                          &
     &                          vmask(i,j)
# endif
# ifdef WET_DRY
                  v(i,j,k,nnew)=v(i,j,k,nnew)*                          &
     &                          vmask_wet(i,j)
# endif
# ifdef WEC
                  v_stokes(i,j,k)=v_stokes(i,j,k)-CFs(i,0)
#  ifdef MASKING
                  v_stokes(i,j,k)=v_stokes(i,j,k)*                      &
     &                            vmask(i,j)
#  endif
#  ifdef WET_DRY
                  v_stokes(i,j,k)=v_stokes(i,j,k)*                      &
     &                            vmask_wet(i,j)
#  endif
# endif
                END DO
              END DO
            END IF
          END IF
!
!  Compute correct mass flux, Hz*v/m.
!
          DO k=N(ng),1,-1
            DO i=IstrT,IendT
              Hvom(i,j,k)=0.5_r8*(Hvom(i,j,k)+v(i,j,k,nnew)*DC(i,k))
# ifdef WEC
              Hvom(i,j,k)=Hvom(i,j,k)+0.5_r8*v_stokes(i,j,k)*DC(i,k)
# endif
              FC(i,0)=FC(i,0)+Hvom(i,j,k)
# ifdef DIAGNOSTICS_UV
              DiaV3wrk(i,j,k,M3rate)=v(i,j,k,nnew)-                     &
     &                               DiaV3wrk(i,j,k,M3rate)
# endif
            END DO
          END DO
          DO i=IstrT,IendT
            FC(i,0)=DC(i,0)*(FC(i,0)-DV_avg2(i,j))      ! recursive
          END DO
          DO k=1,N(ng)
            DO i=IstrT,IendT
              Hvom(i,j,k)=Hvom(i,j,k)-DC(i,k)*FC(i,0)
            END DO
          END DO
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          u(:,:,:,nnew))
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          v(:,:,:,nnew))
#  ifdef WEC
        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          u_stokes(:,:,:))
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          v_stokes(:,:,:))
#  endif
        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          Huon)
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          Hvom)
!
        DO k=1,2
          CALL exchange_u2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            ubar(:,:,k))
          CALL exchange_v2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            vbar(:,:,k))
        END DO
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    u(:,:,:,nnew), v(:,:,:,nnew),                 &
     &                    Huon, Hvom)
#   ifdef WEC
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    u_stokes(:,:,:), v_stokes(:,:,:))
#   endif
!
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    ubar(:,:,1), vbar(:,:,1),                     &
     &                    ubar(:,:,2), vbar(:,:,2))
# endif

      RETURN
      END SUBROUTINE step3d_uv_tile
#endif
      END MODULE step3d_uv_mod
