#include "cppdefs.h"
      MODULE metrics_mod
!
!svn $Id: metrics.F 795 2016-05-11 01:42:43Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine computes various horizontal metric terms.              !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: metrics

      CONTAINS
!
!***********************************************************************
      SUBROUTINE metrics (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
#if defined DIFF_3DCOEF || defined VISC_3DCOEF
      USE mod_mixing
#endif
      USE mod_ocean
#if defined SEDIMENT && defined SED_MORPH && defined SOLVE3D
      USE mod_sedbed
#endif
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
#include "tile.h"
!
      CALL metrics_tile (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp(ng), nnew(ng),                            &
     &                   GRID(ng) % f,                                  &
     &                   GRID(ng) % h,                                  &
     &                   GRID(ng) % pm,                                 &
     &                   GRID(ng) % pn,                                 &
#ifdef MASKING
     &                   GRID(ng) % pmask,                              &
     &                   GRID(ng) % rmask,                              &
#endif
     &                   GRID(ng) % angler,                             &
     &                   GRID(ng) % CosAngler,                          &
     &                   GRID(ng) % SinAngler,                          &
#ifdef SOLVE3D
# ifdef ICESHELF
     &                   GRID(ng) % zice,                               &
# endif
# if defined SEDIMENT && defined SED_MORPH
     &                   SEDBED(ng) % bed_thick,                        &
# endif
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
#endif
#if defined DIFF_GRID || defined DIFF_3DCOEF || \
    defined VISC_GRID || defined VISC_3DCOEF
     &                   GRID(ng) % grdscl,                             &
#endif
#if defined DIFF_3DCOEF
     &                   MIXING(ng) % Hdiffusion,                       &
#endif
#if defined VISC_3DCOEF
     &                   MIXING(ng) % Hviscosity,                       &
#endif
     &                   GRID(ng) % om_p,                               &
     &                   GRID(ng) % om_r,                               &
     &                   GRID(ng) % om_u,                               &
     &                   GRID(ng) % om_v,                               &
     &                   GRID(ng) % on_p,                               &
     &                   GRID(ng) % on_r,                               &
     &                   GRID(ng) % on_u,                               &
     &                   GRID(ng) % on_v,                               &
     &                   GRID(ng) % fomn,                               &
     &                   GRID(ng) % omn,                                &
     &                   GRID(ng) % pnom_p,                             &
     &                   GRID(ng) % pnom_r,                             &
     &                   GRID(ng) % pnom_u,                             &
     &                   GRID(ng) % pnom_v,                             &
     &                   GRID(ng) % pmon_p,                             &
     &                   GRID(ng) % pmon_r,                             &
     &                   GRID(ng) % pmon_u,                             &
     &                   GRID(ng) % pmon_v)
      RETURN
      END SUBROUTINE metrics
!
!***********************************************************************
      SUBROUTINE metrics_tile (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nstp, nnew,                              &
     &                         f, h, pm, pn,                            &
#ifdef MASKING
     &                         pmask, rmask,                            &
#endif
     &                         angler, CosAngler, SinAngler,            &
#ifdef SOLVE3D
# ifdef ICESHELF
     &                         zice,                                    &
# endif
# if defined SEDIMENT && defined SED_MORPH
     &                         bed_thick,                               &
# endif
     &                         Hz, z_r, z_w,                            &
#endif
#if defined DIFF_GRID || defined DIFF_3DCOEF || \
    defined VISC_GRID || defined VISC_3DCOEF
     &                         grdscl,                                  &
#endif
#if defined DIFF_3DCOEF
     &                         Hdiffusion,                              &
#endif
#if defined VISC_3DCOEF
     &                         Hviscosity,                              &
#endif
     &                         om_p, om_r, om_u, om_v,                  &
     &                         on_p, on_r, on_u, on_v,                  &
     &                         fomn, omn,                               &
     &                         pnom_p, pnom_r, pnom_u, pnom_v,          &
     &                         pmon_p, pmon_r, pmon_u, pmon_v)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
#ifdef FOUR_DVAR
      USE mod_fourdvar
#endif
      USE mod_iounits
      USE mod_ncparam
#ifdef NESTING
      USE mod_nesting
#endif
      USE mod_scalars
      USE mod_iounits
!
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_reduce
#endif
      USE exchange_2d_mod
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
#ifdef SOLVE3D
      USE set_depth_mod, ONLY : set_depth_tile
#endif
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew

#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: f(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)

      real(r8), intent(inout) :: h(LBi:,LBj:)
# ifdef MASKING
      real(r8), intent(inout) :: pmask(LBi:,LBj:)
      real(r8), intent(inout) :: rmask(LBi:,LBj:)
# endif
      real(r8), intent(inout) :: angler(LBi:,LBj:)
# ifdef SOLVE3D
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#  endif
#  if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in):: bed_thick(LBi:,LBj:,:)
#  endif
# endif
#if defined DIFF_GRID || defined DIFF_3DCOEF || \
    defined VISC_GRID || defined VISC_3DCOEF
      real(r8), intent(out) :: grdscl(LBi:,LBj:)
# endif
# if defined DIFF_3DCOEF
      real(r8), intent(out) :: Hdiffusion(LBi:,LBj:)
# endif
# if defined VISC_3DCOEF
      real(r8), intent(out) :: Hviscosity(LBi:,LBj:)
# endif
      real(r8), intent(out) :: om_p(LBi:,LBj:)
      real(r8), intent(out) :: om_r(LBi:,LBj:)
      real(r8), intent(out) :: om_u(LBi:,LBj:)
      real(r8), intent(out) :: om_v(LBi:,LBj:)
      real(r8), intent(out) :: on_p(LBi:,LBj:)
      real(r8), intent(out) :: on_r(LBi:,LBj:)
      real(r8), intent(out) :: on_u(LBi:,LBj:)
      real(r8), intent(out) :: on_v(LBi:,LBj:)
      real(r8), intent(out) :: fomn(LBi:,LBj:)
      real(r8), intent(out) :: omn(LBi:,LBj:)
      real(r8), intent(out) :: pnom_p(LBi:,LBj:)
      real(r8), intent(out) :: pnom_r(LBi:,LBj:)
      real(r8), intent(out) :: pnom_u(LBi:,LBj:)
      real(r8), intent(out) :: pnom_v(LBi:,LBj:)
      real(r8), intent(out) :: pmon_p(LBi:,LBj:)
      real(r8), intent(out) :: pmon_r(LBi:,LBj:)
      real(r8), intent(out) :: pmon_u(LBi:,LBj:)
      real(r8), intent(out) :: pmon_v(LBi:,LBj:)
      real(r8), intent(out) :: CosAngler(LBi:,LBj:)
      real(r8), intent(out) :: SinAngler(LBi:,LBj:)
# ifdef SOLVE3D
      real(r8), intent(out) :: Hz(LBi:,LBj:,:)
      real(r8), intent(out) :: z_r(LBi:,LBj:,:)
      real(r8), intent(out) :: z_w(LBi:,LBj:,0:)
# endif
#else
      real(r8), intent(in) :: f(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: h(LBi:UBi,LBj:UBj)
# ifdef MASKING
      real(r8), intent(inout) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(inout) :: angler(LBi:UBi,LBj:UBj)
# ifdef SOLVE3D
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#  endif
#  if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in):: bed_thick(LBi:UBi,LBj:UBj,2)
#  endif
# endif
#if defined DIFF_GRID || defined DIFF_3DCOEF || \
    defined VISC_GRID || defined VISC_3DCOEF
      real(r8), intent(out) :: grdscl(LBi:UBi,LBj:UBj)
# endif
# if defined DIFF_3DCOEF
      real(r8), intent(out) :: Hdiffusion(LBi:UBi,LBj:UBj)
# endif
# if defined VISC_3DCOEF
      real(r8), intent(out) :: Hviscosity(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(out) :: om_p(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: om_r(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: on_p(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: on_r(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: fomn(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: omn(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: pnom_p(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: pnom_r(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: pnom_u(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: pnom_v(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: pmon_p(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: pmon_r(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: pmon_v(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: CosAngler(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: SinAngler(LBi:UBi,LBj:UBj)
# ifdef SOLVE3D
      real(r8), intent(out) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
# endif
#endif
!
!  Local variable declarations.
!
      integer :: NSUB, i, ibry, is, j, k, rec
#ifdef NESTING
      integer :: cr, dg, ig, rg
#endif
      real(r8), parameter :: Large = 1.0E+20_r8
#if defined DIFF_3DCOEF || defined VISC_3DCOEF
      real(r8), parameter :: PecletCoef = 1.0_r8 / 12.0_r8
      real(r8), parameter :: Uscale = 0.1_r8
#endif
      real(r8) :: cff, cff1, cff2
      real(r8) :: my_DXmax, my_DXmin, my_DYmax, my_DYmin
#ifdef SOLVE3D
      real(r8) :: my_DZmax, my_DZmin
#endif
      real(r8) :: my_Cg_Cor, my_Cg_max, my_Cg_min, my_grdmax
#if defined DIFF_3DCOEF
      real(r8) :: my_DiffMax, my_DiffMin
#endif
#if defined VISC_3DCOEF
      real(r8) :: my_ViscMax, my_ViscMin
#endif
#if defined DIFF_3DCOEF || defined VISC_3DCOEF
      character (len=4) :: units
#endif
#if defined FOUR_DVAR
      character (len=50) :: Text
#endif

#ifdef DISTRIBUTE
      real(r8), dimension(14) :: buffer
      character (len=3), dimension(14) :: op_handle
#endif

#ifdef SOLVE3D
      real(r8), dimension(LBi:UBi,LBj:UBj) :: A2d
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute 1/m, 1/n, 1/mn, and f/mn at horizontal RHO-points.
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          om_r(i,j)=1.0_r8/pm(i,j)
          on_r(i,j)=1.0_r8/pn(i,j)
          omn(i,j)=1.0_r8/(pm(i,j)*pn(i,j))
          fomn(i,j)=f(i,j)*omn(i,j)
        END DO
      END DO
!
!  Exchange boundary information.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          om_r)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          on_r)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          omn)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          fomn)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    om_r, on_r, omn, fomn)
#endif
!
!-----------------------------------------------------------------------
!  Compute n/m, and m/n at horizontal RHO-points.
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          pnom_r(i,j)=pn(i,j)/pm(i,j)
          pmon_r(i,j)=pm(i,j)/pn(i,j)
        END DO
      END DO
!
!  Exchange boundary information.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pnom_r)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmon_r)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    pnom_r, pmon_r)
#endif
!
!-----------------------------------------------------------------------
!  Compute m/n, 1/m, and 1/n at horizontal U-points.
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          pmon_u(i,j)=(pm(i-1,j)+pm(i,j))/(pn(i-1,j)+pn(i,j))
          pnom_u(i,j)=(pn(i-1,j)+pn(i,j))/(pm(i-1,j)+pm(i,j))
          om_u(i,j)=2.0_r8/(pm(i-1,j)+pm(i,j))
          on_u(i,j)=2.0_r8/(pn(i-1,j)+pn(i,j))
        END DO
      END DO
!
!  Exchange boundary information.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmon_u)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pnom_u)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          om_u)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          on_u)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    pmon_u, pnom_u, om_u, on_u)
#endif
!
!-----------------------------------------------------------------------
!  Compute n/m, 1/m, and 1/m at horizontal V-points.
!-----------------------------------------------------------------------
!
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          pmon_v(i,j)=(pm(i,j-1)+pm(i,j))/(pn(i,j-1)+pn(i,j))
          pnom_v(i,j)=(pn(i,j-1)+pn(i,j))/(pm(i,j-1)+pm(i,j))
          om_v(i,j)=2.0_r8/(pm(i,j-1)+pm(i,j))
          on_v(i,j)=2.0_r8/(pn(i,j-1)+pn(i,j))
        END DO
      END DO
!
!  Exchange boundary information.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmon_v)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pnom_v)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          om_v)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          on_v)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    pmon_v, pnom_v, om_v, on_v)
#endif
!
!-----------------------------------------------------------------------
!  Compute n/m and m/n at horizontal PSI-points.
!-----------------------------------------------------------------------
!
      DO j=JstrP,JendT
        DO i=IstrP,IendT
          pnom_p(i,j)=(pn(i-1,j-1)+pn(i-1,j)+pn(i,j-1)+pn(i,j))/        &
     &                (pm(i-1,j-1)+pm(i-1,j)+pm(i,j-1)+pm(i,j))
          pmon_p(i,j)=(pm(i-1,j-1)+pm(i-1,j)+pm(i,j-1)+pm(i,j))/        &
     &                (pn(i-1,j-1)+pn(i-1,j)+pn(i,j-1)+pn(i,j))
          om_p(i,j)=4.0_r8/(pm(i-1,j-1)+pm(i-1,j)+pm(i,j-1)+pm(i,j))
          on_p(i,j)=4.0_r8/(pn(i-1,j-1)+pn(i-1,j)+pn(i,j-1)+pn(i,j))
        END DO
      END DO
!
!  Exchange boundary information.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pnom_p)
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmon_p)
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          om_p)
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          on_p)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    pnom_p, pmon_p, om_p, on_p)
#endif

#ifdef MASKING
!
!-----------------------------------------------------------------------
!  Set slipperiness (no-slip) mask at PSI-points.
!-----------------------------------------------------------------------
!
! Set no-slip boundary conditions on land-mask boundaries regardless of
! supplied value of gamma2.
!
      cff1=1.0_r8       ! computation of off-diagonal nonlinear terms
      cff2=2.0_r8
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          IF ((rmask(i-1,j  ).gt.0.5_r8).and.                           &
     &        (rmask(i  ,j  ).gt.0.5_r8).and.                           &
     &        (rmask(i-1,j-1).gt.0.5_r8).and.                           &
     &        (rmask(i  ,j-1).gt.0.5_r8)) THEN
            pmask(i,j)=1.0_r8
          ELSE IF ((rmask(i-1,j  ).lt.0.5_r8).and.                      &
     &             (rmask(i  ,j  ).gt.0.5_r8).and.                      &
     &             (rmask(i-1,j-1).gt.0.5_r8).and.                      &
     &             (rmask(i  ,j-1).gt.0.5_r8)) THEN
            pmask(i,j)=cff1
          ELSE IF ((rmask(i-1,j  ).gt.0.5_r8).and.                      &
     &             (rmask(i  ,j  ).lt.0.5_r8).and.                      &
     &             (rmask(i-1,j-1).gt.0.5_r8).and.                      &
     &             (rmask(i  ,j-1).gt.0.5_r8)) THEN
            pmask(i,j)=cff1
          ELSE IF ((rmask(i-1,j  ).gt.0.5_r8).and.                      &
     &             (rmask(i  ,j  ).gt.0.5_r8).and.                      &
     &             (rmask(i-1,j-1).lt.0.5_r8).and.                      &
     &             (rmask(i  ,j-1).gt.0.5_r8)) THEN
            pmask(i,j)=cff1
          ELSE IF ((rmask(i-1,j  ).gt.0.5_r8).and.                      &
     &             (rmask(i  ,j  ).gt.0.5_r8).and.                      &
     &             (rmask(i-1,j-1).gt.0.5_r8).and.                      &
     &             (rmask(i  ,j-1).lt.0.5_r8)) THEN
            pmask(i,j)=cff1
          ELSE IF ((rmask(i-1,j  ).gt.0.5_r8).and.                      &
     &             (rmask(i  ,j  ).lt.0.5_r8).and.                      &
     &             (rmask(i-1,j-1).gt.0.5_r8).and.                      &
     &             (rmask(i  ,j-1).lt.0.5_r8)) THEN
            pmask(i,j)=cff2
          ELSE IF ((rmask(i-1,j  ).lt.0.5_r8).and.                      &
     &             (rmask(i  ,j  ).gt.0.5_r8).and.                      &
     &             (rmask(i-1,j-1).lt.0.5_r8).and.                      &
     &             (rmask(i  ,j-1).gt.0.5_r8)) THEN
            pmask(i,j)=cff2
          ELSE IF ((rmask(i-1,j  ).gt.0.5_r8).and.                      &
     &             (rmask(i  ,j  ).gt.0.5_r8).and.                      &
     &             (rmask(i-1,j-1).lt.0.5_r8).and.                      &
     &             (rmask(i  ,j-1).lt.0.5_r8)) THEN
            pmask(i,j)=cff2
          ELSE IF ((rmask(i-1,j  ).lt.0.5_r8).and.                      &
     &             (rmask(i  ,j  ).lt.0.5_r8).and.                      &
     &             (rmask(i-1,j-1).gt.0.5_r8).and.                      &
     &             (rmask(i  ,j-1).gt.0.5_r8)) THEN
            pmask(i,j)=cff2
          ELSE
            pmask(i,j)=0.0_r8
          END IF
        END DO
      END DO
!
!  Exchange boundary information.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmask)
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    pmask)
# endif
#endif
!
!-----------------------------------------------------------------------
! Compute cosine and sine of grid rotation angle.
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          CosAngler(i,j)=COS(angler(i,j))
          SinAngler(i,j)=SIN(angler(i,j))
        END DO
      END DO
!
!  Exchange boundary information.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CosAngler)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          SinAngler)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    CosAngler, SinAngler)
#endif

#if defined DIFF_GRID || defined DIFF_3DCOEF || \
    defined VISC_GRID || defined VISC_3DCOEF
!
!-----------------------------------------------------------------------
! Determine the squared root of the area of each grid cell used to
! rescale horizontal mixing by the grid size.
!-----------------------------------------------------------------------
!
      cff=0.0_r8
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          grdscl(i,j)=SQRT(om_r(i,j)*on_r(i,j))
        END DO
      END DO
!
!  Exchange boundary information.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          grdscl)
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    grdscl)
# endif
#endif

#if defined DIFF_3DCOEF || defined VISC_3DCOEF
!
!-----------------------------------------------------------------------
! Compute time invariant, horizontal mixing coefficient using grid
! scale. Following Holland et (1998), Webb et al. (1998), Griffies
! and Hallberg (2000), and Lee et al. (2002), the horizontal mixing
! coefficient can be estimated as:
!
!   Hmixing = 1/12 * Uscale * grdscl           (Harmonic)
!   Bmixing = 1/12 * Uscale * grdscl**3        (Biharmonic)
!-----------------------------------------------------------------------
!
# if defined DIFF_3DCOEF
      my_DiffMin= Large
      my_DiffMax=-Large
# endif
# if defined VISC_3DCOEF
      my_ViscMin= Large
      my_ViscMax=-Large
# endif
      DO j=JstrT,JendT
        DO i=IstrT,IendT
# if defined DIFF_3DCOEF
#  if defined TS_DIF2
          Hdiffusion(i,j)=PecletCoef*Uscale*grdscl(i,j)
#  elif defined TS_DIF4
          Hdiffusion(i,j)=PecletCoef*Uscale*grdscl(i,j)**3
#  endif
          my_DiffMin=MIN(my_DiffMin, Hdiffusion(i,j))
          my_DiffMax=MAX(my_DiffMax, Hdiffusion(i,j))
# endif
# if defined VISC_3DCOEF
#  if defined UV_VIS2
          Hviscosity(i,j)=PecletCoef*Uscale*grdscl(i,j)
#  elif defined UV_VIS4
          Hviscosity(i,j)=PecletCoef*Uscale*grdscl(i,j)**3
#  endif
          my_ViscMin=MIN(my_ViscMin, Hviscosity(i,j))
          my_ViscMax=MAX(my_ViscMax, Hviscosity(i,j))
# endif
        END DO
      END DO
!
!  Exchange boundary information.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
# ifdef DIFF_3DCOEF
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Hdiffusion)
# endif
# ifdef VISC_3DCOEF
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Hviscosity)
# endif
      END IF

# ifdef DISTRIBUTE
#  ifdef DIFF_3DCOEF
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Hdiffusion)
#  endif
#  ifdef VISC_3DCOEF
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Hviscosity)
#  endif
# endif
#endif
!
!-----------------------------------------------------------------------
!  Compute minimum and maximum grid spacing.
!-----------------------------------------------------------------------
#ifdef SOLVE3D
!
!  Compute time invariant depths (use zero free-surface).
!
      DO i=LBi,UBi
        DO j=LBj,UBj
          A2d(i,j)=0.0_r8
        END DO
      END DO

      CALL set_depth_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp, nnew,                                  &
     &                     h,                                           &
# ifdef ICESHELF
     &                     zice,                                        &
# endif
# if defined SEDIMENT && defined SED_MORPH
     &                     bed_thick,                                   &
# endif
     &                     A2d,                                         &
     &                     Hz, z_r, z_w)
#endif
!
!  Compute grid spacing range.
!
      my_DXmin= Large
      my_DXmax=-Large
      my_DYmin= Large
      my_DYmax=-Large
#ifdef SOLVE3D
      my_DZmin= Large
      my_DZmax=-Large
#endif
      my_grdmax=-Large
      DO j=JstrT,JendT
        DO i=IstrT,IendT
#if defined VISC_GRID || defined DIFF_GRID
          cff=grdscl(i,j)
#else
          cff=SQRT(om_r(i,j)*on_r(i,j))
#endif
          my_DXmin=MIN(my_DXmin,om_r(i,j))
          my_DXmax=MAX(my_DXmax,om_r(i,j))
          my_DYmin=MIN(my_DYmin,on_r(i,j))
          my_DYmax=MAX(my_DYmax,on_r(i,j))
          my_grdmax=MAX(my_grdmax,cff)
#ifdef SOLVE3D
          DO k=1,N(ng)
            my_DZmin=MIN(my_DZmin,Hz(i,j,k))
            my_DZmax=MAX(my_DZmax,Hz(i,j,k))
          END DO
#endif
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute gravity waves Courant number.
!-----------------------------------------------------------------------
!
!  The 2D Courant number is defined as:
!
!     Cg = c * dt * SQRT (1/dx^2 + 1/dy^2)
!
!  where c=SQRT(g*h) is gravity wave speed, and dx, dy are grid spacing
!  in each direction.
!
      my_Cg_min= Large
      my_Cg_max=-Large
      my_Cg_Cor=-Large

      DO j=JstrT,JendT
        DO i=IstrT,IendT
#ifdef MASKING
          IF (rmask(i,j).gt.0.0_r8) THEN
#endif
#ifdef ICESHELF
            cff=dtfast(ng)*                                             &
     &          SQRT(g*ABS(h(i,j)-ABS(zice(i,j)))*                      &
     &              (pm(i,j)*pm(i,j)+pn(i,j)*pn(i,j)))
#else
            cff=dtfast(ng)*                                             &
     &          SQRT(g*ABS(h(i,j))*(pm(i,j)*pm(i,j)+pn(i,j)*pn(i,j)))
#endif
            my_Cg_min=MIN(my_Cg_min,cff)
            my_Cg_max=MAX(my_Cg_max,cff)

            cff=dt(ng)*ABS(f(i,j))
            my_Cg_Cor=MAX(my_Cg_Cor,cff)

#ifdef MASKING
          END IF
#endif
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Perform global reductions.
!-----------------------------------------------------------------------
!
#ifdef DISTRIBUTE
      NSUB=1                             ! distributed-memory
#else
      IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                        &
     &    DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
#endif
!$OMP CRITICAL (REDUCTIONS)
      IF (tile_count.eq.0) THEN
        Cg_min=my_Cg_min
        Cg_max=my_Cg_max
        Cg_Cor=my_Cg_Cor
        grdmax(ng)=my_grdmax
        DXmin(ng)=my_DXmin
        DXmax(ng)=my_DXmax
        DYmin(ng)=my_DYmin
        DYmax(ng)=my_DYmax
#ifdef SOLVE3D
        DZmin(ng)=my_DZmin
        DZmax(ng)=my_DZmax
#endif
#ifdef DIFF_3DCOEF
        DiffMin(ng)=my_DiffMin
        DiffMax(ng)=my_DiffMax
#endif
#ifdef VISC_3DCOEF
        ViscMin(ng)=my_ViscMin
        ViscMax(ng)=my_ViscMax
#endif
      ELSE
        Cg_min=MIN(Cg_min,my_Cg_min)
        Cg_max=MAX(Cg_max,my_Cg_max)
        Cg_Cor=MAX(Cg_Cor,my_Cg_Cor)
        grdmax(ng)=MAX(grdmax(ng),my_grdmax)
        DXmin(ng)=MIN(DXmin(ng),my_DXmin)
        DXmax(ng)=MAX(DXmax(ng),my_DXmax)
        DYmin(ng)=MIN(DYmin(ng),my_DYmin)
        DYmax(ng)=MAX(DYmax(ng),my_DYmax)
#ifdef SOLVE3D
        DZmin(ng)=MIN(DZmin(ng),my_DZmin)
        DZmax(ng)=MAX(DZmax(ng),my_DZmax)
#endif
#ifdef DIFF_3DCOEF
        DiffMin(ng)=MIN(DiffMin(ng),my_DiffMin)
        DiffMax(ng)=MAX(DiffMax(ng),my_DiffMax)
#endif
#ifdef VISC_3DCOEF
        ViscMin(ng)=MIN(ViscMin(ng),my_ViscMin)
        ViscMax(ng)=MAX(ViscMax(ng),my_ViscMax)
#endif
      END IF
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
#ifdef DISTRIBUTE
        buffer(1)=Cg_min
        op_handle(1)='MIN'
        buffer(2)=Cg_max
        op_handle(2)='MAX'
        buffer(3)=Cg_Cor
        op_handle(3)='MAX'
        buffer(4)=grdmax(ng)
        op_handle(4)='MAX'
        buffer(5)=DXmin(ng)
        op_handle(5)='MIN'
        buffer(6)=DXmax(ng)
        op_handle(6)='MAX'
        buffer(7)=DYmin(ng)
        op_handle(7)='MIN'
        buffer(8)=DYmax(ng)
        op_handle(8)='MAX'
# ifdef SOLVE3D
        buffer(9)=DZmin(ng)
        op_handle(9)='MIN'
        buffer(10)=DZmax(ng)
        op_handle(10)='MAX'
# else
        buffer(9)=0.0_r8
        op_handle(9)='MIN'
        buffer(10)=0.0_r8
        op_handle(10)='MAX'
# endif
# ifdef DIFF_3DCOEF
        buffer(11)=DiffMin(ng)
        op_handle(11)='MIN'
        buffer(12)=DiffMax(ng)
        op_handle(12)='MAX'
# else
        buffer(11)=0.0_r8
        op_handle(11)='MIN'
        buffer(12)=0.0_r8
        op_handle(12)='MAX'
# endif
# ifdef VISC_3DCOEF
        buffer(13)=ViscMin(ng)
        op_handle(13)='MIN'
        buffer(14)=ViscMax(ng)
        op_handle(14)='MAX'
# else
        buffer(13)=0.0_r8
        op_handle(13)='MIN'
        buffer(14)=0.0_r8
        op_handle(14)='MAX'
# endif
        CALL mp_reduce (ng, model, 14, buffer, op_handle)
        Cg_min=buffer(1)
        Cg_max=buffer(2)
        Cg_Cor=buffer(3)
        grdmax(ng)=buffer(4)
        DXmin(ng)=buffer(5)
        DXmax(ng)=buffer(6)
        DYmin(ng)=buffer(7)
        DYmax(ng)=buffer(8)
# ifdef SOLVE3D
        DZmin(ng)=buffer(9)
        DZmax(ng)=buffer(10)
# endif
# ifdef DIFF_3DCOEF
        DiffMin(ng)=buffer(11)
        DiffMax(ng)=buffer(12)
# endif
# ifdef VISC_3DCOEF
        ViscMin(ng)=buffer(13)
        ViscMax(ng)=buffer(14)
# endif
#endif
        IF (Master.and.LwrtInfo(ng)) THEN
          WRITE(stdout,10) ng,                                          &
     &                     DXmin(ng)/1000.0_r8, DXmax(ng)/1000.0_r8,    &
     &                     DYmin(ng)/1000.0_r8, DYmax(ng)/1000.0_r8
  10      FORMAT (/,' Metrics information for Grid ',i2.2,':',          &
     &            /,' ===============================',/,               &
     &            /,' Minimum X-grid spacing, DXmin = ',1pe15.8,' km',  &
     &            /,' Maximum X-grid spacing, DXmax = ',1pe15.8,' km',  &
     &            /,' Minimum Y-grid spacing, DYmin = ',1pe15.8,' km',  &
     &            /,' Maximum Y-grid spacing, DYmax = ',1pe15.8,' km')
#ifdef SOLVE3D
          WRITE(stdout,20) DZmin(ng), DZmax(ng)
  20      FORMAT (' Minimum Z-grid spacing, DZmin = ',1pe15.8,' m',/,   &
     &            ' Maximum Z-grid spacing, DZmax = ',1pe15.8,' m')
#endif
          WRITE (stdout,30) Cg_min, Cg_max, Cg_Cor
  30      FORMAT (/,' Minimum barotropic Courant Number = ', 1pe15.8,/, &
     &              ' Maximum barotropic Courant Number = ', 1pe15.8,/, &
     &              ' Maximum Coriolis   Courant Number = ', 1pe15.8,/)
# if defined VISC_GRID || defined DIFF_GRID
          WRITE (stdout,40) grdmax(ng)/1000.0_r8
  40      FORMAT (' Horizontal mixing scaled by grid size,',            &
     &            ' GRDMAX = ',1pe15.8,' km')
# endif
# ifdef DIFF_3DCOEF
#  ifdef TS_DIF2
          units='m2/s'
#  elif defined TS_DIF4
          units='m4/s'
#  endif
          WRITE (stdout,50) DiffMin(ng), TRIM(units),                   &
     &                      DiffMax(ng), TRIM(units)
  50      FORMAT (/,' Minimum horizontal diffusion coefficient = ',     &
     &             1pe15.8,1x,a,                                        &
     &            /,' Maximum horizontal diffusion coefficient = ',     &
     &             1pe15.8,1x,a)
# endif
# ifdef VISC_3DCOEF
#  ifdef UV_VIS2
          units='m2/s'
#  elif defined UV_VIS4
          units='m4/s'
#  endif
          WRITE (stdout,60) ViscMin(ng), TRIM(units),                   &
     &                      ViscMax(ng), TRIM(units)
  60      FORMAT (/,' Minimum horizontal viscosity coefficient = ',     &
     &             1pe15.8,1x,a,                                        &
     &            /,' Maximum horizontal viscosity coefficient = ',     &
     &             1pe15.8,1x,a)
# endif
        END IF
      END IF
!$OMP END CRITICAL (REDUCTIONS)

#ifdef NESTING
!
!-----------------------------------------------------------------------
!  If grid refinement, set various important parameters.
!-----------------------------------------------------------------------
!
      IF (ng.eq.Ngrids) THEN
!
!  Set coarser donor grid number to finer grid external contact points.
!  The reciver grid is always finer than the donor grid, DXmax(dg) >
!  DXmax(rg). It is done in terms of the reciver grid.
!
        DO cr=1,Ncontact
          dg=donor_grid(cr)
          rg=receiver_grid(cr)
          IF (RefinedGrid(rg).and.                                        &
     &        (RefineScale(rg).gt.0).and.                                 &
     &         DXmax(dg).gt.DXmax(rg)) THEN
            CoarserDonor(rg)=dg
          END IF
        END DO
!
!  Set donor finer grid number to receiver coarser grid. This is use for
!  two-way feeback between grids. The donor grid is always finer than the
!  receiver grid, DXmax(dg) < DXmax(rg). It is done in terms of the donor
!  grid.
!
        DO cr=1,Ncontact
          dg=donor_grid(cr)
          rg=receiver_grid(cr)
          IF (RefinedGrid(dg).and.                                        &
     &        (RefineScale(dg).gt.0).and.                                 &
     &        DXmax(dg).gt.DXmax(rg)) THEN
            FinerDonor(dg)=rg
          END IF
        END DO
!
!  Set logical indicating which coarser grid is a donor to a finer
!  grid external contact points.
!
        DO dg=1,Ngrids
          DO ig=1,Ngrids
            IF (CoarserDonor(ig).eq.dg) THEN
              DonorToFiner(dg)=.TRUE.
            END IF
          END DO
        END DO
!
!  Set switch indicating which refined grid(s), with RefineScale(dg) > 0,
!  include finer refined grids inside (telescoping refinement).
!
        DO cr=1,Ncontact
          dg=donor_grid(cr)
          rg=receiver_grid(cr)
          IF ((RefineScale(dg).gt.0).and.(CoarserDonor(rg).eq.dg)) THEN
            Telescoping(dg)=.TRUE.
          END IF
        END DO
!
!  Set Number of refined grid time-steps using the ratio between donor
!  and receiver grids.  Optimally, thes values should be based on the
!  spatial refinement ratio for stability. However, the user is allowed
!  to take larger divisible time-step with respect to the donor grid.
!  The user is responsible to set the appropriate refined grid time-step
!  for stability.
!
        DO cr=1,Ncontact
          dg=donor_grid(cr)
          rg=receiver_grid(cr)
          IF (RefinedGrid(rg).and.(CoarserDonor(rg).eq.dg)) THEN
            RefineSteps(rg)=NINT(dt(dg)/dt(rg))
          END IF
        END DO
!
!  Report information.
!
        IF (DOMAIN(ng)%NorthEast_Test(tile)) THEN
          IF (Master.and.ANY(RefinedGrid)) THEN
            WRITE (stdout,70)
  70        FORMAT (/,' Refined Nested Grid(s) Information: ',          &
     &              /,' ==================================',/,/,        &
     &              ' Refined    Donor   Refined   Timestep   Refined', &
     &              /,'  Grid      Grid    Scale      Ratio   ',        &
     &              'Timesteps',/)
            DO ig=1,Ngrids
              IF (RefineScale(ig).gt.0) THEN
                dg=CoarserDonor(ig)
                WRITE (stdout,80) ig, dg, RefineScale(ig),              &
     &                            dt(dg)/dt(ig), RefineSteps(ig)
  80            FORMAT(4x,i2.2,8x,i2.2,7x,i2.2,4x,f9.5,6x,i2.2)
              END IF
            END DO
            WRITE (stdout,90)
  90        FORMAT (/,' WARNING: Usually the number of Refined ',       &
     &              'Timesteps must be the same',/,10x,                 &
     &              'as the Refined Scale for numerical stability.',/)
          END IF
        END IF
!
!  Check refined grid time-step.  The time-step size for refined grids
!  needs to be an exact mode multiple of its coarser donor grid.  In
!  principle, it can be a "RefineScale" factor smaller.  However, other
!  integer smaller or larger factor are allowed such that:
!
!          MOD(dt(dg), dt(rg)) = 0        dg:  donor  coarse grid
!                                         rg:  receiver fine grid
!
        DO ig=1,Ngrids
          IF (RefinedGrid(ig).and.(RefineScale(ig).gt.0)) THEN
            dg=CoarserDonor(ig)
            IF (MOD(dt(dg),dt(ig)).ne.0.0d0) THEN
              IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
                IF (Master) THEN
                  WRITE (stdout,100) ig, dt(ig), dg, dt(dg),            &
     &                               MOD(dt(dg),dt(ig))
  100             FORMAT (/,' METRICS - illegal timestep size for ',    &
     &                    ' finer reciever grid (rg),',/,11x,           &
     &                    'It must be an exact divisible factor from ', &
     &                    'donor grid (dg).',/,/,11x,'rg = ',i2.2,      &
     &                    ', dt(rg) = ',f12.4,/,11x,'dg = ',i2.2,       &
     &                    ', dt(dg) = ',f12.4,4x,                       &
     &                    'MOD[dt(dg), dt(rg)] = ',f10.4,/)
                END IF
              END IF
              exit_flag=5
              RETURN
            END IF
          END IF
        END DO
      END IF
#endif

#ifdef FOUR_DVAR
!
!-----------------------------------------------------------------------
!  Spatial convolution parameters.
!-----------------------------------------------------------------------
!
!  Determine diffusion operator time-step size using FTCS stability
!  criteria:
!
!       Kh DTsizeH / MIN(DXmin,DYmin)^2  < Hgamma / 4
!
!       Kv DTsizeH / DZmin^2  < Vgamma / 2
!
!  where a Hgamma and Vgamma are used to scale the time-step below
!  its theoretical limit for stability and accurary.
!
      cff=MIN(DXmin(ng),DYmin(ng))
      DO rec=1,2
        DO is=1,NstateVar(ng)
#  ifdef SOLVE3D
          IF (is.le.(5+NT(ng))) THEN
#  else
          IF (is.le.3) THEN
#  endif
            DTsizeH(rec,is)=Hgamma(rec)*cff*cff/(4.0_r8*KhMax(ng))
          ELSE              ! use stability factor for surface forcing
            DTsizeH(rec,is)=Hgamma(4)*cff*cff/(4.0_r8*KhMax(ng))
          END IF
# ifdef SOLVE3D
#  ifdef IMPLICIT_VCONV
          DTsizeV(rec,is)=Vgamma(rec)*DZmax(ng)*DZmax(ng)/              &
     &                    (2.0_r8*KvMax(ng))
#  else
          DTsizeV(rec,is)=Vgamma(rec)*DZmin(ng)*DZmin(ng)/              &
     &                    (2.0_r8*KvMax(ng))
#  endif
# endif
        END DO
      END DO

# ifdef ADJUST_BOUNDARY
!
!  Open boundaries convolutions. Recall that the horizontal convolutions
!  are one-dimensional so a factor of 2 is used instead.
!
      DO ibry=1,4
        DO is=1,NstateVar(ng)
          IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
            DTsizeHB(ibry,is)=Hgamma(3)*DYmin(ng)*DYmin(ng)/            &
     &                        (2.0_r8*KhMax(ng))
          ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
            DTsizeHB(ibry,is)=Hgamma(3)*DXmin(ng)*DXmin(ng)/            &
     &                        (2.0_r8*KhMax(ng))
          END IF
#  ifdef SOLVE3D
#   ifdef IMPLICIT_VCONV
          DTsizeVB(ibry,is)=Vgamma(3)*DZmax(ng)*DZmax(ng)/              &
     &                      (2.0_r8*KvMax(ng))
#   else
          DTsizeVB(ibry,is)=Vgamma(3)*DZmin(ng)*DZmin(ng)/              &
     &                      (2.0_r8*KvMax(ng))
#   endif
#  endif
        END DO
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Determine FULL number of integeration steps as function of the
!  spatial decorrelation scale (Hdecay, Vdecay). Notice that in
!  the squared-root filter only half of number of step is used.
!  Thefore, the number of steps are forced to be even.
!-----------------------------------------------------------------------
!
!  Initial conditions, surface forcing and model error covariance.
!
      DO rec=1,2
        DO is=1,NstateVar(ng)
          IF (Hdecay(rec,is,ng).gt.0.0_r8) THEN
            NHsteps(rec,is)=NINT(Hdecay(rec,is,ng)*                     &
     &                           Hdecay(rec,is,ng)/                     &
     &                           (2.0_r8*KhMax(ng)*DTsizeH(rec,is)))
            IF (MOD(NHsteps(rec,is),2).ne.0) THEN
              NHsteps(rec,is)=NHsteps(rec,is)+1
            END IF
          END IF
# ifdef SOLVE3D
          IF (Vdecay(rec,is,ng).gt.0.0_r8) THEN
            NVsteps(rec,is)=NINT(Vdecay(rec,is,ng)*                     &
     &                           Vdecay(rec,is,ng)/                     &
     &                           (2.0_r8*KvMax(ng)*DTsizeV(rec,is)))
#  ifdef IMPLICIT_VCONV
            NVsteps(rec,is)=MAX(1,NVsteps(rec,is))
#  endif
            IF (MOD(NVsteps(rec,is),2).ne.0) THEN
              NVsteps(rec,is)=NVsteps(rec,is)+1
            END IF
# endif
          END IF
        END DO
      END DO

# ifdef ADJUST_BOUNDARY
!
!  Boundary conditions model error covariance.
!
      DO ibry=1,4
        DO is=1,NstateVar(ng)
          IF (HdecayB(is,ibry,ng).gt.0.0_r8) THEN
            NHstepsB(ibry,is)=NINT(HdecayB(is,ibry,ng)*                 &
     &                             HdecayB(is,ibry,ng)/                 &
     &                             (2.0_r8*KhMax(ng)*DTsizeHB(ibry,is)))
            IF (MOD(NHstepsB(ibry,is),2).ne.0) THEN
              NHstepsB(ibry,is)=NHstepsB(ibry,is)+1
            END IF
          END IF
#  ifdef SOLVE3D
          IF (VdecayB(is,ibry,ng).gt.0.0_r8) THEN
            NVstepsB(ibry,is)=NINT(VdecayB(is,ibry,ng)*                 &
     &                             VdecayB(is,ibry,ng)/                 &
     &                             (2.0_r8*KvMax(ng)*DTsizeVB(ibry,is)))
#   ifdef IMPLICIT_VCONV
            NVstepsB(ibry,is)=MAX(1,NVstepsB(ibry,is))
#   endif
            IF (MOD(NVstepsB(ibry,is),2).ne.0) THEN
              NVstepsB(ibry,is)=NVstepsB(ibry,is)+1
            END IF
#  endif
          END IF
        END DO
      END DO
# endif
!
!  Report convolution parameters.
!
      IF (Master.and.LwrtInfo(ng)) THEN
        LwrtInfo(ng)=.FALSE.
        DO rec=1,NSA
          DO is=1,NstateVar(ng)
            IF (rec.eq.1) THEN
              Text='Horizontal convolution, NHstepsI, DTsizeH  ='
            ELSE IF (rec.eq.2) THEN
              Text='Horizontal convolution, NHstepsM, DTsizeH  ='
#  ifdef SOLVE3D
            ELSE IF (is.gt.(5+NT(ng))) THEN
#  else
            ELSE IF (is.gt.3) THEN
#  endif
              Text='Horizontal convolution, NHstepsF, DTsizeH  ='
            END IF
            IF (Hdecay(rec,is,ng).gt.0.0_r8) THEN
              WRITE (stdout,110) TRIM(Text), NHsteps(rec,is),           &
     &                           DTsizeH(rec,is),                       &
     &                           TRIM(Vname(1,idSvar(is)))
            END IF
          END DO
          WRITE (stdout,120)
        END DO
# if defined SOLVE3D && defined VCONVOLUTION
        DO rec=1,NSA
          DO is=1,NstateVar(ng)
            IF (rec.eq.1) THEN
              Text='Vertical   convolution, NVstepsI, DTsizeV  ='
            ELSE IF (rec.eq.2) THEN
              Text='Vertical   convolution, NVstepsM, DTsizeV  ='
#  ifdef SOLVE3D
            ELSE IF (is.gt.(5+NT(ng))) THEN
#  else
            ELSE IF (is.gt.3) THEN
#  endif
              Text='Vertical   convolution, NVstepsF, DTsizeV  ='
            END IF
            IF (Vdecay(rec,is,ng).gt.0.0_r8) THEN
              WRITE (stdout,110) TRIM(Text), NVsteps(rec,is),           &
     &                           DTsizeV(rec,is),                       &
     &                           TRIM(Vname(1,idSvar(is)))
             END IF
          END DO
          WRITE (stdout,120)
        END DO
# endif
# ifdef ADJUST_BOUNDARY
        DO is=1,NstateVar(ng)
          DO ibry=1,4
            IF (Lobc(ibry,is,ng)) THEN
              IF (ibry.eq.iwest) THEN
                Text='West  bry Hconvolution, NHstepsB, DTsizeHB ='
              ELSE IF (ibry.eq.isouth) THEN
                Text='South bry Hconvolution, NHstepsB, DTsizeHB ='
              ELSE IF (ibry.eq.ieast) THEN
                Text='East  bry Hconvolution, NHstepsB, DTsizeHB ='
              ELSE IF (ibry.eq.isouth) THEN
                Text='North bry Hconvolution, NHstepsB, DTsizeHB ='
              END IF
              IF (HdecayB(is,ibry,ng).gt.0.0_r8) THEN
                WRITE (stdout,110) TRIM(Text), NHstepsB(ibry,is),       &
     &                             DTsizeHB(ibry,is),                   &
     &                             TRIM(Vname(1,idSbry(is)))
              END IF
            END IF
          END DO
        END DO
        WRITE (stdout,120)
#  if defined SOLVE3D && defined VCONVOLUTION
        DO is=1,NstateVar(ng)
          DO ibry=1,4
            IF (Lobc(ibry,is,ng)) THEN
              IF (ibry.eq.iwest) THEN
                Text='West  bry Vconvolution, NVstepsB, DTsizeVB ='
              ELSE IF (ibry.eq.isouth) THEN
                Text='South bry Vconvolution, NVstepsB, DTsizeVB ='
              ELSE IF (ibry.eq.ieast) THEN
                Text='East  bry Vconvolution, NVstepsB, DTsizeVB ='
              ELSE IF (ibry.eq.isouth) THEN
                Text='North bry Vconvolution, NVstepsB, DTsizeVB ='
              END IF
              IF (VdecayB(is,ibry,ng).gt.0.0_r8) THEN
                WRITE (stdout,110) TRIM(Text), NVstepsB(ibry,is),       &
     &                             DTsizeVB(ibry,is),                   &
     &                             TRIM(Vname(1,idSbry(is)))
              END IF
            END IF
          END DO
        END DO
#  endif
# endif
      END IF

 110  FORMAT (1X,a,i5,1x,1pe15.8,' s',2x,a)
 120  FORMAT (1x)
#endif

      RETURN
      END SUBROUTINE metrics_tile

      END MODULE metrics_mod
