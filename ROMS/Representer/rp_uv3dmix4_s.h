      SUBROUTINE rp_uv3dmix4 (ng, tile)
!
!svn $Id: rp_uv3dmix4_s.h 795 2016-05-11 01:42:43Z arango $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  This subroutine computes representers tangent linear biharmonic     !
!  mixing  of  momentum,   along  constant  S-surfaces,  from  the     !
!  horizontal  divergence  of  the  stress  tensor.  A  transverse     !
!  isotropy is assumed so the stress tensor is split into vertical     !
!  and horizontal subtensors.                                          !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!      Wajsowicz, R.C, 1993: A consistent formulation of the           !
!         anisotropic stress tensor for use in models of the           !
!         large-scale ocean circulation, JCP, 105, 333-338.            !
!                                                                      !
!      Sadourny, R. and K. Maynard, 1997: Formulations of              !
!         lateral diffusion in geophysical fluid dynamics              !
!         models, In Numerical Methods of Atmospheric and              !
!         Oceanic Modelling. Lin, Laprise, and Ritchie,                !
!         Eds., NRC Research Press, 547-556.                           !
!                                                                      !
!      Griffies, S.M. and R.W. Hallberg, 2000: Biharmonic              !
!         friction with a Smagorinsky-like viscosity for               !
!         use in large-scale eddy-permitting ocean models,             !
!         Monthly Weather Rev., 128, 8, 2935-2946.                     !
!                                                                      !
!  Basic state variables required:  visc4, u, v, Hz.                   !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
#ifdef DIAGNOSTICS_UV
!!    USE mod_diags
#endif
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
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, iRPM, 32)
#endif
      CALL rp_uv3dmix4_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nrhs(ng), nnew(ng),                        &
#ifdef MASKING
     &                       GRID(ng) % pmask,                          &
#endif
     &                       GRID(ng) % Hz,                             &
     &                       GRID(ng) % tl_Hz,                          &
     &                       GRID(ng) % om_p,                           &
     &                       GRID(ng) % om_r,                           &
     &                       GRID(ng) % on_p,                           &
     &                       GRID(ng) % on_r,                           &
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pmon_p,                         &
     &                       GRID(ng) % pmon_r,                         &
     &                       GRID(ng) % pn,                             &
     &                       GRID(ng) % pnom_p,                         &
     &                       GRID(ng) % pnom_r,                         &
     &                       MIXING(ng) % visc4_p,                      &
     &                       MIXING(ng) % visc4_r,                      &
#ifdef DIAGNOSTICS_UV
!!   &                       DIAGS(ng) % DiaRUfrc,                      &
!!   &                       DIAGS(ng) % DiaRVfrc,                      &
!!   &                       DIAGS(ng) % DiaU3wrk,                      &
!!   &                       DIAGS(ng) % DiaV3wrk,                      &
#endif
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v,                             &
     &                       COUPLING(ng) % tl_rufrc,                   &
     &                       COUPLING(ng) % tl_rvfrc,                   &
     &                       OCEAN(ng) % tl_u,                          &
     &                       OCEAN(ng) % tl_v)
#ifdef PROFILE
      CALL wclock_off (ng, iRPM, 32)
#endif

      RETURN
      END SUBROUTINE rp_uv3dmix4

!
!***********************************************************************
      SUBROUTINE rp_uv3dmix4_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nrhs, nnew,                          &
#ifdef MASKING
     &                             pmask,                               &
#endif
     &                             Hz, tl_Hz,                           &
     &                             om_p, om_r, on_p, on_r,              &
     &                             pm, pmon_p, pmon_r,                  &
     &                             pn, pnom_p, pnom_r,                  &
     &                             visc4_p, visc4_r,                    &
#ifdef DIAGNOSTICS_UV
!!   &                             DiaRUfrc, DiaRVfrc,                  &
!!   &                             DiaU3wrk, DiaV3wrk,                  &
#endif
     &                             u, v,                                &
     &                             tl_rufrc, tl_rvfrc, tl_u, tl_v)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, nnew

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: om_p(LBi:,LBj:)
      real(r8), intent(in) :: om_r(LBi:,LBj:)
      real(r8), intent(in) :: on_p(LBi:,LBj:)
      real(r8), intent(in) :: on_r(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pmon_p(LBi:,LBj:)
      real(r8), intent(in) :: pmon_r(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: pnom_p(LBi:,LBj:)
      real(r8), intent(in) :: pnom_r(LBi:,LBj:)
      real(r8), intent(in) :: visc4_p(LBi:,LBj:)
      real(r8), intent(in) :: visc4_r(LBi:,LBj:)

      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)

# ifdef DIAGNOSTICS_UV
!!    real(r8), intent(inout) :: DiaRUfrc(LBi:,LBj:,:,:)
!!    real(r8), intent(inout) :: DiaRVfrc(LBi:,LBj:,:,:)
!!    real(r8), intent(inout) :: DiaU3wrk(LBi:,LBj:,:,:)
!!    real(r8), intent(inout) :: DiaV3wrk(LBi:,LBj:,:,:)
# endif

      real(r8), intent(inout) :: tl_rufrc(LBi:,LBj:)
      real(r8), intent(inout) :: tl_rvfrc(LBi:,LBj:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)

#else

# ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: om_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: visc4_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: visc4_r(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,N(ng),2)

# ifdef DIAGNOSTICS_UV
!!    real(r8), intent(inout) :: DiaRUfrc(LBi:UBi,LBj:UBj,3,NDM2d-1)
!!    real(r8), intent(inout) :: DiaRVfrc(LBi:UBi,LBj:UBj,3,NDM2d-1)
!!    real(r8), intent(inout) :: DiaU3wrk(LBi:UBi,LBj:UBj,N(ng),NDM3d)
!!    real(r8), intent(inout) :: DiaV3wrk(LBi:UBi,LBj:UBj,N(ng),NDM3d)
# endif

      real(r8), intent(inout) :: tl_rufrc(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_rvfrc(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#endif
!
!  Local variable declarations.
!
      integer :: IminU, IminV, ImaxU, ImaxV
      integer :: JminU, JminV, JmaxU, JmaxV
      integer :: i, j, k

      real(r8) :: cff, cff1, cff2
      real(r8) :: tl_cff, tl_cff1, tl_cff2

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: LapU
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: LapV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFx

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_LapU
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_LapV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_UFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_VFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_UFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_VFx

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute horizontal biharmonic viscosity along constant S-surfaces.
!  The biharmonic operator is computed by applying the harmonic
!  operator twice.
!-----------------------------------------------------------------------
!
!  Set local I- and J-ranges.
!
      IF (EWperiodic(ng)) THEN
        IminU=Istr-1
        ImaxU=Iend+1
        IminV=Istr-1
        ImaxV=Iend+1
      ELSE
        IminU=MAX(2,IstrU-1)
        ImaxU=MIN(Iend+1,Lm(ng))
        IminV=MAX(1,Istr-1)
        ImaxV=MIN(Iend+1,Lm(ng))
      END IF
      IF (NSperiodic(ng)) THEN
        JminU=Jstr-1
        JmaxU=Jend+1
        JminV=Jstr-1
        JmaxV=Jend+1
      ELSE
        JminU=MAX(1,Jstr-1)
        JmaxU=MIN(Jend+1,Mm(ng))
        JminV=MAX(2,JstrV-1)
        JmaxV=MIN(Jend+1,Mm(ng))
      END IF
!
!  Compute flux-components of the horizontal divergence of the stress
!  tensor (m4 s^-3/2) in XI- and ETA-directions.  It is assumed here
!  that mixing coefficients are the squared root of the biharmonic
!  viscosity coefficient.  For momentum balance purposes, the
!  thickness "Hz" appears only when computing the second harmonic
!  operator.
!
      K_LOOP : DO k=1,N(ng)
        DO j=JminV-1,JmaxV
          DO i=IminU-1,ImaxU
            cff=visc4_r(i,j)*0.5_r8*                                    &
     &          (pmon_r(i,j)*                                           &
     &           ((pn(i  ,j)+pn(i+1,j))*u(i+1,j,k,nrhs)-                &
     &            (pn(i-1,j)+pn(i  ,j))*u(i  ,j,k,nrhs))-               &
     &           pnom_r(i,j)*                                           &
     &           ((pm(i,j  )+pm(i,j+1))*v(i,j+1,k,nrhs)-                &
     &            (pm(i,j-1)+pm(i,j  ))*v(i,j  ,k,nrhs)))
            tl_cff=visc4_r(i,j)*0.5_r8*                                 &
     &             (pmon_r(i,j)*                                        &
     &              ((pn(i  ,j)+pn(i+1,j))*tl_u(i+1,j,k,nrhs)-          &
     &               (pn(i-1,j)+pn(i  ,j))*tl_u(i  ,j,k,nrhs))-         &
     &              pnom_r(i,j)*                                        &
     &              ((pm(i,j  )+pm(i,j+1))*tl_v(i,j+1,k,nrhs)-          &
     &               (pm(i,j-1)+pm(i,j  ))*tl_v(i,j  ,k,nrhs)))
            UFx(i,j)=on_r(i,j)*on_r(i,j)*cff
            tl_UFx(i,j)=on_r(i,j)*on_r(i,j)*tl_cff
            VFe(i,j)=om_r(i,j)*om_r(i,j)*cff
            tl_VFe(i,j)=om_r(i,j)*om_r(i,j)*tl_cff
          END DO
        END DO
        DO j=JminU,JmaxU+1
          DO i=IminV,ImaxV+1
            cff=visc4_p(i,j)*0.5_r8*                                    &
     &          (pmon_p(i,j)*                                           &
     &           ((pn(i  ,j-1)+pn(i  ,j))*v(i  ,j,k,nrhs)-              &
     &            (pn(i-1,j-1)+pn(i-1,j))*v(i-1,j,k,nrhs))+             &
     &           pnom_p(i,j)*                                           &
     &           ((pm(i-1,j  )+pm(i,j  ))*u(i,j  ,k,nrhs)-              &
     &            (pm(i-1,j-1)+pm(i,j-1))*u(i,j-1,k,nrhs)))
            tl_cff=visc4_p(i,j)*0.5_r8*                                 &
     &             (pmon_p(i,j)*                                        &
     &              ((pn(i  ,j-1)+pn(i  ,j))*tl_v(i  ,j,k,nrhs)-        &
     &               (pn(i-1,j-1)+pn(i-1,j))*tl_v(i-1,j,k,nrhs))+       &
     &              pnom_p(i,j)*                                        &
     &              ((pm(i-1,j  )+pm(i,j  ))*tl_u(i,j  ,k,nrhs)-        &
     &               (pm(i-1,j-1)+pm(i,j-1))*tl_u(i,j-1,k,nrhs)))
#ifdef MASKING
            cff=cff*pmask(i,j)
            tl_cff=tl_cff*pmask(i,j)
#endif
            UFe(i,j)=om_p(i,j)*om_p(i,j)*cff
            tl_UFe(i,j)=om_p(i,j)*om_p(i,j)*tl_cff
            VFx(i,j)=on_p(i,j)*on_p(i,j)*cff
            tl_VFx(i,j)=on_p(i,j)*on_p(i,j)*tl_cff
          END DO
        END DO
!
!  Compute first harmonic operator (m s^-3/2).
!
        DO j=JminU,JmaxU
          DO i=IminU,ImaxU
            LapU(i,j)=0.125_r8*                                         &
     &                (pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))*          &
     &                ((pn(i-1,j)+pn(i,j))*                             &
     &                 (UFx(i,j  )-UFx(i-1,j))+                         &
     &                 (pm(i-1,j)+pm(i,j))*                             &
     &                 (UFe(i,j+1)-UFe(i  ,j)))
            tl_LapU(i,j)=0.125_r8*                                      &
     &                   (pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))*       &
     &                   ((pn(i-1,j)+pn(i,j))*                          &
     &                    (tl_UFx(i,j  )-tl_UFx(i-1,j))+                &
     &                    (pm(i-1,j)+pm(i,j))*                          &
     &                    (tl_UFe(i,j+1)-tl_UFe(i  ,j)))
          END DO
        END DO
        DO j=JminV,JmaxV
          DO i=IminV,ImaxV
            LapV(i,j)=0.125_r8*                                         &
     &                (pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))*          &
     &                ((pn(i,j-1)+pn(i,j))*                             &
     &                 (VFx(i+1,j)-VFx(i,j  ))-                         &
     &                 (pm(i,j-1)+pm(i,j))*                             &
     &                 (VFe(i  ,j)-VFe(i,j-1)))
            tl_LapV(i,j)=0.125_r8*                                      &
     &                   (pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))*       &
     &                   ((pn(i,j-1)+pn(i,j))*                          &
     &                    (tl_VFx(i+1,j)-tl_VFx(i,j  ))-                &
     &                    (pm(i,j-1)+pm(i,j))*                          &
     &                    (tl_VFe(i  ,j)-tl_VFe(i,j-1)))
          END DO
        END DO
!
!  Apply boundary conditions (other than periodic) to the first
!  harmonic operator. These are gradient or closed (free slip or
!  no slip) boundary conditions.
!
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            IF (tl_LBC(iwest,isUvel,ng)%closed) THEN
              DO j=JminU,JmaxU
                LapU(Istr,j)=0.0_r8
                tl_LapU(Istr,j)=0.0_r8
              END DO
            ELSE
              DO j=JminU,JmaxU
                LapU(Istr,j)=LapU(Istr+1,j)
                tl_LapU(Istr,j)=tl_LapU(Istr+1,j)
              END DO
            END IF
            IF (tl_LBC(iwest,isVvel,ng)%closed) THEN
              DO j=JminV,JmaxV
                LapV(Istr-1,j)=gamma2(ng)*LapV(Istr,j)
                tl_LapV(Istr-1,j)=gamma2(ng)*tl_LapV(Istr,j)
              END DO
            ELSE
              DO j=JminV,JmaxV
                LapV(Istr-1,j)=0.0_r8
                tl_LapV(Istr-1,j)=0.0_r8
              END DO
            END IF
          END IF
        END IF
!
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            IF (tl_LBC(ieast,isUvel,ng)%closed) THEN
              DO j=JminU,JmaxU
                LapU(Iend+1,j)=0.0_r8
                tl_LapU(Iend+1,j)=0.0_r8
              END DO
            ELSE
              DO j=JminU,JmaxU
                LapU(Iend+1,j)=LapU(Iend,j)
                tl_LapU(Iend+1,j)=tl_LapU(Iend,j)
              END DO
            END IF
            IF (tl_LBC(ieast,isVvel,ng)%closed) THEN
              DO j=JminV,JmaxV
                LapV(Iend+1,j)=gamma2(ng)*LapV(Iend,j)
                tl_LapV(Iend+1,j)=gamma2(ng)*tl_LapV(Iend,j)
              END DO
            ELSE
              DO j=JminV,JmaxV
                LapV(Iend+1,j)=0.0_r8
                tl_LapV(Iend+1,j)=0.0_r8
              END DO
            END IF
          END IF
        END IF
!
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            IF (tl_LBC(isouth,isUvel,ng)%closed) THEN
              DO i=IminU,ImaxU
                LapU(i,Jstr-1)=gamma2(ng)*LapU(i,Jstr)
                tl_LapU(i,Jstr-1)=gamma2(ng)*tl_LapU(i,Jstr)
              END DO
            ELSE
              DO i=IminU,ImaxU
                LapU(i,Jstr-1)=0.0_r8
                tl_LapU(i,Jstr-1)=0.0_r8
              END DO
            END IF
            IF (tl_LBC(isouth,isVvel,ng)%closed) THEN
              DO i=IminV,ImaxV
                LapV(i,Jstr)=0.0_r8
                tl_LapV(i,Jstr)=0.0_r8
              END DO
            ELSE
              DO i=IminV,ImaxV
                LapV(i,Jstr)=LapV(i,Jstr+1)
                tl_LapV(i,Jstr)=tl_LapV(i,Jstr+1)
              END DO
            END IF
          END IF
        END IF
!
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            IF (tl_LBC(inorth,isUvel,ng)%closed) THEN
              DO i=IminU,ImaxU
                LapU(i,Jend+1)=gamma2(ng)*LapU(i,Jend)
                tl_LapU(i,Jend+1)=gamma2(ng)*tl_LapU(i,Jend)
              END DO
            ELSE
              DO i=IminU,ImaxU
                LapU(i,Jend+1)=0.0_r8
                tl_LapU(i,Jend+1)=0.0_r8
              END DO
            END IF
            IF (tl_LBC(inorth,isVvel,ng)%closed) THEN
              DO i=IminV,ImaxV
                LapV(i,Jend+1)=0.0_r8
                tl_LapV(i,Jend+1)=0.0_r8
              END DO
            ELSE
              DO i=IminV,ImaxV
                LapV(i,Jend+1)=LapV(i,Jend)
                tl_LapV(i,Jend+1)=tl_LapV(i,Jend)
              END DO
            END IF
          END IF
        END IF
!
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(iwest ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            LapU(Istr  ,Jstr-1)=0.5_r8*                                 &
     &                          (LapU(Istr+1,Jstr-1)+                   &
     &                           LapU(Istr  ,Jstr  ))
            tl_LapU(Istr  ,Jstr-1)=0.5_r8*                              &
     &                             (tl_LapU(Istr+1,Jstr-1)+             &
     &                              tl_LapU(Istr  ,Jstr  ))
            LapV(Istr-1,Jstr  )=0.5_r8*                                 &
     &                          (LapV(Istr-1,Jstr+1)+                   &
     &                           LapV(Istr  ,Jstr  ))
            tl_LapV(Istr-1,Jstr  )=0.5_r8*                              &
     &                             (tl_LapV(Istr-1,Jstr+1)+             &
     &                              tl_LapV(Istr  ,Jstr  ))
          END IF
        END IF

        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(ieast ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            LapU(Iend+1,Jstr-1)=0.5_r8*                                 &
     &                          (LapU(Iend  ,Jstr-1)+                   &
     &                           LapU(Iend+1,Jstr  ))
            tl_LapU(Iend+1,Jstr-1)=0.5_r8*                              &
     &                             (tl_LapU(Iend  ,Jstr-1)+             &
     &                              tl_LapU(Iend+1,Jstr  ))
            LapV(Iend+1,Jstr  )=0.5_r8*                                 &
     &                          (LapV(Iend  ,Jstr  )+                   &
     &                           LapV(Iend+1,Jstr+1))
            tl_LapV(Iend+1,Jstr  )=0.5_r8*                              &
     &                             (tl_LapV(Iend  ,Jstr  )+             &
     &                              tl_LapV(Iend+1,Jstr+1))
          END IF
        END IF

        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(iwest ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            LapU(Istr  ,Jend+1)=0.5_r8*                                 &
     &                          (LapU(Istr+1,Jend+1)+                   &
     &                           LapU(Istr  ,Jend  ))
            tl_LapU(Istr  ,Jend+1)=0.5_r8*                              &
     &                             (tl_LapU(Istr+1,Jend+1)+             &
     &                              tl_LapU(Istr  ,Jend  ))
            LapV(Istr-1,Jend+1)=0.5_r8*                                 &
     &                          (LapV(Istr  ,Jend+1)+                   &
     &                           LapV(Istr-1,Jend  ))
            tl_LapV(Istr-1,Jend+1)=0.5_r8*                              &
     &                             (tl_LapV(Istr  ,Jend+1)+             &
     &                              tl_LapV(Istr-1,Jend  ))
          END IF
        END IF

        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(ieast ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            LapU(Iend+1,Jend+1)=0.5_r8*                                 &
     &                          (LapU(Iend  ,Jend+1)+                   &
     &                           LapU(Iend+1,Jend  ))
            tl_LapU(Iend+1,Jend+1)=0.5_r8*                              &
     &                             (tl_LapU(Iend  ,Jend+1)+             &
     &                              tl_LapU(Iend+1,Jend  ))
            LapV(Iend+1,Jend+1)=0.5_r8*                                 &
     &                          (LapV(Iend  ,Jend+1)+                   &
     &                           LapV(Iend+1,Jend  ))
            tl_LapV(Iend+1,Jend+1)=0.5_r8*                              &
     &                             (tl_LapV(Iend  ,Jend+1)+             &
     &                              tl_LapV(Iend+1,Jend  ))
          END IF
        END IF
!
!  Compute flux-components of the horizontal divergence of the
!  harmonic stress tensor (m4/s2) in XI- and ETA-directions.
!
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            cff=visc4_r(i,j)*Hz(i,j,k)*0.5_r8*                          &
     &          (pmon_r(i,j)*                                           &
     &           ((pn(i  ,j)+pn(i+1,j))*LapU(i+1,j)-                    &
     &            (pn(i-1,j)+pn(i  ,j))*LapU(i  ,j))-                   &
     &           pnom_r(i,j)*                                           &
     &           ((pm(i,j  )+pm(i,j+1))*LapV(i,j+1)-                    &
     &            (pm(i,j-1)+pm(i,j  ))*LapV(i,j  )))
            tl_cff=visc4_r(i,j)*0.5_r8*                                 &
     &             (tl_Hz(i,j,k)*                                       &
     &              (pmon_r(i,j)*                                       &
     &               ((pn(i  ,j)+pn(i+1,j))*LapU(i+1,j)-                &
     &                (pn(i-1,j)+pn(i  ,j))*LapU(i  ,j))-               &
     &               pnom_r(i,j)*                                       &
     &               ((pm(i,j  )+pm(i,j+1))*LapV(i,j+1)-                &
     &                (pm(i,j-1)+pm(i,j  ))*LapV(i,j  )))+              &
     &              Hz(i,j,k)*                                          &
     &              (pmon_r(i,j)*                                       &
     &               ((pn(i  ,j)+pn(i+1,j))*tl_LapU(i+1,j)-             &
     &                (pn(i-1,j)+pn(i  ,j))*tl_LapU(i  ,j))-            &
     &               pnom_r(i,j)*                                       &
     &               ((pm(i,j  )+pm(i,j+1))*tl_LapV(i,j+1)-             &
     &                (pm(i,j-1)+pm(i,j  ))*tl_LapV(i,j  ))))-          &
#ifdef TL_IOMS
     &             cff
#endif
!>          UFx(i,j)=on_r(i,j)*on_r(i,j)*cff
!>
            tl_UFx(i,j)=on_r(i,j)*on_r(i,j)*tl_cff
!>          VFe(i,j)=om_r(i,j)*om_r(i,j)*cff
!>
            tl_VFe(i,j)=om_r(i,j)*om_r(i,j)*tl_cff
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend+1
            cff=visc4_p(i,j)*0.125_r8*(Hz(i-1,j  ,k)+Hz(i,j  ,k)+       &
     &                                 Hz(i-1,j-1,k)+Hz(i,j-1,k))*      &
     &          (pmon_p(i,j)*                                           &
     &           ((pn(i  ,j-1)+pn(i  ,j))*LapV(i  ,j)-                  &
     &            (pn(i-1,j-1)+pn(i-1,j))*LapV(i-1,j))+                 &
     &           pnom_p(i,j)*                                           &
     &           ((pm(i-1,j  )+pm(i,j  ))*LapU(i,j  )-                  &
     &            (pm(i-1,j-1)+pm(i,j-1))*LapU(i,j-1)))
            tl_cff=visc4_p(i,j)*0.125_r8*                               &
     &             ((tl_Hz(i-1,j  ,k)+tl_Hz(i,j  ,k)+                   &
     &               tl_Hz(i-1,j-1,k)+tl_Hz(i,j-1,k))*                  &
     &              (pmon_p(i,j)*                                       &
     &               ((pn(i  ,j-1)+pn(i  ,j))*LapV(i  ,j)-              &
     &                (pn(i-1,j-1)+pn(i-1,j))*LapV(i-1,j))+             &
     &               pnom_p(i,j)*                                       &
     &               ((pm(i-1,j  )+pm(i,j  ))*LapU(i,j  )-              &
     &                (pm(i-1,j-1)+pm(i,j-1))*LapU(i,j-1)))+            &
     &              (Hz(i-1,j  ,k)+Hz(i,j  ,k)+                         &
     &               Hz(i-1,j-1,k)+Hz(i,j-1,k))*                        &
     &               (pmon_p(i,j)*                                      &
     &                ((pn(i  ,j-1)+pn(i  ,j))*tl_LapV(i  ,j)-          &
     &                 (pn(i-1,j-1)+pn(i-1,j))*tl_LapV(i-1,j))+         &
     &                pnom_p(i,j)*                                      &
     &                ((pm(i-1,j  )+pm(i,j  ))*tl_LapU(i,j  )-          &
     &                 (pm(i-1,j-1)+pm(i,j-1))*tl_LapU(i,j-1))))-       &
#ifdef TL_IOMS
     &             cff
#endif
#ifdef MASKING
!>          cff=cff*pmask(i,j)
!>
            tl_cff=tl_cff*pmask(i,j)
#endif
!>          UFe(i,j)=om_p(i,j)*om_p(i,j)*cff
!>
            tl_UFe(i,j)=om_p(i,j)*om_p(i,j)*tl_cff
!>          VFx(i,j)=on_p(i,j)*on_p(i,j)*cff
!>
            tl_VFx(i,j)=on_p(i,j)*on_p(i,j)*tl_cff
          END DO
        END DO
!
! Time-step biharmonic, S-surfaces viscosity term.  Notice that
! momentum at this stage is HzU and HzV and has units m2/s.  Add
! contribution for barotropic forcing terms.
!
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=0.25_r8*(pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))
!>          cff1=0.5_r8*((pn(i-1,j)+pn(i,j))*                           &
!>   &                   (UFx(i,j  )-UFx(i-1,j))+                       &
!>   &                   (pm(i-1,j)+pm(i,j))*                           &
!>   &                   (UFe(i,j+1)-UFe(i  ,j)))
!>
            tl_cff1=0.5_r8*((pn(i-1,j)+pn(i,j))*                        &
     &                      (tl_UFx(i,j  )-tl_UFx(i-1,j))+              &
     &                      (pm(i-1,j)+pm(i,j))*                        &
     &                      (tl_UFe(i,j+1)-tl_UFe(i  ,j)))
!>          cff2=dt(ng)*cff*cff1
!>
            tl_cff2=dt(ng)*cff*tl_cff1
!>          rufrc(i,j)=rufrc(i,j)-cff1
!>
            tl_rufrc(i,j)=tl_rufrc(i,j)-tl_cff1
!>          u(i,j,k,nnew)=u(i,j,k,nnew)-cff2
!>
            tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)-tl_cff2
#ifdef DIAGNOSTICS_UV
!!          DiaRUfrc(i,j,3,M2hvis)=DiaRUfrc(i,j,3,M2hvis)-cff1
!!          DiaU3wrk(i,j,k,M3hvis)=-cff2
#endif
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=0.25_r8*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
!>          cff1=0.5_r8*((pn(i,j-1)+pn(i,j))*                           &
!>   &                   (VFx(i+1,j)-VFx(i,j  ))-                       &
!>   &                   (pm(i,j-1)+pm(i,j))*                           &
!>   &                   (VFe(i  ,j)-VFe(i,j-1)))
!>
            tl_cff1=0.5_r8*((pn(i,j-1)+pn(i,j))*                        &
     &                      (tl_VFx(i+1,j)-tl_VFx(i,j  ))-              &
     &                      (pm(i,j-1)+pm(i,j))*                        &
     &                      (tl_VFe(i  ,j)-tl_VFe(i,j-1)))
!>          cff2=dt(ng)*cff*cff1
!>
            tl_cff2=dt(ng)*cff*tl_cff1
!>          rvfrc(i,j)=rvfrc(i,j)-cff1
!>
            tl_rvfrc(i,j)=tl_rvfrc(i,j)-tl_cff1
!>          v(i,j,k,nnew)=v(i,j,k,nnew)-cff2
!>
            tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)-tl_cff2
#ifdef DIAGNOSTICS_UV
!!          DiaRVfrc(i,j,3,M2hvis)=DiaRVfrc(i,j,3,M2hvis)-cff1
!!          DiaV3wrk(i,j,k,M3hvis)=-cff2
#endif
          END DO
        END DO
      END DO K_LOOP

      RETURN
      END SUBROUTINE rp_uv3dmix4_tile
