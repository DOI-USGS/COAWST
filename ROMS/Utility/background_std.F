#include "cppdefs.h"
      MODULE background_std_mod

#ifdef STD_MODEL
!
!!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group            Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  Standard deviation formulation used to model the 4D-Var Background  !
!  Error Covariance matrix, B. It follows the method Mogensen et al.   !
!  (2012) that assumes that the background errors are proportional to  !
!  vertical derivatives of the state. That is, the difference between  !
!  the background state value S_b and true value S_t is due to a       !
!  vertical displacement of the profile:                               !
!                                                                      !
!    S_t(z) ~ S_b(z+deltaz) + [d(S_b)/d(z)]*deltaz                     !
!                                                                      !
!  where deltaz is the displacement.                                   !
!                                                                      !
!  I assumes that the background and true profiles have similar shape  !
!  and the true value of S lies somewhere in the background water      !
!  column. The error is then [d(S_b)/d(z)]*deltaz.                     !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Kara A., P. Rochford, and E. Hulburt, 2000: An optimal definition   !
!      for ocean mixed layer depth, J. Geophys. Res. vol 105, NoC7,    !
!      pp 16, 803-16, 821.                                             !
!                                                                      !
!  Mogensen, K., M.A. Balmaseda, and A.T. Weaver, 2012: The NEMOVAR    !
!    ocean data assimilation system as implemented in the ECMWF ocean  !
!    analysis for system 4. ECMWF Tech. Memorandum 668, 59.            !
!                                                                      !
!  Moore, A., J. Zavala-Garay, H.G. Arango, C.A. Edwards, J. Anderson, !
!    and T. Hoar, 2020: Regional and basin scale applications of       !
!    ensemble adjustement Kalman filter and 4D-Var ocean data          !
!    assimilation systems, Progress in Oceanography, 189, 102450,      !
!    https://doi.org/10.1016/j.pocean.2020.102450.                     !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
      USE mod_param
      USE mod_parallel
# ifdef SOLVE3D
      USE mod_coupling
# endif
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile,                     &
     &                            exchange_u2d_tile, exchange_v2d_tile
# ifdef SOLVE3D
      USE exchange_3d_mod, ONLY : exchange_r3d_tile,                     &
     &                            exchange_u3d_tile, exchange_v3d_tile
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
# endif
# ifdef SOLVE3D
      USE set_depth_mod,   ONLY : set_depth
# endif
!
      implicit none
!
      PRIVATE
      PUBLIC :: background_std
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE background_std (ng, tile, Lbck, Lstd)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Lbck, Lstd
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef SOLVE3D
!
!  Compute background state thickness, depth arrays, thermal expansion,
!  and saline contraction coefficients.
!
      COUPLING(ng) % Zt_avg1 = 0.0_r8

      CALL set_depth (ng, tile, iNLM)
!
# endif
      CALL background_std_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Lbck, Lstd,                             &
# ifdef SOLVE3D
     &                          GRID(ng) % Hz,                          &
     &                          GRID(ng) % z_r,                         &
# endif
# ifdef MASKING
     &                          GRID(ng) % rmask,                       &
     &                          GRID(ng) % umask,                       &
     &                          GRID(ng) % vmask,                       &
# endif
# ifdef SOLVE3D
     &                          OCEAN(ng) % t,                          &
     &                          OCEAN(ng) % u,                          &
     &                          OCEAN(ng) % v,                          &
     &                          OCEAN(ng) % e_t,                        &
     &                          OCEAN(ng) % e_u,                        &
     &                          OCEAN(ng) % e_v,                        &
# endif
     &                          OCEAN(ng) % e_ubar,                     &
     &                          OCEAN(ng) % e_vbar,                     &
     &                          OCEAN(ng) % e_zeta)
!
      IF (Master) WRITE (stdout,10)
 10   FORMAT (/,2x,'BACKGROUND_STD   - computing standard deviation',   &
     &        ' from prior.',/)
!
      RETURN
      END SUBROUTINE background_std
!
!***********************************************************************
      SUBROUTINE background_std_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                Lbck, Lstd,                       &
# ifdef SOLVE3D
     &                                Hz, z_r,                          &
# endif
# ifdef MASKING
     &                                rmask, umask, vmask,              &
# endif
# ifdef SOLVE3D
     &                                t, u, v,                          &
     &                                e_t, e_u, e_v,                    &
# endif
     &                                e_ubar, e_vbar, e_zeta)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lbck, Lstd
!
# ifdef ASSUMED_SHAPE
#  ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
#  endif
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: e_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: e_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: e_v(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(inout) :: e_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: e_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: e_zeta(LBi:,LBj:,:)

# else

#  ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
#  endif
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,2,N(ng))
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,2,N(ng))
      real(r8), intent(inout) :: e_t(LBi:UBi,LBj:UBj,N(ng),NSA,NT(ng))
      real(r8), intent(inout) :: e_u(LBi:UBi,LBj:UBj,N(ng),NSA)
      real(r8), intent(inout) :: e_v(LBi:UBi,LBj:UBj,N(ng),NSA)
#  endif
      real(r8), intent(inout) :: e_ubar(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(inout) :: e_vbar(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(inout) :: e_zeta(LBi:UBi,LBj:UBj,NSA)

# endif
!
!  Local variable declarations.
!
      logical :: base_reached = .FALSE.
      logical :: ml_reached = .FALSE.
!
      integer, parameter :: Norder = 2            ! Shapiro filter order
      integer :: i, j, k, kref, khref, order
!
      real(r8) :: Temp_ref, T_dep, T_high, T_low, T_thvalue
      real(r8) :: sigmabS, sigmabT, sigmabU, SigmabV
      real(r8) :: cff, cff1, cff2, fac
      real(r8) :: href
!
!  Shapiro filter coefficients.
!
      real(r8), dimension(20) ::  filter_coef =                         &
     &   (/ 2.500000E-1_r8,    6.250000E-2_r8,     1.562500E-2_r8,      &
     &      3.906250E-3_r8,    9.765625E-4_r8,     2.44140625E-4_r8,    &
     &      6.103515625E-5_r8, 1.5258789063E-5_r8, 3.814697E-6_r8,      &
     &      9.536743E-7_r8,    2.384186E-7_r8,     5.960464E-8_r8,      &
     &      1.490116E-8_r8,    3.725290E-9_r8,     9.313226E-10_r8,     &
     &      2.328306E-10_r8,   5.820766E-11_r8,    1.455192E-11_r8,     &
     &      3.637979E-12_r8,   9.094947E-13_r8 /)

      real(r8), dimension(N(ng)) :: dSdT, dSdT_filter

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: dTdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: dUdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: dVdz
      real(r8), dimension(N(ng)) :: dTdz_filter
      real(r8), dimension(N(ng)) :: dUdz_filter, dVdz_filter
# ifdef SALINITY
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: dSdz
      real(r8), dimension(N(ng)) :: dSdz_filter
# endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: mld

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
! Compute temperature (dTdz) and salinity (dSdz) shears.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          FC(i,0)=0.0_r8
          dTdz(i,j,0)=0.0_r8
# ifdef SALINITY
          dSdz(i,j,0)=0.0_r8
# endif
        END DO
        DO k=1,N(ng)-1
          DO i=IstrR,IendR
            cff=1.0_r8/(2.0_r8*Hz(i,j,k+1)+                             &
     &          Hz(i,j,k)*(2.0_r8-FC(i,k-1)))
            FC(i,k)=cff*Hz(i,j,k+1)
            dTdz(i,j,k)=cff*(6.0_r8*(t(i,j,k+1,Lbck,itemp)-             &
     &                               t(i,j,k  ,Lbck,itemp))-            &
     &                       Hz(i,j,k)*dTdz(i,j,k-1))
# ifdef SALINITY
            dSdz(i,j,k)=cff*(6.0_r8*(t(i,j,k+1,Lbck,isalt)-             &
     &                               t(i,j,k  ,Lbck,isalt))-            &
     &                       Hz(i,j,k)*dSdz(i,j,k-1))
# endif
          END DO
        END DO
        DO i=IstrR,IendR
          dTdz(i,j,N(ng))=0.0_r8
# ifdef SALINITY
          dSdz(i,j,N(ng))=0.0_r8
# endif
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrR,IendR
            dTdz(i,j,k)=dTdz(i,j,k)-FC(i,k)*dTdz(i,j,k+1)
# ifdef SALINITY
            dSdz(i,j,k)=dSdz(i,j,k)-FC(i,k)*dSdz(i,j,k+1)
# endif
          END DO
        END DO
!
!  Shapiro filter of the vertical derivatives.
!
        DO i=IstrR,IendR
          DO order=1,Norder/2
            IF (order.ne.Norder/2) THEN
# ifdef SALINITY
              dSdz_filter(1)=2.0_r8*(dSdz(i,j,1)-dSdz(i,j,2))
              dSdz_filter(N(ng))=2.0_r8*(dSdz(i,j,N(ng)   )-            &
     &                                   dSdz(i,j,N(ng)-1))
# endif
              dTdz_filter(1)=2.0_r8*(dTdz(i,j,1)-dTdz(i,j,2))
              dTdz_filter(N(ng))=2.0_r8*(dTdz(i,j,N(ng)   )-            &
     &                                   dTdz(i,j,N(ng)-1))
            ELSE
# ifdef SALINITY
              dSdz_filter(1)=0.0_r8
              dSdz_filter(N(ng))=0.0_r8
# endif
              dTdz_filter(1)=0.0_r8
              dTdz_filter(N(ng))=0.0_r8
            END IF
            DO k=2,N(ng)-1
# ifdef SALINITY
              dSdz_filter(k)=2.0_r8*dSdz(i,j,k)-                        &
     &                       dSdz(i,j,k-1)-dSdz(i,j,k+1)
# endif
              dTdz_filter(k)=2.0_r8*dTdz(i,j,k)-                        &
     &                       dTdz(i,j,k-1)-dTdz(i,j,k+1)
            END DO
            DO k=1,N(ng)
# ifdef SALINITY
              dSdz(i,j,k)=dSdz(i,j,k)-                                  &
     &                    filter_coef(Norder/2)*dSdz_filter(k)
# endif
              dTdz(i,j,k)=dTdz(i,j,k)-                                  &
     &                    filter_coef(Norder/2)*dTdz_filter(k)
            END DO
          END DO
        END DO
      END DO
!
!  Compute velocity shears dUdz and dVdz.
!
      DO k=0,N(ng)
        IF ((k.eq.0).or.(k.eq.N(ng))) THEN
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              dUdz(i,j,k)=0.0_r8
            END DO
          END DO
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              dVdz(i,j,k)=0.0_r8
            END DO
          END DO
        ELSE
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              cff=1.0_r8/(0.5_r8*(z_r(i-1,j,k+1)-z_r(i-1,j,k)+          &
     &                            z_r(i  ,j,k+1)-z_r(i  ,j,k)))
              dUdz(i,j,k)=cff*(u(i,j,k+1,Lbck)-                         &
     &                         u(i,j,k  ,Lbck))
            END DO
          END DO
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              cff=1.0_r8/(0.5_r8*(z_r(i,j-1,k+1)-z_r(i,j-1,k)+          &
     &                            z_r(i,j  ,k+1)-z_r(i,j  ,k)))
              dVdz(i,j,k)=cff*(v(i,j,k+1,Lbck)-                         &
     &                         v(i,j,k  ,Lbck))
            END DO
          END DO
        END IF
      END DO
!
!  Shapiro filter of dUdz.
!
      DO j=Jstr-1,Jend+1
        DO i=IstrU-1,Iend+1
          DO order=1,Norder/2
            IF (order.ne.Norder/2) THEN
              dUdz_filter(1)=2.0_r8*(dUdz(i,j,1)-dUdz(i,j,2))
              dUdz_filter(N(ng))=2.0_r8*(dUdz(i,j,N(ng)   )-            &
     &                                   dUdz(i,j,N(ng)-1))
            ELSE
              dUdz_filter(1)=0.0_r8
              dUdz_filter(N(ng))=0.0_r8
            END IF
            DO k=2,N(ng)-1
              dUdz_filter(k)=2.0_r8*dUdz(i,j,k)-                        &
     &                       dUdz(i,j,k-1)-dUdz(i,j,k+1)
            END DO
            DO k=1,N(ng)
              dUdz(i,j,k)=dUdz(i,j,k)-                                  &
     &                    filter_coef(Norder/2)*dUdz_filter(k)
            END DO
          END DO
        END DO
      END DO
!
!  Shapiro filter of dVdz.
!
      DO j=JstrV-1,Jend+1
        DO i=Istr-1,Iend+1
          DO order=1,Norder/2
            IF (order.ne.Norder/2) THEN
              dVdz_filter(1)=2.0_r8*(dVdz(i,j,1)-dVdz(i,j,2))
              dVdz_filter(N(ng))=2.0_r8*(dVdz(i,j,N(ng)   )-            &
     &                                   dVdz(i,j,N(ng)-1))
            ELSE
              dVdz_filter(1)=0.0_r8
              dVdz_filter(N(ng))=0.0_r8
            END IF
            DO k=2,N(ng)-1
              dVdz_filter(k)=2.0_r8*dVdz(i,j,k)-                        &
     &                       dVdz(i,j,k-1)-dVdz(i,j,k+1)
            END DO
            DO k=1,N(ng)
              dVdz(i,j,k)=dVdz(i,j,k)-                                  &
     &                    filter_coef(Norder/2)*dVdz_filter(k)
            END DO
          END DO
        END DO
      END DO

# ifdef COMPUTE_MLD
!
!-----------------------------------------------------------------------
!  Compute mixed layer depth (mld) according to criterion
!  on temperature from the work of Kara et al. (2000).
!-----------------------------------------------------------------------
!
!  Criteria on temperature : 0.8
!
!  Start searching for the mixed layer depth at href.
!  CAUTION: href MUST be greater than the minimum depth in the model
!  domain (usually over land points).
!
      T_thvalue=0.8_r8
!     href = -10.0_r8
      href = -9.0_r8
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
!
!  Compute temperature at reference depth, href.
!
          DO k=N(ng)-1,1,-1
            IF ((z_r(i,j,k+1).gt.href).and.(z_r(i,j,k).le.href)) THEN
              khref = k
            END IF
          END DO
          base_reached =.FALSE.
          kref= khref
!                              Search for an uniform temperature region
          DO k=khref-1,1,-1
            T_dep = ABS(t(i,j,k+1,Lbck,itemp)-t(i,j,k,Lbck,itemp))/     &
     &              ABS(z_r(i,j,k+1)-z_r(i,j,k))
            IF (.not.base_reached.and.(T_dep.lt.0.01_r8*T_thvalue)) THEN
              kref = k
            ELSE
              base_reached = .TRUE.
            END IF
          END DO
!
          IF (.not.base_reached) THEN
            kref = khref
          END IF
!
!  Linearly interpolate temperature to href.
!
          IF (kref.eq.khref) THEN
            DO k=1,N(ng)-1
              IF ((href.ge.z_r(i,j,k)).and.(href.lt.z_r(i,j,k+1))) THEN
                fac=(href-z_r(i,j,k))/(z_r(i,j,k+1)-z_r(i,j,k))
                Temp_ref=fac*t(i,j,k+1,Lbck,itemp)+                     &
                         (1.0_r8-fac)*t(i,j,k,Lbck,itemp)
                mld(i,j)=href
              END IF
            END DO
          ELSE
            Temp_ref=t(i,j,kref,Lbck,itemp)
            mld(i,j)=z_r(i,j,kref)
          END IF
!
          ml_reached=.FALSE.
!
          DO k=kref-1,1,-1
            IF (.not.ml_reached.and.                                    &
     &         (ABS(t(i,j,k,Lbck,itemp)-Temp_ref).gt.T_thvalue)) THEN
              T_high=ABS(t(i,j,k+1,Lbck,itemp)-Temp_ref)
              T_low =ABS(t(i,j,k  ,Lbck,itemp)-Temp_ref)
              mld(i,j)=((T_thvalue-T_high)*z_r(i,j,k)+                   &
     &                  (T_low-T_thvalue)*z_r(i,j,k+1))/(T_low-T_high)
              ml_reached =.TRUE.
            END IF
          END DO
!
          IF (.not.ml_reached) THEN
            mld(i,j)=z_r(i,j,1)
          END IF
        END DO
      END DO
# else
      DO j=JstrR,JendR
        DO i=IstrR,IendR
           mld(i,j)=mld_uniform(ng)
        END DO
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Compute the background (prior) state vector standard deviations.
!-----------------------------------------------------------------------
!
!  Temperature and salinity.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          DO k=1,N(ng)
            cff1=0.5_r8*(dTdz(i,j,k-1)+dTdz(i,j,k))
            sigmabT=MIN(ABS(Sigma_dz(isTvar(itemp),ng)*cff1),           &
     &                  Sigma_max(isTvar(itemp),ng))
            IF (z_r(i,j,k).ge.mld(i,j)) THEN
              sigmabT=MAX(sigmabT, Sigma_ml(isTvar(itemp),ng))
            ELSE
              sigmabT=MAX(sigmabT, Sigma_do(isTvar(itemp),ng))
            END IF
            e_t(i,j,k,Lstd,itemp)=sigmabT
# ifdef MASKING
            e_t(i,j,k,Lstd,itemp)=e_t(i,j,k,Lstd,itemp)*rmask(i,j)
# endif
# ifdef SALINITY
            cff2=0.5_r8*(dSdz(i,j,k-1)+dSdz(i,j,k))
            sigmabS=MIN(ABS(Sigma_dz(isTvar(isalt),ng)*cff2),           &
                        Sigma_max(isTvar(isalt),ng))
            IF (z_r(i,j,k).ge.mld(i,j)) THEN
              sigmabS=MAX(sigmabS, Sigma_ml(isTvar(isalt),ng))
            ELSE
              sigmabS=MAX(sigmabS, Sigma_do(isTvar(isalt),ng))
            END IF
            e_t(i,j,k,Lstd,isalt)=sigmabS
#  ifdef MASKING
            e_t(i,j,k,Lstd,isalt)=e_t(i,j,k,Lstd,isalt)*rmask(i,j)
#  endif
# endif
          END DO
        END DO
      END DO
!
!  U-velocity component.
!
      DO j=Jstr-1,Jend+1
        DO i=IstrU-1,Iend+1
          DO k=1,N(ng)
            cff=0.5_r8*(dUdz(i,j,k-1)+dUdz(i,j,k))
            sigmabU=MIN(ABS(Sigma_dz(isUvel,ng)*cff),                   &
     &                  Sigma_max(isUvel,ng))
            IF (z_r(i,j,k).ge.mld(i,j)) THEN
              sigmabU=MAX(sigmabU, Sigma_ml(isUvel,ng))
            ELSE
              sigmabU=MAX(sigmabU, Sigma_do(isUvel,ng))
            END IF
            e_u(i,j,k,Lstd)=sigmabU
# ifdef MASKING
            e_u(i,j,k,Lstd)=e_u(i,j,k,Lstd)*umask(i,j)
# endif
          END DO
        END DO
      END DO
!
!  V-velocity component.
!
      DO j=JstrV-1,Jend+1
        DO i=Istr-1,Iend+1
          DO k=1,N(ng)
            cff=0.5_r8*(dVdz(i,j,k-1)+dVdz(i,j,k))
            sigmabV=MIN(ABS(Sigma_dz(isVvel,ng)*cff),                   &
     &                  Sigma_max(isVvel,ng))
            IF (z_r(i,j,k).ge.mld(i,j)) THEN
              sigmabV=MAX(sigmabV, Sigma_ml(isVvel,ng))
            ELSE
              sigmabV=MAX(sigmabV, Sigma_do(isVvel,ng))
            END IF
            e_v(i,j,k,Lstd)=sigmabV
# ifdef MASKING
            e_v(i,j,k,Lstd)=e_v(i,j,k,Lstd)*vmask(i,j)
# endif
          END DO
        END DO
      END DO
!
!  Vertically integrated velocity components.  Not used, but needed
!  for I/O manipulations in the split schemes.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          e_ubar(i,j,Lstd)=Sigma_max(isUbar,ng)
# ifdef MASKING
          e_ubar(i,j,Lstd)=e_ubar(i,j,Lstd)*umask(i,j)
# endif
        END DO
      END DO

      DO j=Jstr,JendR
        DO i=IstrR,IendR
          e_vbar(i,j,Lstd)=Sigma_max(isVbar,ng)
# ifdef MASKING
          e_vbar(i,j,Lstd)=e_vbar(i,j,Lstd)*vmask(i,j)
# endif
        END DO
      END DO
!
!  Free surface.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          e_zeta(i,j,Lstd)=Sigma_max(isFsur,ng)
# ifdef MASKING
          e_zeta(i,j,Lstd)=e_zeta(i,j,Lstd)*rmask(i,j)
# endif
        END DO
      END DO
!
!  Exchange boundary information.

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          e_t(:,:,:,Lstd,itemp))
# ifdef SALINITY
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          e_t(:,:,:,Lstd,isalt))
# endif

        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          e_u(:,:,:,Lstd))
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          e_v(:,:,:,Lstd))
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          e_ubar(:,:,Lstd))
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          e_vbar(:,:,Lstd))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          e_zeta(:,:,Lstd))
      END IF

# ifdef DISTRIBUTE
!
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    e_t(:,:,:,Lstd,itemp))
#  ifdef SALINITY
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    e_t(:,:,:,Lstd,isalt))
#  endif
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    e_u(:,:,:,Lstd),                              &
     &                    e_v(:,:,:,Lstd))
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    e_ubar(:,:,Lstd),                             &
     &                    e_vbar(:,:,Lstd),                             &
     &                    e_zeta(:,:,Lstd))
# endif
!
      RETURN
      END SUBROUTINE background_std_tile
#endif
      END MODULE background_std_mod
