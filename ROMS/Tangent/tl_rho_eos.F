#include "cppdefs.h"
      MODULE tl_rho_eos_mod
#if defined TANGENT && defined SOLVE3D
!
!svn $Id: tl_rho_eos.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes  "in situ" density and other associated       !
!  quantitites as a function of potential temperature,  salinity,      !
!  and pressure from a polynomial expression (Jackett & McDougall,     !
!  1992). The polynomial expression was found from fitting to 248      !
!  values  in the  oceanographic  ranges of  salinity,  potential      !
!  temperature,  and pressure.  It  assumes no pressure variation      !
!  along geopotential surfaces, that is, depth (meters; negative)      !
!  and pressure (dbar; assumed negative here) are interchangeable.     !
!                                                                      !
!  Check Values: (T=3 C, S=35.5 PSU, Z=-5000 m)                        !
!                                                                      !
!     alpha = 2.1014611551470d-04 (1/Celsius)                          !
!     beta  = 7.2575037309946d-04 (1/PSU)                              !
!     gamma = 3.9684764511766d-06 (1/Pa)                               !
!     den   = 1050.3639165364     (kg/m3)                              !
!     den1  = 1028.2845117925     (kg/m3)                              !
!     sound = 1548.8815240223     (m/s)                                !
!     bulk  = 23786.056026320     (Pa)                                 !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Jackett, D. R. and T. J. McDougall, 1995, Minimal Adjustment of     !
!    Hydrostatic Profiles to Achieve Static Stability, J. of Atmos.    !
!    and Oceanic Techn., vol. 12, pp. 381-389.                         !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: tl_rho_eos
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_rho_eos (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
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
      CALL wclock_on (ng, iTLM, 14)
# endif
      CALL tl_rho_eos_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng),                                   &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
# ifdef VAR_RHO_2D
     &                     GRID(ng) % Hz,                               &
     &                     GRID(ng) % tl_Hz,                            &
# endif
     &                     GRID(ng) % z_r,                              &
     &                     GRID(ng) % tl_z_r,                           &
     &                     GRID(ng) % z_w,                              &
     &                     GRID(ng) % tl_z_w,                           &
     &                     OCEAN(ng) % t,                               &
     &                     OCEAN(ng) % tl_t,                            &
# ifdef VAR_RHO_2D
     &                     COUPLING(ng) % rhoA,                         &
     &                     COUPLING(ng) % tl_rhoA,                      &
     &                     COUPLING(ng) % rhoS,                         &
     &                     COUPLING(ng) % tl_rhoS,                      &
# endif
# ifdef BV_FREQUENCY_NOT_YET
     &                     MIXING(ng) % tl_bvf,                         &
# endif
# if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
     defined BULK_FLUXES
     &                     MIXING(ng) % alpha,                          &
     &                     MIXING(ng) % tl_alpha,                       &
     &                     MIXING(ng) % beta,                           &
     &                     MIXING(ng) % tl_beta,                        &
#  ifdef LMD_DDMIX_NOT_YET
     &                     MIXING(ng) % tl_alfaobeta,                   &
#  endif
# endif
# if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
     &                     OCEAN(ng) % tl_pden,                         &
# endif
     &                     OCEAN(ng) % rho,                             &
     &                     OCEAN(ng) % tl_rho)
# ifdef PROFILE
      CALL wclock_off (ng, iTLM, 14)
# endif

      RETURN
      END SUBROUTINE tl_rho_eos

# ifdef NONLIN_EOS
!
!***********************************************************************
      SUBROUTINE tl_rho_eos_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nrhs,                                 &
#  ifdef MASKING
     &                            rmask,                                &
#  endif
#  ifdef VAR_RHO_2D
     &                            Hz, tl_Hz,                            &
#  endif
     &                            z_r, tl_z_r,                          &
     &                            z_w, tl_z_w,                          &
     &                            t, tl_t,                              &
#  ifdef VAR_RHO_2D
     &                            rhoA, tl_rhoA,                        &
     &                            rhoS, tl_rhoS,                        &
#  endif
#  ifdef BV_FREQUENCY_NOT_YET
     &                            tl_bvf,                               &
#  endif
#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
      defined BULK_FLUXES
     &                            alpha, tl_alpha,                      &
     &                            beta, tl_beta,                        &
#   ifdef LMD_DDMIX_NOT_YET
     &                            tl_alfaobeta,                         &
#   endif
#  endif
#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
     &                            tl_pden,                              &
#  endif
     &                            rho, tl_rho)
!***********************************************************************
!
      USE mod_param
      USE mod_eoscoef
      USE mod_scalars
#  ifdef SEDIMENT_NOT_YET
      USE mod_sediment
#  endif
!
      USE exchange_2d_mod
      USE exchange_3d_mod
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
#   ifdef VAR_RHO_2D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
#   endif
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
#   ifdef VAR_RHO_2D
      real(r8), intent(in) :: tl_Hz(LBi:,LBj:,:)
#   endif
      real(r8), intent(in) :: tl_z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: tl_t(LBi:,LBj:,:,:,:)
#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
       defined BULK_FLUXES
      real(r8), intent(inout) :: alpha(LBi:,LBj:)
      real(r8), intent(inout) :: beta(LBi:,LBj:)
#   endif
#   ifdef VAR_RHO_2D
      real(r8), intent(out) :: rhoA(LBi:,LBj:)
      real(r8), intent(out) :: rhoS(LBi:,LBj:)
      real(r8), intent(out) :: tl_rhoA(LBi:,LBj:)
      real(r8), intent(out) :: tl_rhoS(LBi:,LBj:)
#   endif
#   ifdef BV_FREQUENCY_NOT_YET
      real(r8), intent(out) :: tl_bvf(LBi:,LBj:,0:)
#   endif
#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
       defined BULK_FLUXES
      real(r8), intent(out) :: tl_alpha(LBi:,LBj:)
      real(r8), intent(out) :: tl_beta(LBi:,LBj:)
#    ifdef LMD_DDMIX_NOT_YET
      real(r8), intent(out) :: tl_alfaobeta(LBi:,LBj:,0:)
#    endif
#   endif
#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
      real(r8), intent(out) :: tl_pden(LBi:,LBj:,:)
#   endif
      real(r8), intent(out) :: rho(LBi:,LBj:,:)
      real(r8), intent(out) :: tl_rho(LBi:,LBj:,:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef VAR_RHO_2D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#   ifdef VAR_RHO_2D
      real(r8), intent(in) :: tl_Hz(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(in) :: tl_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
       defined BULK_FLUXES
      real(r8), intent(inout) :: alpha(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: beta(LBi:UBi,LBj:UBj)
#   endif
#   ifdef VAR_RHO_2D
      real(r8), intent(out) :: rhoA(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rhoS(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tl_rhoA(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tl_rhoS(LBi:UBi,LBj:UBj)
#   endif
#   ifdef BV_FREQUENCY_NOT_YET
      real(r8), intent(out) :: tl_bvf(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
       defined BULK_FLUXES
      real(r8), intent(out) :: tl_alpha(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tl_beta(LBi:UBi,LBj:UBj)
#    ifdef LMD_DDMIX_NOT_YET
      real(r8), intent(out) :: tl_alfaobeta(LBi:UBi,LBj:UBj,0:N(ng))
#    endif
#   endif
#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
      real(r8), intent(out) :: tl_pden(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(out) :: rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: tl_rho(LBi:UBi,LBj:UBj,N(ng))
#  endif
!
!  Local variable declarations.
!
      integer :: i, ised, itrc, j, k

      real(r8) :: SedDen, Tp, Tpr10, Ts, Tt, sqrtTs
      real(r8) :: tl_SedDen, tl_Tp, tl_Tpr10, tl_Ts, tl_Tt, tl_sqrtTs
#  ifdef BV_FREQUENCY_NOT_YET
      real(r8) :: bulk_dn, bulk_up, den_dn, den_up
      real(r8) :: tl_bulk_dn, tl_bulk_up, tl_den_dn, tl_den_up
#  endif
      real(r8) :: cff, cff1, cff2, cff3
      real(r8) :: tl_cff, tl_cff1, tl_cff2, tl_cff3

      real(r8), dimension(0:9) :: C
      real(r8), dimension(0:9) :: tl_C
#  ifdef EOS_TDERIVATIVE
      real(r8), dimension(0:9) :: dCdT(0:9)
      real(r8), dimension(0:9) :: tl_dCdT(0:9)
      real(r8), dimension(0:9) :: d2Cd2T(0:9)

      real(r8), dimension(IminS:ImaxS,N(ng)) :: DbulkDS
      real(r8), dimension(IminS:ImaxS,N(ng)) :: DbulkDT
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Dden1DS
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Dden1DT
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Scof
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Tcof
      real(r8), dimension(IminS:ImaxS,N(ng)) :: wrk

      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_DbulkDS
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_DbulkDT
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_Dden1DS
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_Dden1DT
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_Scof
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_Tcof
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_wrk
#  endif
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk0
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk1
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: den
      real(r8), dimension(IminS:ImaxS,N(ng)) :: den1

      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_bulk
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_bulk0
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_bulk1
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_bulk2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_den
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_den1

#  include "set_bounds.h"
!
!=======================================================================
!  Nonlinear equation of state.  Notice that this equation of state
!  is only valid for potential temperature range of -2C to 40C and
!  a salinity range of 0 PSU to 42 PSU.
!=======================================================================
!
      DO j=JstrT,JendT
        DO k=1,N(ng)
          DO i=IstrT,IendT
!
!  Check temperature and salinity lower values. Assign depth to the
!  pressure.
!
            Tt=MAX(-2.0_r8,t(i,j,k,nrhs,itemp))
            tl_Tt=(0.5_r8-SIGN(0.5_r8,-2.0_r8-                          &
     &                         t(i,j,k,nrhs,itemp)))*                   &
     &            tl_t(i,j,k,nrhs,itemp)
#  ifdef SALINITY
            Ts=MAX(0.0_r8,t(i,j,k,nrhs,isalt))
            tl_Ts=(0.5_r8-SIGN(0.5_r8,-t(i,j,k,nrhs,isalt)))*           &
     &            tl_t(i,j,k,nrhs,isalt)
            sqrtTs=SQRT(Ts)
            IF (Ts.ne.0.0_r8) THEN
              tl_sqrtTs=0.5_r8*tl_Ts/SQRT(Ts)
            ELSE
              tl_sqrtTs=0.0_r8
            END IF
#  else
            Ts=0.0_r8
            tl_Ts=0.0_r8
            sqrtTs=0.0_r8
            tl_sqrtTs=0.0_r8
#  endif
            Tp=z_r(i,j,k)
            tl_Tp=tl_z_r(i,j,k)
            Tpr10=0.1_r8*Tp
            tl_Tpr10=0.1_r8*tl_Tp
!
!-----------------------------------------------------------------------
!  Compute BASIC STATE and tangent linear density (kg/m3) at standard
!  one atmosphere pressure.
!-----------------------------------------------------------------------
!
            C(0)=Q00+Tt*(Q01+Tt*(Q02+Tt*(Q03+Tt*(Q04+Tt*Q05))))
            C(1)=U00+Tt*(U01+Tt*(U02+Tt*(U03+Tt*U04)))
            C(2)=V00+Tt*(V01+Tt*V02)
#  ifdef EOS_TDERIVATIVE
!
            dCdT(0)=Q01+Tt*(2.0_r8*Q02+Tt*(3.0_r8*Q03+Tt*(4.0_r8*Q04+   &
     &                      Tt*5.0_r8*Q05)))
            dCdT(1)=U01+Tt*(2.0_r8*U02+Tt*(3.0_r8*U03+Tt*4.0_r8*U04))
            dCdT(2)=V01+Tt*2.0_r8*V02
#  endif
            tl_C(0)=tl_Tt*dCdT(0)
            tl_C(1)=tl_Tt*dCdT(1)
            tl_C(2)=tl_Tt*dCdT(2)
!
            den1(i,k)=C(0)+Ts*(C(1)+sqrtTs*C(2)+Ts*W00)
            tl_den1(i,k)=tl_C(0)+                                       &
     &                   tl_Ts*(C(1)+sqrtTs*C(2)+Ts*W00)+               &
     &                   Ts*(tl_C(1)+tl_sqrtTs*C(2)+                    &
     &                       sqrtTs*tl_C(2)+tl_Ts*W00)

#  ifdef EOS_TDERIVATIVE
!
!  Compute d(den1)/d(S) and d(den1)/d(T) derivatives used in the
!  computation of thermal expansion and saline contraction
!  coefficients.
!
            d2Cd2T(0)=2.0_r8*Q02+Tt*(6.0_r8*Q03+Tt*(12.0_r8*Q04+        &
     &                               Tt*20.0_r8*Q05))
            d2Cd2T(1)=2.0_r8*U02+Tt*(6.0_r8*U03+Tt*12.0_r8*U04)
            d2Cd2T(2)=2.0_r8*V02
!
            tl_dCdT(0)=tl_Tt*d2Cd2T(0)
            tl_dCdT(1)=tl_Tt*d2Cd2T(1)
            tl_dCdT(2)=tl_Tt*d2Cd2T(2)
!
            Dden1DS(i,k)=C(1)+1.5_r8*C(2)*sqrtTs+2.0_r8*W00*Ts
            Dden1DT(i,k)=dCdT(0)+Ts*(dCdT(1)+sqrtTs*dCdT(2))
!
            tl_Dden1DS(i,k)=tl_C(1)+                                    &
     &                      1.5_r8*(tl_C(2)*sqrtTs+                     &
     &                              C(2)*tl_sqrtTs)+                    &
     &                      2.0_r8*W00*tl_Ts
            tl_Dden1DT(i,k)=tl_dCdT(0)+                                 &
     &                      tl_Ts*(dCdT(1)+sqrtTs*dCdT(2))+             &
     &                      Ts*(tl_dCdT(1)+tl_sqrtTs*dCdT(2)+           &
     &                                     sqrtTs*tl_dCdT(2))
#  endif
!
!-----------------------------------------------------------------------
!  Compute BASIC STATE and tangent linear secant bulk modulus.
!-----------------------------------------------------------------------
!
            C(3)=A00+Tt*(A01+Tt*(A02+Tt*(A03+Tt*A04)))
            C(4)=B00+Tt*(B01+Tt*(B02+Tt*B03))
            C(5)=D00+Tt*(D01+Tt*D02)
            C(6)=E00+Tt*(E01+Tt*(E02+Tt*E03))
            C(7)=F00+Tt*(F01+Tt*F02)
            C(8)=G01+Tt*(G02+Tt*G03)
            C(9)=H00+Tt*(H01+Tt*H02)
#  ifdef EOS_TDERIVATIVE
!
            dCdT(3)=A01+Tt*(2.0_r8*A02+Tt*(3.0_r8*A03+Tt*4.0_r8*A04))
            dCdT(4)=B01+Tt*(2.0_r8*B02+Tt*3.0_r8*B03)
            dCdT(5)=D01+Tt*2.0_r8*D02
            dCdT(6)=E01+Tt*(2.0_r8*E02+Tt*3.0_r8*E03)
            dCdT(7)=F01+Tt*2.0_r8*F02
            dCdT(8)=G02+Tt*2.0_r8*G03
            dCdT(9)=H01+Tt*2.0_r8*H02
#  endif
!
            tl_C(3)=tl_Tt*dCdT(3)
            tl_C(4)=tl_Tt*dCdT(4)
            tl_C(5)=tl_Tt*dCdT(5)
            tl_C(6)=tl_Tt*dCdT(6)
            tl_C(7)=tl_Tt*dCdT(7)
            tl_C(8)=tl_Tt*dCdT(8)
            tl_C(9)=tl_Tt*dCdT(9)
!
            bulk0(i,k)=C(3)+Ts*(C(4)+sqrtTs*C(5))
            bulk1(i,k)=C(6)+Ts*(C(7)+sqrtTs*G00)
            bulk2(i,k)=C(8)+Ts*C(9)
            bulk (i,k)=bulk0(i,k)-Tp*(bulk1(i,k)-Tp*bulk2(i,k))
!
            tl_bulk0(i,k)=tl_C(3)+                                      &
     &                    tl_Ts*(C(4)+sqrtTs*C(5))+                     &
     &                    Ts*(tl_C(4)+tl_sqrtTs*C(5)+                   &
     &                        sqrtTs*tl_C(5))
            tl_bulk1(i,k)=tl_C(6)+                                      &
     &                    tl_Ts*(C(7)+sqrtTs*G00)+                      &
     &                    Ts*(tl_C(7)+tl_sqrtTs*G00)
            tl_bulk2(i,k)=tl_C(8)+tl_Ts*C(9)+Ts*tl_C(9)
            tl_bulk (i,k)=tl_bulk0(i,k)-                                &
     &                     tl_Tp*(bulk1(i,k)-Tp*bulk2(i,k))-            &
     &                     Tp*(tl_bulk1(i,k)-                           &
     &                         tl_Tp*bulk2(i,k)-                        &
     &                         Tp*tl_bulk2(i,k))

#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
      defined BULK_FLUXES
!
!  Compute d(bulk)/d(S) and d(bulk)/d(T) derivatives used
!  in the computation of thermal expansion and saline contraction
!  coefficients.
!
            d2Cd2T(3)=2.0_r8*A02+Tt*(6.0_r8*A03+Tt*12.0_r8*A04)
            d2Cd2T(4)=2.0_r8*B02+Tt*6.0_r8*B03
            d2Cd2T(5)=2.0_r8*D02
            d2Cd2T(6)=2.0_r8*E02+Tt*6.0_r8*E03
            d2Cd2T(7)=2.0_r8*F02
            d2Cd2T(8)=2.0_r8*G03
            d2Cd2T(9)=2.0_r8*H02
!
            tl_dCdT(3)=tl_Tt*d2Cd2T(3)
            tl_dCdT(4)=tl_Tt*d2Cd2T(4)
            tl_dCdT(5)=tl_Tt*d2Cd2T(5)
            tl_dCdT(6)=tl_Tt*d2Cd2T(6)
            tl_dCdT(7)=tl_Tt*d2Cd2T(7)
            tl_dCdT(8)=tl_Tt*d2Cd2T(8)
            tl_dCdT(9)=tl_Tt*d2Cd2T(9)
!
            DbulkDS(i,k)=C(4)+sqrtTs*1.5_r8*C(5)-                       &
     &                   Tp*(C(7)+sqrtTs*1.5_r8*G00-Tp*C(9))
            DbulkDT(i,k)=dCdT(3)+Ts*(dCdT(4)+sqrtTs*dCdT(5))-           &
     &                   Tp*(dCdT(6)+Ts*dCdT(7)-                        &
     &                       Tp*(dCdT(8)+Ts*dCdT(9)))
!
            tl_DbulkDS(i,k)=tl_C(4)+                                    &
     &                      1.5_r8*(tl_sqrtTs*C(5)+                     &
     &                              sqrtTs*tl_C(5))-                    &
     &                      tl_Tp*(C(7)+sqrtTs*1.5_r8*G00-              &
     &                             Tp*C(9))-                            &
     &                      Tp*(tl_C(7)+tl_sqrtTs*1.5_r8*G00-           &
     &                          tl_Tp*C(9)-Tp*tl_C(9))
            tl_DbulkDT(i,k)=tl_dCdT(3)+                                 &
     &                      tl_Ts*(dCdT(4)+sqrtTs*dCdT(5))+             &
     &                      Ts*(tl_dCdT(4)+tl_sqrtTs*dCdT(5)+           &
     &                                     sqrtTs*tl_dCdT(5))-          &
     &                      tl_Tp*(dCdT(6)+Ts*dCdT(7)-                  &
     &                             Tp*(dCdT(8)+Ts*dCdT(9)))-            &
     &                      Tp*(tl_dCdT(6)+tl_Ts*dCdT(7)+Ts*tl_dCdT(7)- &
     &                          tl_Tp*(dCdT(8)+Ts*dCdT(9))-             &
     &                          Tp*(tl_dCdT(8)+tl_Ts*dCdT(9)+           &
     &                                         Ts*tl_dCdT(9)))
#  endif
!
!-----------------------------------------------------------------------
!  Compute local "in situ" density anomaly (kg/m3 - 1000).
!-----------------------------------------------------------------------
!
            cff=1.0_r8/(bulk(i,k)+Tpr10)
            tl_cff=-cff*cff*(tl_bulk(i,k)+tl_Tpr10)
            den(i,k)=den1(i,k)*bulk(i,k)*cff
            tl_den(i,k)=tl_den1(i,k)*bulk(i,k)*cff+                     &
     &                  den1(i,k)*(tl_bulk(i,k)*cff+                    &
     &                             bulk(i,k)*tl_cff)
#  if defined SEDIMENT_NOT_YET && defined SED_DENS_NOT_YET
            SedDen=0.0_r8
            tl_SedDen=0.0_r8
            DO ised=1,NST
              itrc=idsed(ised)
              cff1=1.0_r8/Srho(ised,ng)
              SedDen=SedDen+                                            &
     &               t(i,j,k,nrhs,itrc)*                                &
     &               (Srho(ised,ng)-den(i,k))*cff1
              tl_SedDen=tl_SedDen+                                      &
     &                  (tl_t(i,j,k,nrhs,itrc)*                         &
     &                   (Srho(ised,ng)-den(i,k))-                      &
     &                   t(i,j,k,nrhs,itrc)*                            &
     &                   tl_den(i,k))*cff1
            END DO
            den(i,k)=den(i,k)+SedDen
            tl_den(i,k)=tl_den(i,k)+tl_SedDen
#  endif
            den(i,k)=den(i,k)-1000.0_r8
#  ifdef MASKING
            den(i,k)=den(i,k)*rmask(i,j)
            tl_den(i,k)=tl_den(i,k)*rmask(i,j)
#  endif
          END DO
        END DO

#  ifdef VAR_RHO_2D
!
!-----------------------------------------------------------------------
!  Compute vertical averaged density (rhoA) and density perturbation
!  (rhoS) used in barotropic pressure gradient.
!-----------------------------------------------------------------------
!
        DO i=IstrT,IendT
          cff1=den(i,N(ng))*Hz(i,j,N(ng))
          tl_cff1=tl_den(i,N(ng))*Hz(i,j,N(ng))+                        &
     &            den(i,N(ng))*tl_Hz(i,j,N(ng))
          rhoS(i,j)=0.5_r8*cff1*Hz(i,j,N(ng))
          tl_rhoS(i,j)=0.5_r8*(tl_cff1*Hz(i,j,N(ng))+                   &
     &                         cff1*tl_Hz(i,j,N(ng)))
          rhoA(i,j)=cff1
          tl_rhoA(i,j)=tl_cff1
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrT,IendT
            cff1=den(i,k)*Hz(i,j,k)
            tl_cff1=tl_den(i,k)*Hz(i,j,k)+                              &
     &              den(i,k)*tl_Hz(i,j,k)
            rhoS(i,j)=rhoS(i,j)+Hz(i,j,k)*(rhoA(i,j)+0.5_r8*cff1)
            tl_rhoS(i,j)=tl_rhoS(i,j)+                                  &
     &                   tl_Hz(i,j,k)*(rhoA(i,j)+0.5_r8*cff1)+          &
     &                   Hz(i,j,k)*(tl_rhoA(i,j)+0.5_r8*tl_cff1)
            rhoA(i,j)=rhoA(i,j)+cff1
            tl_rhoA(i,j)=tl_rhoA(i,j)+tl_cff1
          END DO
        END DO
        cff2=1.0_r8/rho0
        DO i=IstrT,IendT
          cff1=1.0_r8/(z_w(i,j,N(ng))-z_w(i,j,0))
          tl_cff1=-cff1*cff1*(tl_z_w(i,j,N(ng))-tl_z_w(i,j,0))
!
!  Here we reverse the order of the NL and TL operations since an
!  intermeridiate value of rhoA and rhoS is needed because they are
!  recursive.
!
          tl_rhoA(i,j)=cff2*(tl_cff1*rhoA(i,j)+cff1*tl_rhoA(i,j))
          rhoA(i,j)=cff2*cff1*rhoA(i,j)
          tl_rhoS(i,j)=2.0_r8*cff2*                                     &
     &                 cff1*(2.0_r8*tl_cff1*rhoS(i,j)+                  &
     &                       cff1*tl_rhoS(i,j))
          rhoS(i,j)=2.0_r8*cff1*cff1*cff2*rhoS(i,j)
        END DO
#  endif

#  if defined BV_FREQUENCY_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute Brunt-Vaisala frequency (1/s2) at horizontal RHO-points
!  and vertical W-points:
!
!                  bvf = - g/rho d(rho)/d(z).
!
!  The density anomaly difference is computed by lowering/rising the
!  water parcel above/below adiabatically at W-point depth "z_w".
!-----------------------------------------------------------------------
!
        DO k=1,N(ng)-1
          DO i=IstrT,IendT
            bulk_up=bulk0(i,k+1)-                                       &
     &              z_w(i,j,k)*(bulk1(i,k+1)-                           &
     &                          bulk2(i,k+1)*z_w(i,j,k))
            tl_bulk_up=tl_bulk0(i,k+1)-                                 &
     &                 tl_z_w(i,j,k)*(bulk1(i,k+1)-                     &
     &                                bulk2(i,k+1)*z_w(i,j,k))-         &
     &                 z_w(i,j,k)*(tl_bulk1(i,k+1)-                     &
     &                             tl_bulk2(i,k+1)*z_w(i,j,k)-          &
     &                             bulk2(i,k+1)*tl_z_w(i,j,k))
            bulk_dn=bulk0(i,k  )-                                       &
     &              z_w(i,j,k)*(bulk1(i,k  )-                           &
     &                          bulk2(i,k  )*z_w(i,j,k))
            tl_bulk_dn=tl_bulk0(i,k  )-                                 &
     &                 tl_z_w(i,j,k)*(bulk1(i,k  )-                     &
     &                                bulk2(i,k  )*z_w(i,j,k))-         &
     &                 z_w(i,j,k)*(tl_bulk1(i,k  )-                     &
     &                             tl_bulk2(i,k  )*z_w(i,j,k)-          &
     &                             bulk2(i,k  )*tl_z_w(i,j,k))
            cff1=1.0_r8/(bulk_up+0.1_r8*z_w(i,j,k))
            cff2=1.0_r8/(bulk_dn+0.1_r8*z_w(i,j,k))
            tl_cff1=-cff1*cff1*(tl_bulk_up+0.1_r8*tl_z_w(i,j,k))
            tl_cff2=-cff2*cff2*(tl_bulk_dn+0.1_r8*tl_z_w(i,j,k))
            den_up=cff1*(den1(i,k+1)*bulk_up)
            den_dn=cff2*(den1(i,k  )*bulk_dn)
            tl_den_up=tl_cff1*(den1(i,k+1)*bulk_up)+                    &
     &                cff1*(tl_den1(i,k+1)*bulk_up+                     &
     &                      den1(i,k+1)*tl_bulk_up)
            tl_den_dn=tl_cff2*(den1(i,k  )*bulk_dn)+                    &
     &                cff2*(tl_den1(i,k  )*bulk_dn+                     &
     &                      den1(i,k  )*tl_bulk_dn)
!>          bvf(i,j,k)=-g*(den_up-den_dn)/                              &
!>   &                 (0.5_r8*(den_up+den_dn)*                         &
!>   &                  (z_r(i,j,k+1)-z_r(i,j,k)))
!>
            cff3=1.0_r8/(0.5_r8*(den_up+den_dn)*                        &
     &                  (z_r(i,j,k+1)-z_r(i,j,k)))
            tl_cff3=-cff3*cff3*                                         &
     &               0.5_r8*((tl_den_up+tl_den_dn)*                     &
     &                       (z_r(i,j,k+1)-z_r(i,j,k))+                 &
     &                       (den_up+den_dn)*                           &
     &                       (tl_z_r(i,j,k+1)-tl_z_r(i,j,k)))
            tl_bvf(i,j,k)=-g*((tl_den_up-tl_den_dn)*cff3+               &
     &                        (den_up-den_dn)*tl_cff3)
          END DO
        END DO
        DO i=IstrT,IendT
!>        bvf(i,j,0)=0.0_r8
!>
          tl_bvf(i,j,0)=0.0_r8
!>        bvf(i,j,N(ng))=0.0_r8
!>
          tl_bvf(i,j,N(ng))=0.0_r8
        END DO
#  endif

#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
      defined BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Compute thermal expansion (1/Celsius) and saline contraction
!  (1/PSU) coefficients.
!-----------------------------------------------------------------------
!
#   ifdef LMD_DDMIX_NOT_YET
        DO k=1,N(ng)
#   else
        DO k=N(ng),N(ng)
#   endif
          DO i=IstrT,IendT
            Tpr10=0.1_r8*z_r(i,j,k)
            tl_Tpr10=0.1_r8*tl_z_r(i,j,k)
!
!  Compute thermal expansion and saline contraction coefficients.
!
            cff=bulk(i,k)+Tpr10
            tl_cff=tl_bulk(i,k)+tl_Tpr10
            cff1=Tpr10*den1(i,k)
            tl_cff1=tl_Tpr10*den1(i,k)+Tpr10*tl_den1(i,k)
            cff2=bulk(i,k)*cff
            tl_cff2=tl_bulk(i,k)*cff+bulk(i,k)*tl_cff
            wrk(i,k)=(den(i,k)+1000.0_r8)*cff*cff
            tl_wrk(i,k)=cff*(cff*tl_den(i,k)+                           &
     &                       2.0_r8*tl_cff*(den(i,k)+1000.0_r8))
            Tcof(i,k)=-(DbulkDT(i,k)*cff1+                              &
     &                  Dden1DT(i,k)*cff2)
            tl_Tcof(i,k)=-(tl_DbulkDT(i,k)*cff1+                        &
     &                     DbulkDT(i,k)*tl_cff1+                        &
     &                     tl_Dden1DT(i,k)*cff2+                        &
     &                     Dden1DT(i,k)*tl_cff2)
            Scof(i,k)= (DbulkDS(i,k)*cff1+                              &
     &                  Dden1DS(i,k)*cff2)
            tl_Scof(i,k)= (tl_DbulkDS(i,k)*cff1+                        &
     &                     DbulkDS(i,k)*tl_cff1+                        &
     &                     tl_Dden1DS(i,k)*cff2+                        &
     &                     Dden1DS(i,k)*tl_cff2)
#   ifdef LMD_DDMIX_NOT_YET
!>          alfaobeta(i,j,k)=Tcof(i,k)/Scof(i,k)
!>
            tl_alfaobeta(i,j,k)=(tl_Tcof(i,k)*Scof(i,k)-                &
     &                           Tcof(i,k)*tl_Scof(i,k))/               &
     &                           (Scof(i,k)*Scof(i,k))
#   endif
          END DO
          IF (k.eq.N(ng)) THEN
            DO i=IstrT,IendT
              cff=1.0_r8/wrk(i,N(ng))
              tl_cff=-cff*cff*tl_wrk(i,N(ng))
              alpha(i,j)=cff*Tcof(i,N(ng))
              tl_alpha(i,j)=tl_cff*Tcof(i,N(ng))+cff*tl_Tcof(i,N(ng))
              beta (i,j)=cff*Scof(i,N(ng))
              tl_beta (i,j)=tl_cff*Scof(i,N(ng))+cff*tl_Scof(i,N(ng))
            END DO
          END IF
        END DO
#  endif
!
!-----------------------------------------------------------------------
!  Load "in situ" density anomaly (kg/m3 - 1000) and potential
!  density anomaly (kg/m3 - 1000) referenced to the surface into global
!  arrays. Notice that this is done in a separate (i,k) DO-loops to
!  facilitate the adjoint.
!-----------------------------------------------------------------------
!
        DO k=1,N(ng)
          DO i=IstrT,IendT
            rho(i,j,k)=den(i,k)
            tl_rho(i,j,k)=tl_den(i,k)
#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
!>          pden(i,j,k)=(den1(i,k)-1000.0_r8)      ! This gives a fatal
!>                                                 ! result in 4D-Var
            tl_pden(i,j,k)=tl_den1(i,k)            ! posterior error...
#   ifdef MASKING
!>          pden(i,j,k)=pden(i,k)*rmask(i,j)
!>
            tl_pden(i,j,k)=tl_pden(i,k)*rmask(i,j)
#   endif
#  endif
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          rho)
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          tl_rho)

#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
!>      CALL exchange_r3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          pden)
!>
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          tl_pden)
#  endif

#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
      defined BULK_FLUXES_NOT_YET
#   ifdef LMD_DDMIX_NOT_YET
!>      CALL exchange_w3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
!>   &                          alfaobeta)
!>
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          tl_alfaobeta)
#   endif
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          alpha)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_alpha)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          beta)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_beta)
#  endif

#  ifdef VAR_RHO_2D
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rhoA)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_rhoA)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rhoS)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_rhoS)
#  endif

#  ifdef BV_FREQUENCY_NOT_YET
!>      CALL exchange_w3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
!>   &                          bvf)
!>
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          tl_bvf)
#  endif
      END IF

#  ifdef DISTRIBUTE
!
      CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rho)
      CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_rho)

#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
!>    CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    pden)
!>
      CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_pden)
#   endif

#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
       defined BULK_FLUXES
#    ifdef LMD_DDMIX_NOT_YET
!>    CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    alfaobeta)
!>
      CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_alfaobeta)
#    endif
      CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    alpha, beta)
      CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_alpha, tl_beta)
#   endif

#   ifdef VAR_RHO_2D
      CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rhoA, rhoS)
      CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_rhoA, tl_rhoS)
#   endif

#   ifdef BV_FREQUENCY_NOT_YET
!>    CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    bvf)
!>
      CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_bvf)
#   endif
#  endif

      RETURN
      END SUBROUTINE tl_rho_eos_tile
# endif

# ifndef NONLIN_EOS
!
!***********************************************************************
      SUBROUTINE tl_rho_eos_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nrhs,                                 &
#  ifdef MASKING
     &                            rmask,                                &
#  endif
#  ifdef VAR_RHO_2D
     &                            Hz, tl_Hz,                            &
#  endif
     &                            z_r, tl_z_r,                          &
     &                            z_w, tl_z_w,                          &
     &                            t, tl_t,                              &
#  ifdef VAR_RHO_2D
     &                            rhoA, tl_rhoA,                        &
     &                            rhoS, tl_rhoS,                        &
#  endif
#  ifdef BV_FREQUENCY_NOT_YET
     &                            tl_bvf,                               &
#  endif
#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
      defined BULK_FLUXES
     &                            alpha, tl_alpha,                      &
     &                            beta,  tl_beta,                       &
#   ifdef LMD_DDMIX_NOT_YET
     &                            tl_alfaobeta,                         &
#   endif
#  endif
#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
     &                            tl_pden,                              &
#  endif
     &                            rho, tl_rho)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
#  ifdef SEDIMENT_NOT_YET
      USE mod_sediment
#  endif
!
      USE exchange_2d_mod
      USE exchange_3d_mod
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
#   ifdef VAR_RHO_2D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
#   endif
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
#   ifdef VAR_RHO_2D
      real(r8), intent(in) :: tl_Hz(LBi:,LBj:,:)
#   endif
      real(r8), intent(in) :: tl_z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: tl_t(LBi:,LBj:,:,:,:)
#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
       defined BULK_FLUXES
      real(r8), intent(inout) :: alpha(LBi:,LBj:)
      real(r8), intent(inout) :: beta(LBi:,LBj:)
#   endif
#   ifdef VAR_RHO_2D
      real(r8), intent(out) :: rhoA(LBi:,LBj:)
      real(r8), intent(out) :: rhoS(LBi:,LBj:)
      real(r8), intent(out) :: tl_rhoA(LBi:,LBj:)
      real(r8), intent(out) :: tl_rhoS(LBi:,LBj:)
#   endif
#   ifdef BV_FREQUENCY_NOT_YET
      real(r8), intent(out) :: tl_bvf(LBi:,LBj:,0:)
#   endif
#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
       defined BULK_FLUXES
      real(r8), intent(out) :: tl_alpha(LBi:,LBj:)
      real(r8), intent(out) :: tl_beta(LBi:,LBj:)
#    ifdef LMD_DDMIX_NOT_YET
      real(r8), intent(out) :: tl_alfaobeta(LBi:,LBj:,0:)
#    endif
#   endif
#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
      real(r8), intent(out) :: tl_pden(LBi:,LBj:,:)
#   endif
      real(r8), intent(out) :: rho(LBi:,LBj:,:)
      real(r8), intent(out) :: tl_rho(LBi:,LBj:,:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef VAR_RHO_2D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#   ifdef VAR_RHO_2D
      real(r8), intent(in) :: tl_Hz(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(in) :: tl_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
       defined BULK_FLUXES
      real(r8), intent(inout) :: alpha(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: beta(LBi:UBi,LBj:UBj)
#   endif
#   ifdef VAR_RHO_2D
      real(r8), intent(out) :: rhoA(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rhoS(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tl_rhoA(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tl_rhoS(LBi:UBi,LBj:UBj)
#   endif
#   ifdef BV_FREQUENCY_NOT_YET
      real(r8), intent(out) :: tl_bvf(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
       defined BULK_FLUXES
      real(r8), intent(out) :: tl_alpha(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tl_beta(LBi:UBi,LBj:UBj)
#    ifdef LMD_DDMIX_NOT_YET
      real(r8), intent(out) :: tl_alfaobeta(LBi:UBi,LBj:UBj,0:N(ng))
#    endif
#   endif
#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
      real(r8), intent(out) :: tl_pden(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(out) :: rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: tl_rho(LBi:UBi,LBj:UBj,N(ng))
#  endif
!
!  Local variable declarations.
!
      integer :: i, ised, itrc, j, k

      real(r8) :: SedDen, cff, cff1, cff2
      real(r8) :: tl_SedDen, tl_cff, tl_cff1

#  include "set_bounds.h"
!
!=======================================================================
!  Tangent linear, linear equation of state.
!=======================================================================
!
!-----------------------------------------------------------------------
!  Compute "in situ" density anomaly (kg/m3 - 1000) using the linear
!  equation of state.
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO k=1,N(ng)
          DO i=IstrT,IendT
            rho(i,j,k)=R0(ng)-                                          &
     &                 R0(ng)*Tcoef(ng)*(t(i,j,k,nrhs,itemp)-T0(ng))
            tl_rho(i,j,k)=-R0(ng)*Tcoef(ng)*tl_t(i,j,k,nrhs,itemp)
#  ifdef SALINITY
            rho(i,j,k)=rho(i,j,k)+                                      &
     &                 R0(ng)*Scoef(ng)*(t(i,j,k,nrhs,isalt)-S0(ng))
            tl_rho(i,j,k)=tl_rho(i,j,k)+                                &
     &                    R0(ng)*Scoef(ng)*tl_t(i,j,k,nrhs,isalt)
#  endif
#  if defined SEDIMENT_NOT_YET && defined SED_DENS_NOT_YET
            SedDen=0.0_r8
            tl_SedDen=0.0_r8
            DO ised=1,NST
              itrc=idsed(ised)
              cff1=1.0_r8/Srho(ised,ng)
              SedDen=SedDen+                                            &
     &               t(i,j,k,nrhs,itrc)*                                &
     &               (Srho(ised,ng)-rho(i,j,k))*cff1
              tl_SedDen=tl_SedDen+                                      &
     &                  (tl_t(i,j,k,nrhs,itrc)*                         &
     &                   (Srho(ised,ng)-rho(i,j,k))-                    &
     &                   t(i,j,k,nrhs,itrc)*                            &
     &                   tl_rho(i,j,k))*cff1
            END DO
            rho(i,j,k)=rho(i,j,k)+SedDen
            tl_rho(i,j,k)=tl_rho(i,j,k)+tl_SedDen
#  endif
            rho(i,j,k)=rho(i,j,k)-1000.0_r8
#  ifdef MASKING
            rho(i,j,k)=rho(i,j,k)*rmask(i,j)
            tl_rho(i,j,k)=tl_rho(i,j,k)*rmask(i,j)
#  endif
#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
!>          pden(i,j,k)=rho(i,j,k)
!>
            tl_pden(i,j,k)=tl_rho(i,j,k)
#  endif
          END DO
        END DO

#  ifdef VAR_RHO_2D
!
!-----------------------------------------------------------------------
!  Compute vertical averaged density (rhoA) and density perturbation
!  used (rhoS) in barotropic pressure gradient.
!-----------------------------------------------------------------------
!
        DO i=IstrT,IendT
          cff1=rho(i,j,N(ng))*Hz(i,j,N(ng))
          tl_cff1=tl_rho(i,j,N(ng))*Hz(i,j,N(ng))+                      &
     &            rho(i,j,N(ng))*tl_Hz(i,j,N(ng))
          rhoS(i,j)=0.5_r8*cff1*Hz(i,j,N(ng))
          tl_rhoS(i,j)=0.5_r8*(tl_cff1*Hz(i,j,N(ng))+                   &
     &                         cff1*tl_Hz(i,j,N(ng)))
          rhoA(i,j)=cff1
          tl_rhoA(i,j)=tl_cff1
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrT,IendT
            cff1=rho(i,j,k)*Hz(i,j,k)
            tl_cff1=tl_rho(i,j,k)*Hz(i,j,k)+                            &
     &              rho(i,j,k)*tl_Hz(i,j,k)
            rhoS(i,j)=rhoS(i,j)+Hz(i,j,k)*(rhoA(i,j)+0.5_r8*cff1)
            tl_rhoS(i,j)=tl_rhoS(i,j)+                                  &
     &                   tl_Hz(i,j,k)*(rhoA(i,j)+0.5_r8*cff1)+          &
     &                   Hz(i,j,k)*(tl_rhoA(i,j)+0.5_r8*tl_cff1)
            rhoA(i,j)=rhoA(i,j)+cff1
            tl_rhoA(i,j)=tl_rhoA(i,j)+tl_cff1
          END DO
        END DO
        cff2=1.0_r8/rho0
        DO i=IstrT,IendT
          cff1=1.0_r8/(z_w(i,j,N(ng))-z_w(i,j,0))
          tl_cff1=-cff1*cff1*(tl_z_w(i,j,N(ng))-tl_z_w(i,j,0))
!
!  Here we reverse the order of the NL and TL operations since an
!  intermeridiate value of rhoA and rhoS is needed because they are
!  recursive.
!
          tl_rhoA(i,j)=cff2*(tl_cff1*rhoA(i,j)+cff1*tl_rhoA(i,j))
          rhoA(i,j)=cff2*cff1*rhoA(i,j)
          tl_rhoS(i,j)=2.0_r8*cff2*                                     &
     &                 cff1*(2.0_r8*tl_cff1*rhoS(i,j)+                  &
     &                       cff1*tl_rhoS(i,j))
          rhoS(i,j)=2.0_r8*cff1*cff1*cff2*rhoS(i,j)
        END DO
#  endif

#  ifdef BV_FREQUENCY_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute Brunt-Vaisala frequency (1/s2) at horizontal RHO-points
!  and vertical W-points.
!-----------------------------------------------------------------------
!
        DO k=1,N(ng)-1
          DO i=IstrT,IendT
!>          bvf(i,j,k)=-gorho0*(rho(i,j,k+1)-rho(i,j,k))/               &
!>   &                         (z_r(i,j,k+1)-z_r(i,j,k))
!>
            cff=1.0_r8/(z_r(i,j,k+1)-z_r(i,j,k))
            tl_cff=-cff*cff*(tl_z_r(i,j,k+1)-tl_z_r(i,j,k))
            tl_bvf(i,j,k)=-gorho0*                                      &
     &                    (cff*(tl_rho(i,j,k+1)-tl_rho(i,j,k))+         &
     &                     tl_cff*(rho(i,j,k+1)-rho(i,j,k)))
          END DO
        END DO
#  endif

#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
      defined BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Compute thermal expansion (1/Celsius) and saline contraction
!  (1/PSU) coefficients.
!-----------------------------------------------------------------------
!
        DO i=IstrT,IendT
          alpha(i,j)=ABS(Tcoef(ng))
          tl_alpha(i,j)=0.0_r8
#   ifdef SALINITY
          beta(i,j)=ABS(Scoef(ng))
          tl_beta(i,j)=0.0_r8
#   else
          beta(i,j)=0.0_r8
          tl_beta(i,j)=0.0_r8
#   endif
        END DO
#   ifdef LMD_DDMIX_NOT_YET
!
!  Compute ratio of thermal expansion and saline contraction
!  coefficients.
!
        IF (Scoef(ng).eq.0.0_r8) THEN
          cff=1.0_r8
        ELSE
          cff=1.0_r8/Scoef(ng)
        END IF
        DO k=1,N(ng)
          DO i=IstrT,IendT
!>          alfaobeta(i,j,k)=cff*Tcoef(ng)
!>
            tl_alfaobeta(i,j,k)=0.0_r8
          END DO
        END DO
#   endif
#  endif
      END DO
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          rho)
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          tl_rho)

#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
!>      CALL exchange_r3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          pden)
!>
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          tl_pden)
#  endif

#  if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
      defined BULK_FLUXES_NOT_YET
#   ifdef LMD_DDMIX_NOT_YET
!>      CALL exchange_w3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
!>   &                          alfaobeta)
!>
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          tl_alfaobeta)
#   endif
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          alpha)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_alpha)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          beta)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_beta)
#  endif

#  ifdef VAR_RHO_2D
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rhoA)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_rhoA)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rhoS)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_rhoS)
#  endif

#  ifdef BV_FREQUENCY_NOT_YET
!>      CALL exchange_w3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
!>   &                          bvf)
!>
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          tl_bvf)
#  endif
      END IF

#  ifdef DISTRIBUTE
!
      CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rho)
      CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_rho)

#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET
!>    CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    pden)
!>
      CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_pden)
#   endif

#   if defined LMD_SKPP_NOT_YET || defined LMD_BKPP_NOT_YET || \
       defined BULK_FLUXES
#    ifdef LMD_DDMIX_NOT_YET
!>    CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    alfaobeta)
!>
      CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_alfaobeta)
#    endif
      CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    alpha, beta)
      CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_alpha, tl_beta)
#   endif

#   ifdef VAR_RHO_2D
      CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rhoA, rhoS)
      CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_rhoA, tl_rhoS)
#   endif

#   ifdef BV_FREQUENCY_NOT_YET
!>    CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    bvf)
!>
      CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_bvf)
#   endif
#  endif

      RETURN
      END SUBROUTINE tl_rho_eos_tile
# endif
#endif
      END MODULE tl_rho_eos_mod
