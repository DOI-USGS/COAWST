#include "cppdefs.h"

#if defined ADJOINT  && defined SOLVE3D       && \
   (defined LMD_SKPP_NOT_YET || defined SOLAR_SOURCE)

      SUBROUTINE ad_lmd_swfrac_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               Zscale, Z, ad_Z, ad_swdk)
!
!svn $Id: ad_lmd_swfrac.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the adjoint fraction of solar shortwave flux  !
!  penetrating to specified depth (times Zscale)  due to  exponential  !
!  decay in Jerlov water type.                                         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Zscale   Scale factor to apply to depth array.                   !
!     Z        Vertical height (meters, negative) for                  !
!                desired solar short-wave fraction.                    !
!     ad_Z     Adjoint vertical height for                             !
!                desired solar short-wave fraction.                    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ad_swdk  Adjoint shortwave (radiation) fractional decay.         !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Paulson, C.A., and J.J. Simpson, 1977: Irradiance meassurements     !
!     in the upper ocean, J. Phys. Oceanogr., 7, 952-956.              !
!                                                                      !
!  This routine was adapted from Bill Large 1995 code.                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_mixing
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS

      real(r8), intent(in) :: Zscale

      real(r8), intent(in) :: Z(IminS:ImaxS,JminS:JmaxS)

      real(r8), intent(inout) :: ad_Z(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(inout) :: ad_swdk(IminS:ImaxS,JminS:JmaxS)
!
!  Local variable declarations.
!
      integer :: Jindex, i, j

      real(r8) :: cff1, cff2
      real(r8) :: ad_cff1, ad_cff2

      real(r8), dimension(IminS:ImaxS) :: fac1, fac2, fac3

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_cff1=0.0_r8
      ad_cff2=0.0_r8
!
!-----------------------------------------------------------------------
!  Use Paulson and Simpson (1977) two wavelength bands solar
!  absorption model.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Jindex=INT(MIXING(ng)%Jwtype(i,j))
          fac1(i)=Zscale/lmd_mu1(Jindex)
          fac2(i)=Zscale/lmd_mu2(Jindex)
          fac3(i)=lmd_r1(Jindex)
        END DO
!!DIR$ VECTOR ALWAYS
        DO i=Istr,Iend
          cff1=EXP(Z(i,j)*fac1(i))
          cff2=EXP(Z(i,j)*fac2(i))
!>        tl_swdk(i,j)=tl_cff1*fac3(i)+                                 &
!>   &                 tl_cff2*(1.0_r8-fac3(i))
!>
          ad_cff1=ad_cff1+fac3(i)*ad_swdk(i,j)
          ad_cff2=ad_cff2+(1.0_r8-fac3(i))*ad_swdk(i,j)
          ad_swdk(i,j)=0.0_r8
!>        tl_cff2=fac2(i)*tl_Z(i,j)*cff2
!>        tl_cff1=fac1(i)*tl_Z(i,j)*cff1
!>
          ad_Z(i,j)=ad_Z(i,j)+                                          &
     &              fac1(i)*cff1*ad_cff1+                               &
     &              fac2(i)*cff2*ad_cff2
          ad_cff2=0.0_r8
          ad_cff1=0.0_r8
        END DO
      END DO
      RETURN
      END SUBROUTINE ad_lmd_swfrac_tile
#else
      SUBROUTINE ad_lmd_swfrac
      RETURN
      END SUBROUTINE ad_lmd_swfrac
#endif

