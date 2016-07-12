#include "cppdefs.h"

#if defined NONLINEAR && (defined LMD_SKPP || defined SOLAR_SOURCE) && \
    defined SOLVE3D

      SUBROUTINE lmd_swfrac_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            Zscale, Z, swdk)
!
!svn $Id: lmd_swfrac.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine computes the  fraction  of  solar shortwave flux    !
!  penetrating to specified depth (times Zscale) due to exponential    !
!  decay in Jerlov water type.                                         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Zscale   scale factor to apply to depth array.                   !
!     Z        vertical height (meters, negative) for                  !
!              desired solar short-wave fraction.                      !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     swdk     shortwave (radiation) fractional decay.                 !
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
      real(r8), intent(out) :: swdk(IminS:ImaxS,JminS:JmaxS)
!
!  Local variable declarations.
!
      integer :: Jindex, i, j

      real(r8), dimension(IminS:ImaxS) :: fac1, fac2, fac3

# include "set_bounds.h"
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
          swdk(i,j)=EXP(Z(i,j)*fac1(i))*fac3(i)+                        &
     &              EXP(Z(i,j)*fac2(i))*(1.0_r8-fac3(i))
        END DO
      END DO
      RETURN
      END SUBROUTINE lmd_swfrac_tile
#else
      SUBROUTINE lmd_swfrac
      RETURN
      END SUBROUTINE lmd_swfrac
#endif

