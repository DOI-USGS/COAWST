#include "cppdefs.h"
#if defined INWAVE_MODEL & defined INWAVE_SWAN_COUPLING
      SUBROUTINE set_inwave_swan_data (ng, tile)
!
!svn $Id: set_inwave_swan_data.F 799 2009-12-08 20:38:55Z jcwarner $
!=======================================================================
!                                                                      !
!  This routine computes the action density for the boundary condition !
!  The action density is derived from the directional                  !
!  wave spectra given by SWAN                                          !
!=======================================================================
!
      USE mod_param
      USE mod_inwave_bound
!
!  Imported variable declarations.

      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
!      CALL wclock_on (ng, iNLM, 4)
# endif
      CALL set_inwave_swan_data_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj)
# ifdef PROFILE
!      CALL wclock_off (ng, iNLM, 4)
# endif
      RETURN
      END SUBROUTINE set_inwave_swan_data
!
!***********************************************************************
      SUBROUTINE set_inwave_swan_data_tile (ng, tile,                   &
     &                                      LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_forces
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_inwave_params
      USE mod_inwave_vars
      USE mod_inwave_bound
      USE mod_inwave_swan
      USE mod_grid
!
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d
#  endif
# endif
      USE set_2dfld_mod
# ifdef SOLVE3D
      USE set_3dfld_mod
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: ILB, IUB, JLB, JUB
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, d, my_tile
      integer :: tidx1, tidx2

      real(r8):: twopi, otwopi, cff1, cff2, cff3
      real(r8):: phiw, cw, dist, toff, alpha
      real(r8):: fac1, fac2, tidx, ratio
      real(r8), parameter :: eps = 0.0001_r8

# include "set_bounds.h"

      twopi=2.0_r8*pi
      otwopi=1.0_r8/(2.0_r8*pi)
!
!     To compute the time index for the BC data, the Ampzeta array
!     is Insteps long (1000sec), but it repeats every 1/df.
!
      ratio=WAVES(ng)%idt/(dt(ng)*WAVES(ng)%dtifac)
      tidx=REAL(MOD(iic(ng),INT((WAVES(ng)%Insteps-1)*ratio)),r8)/ratio+&
     &     1.0_r8
      tidx1=INT(tidx)
      tidx2=tidx+1
      fac2=tidx-REAL(INT(tidx),r8)
      fac1=1.0_r8-fac2

      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        IF (LBC(iwest,isAC3d,ng)%acquire) THEN
          DO d=1,ND
            DO j=Jstr,Jend
              cff2=fac1*WAVES(ng)%Ampzeta(j,tidx1,d)+                   &
     &             fac2*WAVES(ng)%Ampzeta(j,tidx2,d)
              WAVEB(ng)%AC_west(j,d)=cff2*WAVEP(ng)%Tr(Istr-1,j,d)*     &
     &                               otwopi
            END DO
          END DO
        END IF
      END IF
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        IF (LBC(ieast,isAC3d,ng)%acquire) THEN
          DO d=1,ND
            DO j=Jstr,Jend
              cff2=fac1*WAVES(ng)%Ampzeta(j,tidx1,d)+                   &
     &             fac2*WAVES(ng)%Ampzeta(j,tidx2,d)
              WAVEB(ng)%AC_east(j,d)=cff2*WAVEP(ng)%Tr(Iend+1,j,d)*     &
     &                               otwopi
            END DO
          END DO
        END IF
      END IF
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        IF (LBC(inorth,isAC3d,ng)%acquire) THEN
          DO d=1,ND
            DO i=Istr,Iend
              cff2=fac1*WAVES(ng)%Ampzeta(i,tidx1,d)+                   &
     &             fac2*WAVES(ng)%Ampzeta(i,tidx2,d)
              WAVEB(ng)%AC_north(i,d)=cff2*WAVEP(ng)%Tr(i,Jend+1,d)*    &
     &                                otwopi
            END DO
          END DO
        END IF
      END IF
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        IF (LBC(isouth,isAC3d,ng)%acquire) THEN
          DO d=1,ND
            DO i=Istr,Iend
              cff2=fac1*WAVES(ng)%Ampzeta(i,tidx1,d)+                   &
     &             fac2*WAVES(ng)%Ampzeta(i,tidx2,d)
              WAVEB(ng)%AC_south(i,d)=cff2*WAVEP(ng)%Tr(i,Jstr-1,d)*    &
     &                                otwopi
            END DO
          END DO
        END IF
      END IF

      RETURN
      END SUBROUTINE set_inwave_swan_data_tile
#else
      SUBROUTINE set_inwave_swan_data
      RETURN
      END SUBROUTINE set_inwave_swan_data
#endif
