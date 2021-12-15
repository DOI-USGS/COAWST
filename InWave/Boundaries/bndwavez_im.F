#include "cppdefs.h"
      MODULE bndwavez_mod
#ifdef INWAVE_SWAN_COUPLING
# ifdef SOLVE3D
!
!svn $Id: bndwavez_im.F 732 2008-09-07 01:55:51Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!                                                   John C. Warner     !
!                                                                      !
!  This subroutine sets lateral boundary conditions for the            !
!  water levels by adding the bound wave to zeta.                      !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: bndwavez

      CONTAINS
!
!***********************************************************************
      SUBROUTINE bndwavez (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_inwave_swan
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: Insteps

# include "tile.h"

      Insteps=WAVES(ng)%Insteps
!
      CALL bndwavez_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    Insteps)
      RETURN
      END SUBROUTINE bndwavez

!
!***********************************************************************
      SUBROUTINE bndwavez_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Insteps)
!***********************************************************************

      USE mod_forces
      USE mod_grid
      USE mod_ncparam
      USE mod_boundary
      USE mod_scalars
      USE mod_inwave_vars
      USE mod_inwave_swan
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Insteps
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: i, j, d
      integer :: tidx1, tidx2
      real(r8):: cff, cff1, cff2, cff3, ramp, alpha
      real(r8):: fac1, fac2, tidx, ratio
      real(r8):: phiw, cw, dist, toff, bndtot, veltot
      real(r8), parameter :: eps = 0.0001_r8

#  include "set_bounds.h"
#  ifdef RAMP_INWAVE
      ramp=TANH((tdays(ng)-dstart)/0.005_r8)
#  else
      ramp=1.0_r8
#  endif
      ratio=WAVES(ng)%idt/(dt(ng)*WAVES(ng)%dtifac)
      tidx=REAL(MOD(iic(ng),INT((WAVES(ng)%Insteps-1)*ratio)),r8)/ratio+&
     &     1.0_r8
      tidx1=INT(tidx)
      tidx2=tidx+1
      fac2=tidx-REAL(INT(tidx),r8)
      fac1=1.0_r8-fac2
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (LBC(iwest,isFsur,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        IF (LBC(iwest,isAC3d,ng)%acquire) THEN
          DO j=Jstr,Jend
            bndtot=fac1*WAVES(ng)%bndwave_z(j,tidx1)+                   &
     &             fac2*WAVES(ng)%bndwave_z(j,tidx2)
            BOUNDARY(ng)%zeta_west(j)=(BOUNDARY(ng)%zeta_west(j)+       &
     &                                 bndtot)*ramp
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (LBC(ieast,isFsur,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        IF (LBC(ieast,isAC3d,ng)%acquire) THEN
          DO j=Jstr,Jend
            bndtot=fac1*WAVES(ng)%bndwave_z(j,tidx1)+                   &
     &             fac2*WAVES(ng)%bndwave_z(j,tidx2)
            BOUNDARY(ng)%zeta_east(j)=(BOUNDARY(ng)%zeta_east(j)+       &
     &                                 bndtot)*ramp
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (LBC(isouth,isFsur,ng)%acquire.and.                            &
     &    DOMAIN(ng)%Southern_Edge(tile)) THEN
        IF (LBC(isouth,isAC3d,ng)%acquire) THEN
          DO i=Istr,Iend
            bndtot=fac1*WAVES(ng)%bndwave_z(i,tidx1)+                   &
     &             fac2*WAVES(ng)%bndwave_z(i,tidx2)
            BOUNDARY(ng)%zeta_south(i)=(BOUNDARY(ng)%zeta_south(i)+     &
     &                                  bndtot)*ramp
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (LBC(inorth,isFsur,ng)%acquire.and.                            &
     &    DOMAIN(ng)%Northern_Edge(tile)) THEN
        IF (LBC(inorth,isAC3d,ng)%acquire) THEN
          DO i=Istr,Iend
            bndtot=fac1*WAVES(ng)%bndwave_z(i,tidx1)+                   &
     &             fac2*WAVES(ng)%bndwave_z(i,tidx2)
            BOUNDARY(ng)%zeta_north(i)=(BOUNDARY(ng)%zeta_north(i)+     &
     &                                  bndtot)*ramp
          END DO
        END IF
      END IF
!
      RETURN
      END SUBROUTINE bndwavez_tile
# endif
#endif
      END MODULE bndwavez_mod
