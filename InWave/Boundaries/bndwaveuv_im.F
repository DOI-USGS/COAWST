#include "cppdefs.h"
      MODULE bndwaveuv_mod
#ifdef INWAVE_SWAN_COUPLING
# ifdef SOLVE3D
!
!svn $Id: bndwaveuv_im.F 732 2008-09-07 01:55:51Z jcwarner $
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
      PUBLIC  :: bndwaveuv

      CONTAINS
!
!***********************************************************************
      SUBROUTINE bndwaveuv (ng, tile)
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
      CALL bndwaveuv_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     Insteps)
      RETURN
      END SUBROUTINE bndwaveuv

!
!***********************************************************************
      SUBROUTINE bndwaveuv_tile (ng, tile,                              &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Insteps)
!***********************************************************************

      USE mod_forces
      USE mod_grid
      USE mod_ncparam
      USE mod_boundary
      USE mod_scalars
      USE mod_ocean
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
      real(r8):: cff1, ramp
      real(r8):: fac1, fac2, tidx, ratio
      real(r8):: bndtot, veltot
      real(r8), parameter :: eps = 0.0001_r8

#  include "set_bounds.h"
#  ifdef RAMP_INWAVE
      ramp=TANH((tdays(ng)-dstart)/0.005_r8)
#  else
      ramp=1.0_r8
#  endif
      ratio=WAVES(ng)%idt/dt(ng)
      tidx=REAL(MOD(iic(ng),INT((WAVES(ng)%Insteps-1)*ratio)),r8)/ratio+1.0_r8
      tidx1=INT(tidx)
      tidx2=tidx+1
      fac2=tidx-REAL(INT(tidx),r8)
      fac1=1.0_r8-fac2
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
!      IF (LBC(iwest,isFsur,ng)%acquire.and.                             &
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        IF (LBC(iwest,isAC3d,ng)%acquire) THEN
          DO j=Jstr,Jend
# ifdef VARY_ACBC
!           bndtot=WAVES(ng)%bndwave_z(j,tidx)
!           veltot=WAVES(ng)%bndwave_u(j,tidx)
!            bndtot=fac1*WAVES(ng)%bndwave_z(j,tidx1)+                   &
!     &             fac2*WAVES(ng)%bndwave_z(j,tidx2)
            veltot=fac1*WAVES(ng)%bndwave_u(j,tidx1)+                   &
     &             fac2*WAVES(ng)%bndwave_u(j,tidx2)
# else
!           bndtot=WAVES(ng)%bndwave_z(tidx)
!           veltot=WAVES(ng)%bndwave_u(tidx)
!            bndtot=fac1*WAVES(ng)%bndwave_z(tidx1)+                     &
!     &             fac2*WAVES(ng)%bndwave_z(tidx2)
            veltot=fac1*WAVES(ng)%bndwave_u(tidx1)+                     &
     &             fac2*WAVES(ng)%bndwave_u(tidx2)
# endif
!            BOUNDARY(ng)%zeta_west(j)=(BOUNDARY(ng)%zeta_west(j)+       &
!     &                                 bndtot)*ramp
            cff1=1.0_r8/(0.5_r8*(GRID(ng)%h(Istr-1,j)+                  &
     &                           OCEAN(ng)%zeta(Istr-1,j,1)+            &
     &                           GRID(ng)%h(Istr  ,j)+                  &
     &                           OCEAN(ng)%zeta(Istr  ,j,1)))
            IF (LBC(iwest,isUbar,ng)%acquire) THEN
              BOUNDARY(ng)%ubar_west(j)=(BOUNDARY(ng)%ubar_west(j)+     &
     &                                   veltot*cff1)*ramp
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
!      IF (LBC(ieast,isFsur,ng)%acquire.and.                             &
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        IF (LBC(ieast,isAC3d,ng)%acquire) THEN
          DO j=Jstr,Jend
# ifdef VARY_ACBC
!           bndtot=WAVES(ng)%bndwave_z(j,tidx)
!           veltot=WAVES(ng)%bndwave_u(j,tidx)
!            bndtot=fac1*WAVES(ng)%bndwave_z(j,tidx1)+                   &
!     &             fac2*WAVES(ng)%bndwave_z(j,tidx2)
            veltot=fac1*WAVES(ng)%bndwave_u(j,tidx1)+                   &
     &             fac2*WAVES(ng)%bndwave_u(j,tidx2)
# else
!           bndtot=WAVES(ng)%bndwave_z(tidx)
!           veltot=WAVES(ng)%bndwave_u(tidx)
!            bndtot=fac1*WAVES(ng)%bndwave_z(tidx1)+                     &
!     &             fac2*WAVES(ng)%bndwave_z(tidx2)
            veltot=fac1*WAVES(ng)%bndwave_u(tidx1)+                     &
     &             fac2*WAVES(ng)%bndwave_u(tidx2)
# endif
!            BOUNDARY(ng)%zeta_east(j)=(BOUNDARY(ng)%zeta_east(j)+       &
!     &                                 bndtot)*ramp
            cff1=1.0_r8/(0.5_r8*(GRID(ng)%h(Iend  ,j)+                  &
     &                           OCEAN(ng)%zeta(Iend  ,j,1)+            &
     &                           GRID(ng)%h(Iend+1,j)+                  &
     &                           OCEAN(ng)%zeta(Iend+1,j,1)))
            IF (LBC(ieast,isUbar,ng)%acquire) THEN
              BOUNDARY(ng)%ubar_east(j)=(BOUNDARY(ng)%ubar_east(j)+     &
     &                                   veltot*cff1)*ramp
!             BOUNDARY(ng)%ubar_east(j)=(BOUNDARY(ng)%ubar_east(j)+     &
!    &                                   veltot)*ramp
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
!      IF (LBC(isouth,isFsur,ng)%acquire.and.                            &
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        IF (LBC(isouth,isAC3d,ng)%acquire) THEN
          DO i=Istr,Iend
# ifdef VARY_ACBC
!           bndtot=WAVES(ng)%bndwave_z(i,tidx)
!           veltot=WAVES(ng)%bndwave_v(i,tidx)
!            bndtot=fac1*WAVES(ng)%bndwave_z(i,tidx1)+                   &
!     &             fac2*WAVES(ng)%bndwave_z(i,tidx2)
            veltot=fac1*WAVES(ng)%bndwave_v(i,tidx1)+                   &
     &             fac2*WAVES(ng)%bndwave_v(i,tidx2)
# else
!           bndtot=WAVES(ng)%bndwave_z(tidx)
!           veltot=WAVES(ng)%bndwave_v(tidx)
!            bndtot=fac1*WAVES(ng)%bndwave_z(tidx1)+                     &
!     &             fac2*WAVES(ng)%bndwave_z(tidx2)
            veltot=fac1*WAVES(ng)%bndwave_v(tidx1)+                     &
     &             fac2*WAVES(ng)%bndwave_v(tidx2)
# endif
!            BOUNDARY(ng)%zeta_south(i)=(BOUNDARY(ng)%zeta_south(i)+     &
!     &                                  bndtot)*ramp
            cff1=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jstr-1)+                  &
     &                           OCEAN(ng)%zeta(i,Jstr-1,1)+            &
     &                           GRID(ng)%h(i,Jstr  )+                  &
     &                           OCEAN(ng)%zeta(i,Jstr  ,1)))
            IF (LBC(isouth,isVbar,ng)%acquire) THEN
              BOUNDARY(ng)%vbar_south(i)=(BOUNDARY(ng)%vbar_south(i)+   &
     &                                    veltot*cff1)*ramp
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
!      IF (LBC(inorth,isFsur,ng)%acquire.and.                            &
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        IF (LBC(inorth,isAC3d,ng)%acquire) THEN
          DO i=Istr,Iend
# ifdef VARY_ACBC
!           bndtot=WAVES(ng)%bndwave_z(i,tidx)
!           veltot=WAVES(ng)%bndwave_v(i,tidx)
!            bndtot=fac1*WAVES(ng)%bndwave_z(i,tidx1)+                   &
!     &             fac2*WAVES(ng)%bndwave_z(i,tidx2)
            veltot=fac1*WAVES(ng)%bndwave_v(i,tidx1)+                   &
     &             fac2*WAVES(ng)%bndwave_v(i,tidx2)
# else
!           bndtot=WAVES(ng)%bndwave_z(tidx)
!           veltot=WAVES(ng)%bndwave_v(tidx)
!            bndtot=fac1*WAVES(ng)%bndwave_z(tidx1)+                     &
!     &             fac2*WAVES(ng)%bndwave_z(tidx2)
            veltot=fac1*WAVES(ng)%bndwave_v(tidx1)+                     &
     &             fac2*WAVES(ng)%bndwave_v(tidx2)
# endif
!            BOUNDARY(ng)%zeta_north(i)=(BOUNDARY(ng)%zeta_north(i)+     &
!     &                                  bndtot)*ramp
            cff1=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jend  )+                  &
     &                           OCEAN(ng)%zeta(i,Jend  ,1)+            &
     &                           GRID(ng)%h(i,Jend+1)+                  &
     &                           OCEAN(ng)%zeta(i,Jend+1,1)))
            IF (LBC(inorth,isVbar,ng)%acquire) THEN
              BOUNDARY(ng)%vbar_north(i)=(BOUNDARY(ng)%vbar_north(i)+   &
     &                                    veltot*cff1)*ramp
            END IF
          END DO
        END IF
      END IF
!
      RETURN
      END SUBROUTINE bndwaveuv_tile
# endif
#endif
      END MODULE bndwaveuv_mod
