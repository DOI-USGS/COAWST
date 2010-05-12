#include "cppdefs.h"
      MODULE u2dbc_mod
!
!svn $Id: u2dbc_ex.F 429 2009-12-20 17:30:26Z arango $
!=======================================================================
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for vertically     !
!  integrated U-velocity.                                              !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: u2dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE u2dbc (ng, tile, kout)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile, kout

#include "tile.h"
!
      CALL u2dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 krhs(ng), kstp(ng), kout,                        &
     &                 OCEAN(ng) % ubar,                                &
     &                 OCEAN(ng) % vbar,                                &
     &                 OCEAN(ng) % zeta)
      RETURN
      END SUBROUTINE u2dbc
!
!***********************************************************************
      SUBROUTINE u2dbc_tile (ng, Istr, Iend, Jstr, Jend,                &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       krhs, kstp, kout,                          &
     &                       ubar, vbar, zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: krhs, kstp, kout

#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)

      real(r8), intent(inout) :: ubar(LBi:,LBj:,:)
#else
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)

      real(r8), intent(inout) :: ubar(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, know

      real(r8), parameter :: eps = 1.0E-20_r8

      real(r8) :: Ce, Cx, cff, dt2d, dUde, dUdt, dUdx, tau

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: grad

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set time-indices
!-----------------------------------------------------------------------
!
      IF (FIRST_2D_STEP) THEN
        know=krhs
        dt2d=dtfast(ng)
      ELSE IF (PREDICTOR_2D_STEP(ng)) THEN
        know=krhs
        dt2d=2.0_r8*dtfast(ng)
      ELSE
        know=kstp
        dt2d=dtfast(ng)
      END IF

#ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN

# if defined WEST_M2RADIATION
!
!  Western edge, implicit upstream radiation condition.
!
        DO j=Jstr,Jend+1
          grad(1,j)=ubar(Istr  ,j  ,know)-                              &
     &              ubar(Istr  ,j-1,know)
          grad(2,j)=ubar(Istr+1,j  ,know)-                              &
     &              ubar(Istr+1,j-1,know)
        END DO
        DO j=Jstr,Jend
          dUdt=ubar(Istr+1,j,know)-ubar(Istr+1,j,kout)
          dUdx=ubar(Istr+1,j,know)-ubar(Istr+2,j,know)
#  ifdef WEST_M2NUDGING
          tau=M2obc_out(ng,iwest)
          IF ((dUdt*dUdx).lt.0.0_r8) tau=M2obc_in(ng,iwest)
          tau=tau*dt2d
#  endif
          IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
          IF ((dUdt*(grad(Istr+1,j)+grad(Istr+1,j+1))).gt.0.0_r8) THEN
            dUde=grad(Istr+1,j  )
          ELSE
            dUde=grad(Istr+1,j+1)
          END IF
          cff=dUdt/MAX(dUdx*dUdx+dUde*dUde,eps)
          Cx=MIN(1.0_r8,cff*dUdx)
#  ifdef RADIATION_2D
          Ce=MIN(1.0_r8,MAX(-1.0_r8,cff*dUde))
#  else
          Ce=0.0_r8
#  endif
          ubar(Istr,j,kout)=(1.0_r8-Cx)*ubar(Istr,j,know)+              &
     &                      Cx*ubar(Istr+1,j,know)-                     &
     &                      MAX(Ce,0.0_r8)*grad(Istr,j  )-              &
     &                      MIN(Ce,0.0_r8)*grad(Istr,j+1)
#  ifdef WEST_M2NUDGING
         ubar(Istr,j,kout)=ubar(Istr,j,kout)+                           &
     &                     tau*(BOUNDARY(ng)%ubar_west(j)-              &
     &                          ubar(Istr,j,know))
#  endif
#  ifdef MASKING
          ubar(Istr,j,kout)=ubar(Istr,j,kout)*                          &
     &                      GRID(ng)%umask(Istr,j)
#  endif
        END DO

# elif defined WEST_M2FLATHER
!
!  Western edge, Flather boundary condition.
!
        DO j=Jstr,Jend
          ubar(Istr,j,kout)=BOUNDARY(ng)%ubar_west(j)-                  &
     &                      SQRT(g/(0.5_r8*(GRID(ng)%h(Istr-1,j)+       &
     &                                      GRID(ng)%h(Istr  ,j))))*    &
     &                      (0.5_r8*(zeta(Istr-1,j,know)+               &
     &                               zeta(Istr  ,j,know))-              &
     &                      BOUNDARY(ng)%zeta_west(j))
#  ifdef MASKING
          ubar(Istr,j,kout)=ubar(Istr,j,kout)*                          &
     &                      GRID(ng)%umask(Istr,j)
#  endif
        END DO

# elif defined WEST_M2CLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          ubar(Istr,j,kout)=BOUNDARY(ng)%ubar_west(j)
#  ifdef MASKING
          ubar(Istr,j,kout)=ubar(Istr,j,kout)*                          &
     &                      GRID(ng)%umask(Istr,j)
#  endif
        END DO

# elif defined WEST_M2GRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          ubar(Istr,j,kout)=ubar(Istr+1,j,kout)
#  ifdef MASKING
          ubar(Istr,j,kout)=ubar(Istr,j,kout)*                          &
     &                      GRID(ng)%umask(Istr,j)
#  endif
        END DO

# else
!
!  Western edge, closed boundary condition.
!
        DO j=Jstr,Jend
          ubar(Istr,j,kout)=0.0_r8
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN

# if defined EAST_M2RADIATION
!
!  Eastern edge, implicit upstream radiation condition.
!
        DO j=Jstr,Jend+1
          grad(Iend  ,j)=ubar(Iend  ,j  ,know)-                         &
     &                   ubar(Iend  ,j-1,know)
          grad(Iend+1,j)=ubar(Iend+1,j  ,know)-                         &
     &                   ubar(Iend+1,j-1,know)
        END DO
        DO j=Jstr,Jend
          dUdt=ubar(Iend,j,know)-ubar(Iend  ,j,kout)
          dUdx=ubar(Iend,j,know)-ubar(Iend-1,j,know)
#  ifdef EAST_M2NUDGING
          tau=M2obc_out(ng,ieast)
          IF ((dUdt*dUdx).lt.0.0_r8) tau=M2obc_in(ng,ieast)
          tau=tau*dt2d
#  endif
          IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
          IF ((dUdt*(grad(Iend,j)+grad(Iend,j+1))).gt.0.0_r8) THEN
            dUde=grad(Iend,j)
          ELSE
            dUde=grad(Iend,j+1)
          END IF
          cff=dUdt/MAX(dUdx*dUdx+dUde*dUde,eps)
          Cx=MIN(1.0_r8,cff*dUdx)
#  ifdef RADIATION_2D
          Ce=MIN(1.0_r8,MAX(-1.0_r8,cff*dUde))
#  else
          Ce=0.0_r8
#  endif
          ubar(Iend+1,j,kout)=(1.0_r8-Cx)*ubar(Iend+1,j,know)+          &
     &                        Cx*ubar(Iend,j,know)-                     &
     &                        MAX(Ce,0.0_r8)*grad(Iend+1,j  )-          &
     &                        MIN(Ce,0.0_r8)*grad(Iend+1,j+1)
#  ifdef EAST_M2NUDGING
          ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)+                      &
     &                        tau*(BOUNDARY(ng)%ubar_east(j)-           &
     &                             ubar(Iend+1,j,know))
#  endif
#  ifdef MASKING
          ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)*                      &
     &                        GRID(ng)%umask(Iend+1,j)
#  endif
        END DO

# elif defined EAST_M2FLATHER
!
!  Eastern edge, Flather boundary condition.
!
        DO j=Jstr,Jend
          ubar(Iend+1,j,kout)=BOUNDARY(ng)%ubar_east(j)+                &
     &                        SQRT(g/(0.5_r8*(GRID(ng)%h(Iend  ,j)+     &
     &                                        GRID(ng)%h(Iend+1,j))))*  &
     &                        (0.5_r8*(zeta(Iend  ,j,know)+             &
     &                                 zeta(Iend+1,j,know))-            &
     &                        BOUNDARY(ng)%zeta_east(j))
#  ifdef MASKING
          ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)*                      &
     &                        GRID(ng)%umask(Iend+1,j)
#  endif
        END DO

# elif defined EAST_M2CLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          ubar(Iend+1,j,kout)=BOUNDARY(ng)%ubar_east(j)
#  ifdef MASKING
          ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)*                      &
     &                        GRID(ng)%umask(Iend+1,j)
#  endif
        END DO

# elif defined EAST_M2GRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          ubar(Iend+1,j,kout)=ubar(Iend,j,kout)
#  ifdef MASKING
          ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)*                      &
     &                        GRID(ng)%umask(Iend+1,j)
#  endif
        END DO

# else
!
!  Eastern edge, closed boundary condition.
!
        DO j=Jstr,Jend
          ubar(Iend+1,j,kout)=0.0_r8
        END DO
# endif
      END IF
#endif

#ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE) THEN

# if defined SOUTH_M2RADIATION
!
!  Southern edge, implicit upstream radiation condition.
!
        DO i=IstrU-1,Iend
          grad(i,Jstr-1)=ubar(i+1,Jstr-1,know)-                         &
     &                   ubar(i  ,Jstr-1,know)
          grad(i,Jstr  )=ubar(i+1,Jstr  ,know)-                         &
     &                   ubar(i  ,Jstr  ,know)
        END DO
        DO i=IstrU,Iend
          dUdt=ubar(i,Jstr,know)-ubar(i,Jstr  ,kout)
          dUde=ubar(i,Jstr,know)-ubar(i,Jstr+1,know)
#  ifdef SOUTH_M2NUDGING
          tau=M2obc_out(ng,isouth)
          IF ((dUdt*dUde).lt.0.0_r8) tau=M2obc_in(ng,isouth)
          tau=tau*dt2d
#  endif
          IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
          IF ((dUdt*(grad(i-1,Jstr)+grad(i,Jstr))).gt.0.0_r8) THEN
            dUdx=grad(i-1,Jstr)
          ELSE
            dUdx=grad(i  ,Jstr)
          END IF
          cff=dUdt/MAX(dUdx*dUdx+dUde*dUde,eps)
#  ifdef RADIATION_2D
          Cx=MIN(1.0_r8,MAX(-1.0_r8,cff*dUdx))
#  else
          Cx=0.0_r8
#  endif
          Ce=MIN(1.0_r8,cff*dUde)
          ubar(i,Jstr-1,kout)=(1.0_r8-Ce)*ubar(i,Jstr-1,know)+          &
     &                        Ce*ubar(i,Jstr,know)-                     &
     &                        MAX(Cx,0.0_r8)*grad(i-1,Jstr-1)-          &
     &                        MIN(Cx,0.0_r8)*grad(i  ,Jstr-1)
#  ifdef SOUTH_M2NUDGING
          ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)+                      &
     &                        tau*(BOUNDARY(ng)%ubar_south(i)-          &
     &                             ubar(i,Jstr-1,know))
#  endif
#  ifdef MASKING
          ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)*                      &
     &                        GRID(ng)%umask(i,Jstr-1)
#  endif
        END DO

# elif defined SOUTH_M2FLATHER
!
!  Southern edge, Chapman boundary condition.
!
        DO i=IstrU,Iend
          Ce=dt2d*0.5_r8*(GRID(ng)%pn(i-1,Jstr)+                        &
     &                    GRID(ng)%pn(i  ,Jstr))*                       &
     &       SQRT(g*0.5_r8*(GRID(ng)%h(i-1,Jstr)+                       &
     &                      GRID(ng)%h(i  ,Jstr)))
          ubar(i,Jstr-1,kout)=(ubar(i,Jstr-1,know)+                     &
     &                         Ce*ubar(i,Jstr,kout))/                   &
     &                        (1.0_r8+Ce)
#  ifdef MASKING
          ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)*                      &
     &                        GRID(ng)%umask(i,Jstr-1)
#  endif
        END DO

# elif defined SOUTH_M2CLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=IstrU,Iend
          ubar(i,Jstr-1,kout)=BOUNDARY(ng)%ubar_south(i)
#  ifdef MASKING
          ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)*                      &
     &                        GRID(ng)%umask(i,Jstr-1)
#  endif
        END DO

# elif defined SOUTH_M2GRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=IstrU,Iend
          ubar(i,Jstr-1,kout)=ubar(i,Jstr,kout)
#  ifdef MASKING
          ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)*                      &
     &                        GRID(ng)%umask(i,Jstr-1)
#  endif
        END DO

# else
!
!  Southern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
#  ifdef EW_PERIODIC
#   define I_RANGE IstrU,Iend
#  else
#   define I_RANGE Istr,IendR
#  endif
        DO i=I_RANGE
          ubar(i,Jstr-1,kout)=gamma2(ng)*ubar(i,Jstr,kout)
#  ifdef MASKING
          ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)*                      &
     &                        GRID(ng)%umask(i,Jstr-1)
#  endif
        END DO
#  undef I_RANGE
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (NORTHERN_EDGE) THEN

# if defined NORTH_M2RADIATION
!
!  Northern edge, implicit upstream radiation condition.
!
        DO i=IstrU-1,Iend
          grad(i,Jend  )=ubar(i+1,Jend  ,know)-                         &
     &                   ubar(i  ,Jend  ,know)
          grad(i,Jend+1)=ubar(i+1,Jend+1,know)-                         &
     &                   ubar(i  ,Jend+1,know)
        END DO
        DO i=IstrU,Iend
          dUdt=ubar(i,Jend,know)-ubar(i,Jend  ,kout)
          dUde=ubar(i,Jend,know)-ubar(i,Jend-1,know)
#  ifdef NORTH_M2NUDGING
          tau=M2obc_out(ng,inorth)
          IF ((dUdt*dUde).lt.0.0_r8) tau=M2obc_in(ng,inorth)
          tau=tau*dt2d
#  endif
          IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
          IF ((dUdt*(grad(i-1,Jend)+grad(i,Jend))).gt.0.0_r8) THEN
            dUdx=grad(i-1,Jend)
          ELSE
            dUdx=grad(i  ,Jend)
          END IF
          cff=dUdt/MAX(dUdx*dUdx+dUde*dUde,eps)
#  ifdef RADIATION_2D
          Cx=MIN(1.0_r8,MAX(-1.0_r8,cff*dUdx))
#  else
          Cx=0.0_r8
#  endif
          Ce=MIN(1.0_r8,cff*dUde)
          ubar(i,Jend+1,kout)=(1.0_r8-Ce)*ubar(i,Jend+1,know)+          &
     &                        Ce*ubar(i,Jend,know)-                     &
     &                        MAX(Cx,0.0_r8)*grad(i-1,Jend+1)-          &
     &                        MIN(Cx,0.0_r8)*grad(i  ,Jend+1)
#  ifdef NORTH_M2NUDGING
          ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)+                      &
     &                        tau*(BOUNDARY(ng)%ubar_north(i)-          &
     &                             ubar(i,Jend+1,know))
#  endif
#  ifdef MASKING
          ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)*                      &
     &                        GRID(ng)%umask(i,Jend+1)
#  endif
        END DO

# elif defined NORTH_M2FLATHER
!
!  Northern edge, Chapman boundary condition.
!
        DO i=IstrU,Iend
          Ce=dt2d*0.5_r8*(GRID(ng)%pn(i-1,Jend)+                        &
     &                    GRID(ng)%pn(i  ,Jend))*                       &
     &       SQRT(g*0.5_r8*(GRID(ng)%h(i-1,Jend)+                       &
     &                      GRID(ng)%h(i  ,Jend)))
          ubar(i,Jend+1,kout)=(ubar(i,Jend+1,know)+                     &
     &                         Ce*ubar(i,Jend,kout))/                   &
     &                        (1.0_r8+Ce)
#  ifdef MASKING
          ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)*                      &
     &                        GRID(ng)%umask(i,Jend+1)
#  endif
        END DO

# elif defined NORTH_M2CLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=IstrU,Iend
          ubar(i,Jend+1,kout)=BOUNDARY(ng)%ubar_north(i)
#  ifdef MASKING
          ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)*                      &
     &                        GRID(ng)%umask(i,Jend+1)
#  endif
        END DO

# elif defined NORTH_M2GRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=IstrU,Iend
          ubar(i,Jend+1,kout)=ubar(i,Jend,kout)
#  ifdef MASKING
          ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)*                      &
     &                        GRID(ng)%umask(i,Jend+1)
#  endif
        END DO

# else
!
!  Northern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
#  ifdef EW_PERIODIC
#   define I_RANGE IstrU,Iend
#  else
#   define I_RANGE Istr,IendR
#  endif
        DO i=I_RANGE
          ubar(i,Jend+1,kout)=gamma2(ng)*ubar(i,Jend,kout)
#  ifdef MASKING
          ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)*                      &
     &                        GRID(ng)%umask(i,Jend+1)
#  endif
        END DO
#  undef I_RANGE
# endif
      END IF
#endif
#if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        ubar(Istr,Jstr-1,kout)=0.5_r8*(ubar(Istr+1,Jstr-1,kout)+         &
     &                                 ubar(Istr  ,Jstr  ,kout))
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        ubar(Iend+1,Jstr-1,kout)=0.5_r8*(ubar(Iend  ,Jstr-1,kout)+       &
     &                                   ubar(Iend+1,Jstr  ,kout))
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        ubar(Istr,Jend+1,kout)=0.5_r8*(ubar(Istr+1,Jend+1,kout)+         &
     &                                 ubar(Istr  ,Jend  ,kout))
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        ubar(Iend+1,Jend+1,kout)=0.5_r8*(ubar(Iend  ,Jend+1,kout)+       &
     &                                   ubar(Iend+1,Jend  ,kout))
      END IF
#endif

      RETURN
      END SUBROUTINE u2dbc_tile
      END MODULE u2dbc_mod
