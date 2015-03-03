#include "cppdefs.h"
      MODULE us2dbc_mod
#ifdef WEC
!
!svn $Id: u2sdbc_im.F 779 2008-10-04 23:15:47Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for vertically     !
!  integrated Ustokes-velocity.                                              !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: us2dbc, us2dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE us2dbc (ng, tile)
!***********************************************************************
!
      USE mod_param
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
      CALL us2dbc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  OCEAN(ng) % ubar_stokes)
      RETURN
      END SUBROUTINE us2dbc
!
!***********************************************************************
      SUBROUTINE us2dbc_tile (ng, tile,                                 &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       ubar_stokes)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
# ifdef WET_DRY
      USE mod_grid
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ubar_stokes(LBi:,LBj:)
# else
      real(r8), intent(inout) :: ubar_stokes(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, know
      integer :: Imin, Imax

      real(r8), parameter :: eps = 1.0E-20_r8

      real(r8) :: Ce, Cx, cff5
      real(r8) :: bry_pgr, bry_cor, bry_str, bry_val
      real(r8) :: cff, cff1, cff2, dUde, dUdt, dUdx, tau

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!
!  Western edge, implicit upstream radiation condition.
!
        IF (LBC(iwest,isU2Sd,ng)%radiation) THEN
          DO j=Jstr,Jend+1
            grad(Istr  ,j)=ubar_stokes(Istr  ,j  )-                     &
     &                     ubar_stokes(Istr  ,j-1)
            grad(Istr+1,j)=ubar_stokes(Istr+1,j  )-                     &
     &                     ubar_stokes(Istr+1,j-1)
          END DO
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              dUdt=ubar_stokes(Istr+1,j)-ubar_stokes(Istr+1,j)
              dUdx=ubar_stokes(Istr+1,j)-ubar_stokes(Istr+2,j)
              IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
              IF ((dUdt*(grad(Istr+1,j)+grad(Istr+1,j+1))).gt.          &
     &          0.0_r8) THEN
                dUde=grad(Istr+1,j  )
              ELSE
                dUde=grad(Istr+1,j+1)
              END IF
              cff=MAX(dUdx*dUdx+dUde*dUde,eps)
              Cx=dUdt*dUdx
# ifdef RADIATION_2D
              Ce=MIN(cff,MAX(dUdt*dUde,-cff))
# else
              Ce=0.0_r8
# endif
              ubar_stokes(Istr,j)=(cff*ubar_stokes(Istr  ,j)+           &
     &                            Cx *ubar_stokes(Istr+1,j)-            &
     &                            MAX(Ce,0.0_r8)*grad(Istr,j  )-        &
     &                            MIN(Ce,0.0_r8)*grad(Istr,j+1))/       &
     &                           (cff+Cx)
# ifdef MASKING
              ubar_stokes(Istr,j)=ubar_stokes(Istr,j)*                  &
     &                            GRID(ng)%umask(Istr,j)
# endif
# ifdef WET_DRY
              ubar_stokes(Istr,j)=ubar_stokes(Istr,j)*                  &
     &                            GRID(ng)%umask_wet(Istr,j)
# endif
            END IF
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (LBC(iwest,isU2Sd,ng)%clamped) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              ubar_stokes(Istr,j)=BOUNDARY(ng)%ubarstokes_west(j)
# ifdef MASKING
              ubar_stokes(Istr,j)=ubar_stokes(Istr,j)*                  &
     &                            GRID(ng)%umask(Istr,j)
# endif
# ifdef WET_DRY
              ubar_stokes(Istr,j)=ubar_stokes(Istr,j)*                  &
     &                            GRID(ng)%umask_wet(Istr,j)
# endif
            END IF
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (LBC(iwest,isU2Sd,ng)%closed) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              ubar_stokes(Istr,j)=0.0_r8
            END IF
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (LBC(iwest,isU2Sd,ng)%gradient) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              ubar_stokes(Istr,j)=ubar_stokes(Istr+1,j)
# ifdef MASKING
              ubar_stokes(Istr,j)=ubar_stokes(Istr,j)*                  &
     &                            GRID(ng)%umask(Istr,j)
# endif
# ifdef WET_DRY
              ubar_stokes(Istr,j)=ubar_stokes(Istr,j)*                  &
     &                            GRID(ng)%umask_wet(Istr,j)
# endif
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
!
!  Eastern edge, implicit upstream radiation condition.
!
        IF (LBC(ieast,isU2Sd,ng)%radiation) THEN
          DO j=Jstr,Jend+1
            grad(Iend  ,j)=ubar_stokes(Iend  ,j  )-                     &
     &                     ubar_stokes(Iend  ,j-1)
            grad(Iend+1,j)=ubar_stokes(Iend+1,j  )-                     &
     &                     ubar_stokes(Iend+1,j-1)
          END DO
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              dUdt=ubar_stokes(Iend,j)-ubar_stokes(Iend  ,j)
              dUdx=ubar_stokes(Iend,j)-ubar_stokes(Iend-1,j)
              IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
              IF ((dUdt*(grad(Iend,j)+grad(Iend,j+1))).gt.0.0_r8) THEN
                dUde=grad(Iend,j)
              ELSE
                dUde=grad(Iend,j+1)
              END IF
              cff=MAX(dUdx*dUdx+dUde*dUde,eps)
              Cx=dUdt*dUdx
# ifdef RADIATION_2D
              Ce=MIN(cff,MAX(dUdt*dUde,-cff))
# else
              Ce=0.0_r8
# endif
              ubar_stokes(Iend+1,j)=(cff*ubar_stokes(Iend+1,j)+         &
     &                              Cx *ubar_stokes(Iend  ,j)-          &
     &                              MAX(Ce,0.0_r8)*grad(Iend+1,j  )-    &
     &                              MIN(Ce,0.0_r8)*grad(Iend+1,j+1))/   &
     &                             (cff+Cx)
# ifdef MASKING
              ubar_stokes(Iend+1,j)=ubar_stokes(Iend+1,j)*              &
     &                              GRID(ng)%umask(Iend+1,j)
# endif
# ifdef WET_DRY
              ubar_stokes(Iend+1,j)=ubar_stokes(Iend+1,j)*              &
     &                              GRID(ng)%umask_wet(Iend+1,j)
# endif
            END IF
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (LBC(ieast,isU2Sd,ng)%clamped) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              ubar_stokes(Iend+1,j)=BOUNDARY(ng)%ubarstokes_east(j)
# ifdef MASKING
              ubar_stokes(Iend+1,j)=ubar_stokes(Iend+1,j)*              &
     &                              GRID(ng)%umask(Iend+1,j)
# endif
# ifdef WET_DRY
              ubar_stokes(Iend+1,j)=ubar_stokes(Iend+1,j)*              &
     &                              GRID(ng)%umask_wet(Iend+1,j)
# endif
            END IF
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (LBC(ieast,isU2Sd,ng)%closed) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              ubar_stokes(Iend+1,j)=0.0_r8
            END IF
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (LBC(ieast,isU2Sd,ng)%gradient) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              ubar_stokes(Iend+1,j)=ubar_stokes(Iend,j)
# ifdef MASKING
              ubar_stokes(Iend+1,j)=ubar_stokes(Iend+1,j)*              &
     &                              GRID(ng)%umask(Iend+1,j)
# endif
# ifdef WET_DRY
              ubar_stokes(Iend+1,j)=ubar_stokes(Iend+1,j)*              &
     &                              GRID(ng)%umask_wet(Iend+1,j)
# endif
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
!
!  Southern edge, implicit upstream radiation condition.
!
        IF (LBC(isouth,isU2Sd,ng)%radiation) THEN
          DO i=IstrU-1,Iend
            grad(i,Jstr-1)=ubar_stokes(i+1,Jstr-1)-                     &
     &                     ubar_stokes(i  ,Jstr-1)
            grad(i,Jstr  )=ubar_stokes(i+1,Jstr  )-                     &
     &                     ubar_stokes(i  ,Jstr  )
          END DO
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              dUdt=ubar_stokes(i,Jstr)-ubar_stokes(i,Jstr  )
              dUde=ubar_stokes(i,Jstr)-ubar_stokes(i,Jstr+1)
              IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
              IF ((dUdt*(grad(i-1,Jstr)+grad(i,Jstr))).gt.0.0_r8) THEN
                dUdx=grad(i-1,Jstr)
              ELSE
                dUdx=grad(i  ,Jstr)
              END IF
              cff=MAX(dUdx*dUdx+dUde*dUde,eps)
# ifdef RADIATION_2D
              Cx=MIN(cff,MAX(dUdt*dUdx,-cff))
# else
              Cx=0.0_r8
# endif
              Ce=dUdt*dUde
              ubar_stokes(i,Jstr-1)=(cff*ubar_stokes(i,Jstr-1)+         &
     &                              Ce*ubar_stokes(i,Jstr  )-           &
     &                              MAX(Cx,0.0_r8)*grad(i-1,Jstr-1)-    &
     &                              MIN(Cx,0.0_r8)*grad(i  ,Jstr-1))/   &
     &                              (cff+Ce)
# ifdef MASKING
              ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr-1)*              &
     &                              GRID(ng)%umask(i,Jstr-1)
# endif
# ifdef WET_DRY
              ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr-1)*              &
     &                              GRID(ng)%umask_wet(i,Jstr-1)
# endif
            END IF
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (LBC(isouth,isU2Sd,ng)%clamped) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              ubar_stokes(i,Jstr-1)=BOUNDARY(ng)%ubarstokes_south(i)
# ifdef MASKING
              ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr-1)*              &
     &                              GRID(ng)%umask(i,Jstr-1)
# endif
# ifdef WET_DRY
              ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr-1)*              &
     &                              GRID(ng)%umask_wet(i,Jstr-1)
# endif
            END IF
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (LBC(isouth,isU2Sd,ng)%gradient) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr)
# ifdef MASKING
              ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr-1)*              &
     &                              GRID(ng)%umask(i,Jstr-1)
# endif
            END IF
          END DO
!
!  Southern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
        ELSE IF (LBC(isouth,isU2Sd,ng)%closed) THEN
          IF (EWperiodic(ng)) THEN
            Imin=IstrU
            Imax=Iend
          ELSE
            Imin=Istr
            Imax=IendR
          END IF
          DO i=Imin,Imax
            IF (LBC_apply(ng)%south(i)) THEN
              ubar_stokes(i,Jstr-1)=gamma2(ng)*ubar_stokes(i,Jstr)
# ifdef MASKING
              ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr-1)*              &
     &                              GRID(ng)%umask(i,Jstr-1)
# endif
# ifdef WET_DRY
              ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr-1)*              &
     &                              GRID(ng)%umask_wet(i,Jstr-1)
# endif
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
!
!  Northern edge, implicit upstream radiation condition.
!
        IF (LBC(inorth,isU2Sd,ng)%radiation) THEN
          DO i=IstrU-1,Iend
            grad(i,Jend  )=ubar_stokes(i+1,Jend  )-                     &
     &                     ubar_stokes(i  ,Jend  )
            grad(i,Jend+1)=ubar_stokes(i+1,Jend+1)-                     &
     &                     ubar_stokes(i  ,Jend+1)
          END DO
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              dUdt=ubar_stokes(i,Jend)-ubar_stokes(i,Jend  )
              dUde=ubar_stokes(i,Jend)-ubar_stokes(i,Jend-1)
              IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
              IF ((dUdt*(grad(i-1,Jend)+grad(i,Jend))).gt.0.0_r8) THEN
                dUdx=grad(i-1,Jend)
              ELSE
                dUdx=grad(i  ,Jend)
              END IF
              cff=MAX(dUdx*dUdx+dUde*dUde,eps)
# ifdef RADIATION_2D
              Cx=MIN(cff,MAX(dUdt*dUdx,-cff))
# else
              Cx=0.0_r8
# endif
              Ce=dUdt*dUde
              ubar_stokes(i,Jend+1)=(cff*ubar_stokes(i,Jend+1)+         &
     &                              Ce *ubar_stokes(i,Jend  )-          &
     &                              MAX(Cx,0.0_r8)*grad(i-1,Jend+1)-    &
     &                              MIN(Cx,0.0_r8)*grad(i  ,Jend+1))/   &
     &                              (cff+Ce)
# ifdef MASKING
              ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend+1)*              &
     &                              GRID(ng)%umask(i,Jend+1)
# endif
# ifdef WET_DRY
              ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend+1)*              &
     &                              GRID(ng)%umask_wet(i,Jend+1)
# endif
            END IF
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (LBC(inorth,isU2Sd,ng)%clamped) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              ubar_stokes(i,Jend+1)=BOUNDARY(ng)%ubarstokes_north(i)
# ifdef MASKING
              ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend+1)*              &
     &                              GRID(ng)%umask(i,Jend+1)
# endif
# ifdef WET_DRY
              ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend+1)*              &
     &                              GRID(ng)%umask_wet(i,Jend+1)
# endif
            END IF
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (LBC(inorth,isU2Sd,ng)%gradient) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend)
# ifdef MASKING
              ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend+1)*              &
     &                              GRID(ng)%umask(i,Jend+1)
# endif
            END IF
          END DO
!
!  Northern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
        ELSE IF (LBC(inorth,isU2Sd,ng)%closed) THEN
          IF (EWperiodic(ng)) THEN
            Imin=IstrU
            Imax=Iend
          ELSE
            Imin=Istr
            Imax=IendR
          END IF
          DO i=Imin,Imax
            IF (LBC_apply(ng)%north(i)) THEN
              ubar_stokes(i,Jend+1)=gamma2(ng)*ubar_stokes(i,Jend)
# ifdef MASKING
              ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend+1)*              &
     &                              GRID(ng)%umask(i,Jend+1)
# endif
# ifdef WET_DRY
              ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend+1)*              &
     &                              GRID(ng)%umask_wet(i,Jend+1)
# endif
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
              ubar_stokes(Istr,Jstr-1)=0.5_r8*                          &
     &                                 (ubar_stokes(Istr+1,Jstr-1)+     &
     &                                  ubar_stokes(Istr  ,Jstr  ))
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
              ubar_stokes(Iend+1,Jstr-1)=0.5_r8*                        &
     &                                   (ubar_stokes(Iend  ,Jstr-1)+   &
     &                                    ubar_stokes(Iend+1,Jstr  ))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
              ubar_stokes(Istr,Jend+1)=0.5_r8*                          &
     &                                 (ubar_stokes(Istr  ,Jend  )+     &
     &                                  ubar_stokes(Istr+1,Jend+1))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
              ubar_stokes(Iend+1,Jend+1)=0.5_r8*                        &
     &                                   (ubar_stokes(Iend+1,Jend  )+   &
     &                                    ubar_stokes(Iend  ,Jend+1))
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE us2dbc_tile
#endif
      END MODULE us2dbc_mod
