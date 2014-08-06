#include "cppdefs.h"
      MODULE vs3dbc_mod
#if defined WEC && defined SOLVE3D
!
!svn $Id: vs3dbc_im.F 732 2008-09-07 01:55:51Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for total 3D       !
!  Vstokes-velocity.                                                   !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: vs3dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE vs3dbc (ng, tile)
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
      CALL vs3dbc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, N(ng),                      &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  OCEAN(ng) % v_stokes)
      RETURN
      END SUBROUTINE vs3dbc
!
!***********************************************************************
      SUBROUTINE vs3dbc_tile (ng, tile,                                 &
     &                       LBi, UBi, LBj, UBj, UBk,                   &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       v_stokes)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: v_stokes(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: v_stokes(LBi:UBi,LBj:UBj,UBk)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k, Jmin, Jmax

      real(r8), parameter :: eps = 1.0E-20_r8

      real(r8) :: Ce, Cx, cff, dVde, dVdt, dVdx, tau

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
!
!  Southern edge, implicit upstream radiation condition.
!
        IF (LBC(isouth,isV3Sd,ng)%radiation) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend+1
              grad(i,Jstr  )=v_stokes(i  ,Jstr  ,k)-                    &
     &                       v_stokes(i-1,Jstr  ,k)
              grad(i,Jstr+1)=v_stokes(i  ,Jstr+1,k)-                    &
     &                       v_stokes(i-1,Jstr+1,k)
            END DO
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                dVdt=v_stokes(i,Jstr+1,k)-v_stokes(i,Jstr+1,k)
                dVde=v_stokes(i,Jstr+1,k)-v_stokes(i,Jstr+2,k)
                IF ((dVdt*dVde).lt.0.0_r8) dVdt=0.0_r8
                IF ((dVdt*(grad(i,Jstr+1)+grad(i+1,Jstr+1))).gt.        &
     &              0.0_r8) THEN
                  dVdx=grad(i  ,Jstr+1)
                ELSE
                  dVdx=grad(i+1,Jstr+1)
                END IF
                cff=MAX(dVdx*dVdx+dVde*dVde,eps)
# ifdef RADIATION_2D
                Cx=MIN(cff,MAX(dVdt*dVdx,-cff))
# else
                Cx=0.0_r8
# endif
                Ce=dVdt*dVde
                v_stokes(i,Jstr,k)=(cff*v_stokes(i,Jstr  ,k)+           &
     &                             Ce *v_stokes(i,Jstr+1,k)-            &
     &                             MAX(Cx,0.0_r8)*grad(i  ,Jstr)-       &
     &                             MIN(Cx,0.0_r8)*grad(i+1,Jstr))/      &
     &                             (cff+Ce)
# ifdef MASKING
                v_stokes(i,Jstr,k)=v_stokes(i,Jstr,k)*                  &
     &                             GRID(ng)%vmask(i,Jstr)
# endif
# ifdef WET_DRY
                v_stokes(i,Jstr,k)=v_stokes(i,Jstr,k)*                  &
     &                             GRID(ng)%vmask_wet(i,Jstr)
# endif
              END IF
            END DO
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (LBC(isouth,isV3Sd,ng)%clamped) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                v_stokes(i,Jstr,k)=BOUNDARY(ng)%vstokes_south(i,k)
# ifdef MASKING
                v_stokes(i,Jstr,k)=v_stokes(i,Jstr,k)*                  &
     &                             GRID(ng)%vmask(i,Jstr)
# endif
# ifdef WET_DRY
                v_stokes(i,Jstr,k)=v_stokes(i,Jstr,k)*                  &
     &                             GRID(ng)%vmask_wet(i,Jstr)
# endif
              END IF
            END DO
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (LBC(isouth,isV3Sd,ng)%gradient) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                v_stokes(i,Jstr,k)=v_stokes(i,Jstr+1,k)
# ifdef MASKING
                v_stokes(i,Jstr,k)=v_stokes(i,Jstr,k)*                  &
       &                           GRID(ng)%vmask(i,Jstr)
# endif
# ifdef WET_DRY
                v_stokes(i,Jstr,k)=v_stokes(i,Jstr,k)*                  &
     &                             GRID(ng)%vmask_wet(i,Jstr)
# endif
              END IF
            END DO
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE IF (LBC(isouth,isV3Sd,ng)%closed) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                v_stokes(i,Jstr,k)=0.0_r8
              END IF
            END DO
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
        IF (LBC(inorth,isV3Sd,ng)%radiation) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend+1
              grad(i,Jend  )=v_stokes(i  ,Jend  ,k)-                    &
     &                       v_stokes(i-1,Jend  ,k)
              grad(i,Jend+1)=v_stokes(i  ,Jend+1,k)-                    &
     &                       v_stokes(i-1,Jend+1,k)
            END DO
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                dVdt=v_stokes(i,Jend,k)-v_stokes(i,Jend  ,k)
                dVde=v_stokes(i,Jend,k)-v_stokes(i,Jend-1,k)
                IF ((dVdt*dVde).lt.0.0_r8) dVdt=0.0_r8
                IF ((dVdt*(grad(i,Jend)+grad(i+1,Jend))).gt.0.0_r8) THEN
                  dVdx=grad(i  ,Jend)
                ELSE
                  dVdx=grad(i+1,Jend)
                END IF
                cff=MAX(dVdx*dVdx+dVde*dVde,eps)
# ifdef RADIATION_2D
                Cx=MIN(cff,MAX(dVdt*dVdx,-cff))
# else
                Cx=0.0_r8
# endif
                Ce=dVdt*dVde
                v_stokes(i,Jend+1,k)=(cff*v_stokes(i,Jend+1,k)+         &
     &                               Ce *v_stokes(i,Jend  ,k)-          &
     &                               MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-   &
     &                               MIN(Cx,0.0_r8)*grad(i+1,Jend+1))/  &
     &                               (cff+Ce)
# ifdef MASKING
                v_stokes(i,Jend+1,k)=v_stokes(i,Jend+1,k)*              &
     &                               GRID(ng)%vmask(i,Jend+1)
# endif
# ifdef WET_DRY
                v_stokes(i,Jend+1,k)=v_stokes(i,Jend+1,k)*              &
     &                               GRID(ng)%vmask_wet(i,Jend+1)
# endif
              END IF
            END DO
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (LBC(inorth,isV3Sd,ng)%clamped) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                v_stokes(i,Jend+1,k)=BOUNDARY(ng)%v_north(i,k)
# ifdef MASKING
                v_stokes(i,Jend+1,k)=v_stokes(i,Jend+1,k)*              &
     &                               GRID(ng)%vmask(i,Jend+1)
# endif
# ifdef WET_DRY
                v_stokes(i,Jend+1,k)=v_stokes(i,Jend+1,k)*              &
     &                               GRID(ng)%vmask_wet(i,Jend+1)
# endif
              END IF
            END DO
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (LBC(inorth,isV3Sd,ng)%gradient) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                v_stokes(i,Jend+1,k)=v_stokes(i,Jend,k)
# ifdef MASKING
                v_stokes(i,Jend+1,k)=v_stokes(i,Jend+1,k)*              &
     &                                GRID(ng)%vmask(i,Jend+1)
# endif
# ifdef WET_DRY
                v_stokes(i,Jend+1,k)=v_stokes(i,Jend+1,k)*              &
     &                               GRID(ng)%vmask_wet(i,Jend+1)
# endif
            END IF
          END DO
        END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (LBC(inorth,isV3Sd,ng)%closed) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                v_stokes(i,Jend+1,k)=0.0_r8
              END IF
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!
!  Western edge, implicit upstream radiation condition.
!
        IF (LBC(iwest,isV3Sd,ng)%radiation) THEN
          DO k=1,N(ng)
            DO j=JstrV-1,Jend
              grad(Istr-1,j)=v_stokes(Istr-1,j+1,k)-                    &
     &                       v_stokes(Istr-1,j  ,k)
              grad(Istr  ,j)=v_stokes(Istr  ,j+1,k)-                    &
     &                       v_stokes(Istr  ,j  ,k)
            END DO
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                dVdt=v_stokes(Istr,j,k)-v_stokes(Istr  ,j,k)
                dVdx=v_stokes(Istr,j,k)-v_stokes(Istr+1,j,k)
                IF ((dVdt*dVdx).lt.0.0_r8) dVdt=0.0_r8
                IF ((dVdt*(grad(Istr,j-1)+grad(Istr,j))).gt.0.0_r8) THEN
                  dVde=grad(Istr,j-1)
                ELSE
                  dVde=grad(Istr,j  )
                END IF
                cff=MAX(dVdx*dVdx+dVde*dVde,eps)
                Cx=dVdt*dVdx
# ifdef RADIATION_2D
                Ce=MIN(cff,MAX(dVdt*dVde,-cff))
# else
                Ce=0.0_r8
# endif
                v_stokes(Istr-1,j,k)=(cff*v_stokes(Istr-1,j,k)+         &
     &                                Cx *v_stokes(Istr  ,j,k)-         &
     &                                MAX(Ce,0.0_r8)*grad(Istr-1,j-1)-  &
     &                                MIN(Ce,0.0_r8)*grad(Istr-1,j  ))/ &
     &                                (cff+Cx)
# ifdef MASKING
                v_stokes(Istr-1,j,k)=v_stokes(Istr-1,j,k)*              &
     &                               GRID(ng)%vmask(Istr-1,j)
# endif
# ifdef WET_DRY
                v_stokes(Istr-1,j,k)=v_stokes(Istr-1,j,k)*              &
     &                               GRID(ng)%vmask_wet(Istr-1,j)
# endif
              END IF
            END DO
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (LBC(iwest,isV3Sd,ng)%clamped) THEN
          DO k=1,N(ng)
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                v_stokes(Istr-1,j,k)=BOUNDARY(ng)%vstokes_west(j,k)
# ifdef MASKING
                v_stokes(Istr-1,j,k)=v_stokes(Istr-1,j,k)*              &
     &                               GRID(ng)%vmask(Istr-1,j)
# endif
# ifdef WET_DRY
                v_stokes(Istr-1,j,k)=v_stokes(Istr-1,j,k)*              &
     &                               GRID(ng)%vmask_wet(Istr-1,j)
# endif
              END IF
            END DO
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (LBC(iwest,isV3Sd,ng)%gradient) THEN
          DO k=1,N(ng)
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                v_stokes(Istr-1,j,k)=v_stokes(Istr,j,k)
# ifdef MASKING
                v_stokes(Istr-1,j,k)=v_stokes(Istr-1,j,k)*              &
     &                               GRID(ng)%vmask(Istr-1,j)
# endif
# ifdef WET_DRY
                v_stokes(Istr-1,j,k)=v_stokes(Istr-1,j,k)*              &
     &                               GRID(ng)%vmask_wet(Istr-1,j)
# endif
              END IF
            END DO
          END DO
!
!  Western edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
        ELSE IF (LBC(iwest,isV3Sd,ng)%closed) THEN
          IF (NSperiodic(ng)) THEN
            Jmin=JstrV
            Jmax=Jend
          ELSE
            Jmin=Jstr
            Jmax=JendR
          END IF
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              IF (LBC_apply(ng)%west(j)) THEN
                v_stokes(Istr-1,j,k)=gamma2(ng)*v_stokes(Istr,j,k)
# ifdef MASKING
                v_stokes(Istr-1,j,k)=v_stokes(Istr-1,j,k)*              &
     &                               GRID(ng)%vmask(Istr-1,j)
# endif
# ifdef WET_DRY
                v_stokes(Istr-1,j,k)=v_stokes(Istr-1,j,k)*              &
     &                               GRID(ng)%vmask_wet(Istr-1,j)
# endif
              END IF
            END DO
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
        IF (LBC(ieast,isV3Sd,ng)%radiation) THEN
          DO k=1,N(ng)
            DO j=JstrV-1,Jend
              grad(Iend  ,j)=v_stokes(Iend  ,j+1,k)-                    &
     &                       v_stokes(Iend  ,j  ,k)
              grad(Iend+1,j)=v_stokes(Iend+1,j+1,k)-                    &
     &                       v_stokes(Iend+1,j  ,k)
            END DO
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                dVdt=v_stokes(Iend,j,k)-v_stokes(Iend  ,j,k)
                dVdx=v_stokes(Iend,j,k)-v_stokes(Iend-1,j,k)
                IF ((dVdt*dVdx).lt.0.0_r8) dVdt=0.0_r8
                IF ((dVdt*(grad(Iend,j-1)+grad(Iend,j))).gt.0.0_r8) THEN
                  dVde=grad(Iend,j-1)
                ELSE
                  dVde=grad(Iend,j  )
                END IF
                cff=MAX(dVdx*dVdx+dVde*dVde,eps)
                Cx=dVdt*dVdx
# ifdef RADIATION_2D
                Ce=MIN(cff,MAX(dVdt*dVde,-cff))
# else
                Ce=0.0_r8
# endif
                v_stokes(Iend+1,j,k)=(cff*v_stokes(Iend+1,j,k)+         &
     &                               Cx *v_stokes(Iend  ,j,k)-          &
     &                               MAX(Ce,0.0_r8)*grad(Iend+1,j-1)-   &
     &                               MIN(Ce,0.0_r8)*grad(Iend+1,j  ))/  &
     &                               (cff+Cx)
# ifdef MASKING
                v_stokes(Iend+1,j,k)=v_stokes(Iend+1,j,k)*              &
     &                               GRID(ng)%vmask(Iend+1,j)
# endif
# ifdef WET_DRY
                v_stokes(Iend+1,j,k)=v_stokes(Iend+1,j,k)*              &
     &                               GRID(ng)%vmask_wet(Iend+1,j)
# endif
              END IF
            END DO
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (LBC(ieast,isV3Sd,ng)%clamped) THEN
          DO k=1,N(ng)
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                v_stokes(Iend+1,j,k)=BOUNDARY(ng)%vstokes_east(j,k)
# ifdef MASKING
                v_stokes(Iend+1,j,k)=v_stokes(Iend+1,j,k)*              &
     &                               GRID(ng)%vmask(Iend+1,j)
# endif
# ifdef WET_DRY
                v_stokes(Iend+1,j,k)=v_stokes(Iend+1,j,k)*              &
     &                               GRID(ng)%vmask_wet(Iend+1,j)
# endif
              END IF
            END DO
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (LBC(ieast,isV3Sd,ng)%gradient) THEN
          DO k=1,N(ng)
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                v_stokes(Iend+1,j,k)=v_stokes(Iend,j,k)
# ifdef MASKING
                v_stokes(Iend+1,j,k)=v_stokes(Iend+1,j,k)*              &
     &                               GRID(ng)%vmask(Iend+1,j)
# endif
# ifdef WET_DRY
                v_stokes(Iend+1,j,k)=v_stokes(Iend+1,j,k)*              &
     &                               GRID(ng)%vmask_wet(Iend+1,j)
# endif
              END IF
            END DO
          END DO
!
!  Eastern edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
        ELSE IF (LBC(ieast,isV3Sd,ng)%closed) THEN
          IF (NSperiodic(ng)) THEN
            Jmin=JstrV
            Jmax=Jend
          ELSE
            Jmin=Jstr
            Jmax=JendR
          END IF
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              IF (LBC_apply(ng)%east(j)) THEN
                v_stokes(Iend+1,j,k)=gamma2(ng)*v_stokes(Iend,j,k)
# ifdef MASKING
                v_stokes(Iend+1,j,k)=v_stokes(Iend+1,j,k)*              &
     &                               GRID(ng)%vmask(Iend+1,j)
# endif
# ifdef WET_DRY
                v_stokes(Iend+1,j,k)=v_stokes(Iend+1,j,k)*              &
     &                               GRID(ng)%vmask_wet(Iend+1,j)
# endif
              END IF
            END DO
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
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr  )) THEN
            DO k=1,N(ng)
              v_stokes(Istr-1,Jstr,k)=0.5_r8*                           &
     &                                (v_stokes(Istr  ,Jstr  ,k)+       &
     &                                 v_stokes(Istr-1,Jstr+1,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr  )) THEN
            DO k=1,N(ng)
              v_stokes(Iend+1,Jstr,k)=0.5_r8*                           &
     &                                (v_stokes(Iend  ,Jstr  ,k)+       &
     &                                 v_stokes(Iend+1,Jstr+1,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO k=1,N(ng)
              v_stokes(Istr-1,Jend+1,k)=0.5_r8*                         &
     &                                  (v_stokes(Istr-1,Jend  ,k)+     &
     &                                   v_stokes(Istr  ,Jend+1,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO k=1,N(ng)
              v_stokes(Iend+1,Jend+1,k)=0.5_r8*                         &
     &                                  (v_stokes(Iend+1,Jend  ,k)+     &
     &                                   v_stokes(Iend  ,Jend+1,k))
            END DO
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE vs3dbc_tile
#endif
      END MODULE vs3dbc_mod
