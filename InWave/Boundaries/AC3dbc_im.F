#include "cppdefs.h"
      MODULE AC3dbc_mod
#ifdef SOLVE3D
!
!svn $Id: t3dbc_im.F 732 2008-09-07 01:55:51Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for the ITRC-th    !
!  tracer field.                                                       !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: AC3dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE AC3dbc (ng, tile, nout)
!***********************************************************************
!
      USE mod_param
      USE mod_inwave_vars
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, nout
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL AC3dbc_tile (ng, tile,                                       &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 nstp(ng),nout,                                   &
     &                 WAVEP(ng)% AC)
      RETURN
      END SUBROUTINE AC3dbc

!
!***********************************************************************
      SUBROUTINE AC3dbc_tile (ng, tile,                                 &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nstp,nout,                                 &
     &                       AC)
!***********************************************************************

      USE mod_param
      USE mod_inwave_params
      USE mod_inwave_bound
      USE mod_inwave_vars
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nout
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: AC(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: AC(LBi:UBi,LBj:UBj,ND,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k, d_bnd

      real(r8), parameter :: eps =1.0E-20_r8

      real(r8) :: Ce, Cx, cff, dTde, dTdt, dTdx, tau, ramp

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad

# include "set_bounds.h"
# ifdef RAMP_INWAVE
        ramp=TANH((tdays(ng)-dstart)/0.07_r8)
# else
        ramp=1.0_r8
# endif
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!
!  Western edge, implicit upstream radiation condition.
!
        IF (LBC(iwest,isAC3d,ng)%radiation) THEN
          DO k=1,ND
            DO j=Jstr,Jend+1
              grad(Istr-1,j)=AC(Istr-1,j  ,k,nstp)-                     &
     &                       AC(Istr-1,j-1,k,nstp)
# ifdef MASKING
              grad(Istr-1,j)=grad(Istr-1,j)*                            &
     &                       GRID(ng)%vmask(Istr-1,j)
# endif
              grad(Istr  ,j)=AC(Istr  ,j  ,k,nstp)-                     &
     &                       AC(Istr  ,j-1,k,nstp)
# ifdef MASKING
              grad(Istr  ,j)=grad(Istr  ,j)*                            &
     &                       GRID(ng)%vmask(Istr  ,j)
# endif
            END DO
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                dTdt=AC(Istr,j,k,nstp)-AC(Istr  ,j,k,nout)
                dTdx=AC(Istr,j,k,nout)-AC(Istr+1,j,k,nout)

                IF ((dTdt*dTdx).lt.0.0_r8) dTdt=0.0_r8
                IF ((dTdt*(grad(Istr,j)+grad(Istr,j+1))).gt.0.0_r8) THEN
                  dTde=grad(Istr,j  )
                ELSE
                  dTde=grad(Istr,j+1)
                END IF
                cff=MAX(dTdx*dTdx+dTde*dTde,eps)
                Cx=dTdt*dTdx
                Ce=MIN(cff,MAX(dTdt*dTde,-cff))

                AC(Istr-1,j,k,nout)=(cff*AC(Istr-1,j,k,nstp)+             &
     &                                   Cx *AC(Istr  ,j,k,nout)-         &
     &                                   MAX(Ce,0.0_r8)*                  &
     &                                      grad(Istr-1,j  )-             &
     &                                   MIN(Ce,0.0_r8)*                  &
     &                                      grad(Istr-1,j+1))/            &
     &                                  (cff+Cx)

# ifdef MASKING
                AC(Istr-1,j,k,nout)=AC(Istr-1,j,k,nout)*                  &
     &                              GRID(ng)%rmask(Istr-1,j)
# endif
              END IF
            END DO
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (LBC(iwest,isAC3d,ng)%clamped) THEN
          DO k=1,ND
            DO d_bnd=1,WAVEB(ng)%ND_BND
              IF(WAVEB(ng)%WD_BND(d_bnd).EQ.WAVEG(ng)%wd(k))then
                DO j=Jstr,Jend
                  IF (LBC_apply(ng)%west(j)) THEN
                    AC(Istr-1,j,k,nout)=WAVEB(ng)%AC_west(j,d_bnd)
# ifdef MASKING
                    AC(Istr-1,j,k,nout)=AC(Istr-1,j,k,nout)*              &
     &                                  GRID(ng)%rmask(Istr-1,j)
# endif
                  END IF
                END DO
              END IF
            END DO
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (LBC(iwest,isAC3d,ng)%gradient) THEN
          DO k=1,ND
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                AC(Istr-1,j,k,nout)=AC(Istr,j,k,nout)
# ifdef MASKING
                AC(Istr-1,j,k,nout)=AC(Istr-1,j,k,nout)*                  &
     &                              GRID(ng)%rmask(Istr-1,j)
# endif
              END IF
            END DO
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (LBC(iwest,isAC3d,ng)%closed) THEN
          DO k=1,ND
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                AC(Istr-1,j,k,nout)=AC(Istr,j,k,nout)
# ifdef MASKING
                AC(Istr-1,j,k,nout)=AC(Istr-1,j,k,nout)*                  &
     &                              GRID(ng)%rmask(Istr-1,j)
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
        IF (LBC(ieast,isAC3d,ng)%radiation) THEN
          DO k=1,ND
            DO j=Jstr,Jend+1
              grad(Iend  ,j)=AC(Iend  ,j  ,k,nstp)-                     &
     &                       AC(Iend  ,j-1,k,nstp)
# ifdef MASKING
              grad(Iend  ,j)=grad(Iend  ,j)*                            &
     &                       GRID(ng)%vmask(Iend  ,j)
# endif
              grad(Iend+1,j)=AC(Iend+1,j  ,k,nstp)-                     &
     &                       AC(Iend+1,j-1,k,nstp)
# ifdef MASKING
              grad(Iend+1,j)=grad(Iend+1,j)*                            &
     &                       GRID(ng)%vmask(Iend+1,j)
# endif
            END DO
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                dTdt=AC(Iend,j,k,nstp)-AC(Iend  ,j,k,nout)
                dTdx=AC(Iend,j,k,nout)-AC(Iend-1,j,k,nout)
                IF ((dTdt*dTdx).lt.0.0_r8) dTdt=0.0_r8
                IF ((dTdt*(grad(Iend,j)+grad(Iend,j+1))).gt.0.0_r8) THEN
                  dTde=grad(Iend,j  )
                ELSE
                  dTde=grad(Iend,j+1)
                END IF
                cff=MAX(dTdx*dTdx+dTde*dTde,eps)
                Cx=dTdt*dTdx
                Ce=MIN(cff,MAX(dTdt*dTde,-cff))
                AC(Iend+1,j,k,nout)=(cff*AC(Iend+1,j,k,nstp)+           &
     &                                   Cx *AC(Iend  ,j,k,nout)-       &
     &                                   MAX(Ce,0.0_r8)*                &
     &                                      grad(Iend+1,j  )-           &
     &                                   MIN(Ce,0.0_r8)*                &
     &                                      grad(Iend+1,j+1))/          &
     &                                  (cff+Cx)

# ifdef MASKING
                AC(Iend+1,j,k,nout)=AC(Iend+1,j,k,nout)*                &
     &                              GRID(ng)%rmask(Iend+1,j)
# endif
              END IF
            END DO
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (LBC(ieast,isAC3d,ng)%clamped) THEN
          DO d_bnd=1,WAVEB(ng)%ND_BND
            DO k=1,ND
              IF(WAVEB(ng)%WD_BND(d_bnd).EQ.WAVEG(ng)%wd(k)) THEN
                IF (LBC_apply(ng)%east(j)) THEN
                  DO j=Jstr,Jend
                    AC(Iend+1,j,k,nout)=WAVEB(ng)%AC_east(j,d_bnd)
# ifdef MASKING
                    AC(Iend+1,j,k,nout)=AC(Iend+1,j,k,nout)*            &
     &                                  GRID(ng)%rmask(Iend+1,j)
# endif
                  END DO
                ENDIF
              ENDIF
            END DO
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (LBC(ieast,isAC3d,ng)%gradient) THEN
          DO k=1,ND
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                AC(Iend+1,j,k,nout)=AC(Iend,j,k,nout)
# ifdef MASKING
                AC(Iend+1,j,k,nout)=AC(Iend+1,j,k,nout)*                &
     &                              GRID(ng)%rmask(Iend+1,j)
# endif
              END IF
            END DO
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (LBC(ieast,isAC3d,ng)%closed) THEN
          DO k=1,ND
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                AC(Iend+1,j,k,nout)=AC(Iend,j,k,nout)
# ifdef MASKING
                AC(Iend+1,j,k,nout)=AC(Iend+1,j,k,nout)*                &
     &                              GRID(ng)%rmask(Iend+1,j)
# endif
              END IF
            END DO
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
        IF (LBC(isouth,isAC3d,ng)%radiation) THEN
          DO k=1,ND
            DO i=Istr,Iend+1
              grad(i,Jstr  )=AC(i  ,Jstr  ,k,nstp)-                     &
     &                       AC(i-1,Jstr  ,k,nstp)
# ifdef MASKING
              grad(i,Jstr  )=grad(i,Jstr  )*                            &
     &                       GRID(ng)%umask(i,Jstr  )
# endif
              grad(i,Jstr-1)=AC(i  ,Jstr-1,k,nstp)-                     &
     &                       AC(i-1,Jstr-1,k,nstp)
# ifdef MASKING
              grad(i,Jstr-1)=grad(i,Jstr-1)*                            &
     &                       GRID(ng)%umask(i,Jstr-1)
# endif
            END DO
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                dTdt=AC(i,Jstr,k,nstp)-AC(i,Jstr  ,k,nout)
                dTde=AC(i,Jstr,k,nout)-AC(i,Jstr+1,k,nout)
                IF ((dTdt*dTde).lt.0.0_r8) dTdt=0.0_r8
                IF ((dTdt*(grad(i,Jstr)+grad(i+1,Jstr))).gt.0.0_r8) THEN
                  dTdx=grad(i  ,Jstr)
                ELSE
                  dTdx=grad(i+1,Jstr)
                END IF
                cff=MAX(dTdx*dTdx+dTde*dTde,eps)
                Cx=MIN(cff,MAX(dTdt*dTdx,-cff))
                Ce=dTdt*dTde
                AC(i,Jstr-1,k,nout)=(cff*AC(i,Jstr-1,k,nstp)+           &
     &                                   Ce *AC(i,Jstr  ,k,nout)-       &
     &                                   MAX(Cx,0.0_r8)*                &
     &                                      grad(i  ,Jstr-1)-           &
     &                                   MIN(Cx,0.0_r8)*                &
     &                                      grad(i+1,Jstr-1))/          &
     &                                  (cff+Ce)
# ifdef MASKING
                AC(i,Jstr-1,k,nout)=AC(i,Jstr-1,k,nout)*                &
     &                              GRID(ng)%rmask(i,Jstr-1)
# endif
              END IF
            END DO
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (LBC(isouth,isAC3d,ng)%clamped) THEN
          DO d_bnd=1,WAVEB(ng)%ND_BND
            DO k=1,ND
              IF(WAVEB(ng)%WD_BND(d_bnd).EQ.WAVEG(ng)%wd(k))then
                DO i=Istr,Iend
                  IF (LBC_apply(ng)%south(i)) THEN
                    AC(i,Jstr-1,k,nout)=WAVEB(ng)%AC_south(i,d_bnd)*ramp
# ifdef MASKING
                    AC(i,Jstr-1,k,nout)=AC(i,Jstr-1,k,nout)*            &
     &                                  GRID(ng)%rmask(i,Jstr-1)
# endif
                  END IF
                END DO
              END IF
            END DO
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (LBC(isouth,isAC3d,ng)%gradient) THEN
          DO k=1,ND
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                AC(i,Jstr-1,k,nout)=AC(i,Jstr,k,nout)
# ifdef MASKING
                AC(i,Jstr-1,k,nout)=AC(i,Jstr-1,k,nout)*                &
     &                              GRID(ng)%rmask(i,Jstr-1)
# endif
              END IF
            END DO
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE IF (LBC(isouth,isAC3d,ng)%closed) THEN
          DO k=1,ND
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                AC(i,Jstr-1,k,nout)=AC(i,Jstr,k,nout)
# ifdef MASKING
                AC(i,Jstr-1,k,nout)=AC(i,Jstr-1,k,nout)*                &
     &                              GRID(ng)%rmask(i,Jstr-1)
# endif
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
        IF (LBC(inorth,isAC3d,ng)%radiation) THEN
          DO k=1,ND
            DO i=Istr,Iend+1
              grad(i,Jend  )=AC(i  ,Jend  ,k,nstp)-                     &
     &                       AC(i-1,Jend  ,k,nstp)
# ifdef MASKING
              grad(i,Jend  )=grad(i,Jend  )*                            &
     &                       GRID(ng)%umask(i,Jend  )
# endif
              grad(i,Jend+1)=AC(i  ,Jend+1,k,nstp)-                     &
     &                       AC(i-1,Jend+1,k,nstp)
# ifdef MASKING
              grad(i,Jend+1)=grad(i,Jend+1)*                            &
     &                       GRID(ng)%umask(i,Jend+1)
# endif
            END DO
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                dTdt=AC(i,Jend,k,nstp)-AC(i,Jend  ,k,nout)
                dTde=AC(i,Jend,k,nout)-AC(i,Jend-1,k,nout)
                IF ((dTdt*dTde).lt.0.0_r8) dTdt=0.0_r8
                IF ((dTdt*(grad(i,Jend)+grad(i+1,Jend))).gt.0.0_r8) THEN
                  dTdx=grad(i  ,Jend)
                ELSE
                  dTdx=grad(i+1,Jend)
                END IF
                cff=MAX(dTdx*dTdx+dTde*dTde,eps)
                Cx=MIN(cff,MAX(dTdt*dTdx,-cff))
                Ce=dTdt*dTde
                AC(i,Jend+1,k,nout)=(cff*AC(i,Jend+1,k,nstp)+           &
     &                                   Ce *AC(i,Jend  ,k,nout)-       &
     &                                   MAX(Cx,0.0_r8)*                &
     &                                      grad(i  ,Jend+1)-           &
     &                                   MIN(Cx,0.0_r8)*                &
     &                                      grad(i+1,Jend+1))/          &
     &                                  (cff+Ce)
# ifdef MASKING
                AC(i,Jend+1,k,nout)=AC(i,Jend+1,k,nout)*                &
     &                              GRID(ng)%rmask(i,Jend+1)
# endif
              END IF
            END DO
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (LBC(inorth,isAC3d,ng)%clamped) THEN
          DO d_bnd=1,WAVEB(ng)%ND_BND
            DO k=1,ND
              IF(WAVEB(ng)%WD_BND(d_bnd).EQ.WAVEG(ng)%wd(k))then
                DO i=Istr,Iend
                  IF (LBC_apply(ng)%north(i)) THEN
                    AC(i,Jend+1,k,nout)=WAVEB(ng)%AC_north(i,d_bnd)
# ifdef MASKING
                    AC(i,Jend+1,k,nout)=AC(i,Jend+1,k,nout)*            &
     &                                  GRID(ng)%rmask(i,Jend+1)
# endif
                  END IF
                END DO
              END IF
            END DO
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (LBC(inorth,isAC3d,ng)%gradient) THEN
          DO k=1,ND
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                AC(i,Jend+1,k,nout)=AC(i,Jend,k,nout)
# ifdef MASKING
                AC(i,Jend+1,k,nout)=AC(i,Jend+1,k,nout)*                &
     &                              GRID(ng)%rmask(i,Jend+1)
# endif
              END IF
            END DO
          END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (LBC(inorth,isAC3d,ng)%closed) THEN
          DO k=1,ND
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                AC(i,Jend+1,k,nout)=AC(i,Jend,k,nout)
# ifdef MASKING
                AC(i,Jend+1,k,nout)=AC(i,Jend+1,k,nout)*                &
     &                              GRID(ng)%rmask(i,Jend+1)
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
     &        LBC_apply(ng)%west (Jstr-1)) THEN
            DO k=1,ND
              AC(Istr-1,Jstr-1,k,nout)=0.5_r8*                          &
     &                                     (AC(Istr  ,Jstr-1,k,nout)+   &
     &                                      AC(Istr-1,Jstr  ,k,nout))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            DO k=1,ND
              AC(Iend+1,Jstr-1,k,nout)=0.5_r8*                          &
     &                                     (AC(Iend  ,Jstr-1,k,nout)+   &
     &                                      AC(Iend+1,Jstr  ,k,nout))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO k=1,ND
              AC(Istr-1,Jend+1,k,nout)=0.5_r8*                          &
     &                                     (AC(Istr-1,Jend  ,k,nout)+   &
     &                                      AC(Istr  ,Jend+1,k,nout))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO k=1,ND
              AC(Iend+1,Jend+1,k,nout)=0.5_r8*                          &
     &                                     (AC(Iend+1,Jend  ,k,nout)+   &
     &                                      AC(Iend  ,Jend+1,k,nout))
            END DO
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE AC3dbc_tile
#endif
      END MODULE AC3dbc_mod
