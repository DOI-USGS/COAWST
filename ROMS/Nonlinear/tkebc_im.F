#include "cppdefs.h"
      MODULE tkebc_mod
#if defined SOLVE3D && (defined MY25_MIXING || defined GLS_MIXING)
!
!svn $Id: tkebc_im.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine sets lateral boundary conditions for turbulent      !
!  kinetic energy and turbulent length scale variables associated      !
!  with the Mellor and Yamada or GOTM closures.                        !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: tkebc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE tkebc (ng, tile, nout)
!***********************************************************************
!
      USE mod_param
      USE mod_mixing
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
      CALL tkebc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 nout, nstp(ng),                                  &
     &                 MIXING(ng)% gls,                                 &
     &                 MIXING(ng)% tke)
      RETURN
      END SUBROUTINE tkebc
!
!***********************************************************************
      SUBROUTINE tkebc_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, UBk,                   &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nout, nstp,                                &
     &                       gls, tke)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nout, nstp
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: gls(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: tke(LBi:,LBj:,0:,:)
# else
      real(r8), intent(inout) :: gls(LBi:UBi,LBj:UBj,0:UBk,3)
      real(r8), intent(inout) :: tke(LBi:UBi,LBj:UBj,0:UBk,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8), parameter :: eps = 1.0e-20_r8

      real(r8) :: Ce, Cx, cff, dKde, dKdt, dKdx

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: gradL

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!
!  Western edge, implicit upstream radiation condition.
!
        IF (LBC(iwest,isMtke,ng)%radiation) THEN
          DO k=0,N(ng)
            DO j=Jstr,Jend+1
              grad(Istr-1,j)=tke(Istr-1,j  ,k,nstp)-                    &
     &                       tke(Istr-1,j-1,k,nstp)
# ifdef MASKING
              grad(Istr-1,j)=grad(Istr-1,j)*                            &
     &                       GRID(ng)%vmask(Istr-1,j)
# endif
              grad(Istr  ,j)=tke(Istr  ,j  ,k,nstp)-                    &
     &                       tke(Istr  ,j-1,k,nstp)
# ifdef MASKING
              grad(Istr  ,j)=grad(Istr  ,j)*                            &
     &                       GRID(ng)%vmask(Istr  ,j)
# endif
              gradL(Istr-1,j)=gls(Istr-1,j  ,k,nstp)-                   &
     &                        gls(Istr-1,j-1,k,nstp)
# ifdef MASKING
              gradL(Istr-1,j)=gradL(Istr-1,j)*                          &
     &                        GRID(ng)%vmask(Istr-1,j)
# endif
              gradL(Istr  ,j)=gls(Istr  ,j  ,k,nstp)-                   &
     &                        gls(Istr  ,j-1,k,nstp)
# ifdef MASKING
              gradL(Istr  ,j)=gradL(Istr  ,j)*                          &
     &                        GRID(ng)%vmask(Istr  ,j)
# endif
            END DO
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                dKdt=tke(Istr,j,k,nstp)-tke(Istr  ,j,k,nout)
                dKdx=tke(Istr,j,k,nout)-tke(Istr+1,j,k,nout)
                IF ((dKdt*dKdx).lt.0.0_r8) dKdt=0.0_r8
                IF ((dKdt*(grad(Istr,j  )+                              &
     &                     grad(Istr,j+1))).gt.0.0_r8) THEN
                  dKde=grad(Istr,j  )
                ELSE
                  dKde=grad(Istr,j+1)
                END IF
                cff=MAX(dKdx*dKdx+dKde*dKde,eps)
                Cx=dKdt*dKdx
# ifdef RADIATION_2D
                Ce=MIN(cff,MAX(dKdt*dKde,-cff))
# else
                Ce=0.0_r8
# endif
                tke(Istr-1,j,k,nout)=(cff*tke(Istr-1,j,k,nstp)+         &
     &                                Cx *tke(Istr  ,j,k,nout)-         &
     &                                MAX(Ce,0.0_r8)*                   &
     &                                   grad(Istr-1,j  )-              &
     &                                MIN(Ce,0.0_r8)*                   &
     &                                   grad(Istr-1,j+1))/             &
     &                               (cff+Cx)
# ifdef MASKING
                tke(Istr-1,j,k,nout)=tke(Istr-1,j,k,nout)*              &
     &                               GRID(ng)%rmask(Istr-1,j)
# endif
                dKdt=gls(Istr,j,k,nstp)-gls(Istr  ,j,k,nout)
                dKdx=gls(Istr,j,k,nout)-gls(Istr+1,j,k,nout)
                IF ((dKdt*dKdx).lt.0.0_r8) dKdt=0.0_r8
                IF ((dKdt*(gradL(Istr,j  )+                             &
     &                     gradL(Istr,j+1))).gt.0.0_r8) THEN
                  dKde=gradL(Istr,j  )
                ELSE
                  dKde=gradL(Istr,j+1)
                END IF
                cff=MAX(dKdx*dKdx+dKde*dKde,eps)
                Cx=dKdt*dKdx
# ifdef RADIATION_2D
                Ce=MIN(cff,MAX(dKdt*dKde,-cff))
# else
                Ce=0.0_r8
# endif
                gls(Istr-1,j,k,nout)=(cff*gls(Istr-1,j,k,nstp)+         &
     &                                Cx *gls(Istr  ,j,k,nout)-         &
     &                                MAX(Ce,0.0_r8)*                   &
     &                                   gradL(Istr-1,j  )-             &
     &                                MIN(Ce,0.0_r8)*                   &
     &                                   gradL(Istr-1,j+1))/            &
     &                               (cff+Cx)
# ifdef MASKING
                gls(Istr-1,j,k,nout)=gls(Istr-1,j,k,nout)*              &
     &                               GRID(ng)%rmask(Istr-1,j)
# endif
              END IF
            END DO
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (LBC(iwest,isMtke,ng)%gradient) THEN
          DO k=0,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                tke(Istr-1,j,k,nout)=tke(Istr,j,k,nout)
# ifdef MASKING
                tke(Istr-1,j,k,nout)=tke(Istr-1,j,k,nout)*              &
     &                               GRID(ng)%rmask(Istr-1,j)
# endif
                gls(Istr-1,j,k,nout)=gls(Istr,j,k,nout)
# ifdef MASKING
                gls(Istr-1,j,k,nout)=gls(Istr-1,j,k,nout)*              &
     &                               GRID(ng)%rmask(Istr-1,j)
# endif
              END IF
            END DO
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (LBC(iwest,isMtke,ng)%closed) THEN
          DO k=0,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                tke(Istr-1,j,k,nout)=tke(Istr,j,k,nout)
# ifdef MASKING
                tke(Istr-1,j,k,nout)=tke(Istr-1,j,k,nout)*              &
     &                               GRID(ng)%rmask(Istr-1,j)
# endif
                gls(Istr-1,j,k,nout)=gls(Istr,j,k,nout)
# ifdef MASKING
                gls(Istr-1,j,k,nout)=gls(Istr-1,j,k,nout)*              &
     &                               GRID(ng)%rmask(Istr-1,j)
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
        IF (LBC(ieast,isMtke,ng)%radiation) THEN
          DO k=0,N(ng)
            DO j=Jstr,Jend+1
              grad(Iend  ,j)=tke(Iend  ,j  ,k,nstp)-                    &
     &                       tke(Iend  ,j-1,k,nstp)
# ifdef MASKING
              grad(Iend  ,j)=grad(Iend  ,j)*                            &
     &                       GRID(ng)%vmask(Iend  ,j)
# endif
              grad(Iend+1,j)=tke(Iend+1,j  ,k,nstp)-                    &
     &                       tke(Iend+1,j-1,k,nstp)
# ifdef MASKING
              grad(Iend+1,j)=grad(Iend+1,j)*                            &
     &                       GRID(ng)%vmask(Iend+1,j)
# endif
              gradL(Iend  ,j)=gls(Iend  ,j  ,k,nstp)-                   &
     &                        gls(Iend  ,j-1,k,nstp)
# ifdef MASKING
              gradL(Iend  ,j)=gradL(Iend  ,j)*                          &
     &                        GRID(ng)%vmask(Iend  ,j)
# endif
              gradL(Iend+1,j)=gls(Iend+1,j  ,k,nstp)-                   &
     &                        gls(Iend+1,j-1,k,nstp)
# ifdef MASKING
              gradL(Iend+1,j)=gradL(Iend+1,j)*                          &
     &                        GRID(ng)%vmask(Iend+1,j)
# endif
            END DO
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                dKdt=tke(Iend,j,k,nstp)-tke(Iend  ,j,k,nout)
                dKdx=tke(Iend,j,k,nout)-tke(Iend-1,j,k,nout)
                IF ((dKdt*dKdx).lt.0.0_r8) dKdt=0.0_r8
                IF ((dKdt*(grad(Iend,j  )+                              &
     &                     grad(Iend,j+1))).gt.0.0_r8) THEN
                  dKde=grad(Iend,j  )
                ELSE
                  dKde=grad(Iend,j+1)
                END IF
                cff=MAX(dKdx*dKdx+dKde*dKde,eps)
                Cx=dKdt*dKdx
# ifdef RADIATION_2D
                Ce=MIN(cff,MAX(dKdt*dKde,-cff))
# else
                Ce=0.0_r8
# endif
                tke(Iend+1,j,k,nout)=(cff*tke(Iend+1,j,k,nstp)+         &
     &                                Cx *tke(Iend  ,j,k,nout)-         &
     &                                MAX(Ce,0.0_r8)*                   &
     &                                   grad(Iend+1,j  )-              &
     &                                MIN(Ce,0.0_r8)*                   &
     &                                   grad(Iend+1,j+1))/             &
     &                               (cff+Cx)
# ifdef MASKING
                tke(Iend+1,j,k,nout)=tke(Iend+1,j,k,nout)*              &
     &                               GRID(ng)%rmask(Iend+1,j)
# endif
                dKdt=gls(Iend,j,k,nstp)-gls(Iend  ,j,k,nout)
                dKdx=gls(Iend,j,k,nout)-gls(Iend-1,j,k,nout)
                IF ((dKdt*dKdx).lt.0.0_r8) dKdt=0.0_r8
                IF ((dKdt*(gradL(Iend,j  )+                             &
     &                     gradL(Iend,j+1))).gt.0.0_r8) THEN
                  dKde=gradL(Iend,j  )
                ELSE
                  dKde=gradL(Iend,j+1)
                END IF
                cff=MAX(dKdx*dKdx+dKde*dKde,eps)
                Cx=dKdt*dKdx
# ifdef RADIATION_2D
                Ce=MIN(cff,MAX(dKdt*dKde,-cff))
# else
                Ce=0.0_r8
# endif
                gls(Iend+1,j,k,nout)=(cff*gls(Iend+1,j,k,nstp)+         &
     &                                Cx *gls(Iend  ,j,k,nout)-         &
     &                                MAX(Ce,0.0_r8)*                   &
     &                                   gradL(Iend+1,j  )-             &
     &                                MIN(Ce,0.0_r8)*                   &
     &                                   gradL(Iend+1,j+1))/            &
     &                               (cff+Cx)
# ifdef MASKING
                gls(Iend+1,j,k,nout)=gls(Iend+1,j,k,nout)*              &
     &                               GRID(ng)%rmask(Iend+1,j)
# endif
              END IF
            END DO
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (LBC(ieast,isMtke,ng)%gradient) THEN
          DO k=0,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                tke(Iend+1,j,k,nout)=tke(Iend,j,k,nout)
# ifdef MASKING
                tke(Iend+1,j,k,nout)=tke(Iend+1,j,k,nout)*              &
     &                               GRID(ng)%rmask(Iend+1,j)
# endif
                gls(Iend+1,j,k,nout)=gls(Iend,j,k,nout)
# ifdef MASKING
                gls(Iend+1,j,k,nout)=gls(Iend+1,j,k,nout)*              &
     &                               GRID(ng)%rmask(Iend+1,j)
# endif
              END IF
            END DO
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (LBC(ieast,isMtke,ng)%closed) THEN
          DO k=0,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                tke(Iend+1,j,k,nout)=tke(Iend,j,k,nout)
# ifdef MASKING
                tke(Iend+1,j,k,nout)=tke(Iend+1,j,k,nout)*              &
     &                               GRID(ng)%rmask(Iend+1,j)
# endif
                gls(Iend+1,j,k,nout)=gls(Iend,j,k,nout)
# ifdef MASKING
                gls(Iend+1,j,k,nout)=gls(Iend+1,j,k,nout)*              &
     &                               GRID(ng)%rmask(Iend+1,j)
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
        IF (LBC(isouth,isMtke,ng)%radiation) THEN
          DO k=0,N(ng)
            DO i=Istr,Iend+1
              grad(i,Jstr  )=tke(i  ,Jstr  ,k,nstp)-                    &
     &                       tke(i-1,Jstr  ,k,nstp)
# ifdef MASKING
              grad(i,Jstr  )=grad(i,Jstr  )*GRID(ng)%umask(i,Jstr  )
# endif
              grad(i,Jstr-1)=tke(i  ,Jstr-1,k,nstp)-                    &
     &                     tke(i-1,Jstr-1,k,nstp)
# ifdef MASKING
              grad(i,Jstr-1)=grad(i,Jstr-1)*GRID(ng)%umask(i,Jstr-1)
# endif
              gradL(i,Jstr  )=gls(i  ,Jstr  ,k,nstp)-                   &
     &                      gls(i-1,Jstr  ,k,nstp)
# ifdef MASKING
              gradL(i,Jstr  )=gradL(i,Jstr  )*GRID(ng)%umask(i,Jstr  )
# endif
              gradL(i,Jstr-1)=gls(i  ,Jstr-1,k,nstp)-                   &
     &                      gls(i-1,Jstr-1,k,nstp)
# ifdef MASKING
              gradL(i,Jstr-1)=gradL(i,Jstr-1)*GRID(ng)%umask(i,Jstr-1)
# endif
            END DO
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                dKdt=tke(i,Jstr,k,nstp)-tke(i,Jstr  ,k,nout)
                dKde=tke(i,Jstr,k,nout)-tke(i,Jstr+1,k,nout)
                IF ((dKdt*dKde).lt.0.0_r8) dKdt=0.0_r8
                IF ((dKdt*(grad(i  ,Jstr)+                              &
     &                     grad(i+1,Jstr))).gt.0.0_r8) THEN
                  dKdx=grad(i  ,Jstr)
                ELSE
                  dKdx=grad(i+1,Jstr)
                END IF
                cff=MAX(dKdx*dKdx+dKde*dKde, eps)
# ifdef RADIATION_2D
                Cx=MIN(cff,MAX(dKdt*dKdx,-cff))
# else
                Cx=0.0_r8
# endif
                Ce=dKdt*dKde
                tke(i,Jstr-1,k,nout)=(cff*tke(i,Jstr-1,k,nstp)+         &
     &                                Ce *tke(i,Jstr  ,k,nout)-         &
     &                                MAX(Cx,0.0_r8)*                   &
     &                                   grad(i  ,Jstr-1)-              &
     &                                MIN(Cx,0.0_r8)*                   &
     &                                   grad(i+1,Jstr-1))/             &
     &                               (cff+Ce)
# ifdef MASKING
                tke(i,Jstr-1,k,nout)=tke(i,Jstr-1,k,nout)*              &
     &                               GRID(ng)%rmask(i,Jstr-1)
# endif
                dKdt=gls(i,Jstr,k,nstp)-gls(i,Jstr  ,k,nout)
                dKde=gls(i,Jstr,k,nout)-gls(i,Jstr+1,k,nout)
                IF ((dKdt*dKde).lt.0.0_r8) dKdt=0.0_r8
                IF ((dKdt*(gradL(i  ,Jstr)+                             &
     &                     gradL(i+1,Jstr))).gt.0.0_r8) THEN
                  dKdx=gradL(i  ,Jstr)
                ELSE
                  dKdx=gradL(i+1,Jstr)
                END IF
                cff=MAX(dKdx*dKdx+dKde*dKde,eps)
# ifdef RADIATION_2D
                Cx=MIN(cff,MAX(dKdt*dKdx,-cff))
# else
                Cx=0.0_r8
# endif
                Ce=dKdt*dKde
                gls(i,Jstr-1,k,nout)=(cff*gls(i,Jstr-1,k,nstp)+         &
     &                                Ce *gls(i,Jstr  ,k,nout)-         &
     &                                MAX(Cx,0.0_r8)*                   &
     &                                   gradL(i  ,Jstr-1)-             &
     &                                MIN(Cx,0.0_r8)*                   &
     &                                   gradL(i+1,Jstr-1))/            &
     &                               (cff+Ce)
# ifdef MASKING
                gls(i,Jstr-1,k,nout)=gls(i,Jstr-1,k,nout)*              &
     &                               GRID(ng)%rmask(i,Jstr-1)
# endif
              END IF
            END DO
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (LBC(isouth,isMtke,ng)%gradient) THEN
          DO k=0,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                tke(i,Jstr-1,k,nout)=tke(i,Jstr,k,nout)
# ifdef MASKING
                tke(i,Jstr-1,k,nout)=tke(i,Jstr-1,k,nout)*              &
     &                               GRID(ng)%rmask(i,Jstr-1)
# endif
                gls(i,Jstr-1,k,nout)=gls(i,Jstr,k,nout)
# ifdef MASKING
                gls(i,Jstr-1,k,nout)=gls(i,Jstr-1,k,nout)*              &
     &                               GRID(ng)%rmask(i,Jstr-1)
# endif
              END IF
            END DO
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE IF (LBC(isouth,isMtke,ng)%closed) THEN
          DO k=0,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                tke(i,Jstr-1,k,nout)=tke(i,Jstr,k,nout)
# ifdef MASKING
                tke(i,Jstr-1,k,nout)=tke(i,Jstr-1,k,nout)*              &
     &                               GRID(ng)%rmask(i,Jstr-1)
# endif
                gls(i,Jstr-1,k,nout)=gls(i,Jstr,k,nout)
# ifdef MASKING
                gls(i,Jstr-1,k,nout)=gls(i,Jstr-1,k,nout)*              &
     &                               GRID(ng)%rmask(i,Jstr-1)
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
        IF (LBC(inorth,isMtke,ng)%radiation) THEN
          DO k=0,N(ng)
            DO i=Istr,Iend+1
              grad(i,Jend  )=tke(i  ,Jend  ,k,nstp)-                    &
     &                       tke(i-1,Jend  ,k,nstp)
# ifdef MASKING
              grad(i,Jend  )=grad(i,Jend  )*                            &
     &                       GRID(ng)%umask(i,Jend  )
# endif
              grad(i,Jend+1)=tke(i  ,Jend+1,k,nstp)-                    &
     &                       tke(i-1,Jend+1,k,nstp)
# ifdef MASKING
              grad(i,Jend+1)=grad(i,Jend+1)*                            &
     &                       GRID(ng)%umask(i,Jend+1)
# endif
              gradL(i,Jend  )=gls(i  ,Jend  ,k,nstp)-                   &
     &                        gls(i-1,Jend  ,k,nstp)
# ifdef MASKING
              gradL(i,Jend  )=gradL(i,Jend  )*                          &
     &                        GRID(ng)%umask(i,Jend  )
# endif
              gradL(i,Jend+1)=gls(i  ,Jend+1,k,nstp)-                   &
     &                        gls(i-1,Jend+1,k,nstp)
# ifdef MASKING
              gradL(i,Jend+1)=gradL(i,Jend+1)*                          &
     &                        GRID(ng)%umask(i,Jend+1)
# endif
            END DO
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                dKdt=tke(i,Jend,k,nstp)-tke(i,Jend  ,k,nout)
                dKde=tke(i,Jend,k,nout)-tke(i,Jend-1,k,nout)
                IF ((dKdt*dKde).lt.0.0_r8) dKdt=0.0_r8
                IF ((dKdt*(grad(i  ,Jend)+                              &
     &                     grad(i+1,Jend))).gt.0.0_r8) THEN
                  dKdx=grad(i  ,Jend)
                ELSE
                  dKdx=grad(i+1,Jend)
                END IF
                cff=MAX(dKdx*dKdx+dKde*dKde,eps)
# ifdef RADIATION_2D
                Cx=MIN(cff,MAX(dKdt*dKdx,-cff))
# else
                Cx=0.0_r8
# endif
                Ce=dKdt*dKde
                tke(i,Jend+1,k,nout)=(cff*tke(i,Jend+1,k,nstp)+         &
     &                                Ce *tke(i,Jend  ,k,nout)-         &
     &                                MAX(Cx,0.0_r8)*                   &
     &                                   grad(i  ,Jend+1)-              &
     &                                MIN(Cx,0.0_r8)*                   &
     &                                   grad(i+1,Jend+1))/             &
     &                               (cff+Ce)
# ifdef MASKING
                tke(i,Jend+1,k,nout)=tke(i,Jend+1,k,nout)*              &
     &                               GRID(ng)%rmask(i,Jend+1)
# endif
                dKdt=gls(i,Jend,k,nstp)-gls(i,Jend  ,k,nout)
                dKde=gls(i,Jend,k,nout)-gls(i,Jend-1,k,nout)
                IF ((dKdt*dKde).lt.0.0_r8) dKdt=0.0_r8
                IF ((dKdt*(gradL(i  ,Jend)+                             &
     &                     gradL(i+1,Jend))).gt.0.0_r8) THEN
                  dKdx=gradL(i  ,Jend)
                ELSE
                  dKdx=gradL(i+1,Jend)
                END IF
                cff=MAX(dKdx*dKdx+dKde*dKde,eps)
# ifdef RADIATION_2D
                Cx=MIN(cff,MAX(dKdt*dKdx,-cff))
# else
                Cx=0.0_r8
# endif
                Ce=dKdt*dKde
                gls(i,Jend+1,k,nout)=(cff*gls(i,Jend+1,k,nstp)+         &
     &                                Ce *gls(i,Jend  ,k,nout)-         &
     &                                MAX(Cx,0.0_r8)*                   &
     &                                   gradL(i  ,Jend+1)-             &
     &                                MIN(Cx,0.0_r8)*                   &
     &                                   gradL(i+1,Jend+1))/            &
     &                               (cff+Ce)
# ifdef MASKING
                gls(i,Jend+1,k,nout)=gls(i,Jend+1,k,nout)*              &
     &                               GRID(ng)%rmask(i,Jend+1)
# endif
              END IF
            END DO
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (LBC(inorth,isMtke,ng)%gradient) THEN
          DO k=0,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                tke(i,Jend+1,k,nout)=tke(i,Jend,k,nout)
# ifdef MASKING
                tke(i,Jend+1,k,nout)=tke(i,Jend+1,k,nout)*              &
     &                               GRID(ng)%rmask(i,Jend+1)
# endif
                gls(i,Jend+1,k,nout)=gls(i,Jend,k,nout)
# ifdef MASKING
                gls(i,Jend+1,k,nout)=gls(i,Jend+1,k,nout)*              &
     &                               GRID(ng)%rmask(i,Jend+1)
# endif
              END IF
            END DO
          END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (LBC(inorth,isMtke,ng)%closed) THEN
          DO k=0,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                tke(i,Jend+1,k,nout)=tke(i,Jend,k,nout)
# ifdef MASKING
                tke(i,Jend+1,k,nout)=tke(i,Jend+1,k,nout)*              &
     &                               GRID(ng)%rmask(i,Jend+1)
# endif
                gls(i,Jend+1,k,nout)=gls(i,Jend,k,nout)
# ifdef MASKING
                gls(i,Jend+1,k,nout)=gls(i,Jend+1,k,nout)*              &
     &                               GRID(ng)%rmask(i,Jend+1)
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
            DO k=0,N(ng)
              tke(Istr-1,Jstr-1,k,nout)=0.5_r8*                         &
     &                                  (tke(Istr  ,Jstr-1,k,nout)+     &
     &                                   tke(Istr-1,Jstr  ,k,nout))
              gls(Istr-1,Jstr-1,k,nout)=0.5_r8*                         &
     &                                  (gls(Istr  ,Jstr-1,k,nout)+     &
     &                                   gls(Istr-1,Jstr  ,k,nout))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            DO k=0,N(ng)
              tke(Iend+1,Jstr-1,k,nout)=0.5_r8*                         &
     &                                  (tke(Iend  ,Jstr-1,k,nout)+     &
     &                                   tke(Iend+1,Jstr  ,k,nout))
               gls(Iend+1,Jstr-1,k,nout)=0.5_r8*                        &
     &                                   (gls(Iend  ,Jstr-1,k,nout)+    &
     &                                    gls(Iend+1,Jstr  ,k,nout))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO k=0,N(ng)
              tke(Istr-1,Jend+1,k,nout)=0.5_r8*                         &
     &                                  (tke(Istr  ,Jend+1,k,nout)+     &
     &                                   tke(Istr-1,Jend  ,k,nout))
              gls(Istr-1,Jend+1,k,nout)=0.5_r8*                         &
     &                                  (gls(Istr  ,Jend+1,k,nout)+     &
     &                                   gls(Istr-1,Jend  ,k,nout))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO k=0,N(ng)
              tke(Iend+1,Jend+1,k,nout)=0.5_r8*                         &
     &                                  (tke(Iend  ,Jend+1,k,nout)+     &
     &                                   tke(Iend+1,Jend  ,k,nout))
              gls(Iend+1,Jend+1,k,nout)=0.5_r8*                         &
     &                                  (gls(Iend  ,Jend+1,k,nout)+     &
     &                                   gls(Iend+1,Jend  ,k,nout))
            END DO
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE tkebc_tile
#endif
      END MODULE tkebc_mod
