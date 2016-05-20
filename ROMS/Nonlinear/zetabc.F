#include "cppdefs.h"
      MODULE zetabc_mod
!
!svn $Id: zetabc.F 771 2008-09-30 02:30:00Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets lateral boundary conditions for free-surface.     !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: zetabc_tile
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE zetabc (ng, tile, kout)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, kout
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL zetabc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  krhs(ng), kstp(ng), kout,                       &
     &                  OCEAN(ng) % zeta)

      RETURN
      END SUBROUTINE zetabc
!
!***********************************************************************
      SUBROUTINE zetabc_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        krhs, kstp, kout,                         &
     &                        zeta)
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
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: krhs, kstp, kout
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: zeta(LBi:,LBj:,:)
#else
      real(r8), intent(inout) :: zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: i, j, know

      real(r8), parameter :: eps =1.0E-20_r8

      real(r8) :: Ce, Cx
      real(r8) :: cff, cff1, cff2, dt2d, dZde, dZdt, dZdx, tau

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad

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
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!
!  Western edge, implicit upstream radiation condition.
!
        IF (LBC(iwest,isFsur,ng)%radiation) THEN
          DO j=Jstr,Jend+1
            grad(Istr-1,j)=zeta(Istr-1,j  ,know)-                       &
     &                     zeta(Istr-1,j-1,know)
#ifdef MASKING
            grad(Istr-1,j)=grad(Istr-1,j)*GRID(ng)%vmask(Istr-1,j)
#endif
            grad(Istr,j)=zeta(Istr,j  ,know)-                           &
     &                   zeta(Istr,j-1,know)
#ifdef MASKING
            grad(Istr,j)=grad(Istr,j)*GRID(ng)%vmask(Istr,j)
#endif
          END DO
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              dZdt=zeta(Istr,j,know)-zeta(Istr  ,j,kout)
              dZdx=zeta(Istr,j,kout)-zeta(Istr+1,j,kout)

              IF (LBC(iwest,isFsur,ng)%nudging) THEN
                IF ((dZdt*dZdx).lt.0.0_r8) THEN
                  tau=FSobc_in(ng,iwest)
                ELSE
                  tau=FSobc_out(ng,iwest)
                END IF
                tau=tau*dt2d
              END IF

              IF ((dZdt*dZdx).lt.0.0_r8) dZdt=0.0_r8
              IF ((dZdt*(grad(Istr,j)+grad(Istr,j+1))).gt.0.0_r8) THEN
                dZde=grad(Istr,j  )
              ELSE
                dZde=grad(Istr,j+1)
              END IF
              cff=MAX(dZdx*dZdx+dZde*dZde,eps)
              Cx=dZdt*dZdx
#ifdef RADIATION_2D
              Ce=MIN(cff,MAX(dZdt*dZde,-cff))
#else
              Ce=0.0_r8
#endif
#if defined CELERITY_WRITE && defined FORWARD_WRITE
              BOUNDARY(ng)%zeta_west_Cx(j)=Cx
              BOUNDARY(ng)%zeta_west_Ce(j)=Ce
              BOUNDARY(ng)%zeta_west_C2(j)=cff
#endif
              zeta(Istr-1,j,kout)=(cff*zeta(Istr-1,j,know)+             &
     &                             Cx *zeta(Istr  ,j,kout)-             &
     &                             MAX(Ce,0.0_r8)*grad(Istr-1,j  )-     &
     &                             MIN(Ce,0.0_r8)*grad(Istr-1,j+1))/    &
     &                            (cff+Cx)

              IF (LBC(iwest,isFsur,ng)%nudging) THEN
                zeta(Istr-1,j,kout)=zeta(Istr-1,j,kout)+                &
     &                              tau*(BOUNDARY(ng)%zeta_west(j)-     &
     &                                   zeta(Istr-1,j,know))
              END IF
#ifdef MASKING
              zeta(Istr-1,j,kout)=zeta(Istr-1,j,kout)*                  &
     &                            GRID(ng)%rmask(Istr-1,j)
#endif
            END IF
          END DO
!
!  Western edge, explicit Chapman boundary condition.
!
        ELSE IF (LBC(iwest,isFsur,ng)%Chapman_explicit) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              cff=dt2d*GRID(ng)%pm(Istr,j)
              cff1=SQRT(g*(GRID(ng)%h(Istr,j)+                          &
     &                     zeta(Istr,j,know)))
              Cx=cff*cff1
              zeta(Istr-1,j,kout)=(1.0_r8-Cx)*zeta(Istr-1,j,know)+      &
     &                            Cx*zeta(Istr,j,know)
#ifdef MASKING
              zeta(Istr-1,j,kout)=zeta(Istr-1,j,kout)*                  &
     &                            GRID(ng)%rmask(Istr-1,j)
#endif
            END IF
          END DO
!
!  Western edge, implicit Chapman boundary condition.
!
        ELSE IF (LBC(iwest,isFsur,ng)%Chapman_implicit) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              cff=dt2d*GRID(ng)%pm(Istr,j)
#  ifdef ICESHELF
              cff1=SQRT(g*(GRID(ng)%h(Istr,j)-                          &
     &                     ABS(GRID(ng)%zice(Istr,j))+                  &
     &                     zeta(Istr,j,know)))
#  else
              cff1=SQRT(g*(GRID(ng)%h(Istr,j)+                          &
     &                     zeta(Istr,j,know)))
#  endif
              Cx=cff*cff1
              cff2=1.0_r8/(1.0_r8+Cx)
              zeta(Istr-1,j,kout)=cff2*(zeta(Istr-1,j,know)+            &
     &                                  Cx*zeta(Istr,j,kout))
#ifdef MASKING
              zeta(Istr-1,j,kout)=zeta(Istr-1,j,kout)*                  &
     &                            GRID(ng)%rmask(Istr-1,j)
#endif
            END IF
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (LBC(iwest,isFsur,ng)%clamped) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              zeta(Istr-1,j,kout)=BOUNDARY(ng)%zeta_west(j)
#ifdef MASKING
              zeta(Istr-1,j,kout)=zeta(Istr-1,j,kout)*                  &
     &                            GRID(ng)%rmask(Istr-1,j)
#endif
            END IF
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (LBC(iwest,isFsur,ng)%gradient) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              zeta(Istr-1,j,kout)=zeta(Istr,j,kout)
#ifdef MASKING
              zeta(Istr-1,j,kout)=zeta(Istr-1,j,kout)*                  &
     &                            GRID(ng)%rmask(Istr-1,j)
#endif
            END IF
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (LBC(iwest,isFsur,ng)%closed) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              zeta(Istr-1,j,kout)=zeta(Istr,j,kout)
#ifdef MASKING
              zeta(Istr-1,j,kout)=zeta(Istr-1,j,kout)*                  &
     &                            GRID(ng)%rmask(Istr-1,j)
#endif
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
        IF (LBC(ieast,isFsur,ng)%radiation) THEN
          DO j=Jstr,Jend+1
            grad(Iend  ,j)=zeta(Iend  ,j  ,know)-                       &
     &                     zeta(Iend  ,j-1,know)
#ifdef MASKING
            grad(Iend  ,j)=grad(Iend  ,j)*GRID(ng)%vmask(Iend  ,j)
#endif
            grad(Iend+1,j)=zeta(Iend+1,j  ,know)-                       &
     &                     zeta(Iend+1,j-1,know)
#ifdef MASKING
            grad(Iend+1,j)=grad(Iend+1,j)*GRID(ng)%vmask(Iend+1,j)
#endif
          END DO
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              dZdt=zeta(Iend,j,know)-zeta(Iend  ,j,kout)
              dZdx=zeta(Iend,j,kout)-zeta(Iend-1,j,kout)

              IF (LBC(ieast,isFsur,ng)%nudging) THEN
                IF ((dZdt*dZdx).lt.0.0_r8) THEN
                  tau=FSobc_in(ng,ieast)
                ELSE
                  tau=FSobc_out(ng,ieast)
                END IF
                tau=tau*dt2d
              END IF

              IF ((dZdt*dZdx).lt.0.0_r8) dZdt=0.0_r8
              IF ((dZdt*(grad(Iend,j)+grad(Iend,j+1))).gt.0.0_r8) THEN
                dZde=grad(Iend,j  )
              ELSE
                dZde=grad(Iend,j+1)
              END IF
              cff=MAX(dZdx*dZdx+dZde*dZde,eps)
              Cx=dZdt*dZdx
#ifdef RADIATION_2D
              Ce=MIN(cff,MAX(dZdt*dZde,-cff))
#else
              Ce=0.0_r8
#endif
#if defined CELERITY_WRITE && defined FORWARD_WRITE
              BOUNDARY(ng)%zeta_east_Cx(j)=Cx
              BOUNDARY(ng)%zeta_east_Ce(j)=Ce
              BOUNDARY(ng)%zeta_east_C2(j)=cff
#endif
              zeta(Iend+1,j,kout)=(cff*zeta(Iend+1,j,know)+             &
     &                             Cx *zeta(Iend  ,j,kout)-             &
     &                             MAX(Ce,0.0_r8)*grad(Iend+1,j  )-     &
     &                             MIN(Ce,0.0_r8)*grad(Iend+1,j+1))/    &
     &                            (cff+Cx)

              IF (LBC(ieast,isFsur,ng)%nudging) THEN
                zeta(Iend+1,j,kout)=zeta(Iend+1,j,kout)+                &
     &                              tau*(BOUNDARY(ng)%zeta_east(j)-     &
     &                                   zeta(Iend+1,j,know))
              END IF
#ifdef MASKING
              zeta(Iend+1,j,kout)=zeta(Iend+1,j,kout)*                  &
     &                            GRID(ng)%rmask(Iend+1,j)
#endif
            END IF
          END DO
!
!  Eastern edge, explicit Chapman boundary condition.
!
        ELSE IF (LBC(ieast,isFsur,ng)%Chapman_explicit) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              cff=dt2d*GRID(ng)%pm(Iend,j)
              cff1=SQRT(g*(GRID(ng)%h(Iend,j)+                          &
     &                     zeta(Iend,j,know)))
              Cx=cff*cff1
              zeta(Iend+1,j,kout)=(1.0_r8-Cx)*zeta(Iend+1,j,know)+      &
     &                            Cx*zeta(Iend,j,know)
#ifdef MASKING
              zeta(Iend+1,j,kout)=zeta(Iend+1,j,kout)*                  &
     &                            GRID(ng)%rmask(Iend+1,j)
#endif
            END IF
          END DO
!
!  Eastern edge, implicit Chapman boundary condition.
!
        ELSE IF (LBC(ieast,isFsur,ng)%Chapman_implicit) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              cff=dt2d*GRID(ng)%pm(Iend,j)
#  ifdef ICESHELF
              cff1=SQRT(g*(GRID(ng)%h(Iend,j)-                          &
     &                     ABS(GRID(ng)%zice(Iend,j))+                  &
     &                     zeta(Iend,j,know)))
#  else
              cff1=SQRT(g*(GRID(ng)%h(Iend,j)+                          &
     &                     zeta(Iend,j,know)))
#  endif
              Cx=cff*cff1
              cff2=1.0_r8/(1.0_r8+Cx)
              zeta(Iend+1,j,kout)=cff2*(zeta(Iend+1,j,know)+            &
     &                                  Cx*zeta(Iend,j,kout))
#ifdef MASKING
              zeta(Iend+1,j,kout)=zeta(Iend+1,j,kout)*                  &
     &                            GRID(ng)%rmask(Iend+1,j)
#endif
            END IF
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (LBC(ieast,isFsur,ng)%clamped) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              zeta(Iend+1,j,kout)=BOUNDARY(ng)%zeta_east(j)
#ifdef MASKING
              zeta(Iend+1,j,kout)=zeta(Iend+1,j,kout)*                  &
     &                            GRID(ng)%rmask(Iend+1,j)
#endif
            END IF
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (LBC(ieast,isFsur,ng)%gradient) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              zeta(Iend+1,j,kout)=zeta(Iend,j,kout)
#ifdef MASKING
              zeta(Iend+1,j,kout)=zeta(Iend+1,j,kout)*                  &
     &                            GRID(ng)%rmask(Iend+1,j)
#endif
            END IF
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (LBC(ieast,isFsur,ng)%closed) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              zeta(Iend+1,j,kout)=zeta(Iend,j,kout)
#ifdef MASKING
              zeta(Iend+1,j,kout)=zeta(Iend+1,j,kout)*                  &
     &                            GRID(ng)%rmask(Iend+1,j)
#endif
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
        IF (LBC(isouth,isFsur,ng)%radiation) THEN
          DO i=Istr,Iend+1
            grad(i,Jstr  )=zeta(i  ,Jstr,know)-                         &
     &                     zeta(i-1,Jstr,know)
#ifdef MASKING
            grad(i,Jstr  )=grad(i,Jstr  )*GRID(ng)%umask(i,Jstr  )
#endif
            grad(i,Jstr-1)=zeta(i  ,Jstr-1,know)-                       &
     &                     zeta(i-1,Jstr-1,know)
#ifdef MASKING
            grad(i,Jstr-1)=grad(i,Jstr-1)*GRID(ng)%umask(i,Jstr-1)
#endif
          END DO
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              dZdt=zeta(i,Jstr,know)-zeta(i,Jstr  ,kout)
              dZde=zeta(i,Jstr,kout)-zeta(i,Jstr-1,kout)

              IF (LBC(isouth,isFsur,ng)%nudging) THEN
                IF ((dZdt*dZde).lt.0.0_r8) THEN
                  tau=FSobc_in(ng,isouth)
                ELSE
                  tau=FSobc_out(ng,isouth)
                END IF
                tau=tau*dt2d
              END IF

              IF ((dZdt*dZde).lt.0.0_r8) dZdt=0.0_r8
              IF ((dZdt*(grad(i,Jstr)+grad(i+1,Jstr))).gt.0.0_r8) THEN
                dZdx=grad(i  ,Jstr)
              ELSE
                dZdx=grad(i+1,Jstr)
              END IF
              cff=MAX(dZdx*dZdx+dZde*dZde,eps)
#ifdef RADIATION_2D
              Cx=MIN(cff,MAX(dZdt*dZdx,-cff))
#else
              Cx=0.0_r8
#endif
              Ce=dZdt*dZde
#if defined CELERITY_WRITE && defined FORWARD_WRITE
              BOUNDARY(ng)%zeta_south_Cx(i)=Cx
              BOUNDARY(ng)%zeta_south_Ce(i)=Ce
              BOUNDARY(ng)%zeta_south_C2(i)=cff
#endif
              zeta(i,Jstr-1,kout)=(cff*zeta(i,Jstr-1,know)+             &
     &                             Ce *zeta(i,Jstr  ,kout)-             &
     &                             MAX(Cx,0.0_r8)*grad(i  ,Jstr)-       &
     &                             MIN(Cx,0.0_r8)*grad(i+1,Jstr))/      &
     &                            (cff+Ce)

              IF (LBC(isouth,isFsur,ng)%nudging) THEN
                zeta(i,Jstr-1,kout)=zeta(i,Jstr-1,kout)+                &
     &                              tau*(BOUNDARY(ng)%zeta_south(i)-    &
     &                                   zeta(i,Jstr-1,know))
              END IF
#ifdef MASKING
              zeta(i,Jstr-1,kout)=zeta(i,Jstr-1,kout)*                  &
     &                            GRID(ng)%rmask(i,Jstr-1)
#endif
            END IF
          END DO
!
!  Southern edge, explicit Chapman boundary condition.
!
        ELSE IF (LBC(isouth,isFsur,ng)%Chapman_explicit) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              cff=dt2d*GRID(ng)%pn(i,Jstr)
              cff1=SQRT(g*(GRID(ng)%h(i,Jstr)+                          &
     &                     zeta(i,Jstr,know)))
              Ce=cff*cff1
              zeta(i,Jstr-1,kout)=(1.0_r8-Ce)*zeta(i,Jstr-1,know)+      &
     &                            Ce*zeta(i,Jstr,know)
#ifdef MASKING
              zeta(i,Jstr-1,kout)=zeta(i,Jstr-1,kout)*                  &
     &                            GRID(ng)%rmask(i,Jstr-1)
#endif
            END IF
          END DO
!
!  Southern edge, implicit Chapman boundary condition.
!
        ELSE IF (LBC(isouth,isFsur,ng)%Chapman_implicit) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              cff=dt2d*GRID(ng)%pn(i,Jstr)
#  ifdef ICESHELF
              cff1=SQRT(g*(GRID(ng)%h(i,Jstr)-                          &
     &                     ABS(GRID(ng)%zice(i,Jstr))+                  &
     &                     zeta(i,Jstr,know)))
#  else
              cff1=SQRT(g*(GRID(ng)%h(i,Jstr)+                          &
     &                     zeta(i,Jstr,know)))
#  endif
              Ce=cff*cff1
              cff2=1.0_r8/(1.0_r8+Ce)
              zeta(i,Jstr-1,kout)=cff2*(zeta(i,Jstr-1,know)+            &
     &                                  Ce*zeta(i,Jstr,kout))
#ifdef MASKING
              zeta(i,Jstr-1,kout)=zeta(i,Jstr-1,kout)*                  &
     &                            GRID(ng)%rmask(i,Jstr-1)
#endif
            END IF
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (LBC(isouth,isFsur,ng)%clamped) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              zeta(i,Jstr-1,kout)=BOUNDARY(ng)%zeta_south(i)
#ifdef MASKING
              zeta(i,Jstr-1,kout)=zeta(i,Jstr-1,kout)*                  &
     &                            GRID(ng)%rmask(i,Jstr-1)
#endif
            END IF
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (LBC(isouth,isFsur,ng)%gradient) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              zeta(i,Jstr-1,kout)=zeta(i,Jstr,kout)
#ifdef MASKING
              zeta(i,Jstr-1,kout)=zeta(i,Jstr-1,kout)*                  &
     &                            GRID(ng)%rmask(i,Jstr-1)
#endif
            END IF
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE IF (LBC(isouth,isFsur,ng)%closed) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              zeta(i,Jstr-1,kout)=zeta(i,Jstr,kout)
#ifdef MASKING
              zeta(i,Jstr-1,kout)=zeta(i,Jstr-1,kout)*                  &
     &                            GRID(ng)%rmask(i,Jstr-1)
#endif
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
        IF (LBC(inorth,isFsur,ng)%radiation) THEN
          DO i=Istr,Iend+1
            grad(i,Jend  )=zeta(i  ,Jend  ,know)-                       &
     &                     zeta(i-1,Jend  ,know)
#ifdef MASKING
            grad(i,Jend  )=grad(i,Jend  )*GRID(ng)%umask(i,Jend  )
#endif
            grad(i,Jend+1)=zeta(i  ,Jend+1,know)-                       &
     &                     zeta(i-1,Jend+1,know)
#ifdef MASKING
            grad(i,Jend+1)=grad(i,Jend+1)*GRID(ng)%umask(i,Jend+1)
#endif
          END DO
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              dZdt=zeta(i,Jend,know)-zeta(i,Jend  ,kout)
              dZde=zeta(i,Jend,kout)-zeta(i,Jend-1,kout)

              IF (LBC(inorth,isFsur,ng)%nudging) THEN
                IF ((dZdt*dZde).lt.0.0_r8) THEN
                  tau=FSobc_in(ng,inorth)
                ELSE
                  tau=FSobc_out(ng,inorth)
                END IF
                tau=tau*dt2d
              END IF

              IF ((dZdt*dZde).lt.0.0_r8) dZdt=0.0_r8
              IF ((dZdt*(grad(i,Jend)+grad(i+1,Jend))).gt.0.0_r8) THEN
                dZdx=grad(i  ,Jend)
              ELSE
                dZdx=grad(i+1,Jend)
              END IF
              cff=MAX(dZdx*dZdx+dZde*dZde,eps)
#ifdef RADIATION_2D
              Cx=MIN(cff,MAX(dZdt*dZdx,-cff))
#else
              Cx=0.0_r8
#endif
              Ce=dZdt*dZde
#if defined CELERITY_WRITE && defined FORWARD_WRITE
              BOUNDARY(ng)%zeta_north_Cx(i)=Cx
              BOUNDARY(ng)%zeta_north_Ce(i)=Ce
              BOUNDARY(ng)%zeta_north_C2(i)=cff
#endif
              zeta(i,Jend+1,kout)=(cff*zeta(i,Jend+1,know)+             &
     &                             Ce *zeta(i,Jend  ,kout)-             &
     &                             MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-     &
     &                             MIN(Cx,0.0_r8)*grad(i+1,Jend+1))/    &
     &                            (cff+Ce)

              IF (LBC(inorth,isFsur,ng)%nudging) THEN
                zeta(i,Jend+1,kout)=zeta(i,Jend+1,kout)+                &
     &                              tau*(BOUNDARY(ng)%zeta_north(i)-    &
     &                                   zeta(i,Jend+1,know))
              END IF
#ifdef MASKING
              zeta(i,Jend+1,kout)=zeta(i,Jend+1,kout)*                  &
     &                            GRID(ng)%rmask(i,Jend+1)
#endif
            END IF
          END DO
!
!  Northern edge, explicit Chapman boundary condition.
!
        ELSE IF (LBC(inorth,isFsur,ng)%Chapman_explicit) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              cff=dt2d*GRID(ng)%pn(i,Jend)
              cff1=SQRT(g*(GRID(ng)%h(i,Jend)+                          &
     &                     zeta(i,Jend,know)))
              Ce=cff*cff1
              zeta(i,Jend+1,kout)=(1.0_r8-Ce)*zeta(i,Jend+1,know)+      &
     &                            Ce*zeta(i,Jend,know)
#ifdef MASKING
              zeta(i,Jend+1,kout)=zeta(i,Jend+1,kout)*                  &
     &                            GRID(ng)%rmask(i,Jend+1)
#endif
            END IF
          END DO
!
!  Northern edge, implicit Chapman boundary condition.
!
        ELSE IF (LBC(inorth,isFsur,ng)%Chapman_implicit) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              cff=dt2d*GRID(ng)%pn(i,Jend)
#  ifdef ICESHELF
              cff1=SQRT(g*(GRID(ng)%h(i,Jend)-                          &
     &                     ABS(GRID(ng)%zice(i,Jend))+                  &
     &                     zeta(i,Jend,know)))
#  else
              cff1=SQRT(g*(GRID(ng)%h(i,Jend)+                          &
     &                     zeta(i,Jend,know)))
#  endif
              Ce=cff*cff1
              cff2=1.0_r8/(1.0_r8+Ce)
              zeta(i,Jend+1,kout)=cff2*(zeta(i,Jend+1,know)+            &
     &                                  Ce*zeta(i,Jend,kout))
#ifdef MASKING
              zeta(i,Jend+1,kout)=zeta(i,Jend+1,kout)*                  &
     &                            GRID(ng)%rmask(i,Jend+1)
#endif
            END IF
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (LBC(inorth,isFsur,ng)%clamped) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              zeta(i,Jend+1,kout)=BOUNDARY(ng)%zeta_north(i)
#ifdef MASKING
              zeta(i,Jend+1,kout)=zeta(i,Jend+1,kout)*                  &
     &                            GRID(ng)%rmask(i,Jend+1)
#endif
            END IF
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (LBC(inorth,isFsur,ng)%gradient) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              zeta(i,Jend+1,kout)=zeta(i,Jend,kout)
#ifdef MASKING
              zeta(i,Jend+1,kout)=zeta(i,Jend+1,kout)*                  &
     &                            GRID(ng)%rmask(i,Jend+1)
#endif
            END IF
          END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (LBC(inorth,isFsur,ng)%closed) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              zeta(i,Jend+1,kout)=zeta(i,Jend,kout)
#ifdef MASKING
              zeta(i,Jend+1,kout)=zeta(i,Jend+1,kout)*                  &
     &                            GRID(ng)%rmask(i,Jend+1)
#endif
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
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
            zeta(Istr-1,Jstr-1,kout)=0.5_r8*(zeta(Istr  ,Jstr-1,kout)+  &
     &                                       zeta(Istr-1,Jstr  ,kout))
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            zeta(Iend+1,Jstr-1,kout)=0.5_r8*(zeta(Iend  ,Jstr-1,kout)+  &
     &                                       zeta(Iend+1,Jstr  ,kout))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            zeta(Istr-1,Jend+1,kout)=0.5_r8*(zeta(Istr-1,Jend  ,kout)+  &
     &                                       zeta(Istr  ,Jend+1,kout))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            zeta(Iend+1,Jend+1,kout)=0.5_r8*(zeta(Iend+1,Jend  ,kout)+  &
     &                                       zeta(Iend  ,Jend+1,kout))
          END IF
        END IF
      END IF

#if defined WET_DRY
!
!-----------------------------------------------------------------------
! Ensure that water level on boundary cells is above bed elevation.
!-----------------------------------------------------------------------
!
      cff=Dcrit(ng)-eps
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              IF (zeta(Istr-1,j,kout).le.                               &
     &            (Dcrit(ng)-GRID(ng)%h(Istr-1,j))) THEN
                zeta(Istr-1,j,kout)=cff-GRID(ng)%h(Istr-1,j)
              END IF
            END IF
          END DO
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              IF (zeta(Iend+1,j,kout).le.                               &
     &            (Dcrit(ng)-GRID(ng)%h(Iend+1,j))) THEN
                zeta(Iend+1,j,kout)=cff-GRID(ng)%h(Iend+1,j)
              END IF
            END IF
          END DO
        END IF
      END IF

      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              IF (zeta(i,Jstr-1,kout).le.                               &
     &            (Dcrit(ng)-GRID(ng)%h(i,Jstr-1))) THEN
                zeta(i,Jstr-1,kout)=cff-GRID(ng)%h(i,Jstr-1)
              END IF
            END IF
          END DO
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              IF (zeta(i,Jend+1,kout).le.                               &
     &            (Dcrit(ng)-GRID(ng)%h(i,Jend+1))) THEN
                zeta(i,Jend+1,kout)=cff-GRID(ng)%h(i,Jend+1)
              END IF
            END IF
          END DO
        END IF
      END IF

      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
            IF (zeta(Istr-1,Jstr-1,kout).le.                            &
     &          (Dcrit(ng)-GRID(ng)%h(Istr-1,Jstr-1))) THEN
              zeta(Istr-1,Jstr-1,kout)=cff-GRID(ng)%h(Istr-1,Jstr-1)
            END IF
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            IF (zeta(Iend+1,Jstr-1,kout).le.                            &
     &          (Dcrit(ng)-GRID(ng)%h(Iend+1,Jstr-1))) THEN
              zeta(Iend+1,Jstr-1,kout)=cff-GRID(ng)%h(Iend+1,Jstr-1)
            END IF
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            IF (zeta(Istr-1,Jend+1,kout).le.                            &
     &          (Dcrit(ng)-GRID(ng)%h(Istr-1,Jend+1))) THEN
              zeta(Istr-1,Jend+1,kout)=cff-GRID(ng)%h(Istr-1,Jend+1)
            END IF
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            IF (zeta(Iend+1,Jend+1,kout).le.                            &
     &          (Dcrit(ng)-GRID(ng)%h(Iend+1,Jend+1))) THEN
              zeta(Iend+1,Jend+1,kout)=cff-GRID(ng)%h(Iend+1,Jend+1)
            END IF
          END IF
        END IF
      END IF
#endif

      RETURN

      END SUBROUTINE zetabc_tile
      END MODULE zetabc_mod
