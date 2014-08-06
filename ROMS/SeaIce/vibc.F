#include "cppdefs.h"
      MODULE vibc_mod
#ifdef ICE_MODEL
!***********************************************************************
!  Compute the lateral boundary conditions on the ice V-velocity.
!***********************************************************************

      implicit none

      PRIVATE
      PUBLIC vibc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE vibc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ice
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

#include "tile.h"
!
      CALL vibc_tile (ng, tile,                                         &
     &                LBi, UBi, LBj, UBj,                               &
     &                IminS, ImaxS, JminS, JmaxS,                       &
     &                liuol(ng), liunw(ng),                             &
     &                ICE(ng) % vi)
      RETURN
      END SUBROUTINE vibc
!
!***********************************************************************
      SUBROUTINE vibc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liuol, liunw,                               &
     &                      vi)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_boundary
      USE mod_grid
      USE mod_scalars

      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: liuol, liunw

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: vi(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: vi(LBi:UBi,LBj:UBj,2)
# endif
!
!  Local variable declarations.
!
      integer :: i, Jmax, Jmin, j, know
      real(r8), parameter :: eps =1.0E-20_r8
      real(r8) :: Ce, Cx
      real(r8):: cff, dVde, dVdt, dVdx, tau

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set time-indices
!-----------------------------------------------------------------------
!
        know=liuol
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        IF (LBC(iwest,isVice,ng)%radiation) THEN
!
!  Western edge, implicit upstream radiation condition.
!
          DO j=JstrP-1,Jend
            grad(Istr-1,j)=vi(Istr-1,j+1,know)-                         &
     &                     vi(Istr-1,j  ,know)
            grad(Istr  ,j)=vi(Istr  ,j+1,know)-                         &
     &                     vi(Istr  ,j  ,know)
          END DO
          DO j=JstrP,Jend
            dVdt=vi(Istr,j,know)-vi(Istr  ,j,liunw)
            dVdx=vi(Istr,j,liunw)-vi(Istr+1,j,liunw)
            IF (LBC(iwest,isVice,ng)%nudging) THEN
              IF ((dVdt*dVdx).lt.0.0_r8) THEN
                tau=M2obc_in(ng,iwest)
              ELSE
                tau=M2obc_out(ng,iwest)
              END IF
              tau=tau*dt(ng)
            END IF
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
            vi(Istr-1,j,liunw)=(cff*vi(Istr-1,j,know)+                  &
     &                         Cx *vi(Istr  ,j,liunw)-                  &
     &                         MAX(Ce,0.0_r8)*grad(Istr-1,j-1)-         &
     &                         MIN(Ce,0.0_r8)*grad(Istr-1,j  ))/        &
     &                        (cff+Cx)
            IF (LBC(iwest,isVice,ng)%nudging) THEN
              vi(Istr-1,j,liunw)=vi(Istr-1,j,liunw)+                    &
     &                        tau*(BOUNDARY(ng)%vi_west(j)-             &
     &                             vi(Istr-1,j,know))
            END IF
# ifdef MASKING
            vi(Istr-1,j,liunw)=vi(Istr-1,j,liunw)*                      &
     &                        GRID(ng)%vmask(Istr-1,j)
# endif
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (LBC(iwest,isVice,ng)%clamped) THEN
          DO j=JstrP,Jend
            vi(0,j,liunw)=BOUNDARY(ng)%vi_west(j)
# ifdef MASKING
            vi(0,j,liunw)=vi(0,j,liunw)*                                &
     &                   GRID(ng)%vmask(0,j)
# endif
# ifdef WET_DRY
            vi(0,j,liunw)=vi(0,j,liunw)*                                &
     &                   GRID(ng)%vmask_wet(0,j)
# endif
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (LBC(iwest,isVice,ng)%gradient) THEN
          DO j=JstrP,Jend
            vi(0,j,liunw)=vi(1,j,liunw)
# ifdef MASKING
            vi(0,j,liunw)=vi(0,j,liunw)*                                &
     &                   GRID(ng)%vmask(0,j)
# endif
# ifdef WET_DRY
            vi(0,j,liunw)=vi(0,j,liunw)*                                &
     &                   GRID(ng)%vmask_wet(0,j)
# endif
          END DO
!
!  Western edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
        ELSE IF (LBC(iwest,isVice,ng)%closed) THEN
          IF (NSperiodic(ng)) THEN
            Jmin=JstrP
            Jmax=Jend
          ELSE
            Jmin=Jstr
            Jmax=JendT
          END IF
          DO j=Jmin,Jmax
            vi(0,j,liunw)=gamma2(ng)*vi(0,j,liunw)
# ifdef MASKING
            vi(0,j,liunw)=vi(0,j,liunw)*                                &
     &                   GRID(ng)%vmask(0,j)
# endif
# ifdef WET_DRY
            vi(0,j,liunw)=vi(0,j,liunw)*                                &
     &                   GRID(ng)%vmask_wet(0,j)
# endif
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        IF (LBC(ieast,isVice,ng)%radiation) THEN
!
!  Eastern edge, implicit upstream radiation condition.
!
          DO j=JstrP-1,Jend
            grad(Iend  ,j)=vi(Iend  ,j+1,know)-                         &
     &                     vi(Iend  ,j  ,know)
            grad(Iend+1,j)=vi(Iend+1,j+1,know)-                         &
     &                   vi(Iend+1,j  ,know)
          END DO
          DO j=JstrP,Jend
            dVdt=vi(Iend,j,know)-vi(Iend  ,j,liunw)
            dVdx=vi(Iend,j,liunw)-vi(Iend-1,j,liunw)
            IF (LBC(ieast,isVice,ng)%nudging) THEN
              IF ((dVdt*dVdx).lt.0.0_r8) THEN
                tau=M2obc_in(ng,ieast)
              ELSE
                tau=M2obc_out(ng,ieast)
              END IF
              tau=tau*dt(ng)
            END IF
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
            vi(Iend+1,j,liunw)=(cff*vi(Iend+1,j,know)+                  &
     &                         Cx *vi(Iend  ,j,liunw)-                  &
     &                         MAX(Ce,0.0_r8)*grad(Iend+1,j-1)-         &
     &                         MIN(Ce,0.0_r8)*grad(Iend+1,j  ))/        &
     &                        (cff+Cx)
            IF (LBC(ieast,isVice,ng)%nudging) THEN
              vi(Iend+1,j,liunw)=vi(Iend+1,j,liunw)+                    &
     &                        tau*(BOUNDARY(ng)%vi_east(j)-             &
     &                             vi(Iend+1,j,know))
            END IF
# ifdef MASKING
            vi(Iend+1,j,liunw)=vi(Iend+1,j,liunw)*                      &
     &                        GRID(ng)%vmask(Iend+1,j)
# endif
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (LBC(ieast,isVice,ng)%clamped) THEN
          DO j=JstrP,Jend
            vi(Lm(ng)+1,j,liunw)=BOUNDARY(ng)%vi_east(j)
# ifdef MASKING
            vi(Lm(ng)+1,j,liunw)=vi(Lm(ng)+1,j,liunw)*                  &
     &                          GRID(ng)%vmask(Lm(ng)+1,j)
# endif
# ifdef WET_DRY
            vi(Lm(ng)+1,j,liunw)=vi(Lm(ng)+1,j,liunw)*                  &
     &                          GRID(ng)%vmask_wet(Lm(ng)+1,j)
# endif
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (LBC(ieast,isVice,ng)%gradient) THEN
          DO j=JstrP,Jend
            vi(Lm(ng)+1,j,liunw)=vi(Lm(ng),j,liunw)
# ifdef MASKING
            vi(Lm(ng)+1,j,liunw)=vi(Lm(ng)+1,j,liunw)*                  &
     &                       GRID(ng)%vmask(Lm(ng)+1,j)
# endif
# ifdef WET_DRY
            vi(Lm(ng)+1,j,liunw)=vi(Lm(ng)+1,j,liunw)*                  &
     &                       GRID(ng)%vmask_wet(Lm(ng)+1,j)
# endif
          END DO
!
!  Eastern edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
        ELSE IF (LBC(ieast,isVice,ng)%closed) THEN
          IF (NSperiodic(ng)) THEN
            Jmin=JstrP
            Jmax=Jend
          ELSE
            Jmin=Jstr
            Jmax=JendT
          END IF
          DO j=Jmin,Jmax
            vi(Lm(ng)+1,j,liunw)=gamma2(ng)*vi(Lm(ng),j,liunw)
# ifdef MASKING
            vi(Lm(ng)+1,j,liunw)=vi(Lm(ng)+1,j,liunw)*                  &
     &                          GRID(ng)%vmask(Lm(ng)+1,j)
# endif
# ifdef WET_DRY
            vi(Lm(ng)+1,j,liunw)=vi(Lm(ng)+1,j,liunw)*                  &
     &                          GRID(ng)%vmask_wet(Lm(ng)+1,j)
# endif
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        IF (LBC(isouth,isVice,ng)%radiation) THEN
!
!  Southern edge, implicit upstream radiation condition.
!
          DO i=Istr,Iend+1
            grad(i,Jstr  )=vi(i  ,Jstr  ,know)-                         &
     &                     vi(i-1,Jstr  ,know)
            grad(i,Jstr+1)=vi(i  ,Jstr+1,know)-                         &
     &                     vi(i-1,Jstr+1,know)
          END DO
          DO i=Istr,Iend
            dVdt=vi(i,Jstr+1,know)-vi(i,Jstr+1,liunw)
            dVde=vi(i,Jstr+1,liunw)-vi(i,Jstr+2,liunw)
            IF (LBC(isouth,isVice,ng)%nudging) THEN
              IF ((dVdt*dVde).lt.0.0_r8) THEN
                tau=M2obc_in(ng,isouth)
              ELSE
                tau=M2obc_out(ng,isouth)
              END IF
              tau=tau*dt(ng)
            END IF
            IF ((dVdt*dVde).lt.0.0_r8) dVdt=0.0_r8
            IF ((dVdt*(grad(i,Jstr+1)+grad(i+1,Jstr+1))).gt.0.0_r8) THEN
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
            vi(i,Jstr,liunw)=(cff*vi(i,Jstr  ,know)+                    &
     &                       Ce *vi(i,Jstr+1,liunw)-                    &
     &                       MAX(Cx,0.0_r8)*grad(i  ,Jstr)-             &
     &                       MIN(Cx,0.0_r8)*grad(i+1,Jstr))/            &
     &                      (cff+Ce)
            IF (LBC(isouth,isVice,ng)%nudging) THEN
              vi(i,Jstr,liunw)=vi(i,Jstr,liunw)+                        &
     &                      tau*(BOUNDARY(ng)%vi_south(i)-              &
     &                           vi(i,Jstr,know))
            END IF
# ifdef MASKING
            vi(i,Jstr,liunw)=vi(i,Jstr,liunw)*                          &
     &                      GRID(ng)%vmask(i,Jstr)
# endif
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (LBC(isouth,isVice,ng)%clamped) THEN
          DO i=Istr,Iend
            vi(i,1,liunw)=BOUNDARY(ng)%vi_south(i)
# ifdef MASKING
            vi(i,1,liunw)=vi(i,1,liunw)*                                &
     &                   GRID(ng)%vmask(i,1)
# endif
# ifdef WET_DRY
            vi(i,1,liunw)=vi(i,1,liunw)*                                &
     &                   GRID(ng)%vmask_wet(i,1)
# endif
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (LBC(isouth,isVice,ng)%gradient) THEN
          DO i=Istr,Iend
            vi(i,1,liunw)=vi(i,2,liunw)
# ifdef MASKING
            vi(i,1,liunw)=vi(i,1,liunw)*                                &
     &                   GRID(ng)%vmask(i,1)
# endif
# ifdef WET_DRY
            vi(i,1,liunw)=vi(i,1,liunw)*                                &
     &                   GRID(ng)%vmask_wet(i,1)
# endif
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE IF (LBC(isouth,isVice,ng)%closed) THEN
          DO i=Istr,Iend
            vi(i,1,liunw)=0.0_r8
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        IF (LBC(inorth,isVice,ng)%radiation) THEN
!
!  Northern edge, implicit upstream radiation condition.
!
          DO i=Istr,Iend+1
            grad(i,Jend  )=vi(i  ,Jend  ,know)-                         &
     &                     vi(i-1,Jend  ,know)
            grad(i,Jend+1)=vi(i  ,Jend+1,know)-                         &
     &                     vi(i-1,Jend+1,know)
          END DO
          DO i=Istr,Iend
            dVdt=vi(i,Jend,know)-vi(i,Jend  ,liunw)
            dVde=vi(i,Jend,liunw)-vi(i,Jend-1,liunw)
            IF (LBC(inorth,isVice,ng)%nudging) THEN
              IF ((dVdt*dVde).lt.0.0_r8) THEN
                tau=M2obc_in(ng,inorth)
              ELSE
                tau=M2obc_out(ng,inorth)
              END IF
              tau=tau*dt(ng)
            END IF
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
            vi(i,Jend+1,liunw)=(cff*vi(i,Jend+1,know)+                  &
     &                         Ce *vi(i,Jend  ,liunw)-                  &
     &                         MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-         &
     &                         MIN(Cx,0.0_r8)*grad(i+1,Jend+1))/        &
     &                        (cff+Ce)
            IF (LBC(inorth,isVice,ng)%nudging) THEN
              vi(i,Jend+1,liunw)=vi(i,Jend+1,liunw)+                    &
     &                         tau*(BOUNDARY(ng)%vi_north(i)-           &
     &                              vi(i,Jend+1,know))
            END IF
# ifdef MASKING
            vi(i,Jend+1,liunw)=vi(i,Jend+1,liunw)* &
     &                        GRID(ng)%vmask(i,Jend+1)
# endif
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (LBC(inorth,isVice,ng)%clamped) THEN
          DO i=Istr,Iend
            vi(i,Mm(ng)+1,liunw)=BOUNDARY(ng)%vi_north(i)
# ifdef MASKING
            vi(i,Mm(ng)+1,liunw)=vi(i,Mm(ng)+1,liunw)*                  &
     &                          GRID(ng)%vmask(i,Mm(ng))
# endif
# ifdef WET_DRY
            vi(i,Mm(ng)+1,liunw)=vi(i,Mm(ng)+1,liunw)*                  &
     &                          GRID(ng)%vmask_wet(i,Mm(ng))
# endif
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (LBC(inorth,isVice,ng)%gradient) THEN
          DO i=Istr,Iend
            vi(i,Mm(ng)+1,liunw)=vi(i,Mm(ng),liunw)
# ifdef MASKING
            vi(i,Mm(ng)+1,liunw)=vi(i,Mm(ng)+1,liunw)*                  &
     &                          GRID(ng)%vmask(i,Mm(ng))
# endif
# ifdef WET_DRY
            vi(i,Mm(ng)+1,liunw)=vi(i,Mm(ng)+1,liunw)*                  &
     &                          GRID(ng)%vmask_wet(i,Mm(ng))
# endif
          END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (LBC(inorth,isVice,ng)%closed) THEN
          DO i=Istr,Iend
            vi(i,Mm(ng)+1,liunw)=0.0_r8
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
          vi(0,1,liunw)=0.5_r8*(vi(0,2,liunw)+                          &
     &                         vi(1,1,liunw))
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          vi(Lm(ng)+1,1,liunw)=0.5_r8*(vi(Lm(ng)  ,1,liunw)+            &
     &                                vi(Lm(ng)+1,2,liunw))
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          vi(0,Mm(ng)+1,liunw)=0.5_r8*(vi(0,Mm(ng)  ,liunw)+            &
     &                                vi(1,Mm(ng)+1,liunw))
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          vi(Lm(ng)+1,Mm(ng)+1,liunw)=0.5_r8*                           &
     &                               (vi(Lm(ng)+1,Mm(ng)  ,liunw)+      &
     &                                vi(Lm(ng)  ,Mm(ng)+1,liunw))
        END IF
      END IF
      RETURN
      END SUBROUTINE vibc_tile
#endif

      END MODULE vibc_mod
