#include "cppdefs.h"
      MODULE ice_uibc_mod

#ifdef ICE_MODEL
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license           W. Paul Budgell    !
!    See License_ROMS.md                            Katherine Hedstrom !
!================================================== Hernan G. Arango ===
!                                                                      !
!  Sets the lateral boundary conditions on the the ice U-velocity.     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ice
      USE mod_scalars
      USE mod_boundary
!
      implicit none
!
      PRIVATE
      PUBLIC  ice_uibc_tile
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_uibc (ng, tile, model)
!***********************************************************************
!
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
#include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 42, __LINE__, MyFile)
# endif
      CALL ice_uibc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liuol(ng), liunw(ng),                         &
     &                    ICE(ng) % Si(:,:,:,isUice))
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE ice_uibc
!
!***********************************************************************
      SUBROUTINE ice_uibc_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          liuol, liunw,                           &
     &                          ui)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: liuol, liunw
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ui(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: ui(LBi:UBi,LBj:UBj,2)
# endif
!
!  Local variable declarations.
!
      integer :: i, Imax, Imin, j, know
!
      real(r8), parameter :: eps =1.0E-20_r8
      real(r8) :: Ce, Cx, cff, dUde, dUdt, dUdx, tau

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
        IF (LBC(iwest,ibICE(isUice),ng)%radiation) THEN
!
!  Western edge, implicit upstream radiation condition.
!
          DO j=Jstr,Jend+1
            grad(Istr  ,j)=ui(Istr  ,j  ,know)-                         &
     &                     ui(Istr  ,j-1,know)
            grad(Istr+1,j)=ui(Istr+1,j  ,know)-                         &
     &                     ui(Istr+1,j-1,know)
          END DO
          DO j=Jstr,Jend
            dUdt=ui(Istr+1,j,know )-ui(Istr+1,j,liunw)
            dUdx=ui(Istr+1,j,liunw)-ui(Istr+2,j,liunw)
            IF (LBC(iwest,ibICE(isUice),ng)%nudging) THEN
              IF ((dUdt*dUdx).lt.0.0_r8) THEN
                tau=M2obc_in(ng,iwest)
              ELSE
                tau=M2obc_out(ng,iwest)
              END IF
              tau=tau*dt(ng)
            END IF
            IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
            IF ((dUdt*(grad(Istr+1,j)+grad(Istr+1,j+1))).gt.0.0_r8) THEN
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
            ui(Istr,j,liunw)=(cff*ui(Istr  ,j,know)+                    &
     &                        Cx *ui(Istr+1,j,liunw)-                   &
     &                        MAX(Ce,0.0_r8)*grad(Istr,j  )-            &
     &                        MIN(Ce,0.0_r8)*grad(Istr,j+1))/           &
     &                       (cff+Cx)
            IF (LBC(iwest,ibICE(isUice),ng)%nudging) THEN
              ui(Istr,j,liunw)=ui(Istr,j,liunw)+                        &
     &                         tau*(BOUNDARY(ng)%ice_west(isUice,j)-    &
     &                              ui(Istr,j,know))
            END IF
# ifdef MASKING
            ui(Istr,j,liunw)=ui(Istr,j,liunw)*                          &
     &                       GRID(ng)%umask(Istr,j)
# endif
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (LBC(iwest,ibICE(isUice),ng)%clamped) THEN
          DO j=Jstr,Jend
            ui(Istr,j,liunw)=BOUNDARY(ng)%ice_west(isUice,j)
# ifdef MASKING
            ui(Istr,j,liunw)=ui(Istr,j,liunw)*                          &
     &                       GRID(ng)%umask(Istr,j)
# endif
# ifdef WET_DRY
            ui(Istr,j,liunw)=ui(Istr,j,liunw)*                          &
     &                       GRID(ng)%umask_wet(Istr,j)
# endif
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (LBC(iwest,ibICE(isUice),ng)%gradient) THEN
          DO j=Jstr,Jend
            ui(Istr,j,liunw)=ui(Istr+1,j,liunw)
# ifdef MASKING
            ui(Istr,j,liunw)=ui(Istr,j,liunw)*                          &
     &                       GRID(ng)%umask(Istr,j)
# endif
# ifdef WET_DRY
            ui(Istr,j,liunw)=ui(Istr,j,liunw)*                          &
     &                       GRID(ng)%umask_wet(Istr,j)
# endif
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (LBC(iwest,ibICE(isUice),ng)%closed) THEN
          DO j=Jstr,Jend
            ui(Istr,j,liunw)=0.0_r8
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        IF (LBC(ieast,ibICE(isUice),ng)%radiation) THEN
!
!  Eastern edge, implicit upstream radiation condition.
!
          DO j=Jstr,Jend+1
            grad(Iend  ,j)=ui(Iend  ,j  ,know)-                         &
     &                     ui(Iend  ,j-1,know)
            grad(Iend+1,j)=ui(Iend+1,j  ,know)-                         &
     &                     ui(Iend+1,j-1,know)
          END DO
          DO j=Jstr,Jend
            dUdt=ui(Iend,j,know )-ui(Iend  ,j,liunw)
            dUdx=ui(Iend,j,liunw)-ui(Iend-1,j,liunw)
            IF (LBC(ieast,ibICE(isUice),ng)%nudging) THEN
              IF ((dUdt*dUdx).lt.0.0_r8) THEN
                tau=M2obc_in(ng,ieast)
              ELSE
              tau=M2obc_out(ng,ieast)
              END IF
              tau=tau*dt(ng)
            END IF
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
            ui(Iend+1,j,liunw)=(cff*ui(Iend+1,j,know)+                  &
     &                          Cx *ui(Iend  ,j,liunw)-                 &
     &                          MAX(Ce,0.0_r8)*grad(Iend+1,j  )-        &
     &                          MIN(Ce,0.0_r8)*grad(Iend+1,j+1))/       &
     &                         (cff+Cx)
            IF (LBC(ieast,ibICE(isUice),ng)%nudging) THEN
              ui(Iend+1,j,liunw)=ui(Iend+1,j,liunw)+ &
     &                           tau*(BOUNDARY(ng)%ice_east(isUice,j)-  &
     &                                ui(Iend+1,j,know))
            END IF
# ifdef MASKING
            ui(Iend+1,j,liunw)=ui(Iend+1,j,liunw)*                      &
     &                         GRID(ng)%umask(Iend+1,j)
# endif
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (LBC(ieast,ibICE(isUice),ng)%clamped) THEN
          DO j=Jstr,Jend
            ui(Iend+1,j,liunw)=BOUNDARY(ng)%ice_east(isUice,j)
# ifdef MASKING
            ui(Iend+1,j,liunw)=ui(Iend+1,j,liunw)*                      &
     &                         GRID(ng)%umask(Iend+1,j)
# endif
# ifdef WET_DRY
            ui(Iend+1,j,liunw)=ui(Iend+1,j,liunw)*                      &
     &                          GRID(ng)%umask_wet(Iend+1,j)
# endif
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (LBC(ieast,ibICE(isUice),ng)%gradient) THEN
          DO j=Jstr,Jend
            ui(Iend+1,j,liunw)=ui(Iend,j,liunw)
# ifdef MASKING
            ui(Iend+1,j,liunw)=ui(Iend+1,j,liunw)*                      &
     &                         GRID(ng)%umask(Iend+1,j)
# endif
# ifdef WET_DRY
            ui(Iend+1,j,liunw)=ui(Iend+1,j,liunw)*                      &
     &                         GRID(ng)%umask_wet(Iend+1,j)
# endif
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (LBC(ieast,ibICE(isUice),ng)%closed) THEN
          DO j=Jstr,Jend
            ui(Iend+1,j,liunw)=0.0_r8
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        IF (LBC(isouth,ibICE(isUice),ng)%radiation) THEN
!
!  Southern edge, implicit upstream radiation condition.
!
          DO i=IstrP-1,Iend
            grad(i,Jstr-1)=ui(i+1,Jstr-1,know)-                         &
     &                     ui(i  ,Jstr-1,know)
            grad(i,Jstr  )=ui(i+1,Jstr  ,know)-                         &
     &                     ui(i  ,Jstr  ,know)
          END DO
          DO i=IstrP,Iend
            dUdt=ui(i,Jstr,know )-ui(i,Jstr  ,liunw)
            dUde=ui(i,Jstr,liunw)-ui(i,Jstr+1,liunw)
            IF (LBC(isouth,ibICE(isUice),ng)%nudging) THEN
              IF ((dUdt*dUde).lt.0.0_r8) THEN
                tau=M2obc_in(ng,isouth)
              ELSE
                tau=M2obc_out(ng,isouth)
              END IF
              tau=tau*dt(ng)
            END IF
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
            ui(i,Jstr-1,liunw)=(cff*ui(i,Jstr-1,know)+                  &
     &                          Ce *ui(i,Jstr  ,liunw)-                 &
     &                          MAX(Cx,0.0_r8)*grad(i-1,Jstr-1)-        &
     &                          MIN(Cx,0.0_r8)*grad(i  ,Jstr-1))/       &
     &                         (cff+Ce)
            IF (LBC(isouth,ibICE(isUice),ng)%nudging) THEN
              ui(i,Jstr-1,liunw)=ui(i,Jstr-1,liunw)+                    &
     &                           tau*(BOUNDARY(ng)%ice_south(isUice,i)- &
     &                                ui(i,Jstr-1,know))
            END IF
# ifdef MASKING
            ui(i,Jstr-1,liunw)=ui(i,Jstr-1,liunw)*                      &
     &                         GRID(ng)%umask(i,Jstr-1)
# endif
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (LBC(isouth,ibICE(isUice),ng)%clamped) THEN
          DO i=IstrP,Iend
            ui(i,Jstr-1,liunw)=BOUNDARY(ng)%ice_south(isUice,i)
# ifdef MASKING
            ui(i,Jstr-1,liunw)=ui(i,Jstr-1,liunw)*                      &
     &                         GRID(ng)%umask(i,Jstr-1)
# endif
# ifdef WET_DRY
            ui(i,Jstr-1,liunw)=ui(i,Jstr-1,liunw)*                      &
     &                         GRID(ng)%umask_wet(i,Jstr-1)
# endif
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (LBC(isouth,ibICE(isUice),ng)%gradient) THEN
          DO i=IstrP,Iend
            ui(i,Jstr-1,liunw)=ui(i,1,liunw)
# ifdef MASKING
            ui(i,Jstr-1,liunw)=ui(i,Jstr-1,liunw)*                      &
     &                         GRID(ng)%umask(i,Jstr-1)
# endif
# ifdef WET_DRY
            ui(i,Jstr-1,liunw)=ui(i,Jstr-1,liunw)*                      &
     &                         GRID(ng)%umask_wet(i,Jstr-1)
# endif
          END DO
!
!  Southern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
        ELSE IF (LBC(isouth,ibICE(isUice),ng)%closed) THEN
          IF (EWperiodic(ng)) THEN
            Imin=IstrP
            Imax=Iend
          ELSE
            Imin=Istr
            Imax=IendT
          END IF
          DO i=Imin,Imax
            ui(i,Jstr-1,liunw)=gamma2(ng)*ui(i,1,liunw)
# ifdef MASKING
            ui(i,Jstr-1,liunw)=ui(i,Jstr-1,liunw)*                      &
     &                         GRID(ng)%umask(i,Jstr-1)
# endif
# ifdef WET_DRY
            ui(i,Jstr-1,liunw)=ui(i,Jstr-1,liunw)*                      &
     &                         GRID(ng)%umask_wet(i,Jstr-1)
# endif
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        IF (LBC(inorth,ibICE(isUice),ng)%radiation) THEN
!
!  Northern edge, implicit upstream radiation condition.
!
          DO i=IstrP-1,Iend
            grad(i,Jend  )=ui(i+1,Jend  ,know)-                         &
     &                     ui(i  ,Jend  ,know)
            grad(i,Jend+1)=ui(i+1,Jend+1,know)-                         &
     &                     ui(i  ,Jend+1,know)
          END DO
          DO i=IstrP,Iend
            dUdt=ui(i,Jend,know )-ui(i,Jend  ,liunw)
            dUde=ui(i,Jend,liunw)-ui(i,Jend-1,liunw)
            IF (LBC(inorth,ibICE(isUice),ng)%nudging) THEN
              IF ((dUdt*dUde).lt.0.0_r8) THEN
                tau=M2obc_in(ng,inorth)
              ELSE
                tau=M2obc_out(ng,inorth)
              END IF
              tau=tau*dt(ng)
            END IF
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
            ui(i,Jend+1,liunw)=(cff*ui(i,Jend+1,know)+                  &
     &                          Ce *ui(i,Jend  ,liunw)-                 &
     &                          MAX(Cx,0.0_r8)*grad(i-1,Jend+1)-        &
     &                          MIN(Cx,0.0_r8)*grad(i  ,Jend+1))/       &
     &                         (cff+Ce)
# ifdef NORTH_MINUDGING
            IF (LBC(inorth,ibICE(isUice),ng)%nudging) THEN
              ui(i,Jend+1,liunw)=ui(i,Jend+1,liunw)+                    &
     &                           tau*(BOUNDARY(ng)%ice_north(isUice,i)- &
     &                                ui(i,Jend+1,know))
# endif
# ifdef MASKING
            ui(i,Jend+1,liunw)=ui(i,Jend+1,liunw)*                      &
     &                         GRID(ng)%umask(i,Jend+1)
# endif
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (LBC(inorth,ibICE(isUice),ng)%clamped) THEN
          DO i=IstrP,Iend
            ui(i,Jend+1,liunw)=BOUNDARY(ng)%ice_north(isUice,i)
# ifdef MASKING
            ui(i,Jend+1,liunw)=ui(i,Jend+1,liunw)*                      &
     &                         GRID(ng)%umask(i,Jend+1)
# endif
# ifdef WET_DRY
            ui(i,Jend+1,liunw)=ui(i,Jend+1,liunw)*                      &
     &                         GRID(ng)%umask_wet(i,Jend+1)
# endif
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (LBC(inorth,ibICE(isUice),ng)%gradient) THEN
          DO i=IstrP,Iend
            ui(i,Jend+1,liunw)=ui(i,Jend,liunw)
# ifdef MASKING
            ui(i,Jend+1,liunw)=ui(i,Jend+1,liunw)*                      &
     &                         GRID(ng)%umask(i,Jend+1)
# endif
# ifdef WET_DRY
            ui(i,Jend+1,liunw)=ui(i,Jend+1,liunw)*                      &
     &                         GRID(ng)%umask_wet(i,Jend+1)
# endif
          END DO
!
!  Northern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
        ELSE IF (LBC(inorth,ibICE(isUice),ng)%closed) THEN
          IF (EWperiodic(ng)) THEN
            Imin=IstrP
            Imax=Iend
          ELSE
            Imin=Istr
            Imax=IendT
          END IF
          DO i=Imin,Imax
            ui(i,Jend+1,liunw)=gamma2(ng)*ui(i,Jend,liunw)
# ifdef MASKING
            ui(i,Jend+1,liunw)=ui(i,Jend+1,liunw)*                      &
     &                         GRID(ng)%umask(i,Jend+1)
# endif
# ifdef WET_DRY
            ui(i,Jend+1,liunw)=ui(i,Jend+1,liunw)*                      &
     &                         GRID(ng)%umask_wet(i,Jend+1)
# endif
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
          ui(Istr,Jstr-1,liunw)=0.5_r8*(ui(Istr+1,Jstr-1,liunw)+        &
     &                                  ui(Istr  ,Jstr  ,liunw))
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          ui(Iend+1,Jstr-1,liunw)=0.5_r8*(ui(Iend  ,Jstr-1,liunw)+      &
     &                                    ui(Iend+1,Jstr  ,liunw))
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          ui(Istr,Jend+1,liunw)=0.5_r8*(ui(Istr+1,Jend+1,liunw)+        &
     &                                  ui(Istr  ,Jend  ,liunw))
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          ui(Iend+1,Jend+1,liunw)=0.5_r8*(ui(Iend  ,Jend+1,liunw)+      &
     &                                    ui(Iend+1,Jend  ,liunw))
        END IF
      END IF
!
      RETURN
      END SUBROUTINE ice_uibc_tile
#endif
      END MODULE ice_uibc_mod
