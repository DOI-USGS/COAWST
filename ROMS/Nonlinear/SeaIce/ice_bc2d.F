#include "cppdefs.h"
      MODULE ice_bc2d_mod

#ifdef ICE_MODEL
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license           W. Paul Budgell    !
!    See License_ROMS.md                            Katherine Hedstrom !
!================================================== Hernan G. Arango ===
!                                                                      !
!  Sets lateral boundary conditions for any 2D ice field with array    !
!  index, ifld.
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
      PUBLIC  ice_bc2d_tile
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_bc2d_tile (ng, tile, model, ifld,                  &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          liold, linew,                           &
     &                          ui, vi, field, S)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, ifld
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: liold, linew
!
      TYPE(T_LBC), intent(in) :: S(4)
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)    :: ui(LBi:,LBj:,:)
      real(r8), intent(in)    :: vi(LBi:,LBj:,:)
      real(r8), intent(inout) :: field(LBi:,LBj:,:)
# else
      real(r8), intent(in)    :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: field(LBi:UBi,LBj:UBj,2)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, know
!
      real(r8), parameter :: eps =1.0E-20_r8
      real(r8) :: Ce, Cx, cff, dTde, dTdt, dTdx, tau

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set time-indices.
!-----------------------------------------------------------------------
!
        know=liold
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!
!  Western edge, implicit upstream radiation condition.
!
        IF (S(iwest)%radiation) THEN
          DO j=Jstr,Jend+1
            grad(Istr-1,j)=field(Istr-1,j,know)-field(Istr-1,j-1,know)
# ifdef MASKING
            grad(Istr-1,j)=grad(Istr-1,j)*GRID(ng)%vmask(Istr-1,j)
# endif
            grad(Istr,j)=field(Istr,j,know)-field(Istr,j-1,know)
# ifdef MASKING
            grad(Istr,j)=grad(Istr,j)*GRID(ng)%vmask(Istr,j)
# endif
          END DO
          DO j=Jstr,Jend
            dTdt=field(Istr,j,know )-field(Istr  ,j,linew)
            dTdx=field(Istr,j,linew)-field(Istr+1,j,linew)
            IF (S(iwest)%nudging) THEN
              tau=Tobc_out(1,ng,iwest)
              IF ((dTdt*dTdx).lt.0.0_r8) tau=Tobc_in(1,ng,iwest)
              tau=tau*dt(ng)
            END IF
            IF ((dTdt*dTdx).lt.0.0_r8) dTdt=0.0_r8
            IF ((dTdt*(grad(Istr,j)+grad(Istr,j+1))).gt.0.0_r8) THEN
              dTde=grad(Istr,j  )
            ELSE
              dTde=grad(Istr,j+1)
            END IF
            cff=MAX(dTdx*dTdx+dTde*dTde,eps)
            Cx=dTdt*dTdx
# ifdef RADIATION_2D
            Ce=MIN(cff,MAX(dTdt*dTde,-cff))
# else
            Ce=0.0_r8
# endif
            field(Istr-1,j,linew)=(cff*field(Istr-1,j,know)+            &
     &                             Cx *field(Istr  ,j,linew)-           &
     &                             MAX(Ce,0.0_r8)*grad(Istr-1,j  )-     &
     &                             MIN(Ce,0.0_r8)*grad(Istr-1,j+1))/    &
     &                            (cff+Cx)
            IF (S(iwest)%nudging) THEN
              field(Istr-1,j,linew)=field(Istr-1,j,linew)+              &
     &                            tau*(BOUNDARY(ng)%ice_west(ifld,j)-   &
     &                                 field(Istr-1,j,know))
            END IF
# ifdef MASKING
            field(Istr-1,j,linew)=field(Istr-1,j,linew)*                &
     &                            GRID(ng)%rmask(Istr-1,j)
# endif
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (S(iwest)%clamped) THEN
          DO j=Jstr,Jend
            field(Istr-1,j,linew)=BOUNDARY(ng)%ice_west(ifld,j)
# ifdef MASKING
            field(Istr-1,j,linew)=field(Istr-1,j,linew)*                &
     &                            GRID(ng)%rmask(Istr-1,j)
# endif
# ifdef WET_DRY
            field(Istr-1,j,linew)=field(Istr-1,j,linew)*                &
     &                            GRID(ng)%rmask_wet(Istr-1,j)
# endif
          END DO
!
!  Western edge, clamped on inflow, gradient on outflow.
!
        ELSE IF (S(iwest)%mixed) THEN
          DO j=Jstr,Jend
            IF (ui(1,j,linew).ge.0.0_r8) THEN
              field(Istr-1,j,linew)=BOUNDARY(ng)%ice_west(ifld,j)
# ifdef MASKING
              field(Istr-1,j,linew)=field(Istr-1,j,linew)*              &
     &                              GRID(ng)%rmask(Istr-1,j)
# endif
# ifdef WET_DRY
              field(Istr-1,j,linew)=field(Istr-1,j,linew)*              &
     &                              GRID(ng)%rmask_wet(Istr-1,j)
# endif
            ELSE
              field(Istr-1,j,linew)=field(Istr,j,liold)
# ifdef MASKING
              field(Istr-1,j,linew)=field(Istr-1,j,linew)*              &
     &                              GRID(ng)%rmask(Istr-1,j)
# endif
# ifdef WET_DRY
              field(Istr-1,j,linew)=field(Istr-1,j,linew)*              &
     &                              GRID(ng)%rmask_wet(Istr-1,j)
# endif
            END IF
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (S(iwest)%closed) THEN
          DO j=Jstr,Jend
            field(Istr-1,j,linew)=field(Istr,j,linew)
# ifdef MASKING
            field(Istr-1,j,linew)=field(Istr-1,j,linew)*                &
     &                            GRID(ng)%rmask(Istr-1,j)
# endif
# ifdef WET_DRY
            field(Istr-1,j,linew)=field(Istr-1,j,linew)*                &
     &                            GRID(ng)%rmask_wet(Istr-1,j)
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
!
!  Eastern edge, implicit upstream radiation condition.
!
        IF (S(ieast)%radiation) THEN
          DO j=Jstr,Jend+1
            grad(Iend,j)=field(Iend,j,know)-field(Iend,j-1,know)
# ifdef MASKING
            grad(Iend,j)=grad(Iend,j)*GRID(ng)%vmask(Iend  ,j)
# endif
            grad(Iend+1,j)=field(Iend+1,j,know)-field(Iend+1,j-1,know)
# ifdef MASKING
            grad(Iend+1,j)=grad(Iend+1,j)*GRID(ng)%vmask(Iend+1,j)
# endif
          END DO
          DO j=Jstr,Jend
            dTdt=field(Iend,j,know )-field(Iend  ,j,linew)
            dTdx=field(Iend,j,linew)-field(Iend-1,j,linew)
            IF (S(ieast)%nudging) THEN
              tau=Tobc_out(1,ng,ieast)
              IF ((dTdt*dTdx).lt.0.0_r8) tau=Tobc_in(1,ng,ieast)
              tau=tau*dt(ng)
            END IF
            IF ((dTdt*dTdx).lt.0.0_r8) dTdt=0.0_r8
              IF ((dTdt*(grad(Iend,j)+grad(Iend,j+1))).gt.0.0_r8) THEN
              dTde=grad(Iend,j  )
            ELSE
              dTde=grad(Iend,j+1)
            END IF
            cff=MAX(dTdx*dTdx+dTde*dTde,eps)
            Cx=dTdt*dTdx
# ifdef RADIATION_2D
            Ce=MIN(cff,MAX(dTdt*dTde,-cff))
# else
            Ce=0.0_r8
# endif
            field(Iend+1,j,linew)=(cff*field(Iend+1,j,know)+            &
     &                             Cx *field(Iend  ,j,linew)-           &
     &                             MAX(Ce,0.0_r8)*grad(Iend+1,j  )-     &
     &                             MIN(Ce,0.0_r8)*grad(Iend+1,j+1))/    &
     &                            (cff+Cx)
            IF (S(ieast)%nudging) THEN
              field(Iend+1,j,linew)=field(Iend+1,j,linew)+              &
     &                            tau*(BOUNDARY(ng)%ice_east(ifld,j)-   &
     &                                 field(Iend+1,j,know))
            END IF
# ifdef MASKING
            field(Iend+1,j,linew)=field(Iend+1,j,linew)*                &
     &                            GRID(ng)%rmask(Iend+1,j)
# endif
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (S(ieast)%clamped) THEN
          DO j=Jstr,Jend
            field(Iend+1,j,linew)=BOUNDARY(ng)%ice_east(ifld,j)
# ifdef MASKING
            field(Iend+1,j,linew)=field(Iend+1,j,linew)*                &
     &                            GRID(ng)%rmask(Iend+1,j)
# endif
# ifdef WET_DRY
            field(Iend+1,j,linew)=field(Iend+1,j,linew)*                &
     &                            GRID(ng)%rmask_wet(Iend+1,j)
# endif
          END DO
!
!  Eastern edge, clamped on inflow, gradient on outflow.
!
        ELSE IF (S(ieast)%mixed) THEN
          DO j=Jstr,Jend
            IF (ui(Iend+1,j,linew).le.0.0_r8) THEN
              field(Iend+1,j,linew)=BOUNDARY(ng)%ice_east(ifld,j)
# ifdef MASKING
              field(Iend+1,j,linew)=field(Iend+1,j,linew)*              &
     &                              GRID(ng)%rmask(Iend+1,j)
# endif
# ifdef WET_DRY
              field(Iend+1,j,linew)=field(Iend+1,j,linew)*              &
     &                              GRID(ng)%rmask_wet(Iend+1,j)
# endif
            ELSE
              field(Iend+1,j,linew)=field(Iend,j,liold)
# ifdef MASKING
              field(Iend+1,j,linew)=field(Iend+1,j,linew)*              &
     &                              GRID(ng)%rmask(Iend+1,j)
# endif
# ifdef WET_DRY
              field(Iend+1,j,linew)=field(Iend+1,j,linew)*              &
     &                              GRID(ng)%rmask_wet(Iend+1,j)
# endif
            END IF
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (S(ieast)%closed) THEN
          DO j=Jstr,Jend
            field(Iend+1,j,linew)=field(Iend,j,linew)
# ifdef MASKING
            field(Iend+1,j,linew)=field(Iend+1,j,linew)*                &
     &                            GRID(ng)%rmask(Iend+1,j)
# endif
# ifdef WET_DRY
            field(Iend+1,j,linew)=field(Iend+1,j,linew)*                &
     &                            GRID(ng)%rmask_wet(Iend+1,j)
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
!
!  Southern edge, implicit upstream radiation condition.
!
        IF (S(isouth)%radiation) THEN
          DO i=Istr,Iend+1
            grad(i,Jstr)=field(i,Jstr,know)-field(i-1,Jstr,know)
# ifdef MASKING
            grad(i,Jstr)=grad(i,Jstr)*GRID(ng)%umask(i,Jstr)
# endif
            grad(i,Jstr-1)=field(i,Jstr-1,know)-field(i-1,Jstr-1,know)
# ifdef MASKING
            grad(i,Jstr-1)=grad(i,Jstr-1)*GRID(ng)%umask(i,Jstr-1)
# endif
          END DO
          DO i=Istr,Iend
            dTdt=field(i,Jstr,know )-field(i,Jstr  ,linew)
            dTde=field(i,Jstr,linew)-field(i,Jstr+1,linew)
            IF (S(isouth)%nudging) THEN
              tau=Tobc_out(1,ng,isouth)
              IF ((dTdt*dTde).lt.0.0_r8) tau=Tobc_in(1,ng,isouth)
              tau=tau*dt(ng)
            END IF
            IF ((dTdt*dTde).lt.0.0_r8) dTdt=0.0_r8
            IF ((dTdt*(grad(i,Jstr)+grad(i+1,Jstr))).gt.0.0_r8) THEN
              dTdx=grad(i  ,Jstr)
            ELSE
              dTdx=grad(i+1,Jstr)
            END IF
            cff=MAX(dTdx*dTdx+dTde*dTde,eps)
# ifdef RADIATION_2D
            Cx=MIN(cff,MAX(dTdt*dTdx,-cff))
# else
            Cx=0.0_r8
# endif
            Ce=dTdt*dTde
            field(i,Jstr-1,linew)=(cff*field(i,Jstr-1,know)+            &
     &                             Ce *field(i,Jstr  ,linew)-           &
     &                             MAX(Cx,0.0_r8)*grad(i  ,Jstr-1)-     &
     &                             MIN(Cx,0.0_r8)*grad(i+1,Jstr-1))/    &
     &                            (cff+Ce)
            IF (S(isouth)%nudging) THEN
              field(i,Jstr-1,linew)=field(i,Jstr-1,linew)+              &
     &                            tau*(BOUNDARY(ng)%ice_south(ifld,i)-  &
     &                                 field(i,Jstr-1,know))
            END IF
# ifdef MASKING
            field(i,Jstr-1,linew)=field(i,Jstr-1,linew)*                &
     &                            GRID(ng)%rmask(i,Jstr-1)
# endif
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (S(isouth)%clamped) THEN
          DO i=Istr,Iend
            field(i,Jstr-1,linew)=BOUNDARY(ng)%ice_south(ifld,i)
# ifdef MASKING
            field(i,Jstr-1,linew)=field(i,Jstr-1,linew)*                &
     &                            GRID(ng)%rmask(i,Jstr-1)
# endif
# ifdef WET_DRY
            field(i,Jstr-1,linew)=field(i,Jstr-1,linew)*                &
     &                           GRID(ng)%rmask_wet(i,Jstr-1)
# endif
          END DO
!
!  Southern edge, clamped on inflow, gradient on outflow.
!
        ELSE IF (S(isouth)%mixed) THEN
          DO i=Istr,Iend
            IF (vi(i,1,linew).ge.0._r8) THEN
              field(i,Jstr-1,linew)=BOUNDARY(ng)%ice_south(ifld,i)
# ifdef MASKING
              field(i,Jstr-1,linew)=field(i,Jstr-1,linew)*              &
     &                              GRID(ng)%rmask(i,Jstr-1)
# endif
# ifdef WET_DRY
              field(i,Jstr-1,linew)=field(i,Jstr-1,linew)*              &
     &                              GRID(ng)%rmask_wet(i,Jstr-1)
# endif
            ELSE
              field(i,Jstr-1,linew)=field(i,Jstr,liold)
# ifdef MASKING
              field(i,Jstr-1,linew)=field(i,Jstr-1,linew)*              &
     &                              GRID(ng)%rmask(i,Jstr-1)
# endif
# ifdef WET_DRY
              field(i,Jstr-1,linew)=field(i,Jstr-1,linew)*              &
     &                              GRID(ng)%rmask_wet(i,Jstr-1)
# endif
            END IF
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE IF (S(isouth)%closed) THEN
          DO i=Istr,Iend
            field(i,Jstr-1,linew)=field(i,Jstr,linew)
# ifdef MASKING
            field(i,Jstr-1,linew)=field(i,Jstr-1,linew)*                &
     &                            GRID(ng)%rmask(i,Jstr-1)
# endif
# ifdef WET_DRY
            field(i,Jstr-1,linew)=field(i,Jstr-1,linew)*                &
     &                            GRID(ng)%rmask_wet(i,Jstr-1)
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
!
!  Northern edge, implicit upstream radiation condition.
!
        IF (S(inorth)%radiation) THEN
          DO i=Istr,Iend+1
            grad(i,Jend)=field(i,Jend,know)-field(i-1,Jend,know)
# ifdef MASKING
            grad(i,Jend)=grad(i,Jend)*GRID(ng)%umask(i,Jend)
# endif
            grad(i,Jend+1)=field(i,Jend+1,know)-field(i-1,Jend+1,know)
# ifdef MASKING
            grad(i,Jend+1)=grad(i,Jend+1)*GRID(ng)%umask(i,Jend+1)
# endif
          END DO
          DO i=Istr,Iend
            dTdt=field(i,Jend,know )-field(i,Jend  ,linew)
            dTde=field(i,Jend,linew)-field(i,Jend-1,linew)
            IF (S(inorth)%nudging) THEN
              tau=Tobc_out(1,ng,inorth)
              IF ((dTdt*dTde).lt.0.0_r8) tau=Tobc_in(1,ng,inorth)
              tau=tau*dt(ng)
            END IF
            IF ((dTdt*dTde).lt.0.0_r8) dTdt=0.0_r8
            IF ((dTdt*(grad(i,Jend)+grad(i+1,Jend))).gt.0.0_r8) THEN
              dTdx=grad(i  ,Jend)
            ELSE
              dTdx=grad(i+1,Jend)
            END IF
            cff=MAX(dTdx*dTdx+dTde*dTde,eps)
# ifdef RADIATION_2D
            Cx=MIN(cff,MAX(dTdt*dTdx,-cff))
# else
            Cx=0.0_r8
# endif
            Ce=dTdt*dTde
            field(i,Jend+1,linew)=(cff*field(i,Jend+1,know)+            &
     &                             Ce *field(i,Jend  ,linew)-           &
     &                             MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-     &
     &                             MIN(Cx,0.0_r8)*grad(i+1,Jend+1))/    &
     &                            (cff+Ce)
            IF (S(inorth)%nudging) THEN
              field(i,Jend+1,linew)=field(i,Jend+1,linew)+              &
     &                              tau*(BOUNDARY(ng)%ice_north(ifld,i)-&
     &                                   field(i,Jend+1,know))
            END IF
# ifdef MASKING
            field(i,Jend+1,linew)=field(i,Jend+1,linew)*                &
     &                            GRID(ng)%rmask(i,Jend+1)
# endif
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (S(inorth)%clamped) THEN
          DO i=Istr,Iend
            field(i,Jend+1,linew)=BOUNDARY(ng)%ice_north(ifld,i)
# ifdef MASKING
            field(i,Jend+1,linew)=field(i,Jend+1,linew)*                &
     &                            GRID(ng)%rmask(i,Jend+1)
# endif
# ifdef WET_DRY
            field(i,Jend+1,linew)=field(i,Jend+1,linew)*                &
     &                            GRID(ng)%rmask_wet(i,Jend+1)
# endif
          END DO
!
!  Northern edge, clamped on inflow, gradient on outflow.
!
        ELSE IF (S(inorth)%mixed) THEN
          DO i=Istr,Iend
            IF (vi(i,Jend+1,linew).le.0.0_r8) THEN
              field(i,Jend+1,linew)=BOUNDARY(ng)%ice_north(ifld,i)
# ifdef MASKING
              field(i,Jend+1,linew)=field(i,Mm(ng)+1,linew)*            &
     &                              GRID(ng)%rmask(i,Jend+1)
# endif
# ifdef WET_DRY
              field(i,Jend+1,linew)=field(i,Jend+1,linew)*              &
     &                              GRID(ng)%rmask_wet(i,Jend+1)
# endif
            ELSE
              field(i,Jend+1,linew)=field(i,Jend,liold)
# ifdef MASKING
              field(i,Jend+1,linew)=field(i,Jend+1,linew)*              &
     &                              GRID(ng)%rmask(i,Jend+1)
# endif
# ifdef WET_DRY
              field(i,Jend+1,linew)=field(i,Jend+1,linew)*              &
     &                              GRID(ng)%rmask_wet(i,Jend+1)
# endif
            END IF
          END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (S(inorth)%closed) THEN
          DO i=Istr,Iend
            field(i,Jend+1,linew)=field(i,Jend,linew)
# ifdef MASKING
            field(i,Jend+1,linew)=field(i,Jend+1,linew)*                &
     &                            GRID(ng)%rmask(i,Jend+1)
# endif
# ifdef WET_DRY
            field(i,Jend+1,linew)=field(i,Jend+1,linew)*                &
     &                            GRID(ng)%rmask_wet(i,Jend+1)
# endif
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng).and. .not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          field(Istr-1,Jstr-1,linew)=0.5_r8*                            &
     &                               (field(Istr  ,Jstr-1,linew)+       &
     &                                field(Istr-1,Jstr  ,linew))
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          field(Iend+1,Jstr-1,linew)=0.5_r8*                            &
     &                               (field(Iend+1,Jstr  ,linew)+       &
     &                                field(Iend  ,Jstr-1,linew))
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          field(Istr-1,Jend+1,linew)=0.5_r8*                            &
     &                               (field(Istr-1,Jend  ,linew)+       &
     &                                field(Istr  ,Jend+1,linew))
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          field(Iend+1,Jend+1,linew)=0.5_r8*                            &
     &                               (field(Iend+1,Jend  ,linew)+       &
     &                                field(Iend  ,Jend+1,linew))
        END IF
      END IF
      RETURN
      END SUBROUTINE ice_bc2d_tile
#endif

      END MODULE ice_bc2d_mod
