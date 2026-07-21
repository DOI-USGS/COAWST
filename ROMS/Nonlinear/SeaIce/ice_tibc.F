#include "cppdefs.h"
      MODULE ice_tibc_mod

#ifdef ICE_MODEL
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license           W. Paul Budgell    !
!    See License_ROMS.md                            Katherine Hedstrom !
!================================================== Hernan G. Arango ===
!                                                                      !
!  Sets the lateral boundary conditions on the internal ice            !
!  temperature, ti.                                                    !
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
      PUBLIC  ice_tibc_tile
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_tibc (ng, tile, model)
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
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 42, __LINE__, MyFile)
# endif
      CALL ice_tibc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    liold(ng), linew(ng),                         &
     &                    ICE(ng) % Si(:,:,:,isUice),                   &
     &                    ICE(ng) % Si(:,:,:,isVice),                   &
     &                    ICE(ng) % Si(:,:,:,isHice),                   &
     &                    ICE(ng) % Si(:,:,:,isTice),                   &
     &                    ICE(ng) % Si(:,:,:,isEnth))
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE ice_tibc
!
!***********************************************************************
      SUBROUTINE ice_tibc_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          liold, linew,                           &
     &                          ui, vi, hi, ti, enthalpy)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: liold, linew
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)    :: ui(LBi:,LBj:,:)
      real(r8), intent(in)    :: vi(LBi:,LBj:,:)
      real(r8), intent(in)    :: hi(LBi:,LBj:,:)
      real(r8), intent(inout) :: ti(LBi:,LBj:,:)
      real(r8), intent(inout) :: enthalpy(LBi:,LBj:,:)
# else
      real(r8), intent(in)    :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: hi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ti(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: enthalpy(LBi:UBi,LBj:UBj,2)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, know
!
      real(r8), parameter :: eps=1.0E-6_r8

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set time-indices
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
!  Western edge, clamped boundary condition.
!
        IF (LBC(iwest,ibICE(isTice),ng)%clamped) THEN
          DO j=Jstr,Jend
            enthalpy(Istr-1,j,linew)=BOUNDARY(ng)%ice_west(isHice,j)*   &
     &                               BOUNDARY(ng)%ice_west(isTice,j)
# ifdef MASKING
            enthalpy(Istr-1,j,linew)=enthalpy(Istr-1,j,linew)*          &
     &                               GRID(ng)%rmask(Istr-1,j)
# endif
# ifdef WET_DRY
            enthalpy(Istr-1,j,linew)=enthalpy(Istr-1,j,linew)*          &
     &                               GRID(ng)%rmask_wet(Istr-1,j)
# endif
            ti(Istr-1,j,linew)=enthalpy(Istr-1,j,linew)/                &
     &                         MAX(hi(Istr-1,j,linew),eps)
            IF (hi(Istr-1,j,linew).le.min_hi(ng)) THEN
              enthalpy(Istr-1,j,linew)=0.0_r8
              ti(Istr-1,j,linew)=0.0_r8
            END IF
          END DO
!
!  Western edge, clamped on inflow, gradient on outflow.
!
        ELSE IF (LBC(iwest,ibICE(isTice),ng)%mixed) THEN
          DO j=Jstr,Jend
            IF (ui(Istr,j,linew).ge.0.0_r8) THEN
              enthalpy(Istr-1,j,linew)=BOUNDARY(ng)%ice_west(isHice,j)* &
     &                                 BOUNDARY(ng)%ice_west(isTice,j)
# ifdef MASKING
              enthalpy(Istr-1,j,linew)=enthalpy(Istr-1,j,linew)*        &
     &                                 GRID(ng)%rmask(Istr-1,j)
# endif
# ifdef WET_DRY
              enthalpy(Istr-1,j,linew)=enthalpy(Istr-1,j,linew)*        &
     &                                 GRID(ng)%rmask_wet(Istr-1,j)
# endif
            ELSE
              enthalpy(Istr-1,j,linew)=enthalpy(Istr,j,liold)
# ifdef MASKING
              enthalpy(Istr-1,j,linew)=enthalpy(Istr-1,j,linew)*        &
     &                                 GRID(ng)%rmask(Istr-1,j)
# endif
# ifdef WET_DRY
              enthalpy(Istr-1,j,linew)=enthalpy(Istr-1,j,linew)*        &
     &                                 GRID(ng)%rmask_wet(Istr-1,j)
# endif
            END IF
            ti(Istr-1,j,linew)=enthalpy(Istr-1,j,linew)/                &
     &                         MAX(hi(Istr-1,j,linew),eps)
            IF (hi(Istr-1,j,linew).le.min_hi(ng)) THEN
              enthalpy(Istr-1,j,linew)=0.0_r8
              ti(Istr-1,j,linew)=0.0_r8
            END IF
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (LBC(iwest,ibICE(isTice),ng)%gradient) THEN
          DO j=Jstr,Jend
            enthalpy(Istr-1,j,linew)=hi(Istr,j,linew)*                  &
     &                               ti(Istr,j,linew)
# ifdef MASKING
            enthalpy(Istr-1,j,linew)=enthalpy(Istr-1,j,linew)*          &
     &                               GRID(ng)%rmask(Istr-1,j)
# endif
# ifdef WET_DRY
            enthalpy(Istr-1,j,linew)=enthalpy(Istr-1,j,linew)*          &
     &                               GRID(ng)%rmask_wet(Istr-1,j)
# endif
            ti(Istr-1,j,linew)=enthalpy(Istr-1,j,linew)/                &
     &                         MAX(hi(Istr-1,j,linew),eps)
            IF (hi(Istr-1,j,linew).le.min_hi(ng)) THEN
              enthalpy(Istr-1,j,linew)=0.0_r8
              ti(Istr-1,j,linew)=0.0_r8
            END IF
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (LBC(iwest,ibICE(isTice),ng)%closed) THEN
          DO j=Jstr,Jend
            enthalpy(Istr-1,j,linew)=hi(Istr,j,linew)*                  &
     &                               ti(Istr,j,linew)
!           ti(Istr-1,j,linew)=enthalpy(Istr-1,j,linew)/                &
!    &                         MAX(hi(Istr-1,j,linew),eps)
             ti(Istr-1,j,linew)=ti(Istr,j,linew)
# ifdef MASKING
            enthalpy(Istr-1,j,linew)=enthalpy(Istr-1,j,linew)*          &
     &                               GRID(ng)%rmask(Istr-1,j)
            ti(Istr-1,j,linew)=ti(Istr-1,j,linew)*                      &
     &                         GRID(ng)%rmask(Istr-1,j)
# endif
# ifdef WET_DRY
            enthalpy(Istr-1,j,linew)=enthalpy(Istr-1,j,linew)*          &
     &                               GRID(ng)%rmask_wet(Istr-1,j)
            ti(Istr-1,j,linew)=ti(Istr-1,j,linew)*                      &
     &                         GRID(ng)%rmask_wet(Istr-1,j)
# endif
            IF (hi(Istr-1,j,linew).le.min_hi(ng)) THEN
              enthalpy(Istr-1,j,linew)=0.0_r8
              ti(Istr-1,j,linew)=0.0_r8
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
        IF (LBC(ieast,ibICE(isTice),ng)%clamped) THEN
!
!  Eastern edge, clamped boundary condition.
!
          DO j=Jstr,Jend
            enthalpy(Iend+1,j,linew)=BOUNDARY(ng)%ice_east(isHice,j)*   &
     &                               BOUNDARY(ng)%ice_east(isTice,j)
# ifdef MASKING
            enthalpy(Iend+1,j,linew)=enthalpy(Iend+1,j,linew)*          &
     &                               GRID(ng)%rmask(Iend+1,j)
# endif
# ifdef WET_DRY
            enthalpy(Iend+1,j,linew)=enthalpy(Iend+1,j,linew)*          &
     &                               GRID(ng)%rmask_wet(Iend+1,j)
# endif
            ti(Iend+1,j,linew)=enthalpy(Iend+1,j,linew)/                &
     &                         MAX(hi(Iend+1,j,linew),eps)
            IF (hi(Iend+1,j,linew).le.min_hi(ng)) THEN
              enthalpy(Iend+1,j,linew)=0.0_r8
              ti(Iend+1,j,linew)=0.0_r8
            END IF
          END DO
!
!  Eastern edge, clamped on inflow, gradient on outflow.
!
        ELSE IF (LBC(ieast,ibICE(isTice),ng)%mixed) THEN
          DO j=Jstr,Jend
            IF (ui(Iend+1,j,linew).le.0.0_r8) THEN
              enthalpy(Iend+1,j,linew)=BOUNDARY(ng)%ice_east(isHice,j)* &
     &                                 BOUNDARY(ng)%ice_east(isTice,j)
# ifdef MASKING
              enthalpy(Iend+1,j,linew)=enthalpy(Iend+1,j,linew)*        &
     &                                 GRID(ng)%rmask(Iend+1,j)
# endif
# ifdef WET_DRY
              enthalpy(Iend+1,j,linew)=enthalpy(Iend+1,j,linew)*        &
     &                                 GRID(ng)%rmask_wet(Iend+1,j)
# endif
            ELSE
              enthalpy(Iend+1,j,linew)=hi(Iend,j,liold)*                &
     &                                 ti(Iend,j,liold)
# ifdef MASKING
              enthalpy(Iend+1,j,linew)=enthalpy(Iend+1,j,linew)*        &
     &                                 GRID(ng)%rmask(Iend+1,j)
# endif
# ifdef WET_DRY
              enthalpy(Iend+1,j,linew)=enthalpy(Iend+1,j,linew)*        &
     &                                 GRID(ng)%rmask_wet(Iend+1,j)
# endif
              ti(Iend+1,j,linew)=enthalpy(Iend+1,j,linew)/              &
     &                           MAX(hi(Iend+1,j,linew),eps)
              IF (hi(Iend+1,j,linew).le.min_hi(ng)) THEN
                enthalpy(Iend+1,j,linew)=0.0_r8
                ti(Iend+1,j,linew)=0.0_r8
              END IF
            END IF
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (LBC(ieast,ibICE(isTice),ng)%gradient) THEN
          DO j=Jstr,Jend
            enthalpy(Iend+1,j,linew)=hi(Iend,j,linew)*                  &
     &                               ti(Iend,j,linew)
# ifdef MASKING
            enthalpy(Iend+1,j,linew)=enthalpy(Iend+1,j,linew)*          &
     &                               GRID(ng)%rmask(Iend+1,j)
# endif
# ifdef WET_DRY
            enthalpy(Iend+1,j,linew)=enthalpy(Iend+1,j,linew)*          &
     &                               GRID(ng)%rmask_wet(Iend+1,j)
# endif
            ti(Iend+1,j,linew)=enthalpy(Iend+1,j,linew)/                &
     &                         MAX(hi(Iend+1,j,linew),eps)
            IF (hi(Iend+1,j,linew).le.min_hi(ng)) THEN
              enthalpy(Iend+1,j,linew)=0.0_r8
              ti(Iend+1,j,linew)=0.0_r8
            END IF
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (LBC(ieast,ibICE(isTice),ng)%closed) THEN
          DO j=Jstr,Jend
            enthalpy(Iend+1,j,linew)=hi(Iend,j,linew)*                  &
     &                               ti(Iend,j,linew)
!           ti(Iend+1,j,linew)=enthalpy(Iend+1,j,linew)/                &
!    &                         MAX(hi(Iend+1,j,linew),eps)
            ti(Iend+1,j,linew)=ti(Iend,j,linew)
# ifdef MASKING
            enthalpy(Iend+1,j,linew)=enthalpy(Iend+1,j,linew)*          &
     &                               GRID(ng)%rmask(Iend+1,j)
            ti(Iend+1,j,linew)=ti(Iend+1,j,linew)*                      &
     &                         GRID(ng)%rmask(Iend+1,j)
# endif
# ifdef WET_DRY
            enthalpy(Iend+1,j,linew)=enthalpy(Iend+1,j,linew)*          &
     &                               GRID(ng)%rmask_wet(Iend+1,j)
            ti(Iend+1,j,linew)=ti(Iend+1,j,linew)*                      &
     &                         GRID(ng)%rmask_wet(Iend+1,j)
# endif
            IF (hi(Iend+1,j,linew).le.min_hi(ng)) THEN
              enthalpy(Iend+1,j,linew)=0.0_r8
              ti(Iend+1,j,linew)=0.0_r8
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
        IF (LBC(isouth,ibICE(isTice),ng)%clamped) THEN
!
!  Southern edge, clamped boundary condition.
!
          DO i=Istr,Iend
            enthalpy(i,Jstr-1,linew)=BOUNDARY(ng)%ice_south(isHice,i)*  &
     &                               BOUNDARY(ng)%ice_south(isTice,i)
# ifdef MASKING
            enthalpy(i,Jstr-1,linew)=enthalpy(i,Jstr-1,linew)*          &
     &                               GRID(ng)%rmask(i,Jstr-1)
# endif
# ifdef WET_DRY
            enthalpy(i,Jstr-1,linew)=enthalpy(i,Jstr-1,linew)*          &
     &                               GRID(ng)%rmask_wet(i,Jstr-1)
# endif
            ti(i,Jstr-1,linew)=enthalpy(i,Jstr-1,linew)/                &
     &                         MAX(hi(i,Jstr-1,linew),eps)
            IF (hi(i,Jstr-1,linew).le.min_hi(ng)) THEN
              enthalpy(i,Jstr-1,linew)=0.0_r8
              ti(i,Jstr-1,linew)=0.0_r8
            END IF
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (LBC(isouth,ibICE(isTice),ng)%mixed) THEN
          DO i=Istr,Iend
            IF (vi(i,1,linew).ge.0.0_r8) THEN
              enthalpy(i,Jstr-1,linew)=BOUNDARY(ng)%ice_south(isHice,i)*&
     &                                 BOUNDARY(ng)%ice_south(isTice,i)
# ifdef MASKING
              enthalpy(i,Jstr-1,linew)=enthalpy(i,Jstr-1,linew)*        &
     &                                 GRID(ng)%rmask(i,Jstr-1)
# endif
# ifdef WET_DRY
              enthalpy(i,Jstr-1,linew)=enthalpy(i,Jstr-1,linew)*        &
     &                                 GRID(ng)%rmask_wet(i,Jstr-1)
# endif
            ELSE
              enthalpy(i,Jstr-1,linew)=enthalpy(i,Jstr,liold)
# ifdef MASKING
              enthalpy(i,Jstr-1,linew)=enthalpy(i,Jstr-1,linew)*        &
     &                                 GRID(ng)%rmask(i,Jstr-1)
# endif
# ifdef WET_DRY
              enthalpy(i,Jstr-1,linew)=enthalpy(i,Jstr-1,linew)*        &
     &                                 GRID(ng)%rmask_wet(i,Jstr-1)
# endif
              ti(i,Jstr-1,linew)=enthalpy(i,Jstr-1,linew)/              &
     &                           MAX(hi(i,Jstr-1,linew),eps)
              IF (hi(i,Jstr-1,linew).le.min_hi(ng)) THEN
                enthalpy(i,Jstr-1,linew)=0.0_r8
                ti(i,Jstr-1,linew)=0.0_r8
              END IF
            ENDIF
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (LBC(isouth,ibICE(isTice),ng)%gradient) THEN
          DO i=Istr,Iend
            enthalpy(i,Jstr-1,linew)=hi(i,Jstr,linew)*                  &
     &                               ti(i,Jstr,linew)
# ifdef MASKING
            enthalpy(i,Jstr-1,linew)=enthalpy(i,Jstr-1,linew)*          &
     &                               GRID(ng)%rmask(i,Jstr-1)
# endif
# ifdef WET_DRY
            enthalpy(i,Jstr-1,linew)=enthalpy(i,Jstr-1,linew)*          &
     &                               GRID(ng)%rmask_wet(i,Jstr-1)
# endif
            ti(i,Jstr-1,linew)=enthalpy(i,Jstr-1,linew)/                &
     &                         MAX(hi(i,Jstr-1,linew),eps)
            IF (hi(i,Jstr-1,linew).le.min_hi(ng)) THEN
              enthalpy(i,Jstr-1,linew)=0.0_r8
              ti(i,Jstr-1,linew)=0.0_r8
            END IF
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE IF (LBC(isouth,ibICE(isTice),ng)%closed) THEN
          DO i=Istr,Iend
            enthalpy(i,Jstr-1,linew)=enthalpy(i,Jstr,linew)
!           ti(i,Jstr-1,linew)=enthalpy(i,Jstr-1,linew)/                &
!     &                         MAX(hi(i,Jstr-1,linew),eps)
            ti(i,Jstr-1,linew)=ti(i,Jstr,linew)
# ifdef MASKING
            enthalpy(i,Jstr-1,linew)=enthalpy(i,Jstr-1,linew)*          &
     &                               GRID(ng)%rmask(i,Jstr-1)
            ti(i,Jstr-1,linew)=ti(i,Jstr-1,linew)*                      &
     &                         GRID(ng)%rmask(i,Jstr-1)
# endif
# ifdef WET_DRY
            enthalpy(i,Jstr-1,linew)=enthalpy(i,Jstr-1,linew)*          &
     &                               GRID(ng)%rmask_wet(i,Jstr-1)
            ti(i,Jstr-1,linew)=ti(i,Jstr-1,linew)*                      &
     &                         GRID(ng)%rmask_wet(i,Jstr-1)
# endif
            IF (hi(i,Jstr-1,linew).le.min_hi(ng)) THEN
              enthalpy(i,Jstr-1,linew)=0.0_r8
              ti(i,Jstr-1,linew)=0.0_r8
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
        IF (LBC(inorth,ibICE(isTice),ng)%clamped) THEN
!
!  Northern edge, clamped boundary condition.
!
          DO i=Istr,Iend
            enthalpy(i,Jend+1,linew)=BOUNDARY(ng)%ice_north(isHice,i)*  &
     &                               BOUNDARY(ng)%ice_north(isTice,i)
# ifdef MASKING
            enthalpy(i,Jend+1,linew)=enthalpy(i,Jend+1,linew)*          &
     &                               GRID(ng)%rmask(i,Jend+1)
# endif
# ifdef WET_DRY
            enthalpy(i,Jend+1,linew)=enthalpy(i,Jend+1,linew)*          &
     &                               GRID(ng)%rmask_wet(i,Jend+1)
# endif
            ti(i,Jend+1,linew)=enthalpy(i,Jend+1,linew)/                &
     &                         MAX(hi(i,Jend+1,linew),eps)
            IF (hi(i,Jend+1,linew).le.min_hi(ng)) THEN
              enthalpy(i,Jend+1,linew)=0.0_r8
              ti(i,Jend+1,linew)=0.0_r8
            END IF
          END DO
!
!  Northern edge, clamped on inflow, gradient on outflow.
!
        ELSE IF (LBC(inorth,ibICE(isTice),ng)%mixed) THEN
          DO i=Istr,Iend
            IF (vi(i,Jend+1,linew).le.0.0_r8) THEN
              enthalpy(i,Jend+1,linew)=BOUNDARY(ng)%ice_north(isHice,i)*&
     &                                 BOUNDARY(ng)%ice_north(isTice,i)
# ifdef MASKING
              enthalpy(i,Jend+1,linew)=enthalpy(i,Jend+1,linew)*        &
     &                                 GRID(ng)%rmask(i,Jend+1)
# endif
# ifdef WET_DRY
              enthalpy(i,Jend+1,linew)=enthalpy(i,Jend+1,linew)*        &
     &                                 GRID(ng)%rmask_wet(i,Jend+1)
# endif
            ELSE
              enthalpy(i,Jend+1,linew)=enthalpy(i,Jend,liold)
# ifdef MASKING
              enthalpy(i,Jend+1,linew)=enthalpy(i,Jend+1,linew)*        &
     &                                 GRID(ng)%rmask(i,Jend+1)
# endif
# ifdef WET_DRY
              enthalpy(i,Jend+1,linew)=enthalpy(i,Jend+1,linew)*        &
     &                                 GRID(ng)%rmask_wet(i,Jend+1)
# endif
            ENDIF
            ti(i,Jend+1,linew)=enthalpy(i,Jend+1,linew)/                &
     &                         MAX(hi(i,Jend+1,linew),eps)
            IF (hi(i,Jend+1,linew).le.min_hi(ng)) THEN
              enthalpy(i,Jend+1,linew)=0.0_r8
              ti(i,Jend+1,linew)=0.0_r8
            END IF
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (LBC(inorth,ibICE(isTice),ng)%gradient) THEN
          DO i=Istr,Iend
            enthalpy(i,Jend+1,linew)=hi(i,Jend,linew)*                  &
     &                               ti(i,Jend,linew)
# ifdef MASKING
            enthalpy(i,Jend+1,linew)=enthalpy(i,Jend+1,linew)*          &
     &                               GRID(ng)%rmask(i,Jend+1)
# endif
# ifdef WET_DRY
            enthalpy(i,Jend+1,linew)=enthalpy(i,Jend+1,linew)*          &
     &                               GRID(ng)%rmask_wet(i,Jend+1)
# endif
# ifdef SOGLOBEC
            ti(i,Jend+1,linew)=enthalpy(i,Jend+1,linew)/                &
     &                         MAX(hi(i,Jend,linew),eps)
# else
            ti(i,Jend+1,linew)=enthalpy(i,Jend+1,linew)/                &
     &                         MAX(hi(i,Jend+1,linew),eps)
# endif
            IF (hi(i,Jend+1,linew).le.min_hi(ng)) THEN
              enthalpy(i,Jend+1,linew)=0.0_r8
              ti(i,Jend+1,linew)=0.0_r8
            END IF
          END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (LBC(inorth,ibICE(isTice),ng)%closed) THEN
          DO i=Istr,Iend
            enthalpy(i,Jend+1,linew)=hi(i,Jend,linew)*                  &
     &                               ti(i,Jend,linew)
!           ti(i,Jend+1,linew)=enthalpy(i,Jend+1,linew)/                &
!     &                        MAX(hi(i,Jend+1,linew),eps)
            ti(i,Jend+1,linew)=ti(i,Jend,linew)
# ifdef MASKING
            enthalpy(i,Jend+1,linew)=enthalpy(i,Jend+1,linew)*          &
     &                               GRID(ng)%rmask(i,Jend+1)
            ti(i,Jend+1,linew)=ti(i,Jend+1,linew)*                      &
     &                         GRID(ng)%rmask(i,Jend+1)
# endif
# ifdef WET_DRY
            enthalpy(i,Jend+1,linew)=enthalpy(i,Jend+1,linew)*          &
     &                               GRID(ng)%rmask_wet(i,Jend+1)
            ti(i,Jend+1,linew)=ti(i,Jend+1,linew)*                      &
     &                         GRID(ng)%rmask_wet(i,Jend+1)
# endif
            IF (hi(i,Jend+1,linew).le.min_hi(ng)) THEN
              enthalpy(i,Jend+1,linew)=0.0_r8
              ti(i,Jend+1,linew)=0.0_r8
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
          enthalpy(Istr-1,Jstr-1,linew)=0.5_r8*                         &
     &                                  (enthalpy(Istr  ,Jstr-1,linew)+ &
     &                                   enthalpy(Istr-1,Jstr  ,linew))
!         ti(Istr-1,Jstr-1,linew)=enthalpy(Istr-1,Jstr-1,linew)/        &
!    &                            MAX(hi(Istr-1,Jstr-1,linew),eps)
          ti(Istr-1,Jstr-1,linew)=0.5_r8*(ti(Istr  ,Jstr-1,linew)+      &
     &                                    ti(Istr-1,Jstr  ,linew))
# ifdef MASKING
          enthalpy(Istr-1,Jstr-1,linew)=enthalpy(Istr-1,Jstr-1,linew)*  &
     &                                  GRID(ng)%rmask(Istr-1,Jstr-1)
          ti(Istr-1,Jstr-1,linew)=ti(Istr-1,Jstr-1,linew)*              &
     &                            GRID(ng)%rmask(Istr-1,Jstr-1)
# endif
# ifdef WET_DRY
          enthalpy(Istr-1,Jstr-1,linew)=enthalpy(Istr-1,Jstr-1,linew)*  &
     &                                 GRID(ng)%rmask_wet(Istr-1,Jstr-1)
          ti(Istr-1,Jstr-1,linew)=ti(Istr-1,Jstr-1,linew)*              &
     &                            GRID(ng)%rmask_wet(Istr-1,Jstr-1)
# endif
          IF (hi(Istr-1,Jstr-1,linew).le.min_hi(ng)) THEN
            enthalpy(Istr-1,Jstr-1,linew)=0.0_r8
              ti(Istr-1,Jstr-1,linew)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          enthalpy(Iend+1,Jstr-1,linew)=0.5_r8*                         &
     &                                  (enthalpy(Iend+1,Jstr  ,linew)+ &
     &                                   enthalpy(Iend  ,Jstr-1,linew))
!         ti(Iend+1,Jstr-1,linew)=enthalpy(Iend+1,Jstr-1,linew)/        &
!    &                            MAX(hi(Iend+1,Jstr-1,linew),eps)
          ti(Iend+1,Jstr-1,linew)=0.5_r8*(ti(Iend+1,Jstr  ,linew)+      &
     &                                    ti(Iend  ,Jstr-1,linew))
# ifdef MASKING
          enthalpy(Iend+1,Jstr-1,linew)=enthalpy(Iend+1,Jstr-1,linew)*  &
     &                                  GRID(ng)%rmask(Iend+1,Jstr-1)
          ti(Iend+1,Jstr-1,linew)=ti(Iend+1,Jstr-1,linew)*              &
     &                            GRID(ng)%rmask(Iend+1,Jstr-1)
# endif
# ifdef WET_DRY
          enthalpy(Iend+1,Jstr-1,linew)=enthalpy(Iend+1,Jstr-1,linew)*  &
     &                                 GRID(ng)%rmask_wet(Iend+1,Jstr-1)
          ti(Iend+1,Jstr-1,linew)=ti(Iend+1,Jstr-1,linew)*              &
     &                            GRID(ng)%rmask_wet(Iend+1,Jstr-1)
# endif
          IF (hi(Iend+1,Jstr-1,linew).LE.min_hi(ng)) THEN
            enthalpy(Iend+1,Jstr-1,linew)=0.0_r8
            ti(Iend+1,Jstr-1,linew)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          enthalpy(Istr-1,Jend+1,linew)=0.5_r8*                         &
     &                                  (enthalpy(Istr-1,Jend  ,linew)+ &
     &                                   enthalpy(Istr  ,Jend+1,linew))
!         ti(Istr-1,Jend+1,linew)=enthalpy(Istr-1,Jend+1,linew)/        &
!    &                            MAX(hi(Istr-1,Jend+1,linew),eps)
          ti(Istr-1,Jend+1,linew)=0.5_r8*(ti(Istr-1,Jend  ,linew)+      &
     &                                    ti(Istr  ,Jend+1,linew))
# ifdef MASKING
          enthalpy(Istr-1,Jend+1,linew)=enthalpy(Istr-1,Jend+1,linew)*  &
     &                                  GRID(ng)%rmask(Istr-1,Jend+1)
          ti(Istr-1,Jend+1,linew)=ti(Istr-1,Jend+1,linew)*              &
     &                            GRID(ng)%rmask(Istr-1,Jend+1)
# endif
# ifdef WET_DRY
          enthalpy(Istr-1,Jend+1,linew)=enthalpy(Istr-1,Jend+1,linew)*  &
     &                                 GRID(ng)%rmask_wet(Istr-1,Jend+1)
          ti(Istr-1,Jend+1,linew)=ti(Istr-1,Jend+1,linew)*              &
     &                            GRID(ng)%rmask_wet(Istr-1,Jend+1)
# endif
          IF (hi(Istr-1,Jend+1,linew).le.min_hi(ng)) THEN
            enthalpy(Istr-1,Jend+1,linew)=0.0_r8
            ti(Istr-1,Jend+1,linew)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          enthalpy(Iend+1,Jend+1,linew)=0.5_r8*                         &
     &                                  (enthalpy(Iend+1,Jend  ,linew)+ &
     &                                   enthalpy(Iend  ,Jend+1,linew))
!         ti(Iend+1,Jend+1,linew)=enthalpy(Iend+1,Jend+1,linew)/        &
!    &                            MAX(hi(Iend+1,Jend+1,linew),eps)
          ti(Iend+1,Jend+1,linew)=0.5_r8*(ti(Iend+1,Jend  ,linew)+      &
     &                                    ti(Iend  ,Jend+1,linew))
# ifdef MASKING
          enthalpy(Iend+1,Jend+1,linew)=enthalpy(Iend+1,Jend+1,linew)*  &
     &                                  GRID(ng)%rmask(Iend+1,Jend+1)
          ti(Iend+1,Jend+1,linew)=ti(Iend+1,Jend+1,linew)*              &
     &                            GRID(ng)%rmask(Iend+1,Jend+1)
# endif
# ifdef WET_DRY
          enthalpy(Iend+1,Jend+1,linew)=enthalpy(Iend+1,Jend+1,linew)*  &
     &                                 GRID(ng)%rmask_wet(Iend+1,Jend+1)
          ti(Iend+1,Jend+1,linew)=ti(Iend+1,Jend+1,linew)*              &
     &                            GRID(ng)%rmask_wet(Iend+1,Jend+1)
# endif
          IF (hi(Iend+1,Jend+1,linew).le.min_hi(ng)) THEN
            enthalpy(Iend+1,Jend+1,linew)=0.0_r8
            ti(Iend+1,Jend+1,linew)=0.0_r8
          END IF
        END IF
      END IF
!
      RETURN
      END SUBROUTINE ice_tibc_tile
#endif
      END MODULE ice_tibc_mod
