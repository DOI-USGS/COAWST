#include "cppdefs.h"
      MODULE v3dbc_mod
#ifdef SOLVE3D
!
!svn $Id: v3dbc_im.F 732 2008-09-07 01:55:51Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for total 3D       !
!  V-velocity.                                                         !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: v3dbc_tile
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE v3dbc (ng, tile, nout)
!***********************************************************************
!
      USE mod_param
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
      CALL v3dbc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 nstp(ng), nout,                                  &
     &                 OCEAN(ng) % v)

      RETURN
      END SUBROUTINE v3dbc
!
!***********************************************************************
      SUBROUTINE v3dbc_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, UBk,                   &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nstp, nout,                                &
     &                       v)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_clima
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nout
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,UBk,2)
# endif
!
!  Local variable declarations.
!
      integer :: Jmin, Jmax
      integer :: i, j, k

      real(r8), parameter :: eps = 1.0E-20_r8

      real(r8) :: Ce, Cx, cff, dVde, dVdt, dVdx
      real(r8) :: obc_in, obc_out, tau

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
        IF (LBC(isouth,isVvel,ng)%radiation) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend+1
              grad(i,Jstr  )=v(i  ,Jstr  ,k,nstp)-                      &
     &                       v(i-1,Jstr  ,k,nstp)
              grad(i,Jstr+1)=v(i  ,Jstr+1,k,nstp)-                      &
     &                       v(i-1,Jstr+1,k,nstp)
            END DO
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                dVdt=v(i,Jstr+1,k,nstp)-v(i,Jstr+1,k,nout)
                dVde=v(i,Jstr+1,k,nout)-v(i,Jstr+2,k,nout)

                IF (LBC(isouth,isVvel,ng)%nudging) THEN
                  IF (LnudgeM3CLM(ng)) THEN
                    obc_out=0.5_r8*                                     &
     &                      (CLIMA(ng)%M3nudgcof(i,Jstr-1,k)+           &
     &                       CLIMA(ng)%M3nudgcof(i,Jstr  ,k))
                    obc_in =obcfac(ng)*obc_out
                  ELSE
                    obc_out=M3obc_out(ng,isouth)
                    obc_in =M3obc_in (ng,isouth)
                  END IF
                  IF ((dVdt*dVde).lt.0.0_r8) THEN
                    tau=obc_in
                  ELSE
                    tau=obc_out
                  END IF
                  tau=tau*dt(ng)
                END IF

                IF ((dVdt*dVde).lt.0.0_r8) dVdt=0.0_r8
                IF ((dVdt*(grad(i  ,Jstr+1)+                            &
     &                     grad(i+1,Jstr+1))).gt.0.0_r8) THEN
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
# if defined CELERITY_WRITE && defined FORWARD_WRITE
                BOUNDARY(ng)%v_south_Cx(i,k)=Cx
                BOUNDARY(ng)%v_south_Ce(i,k)=Ce
                BOUNDARY(ng)%v_south_C2(i,k)=cff
# endif
                v(i,Jstr,k,nout)=(cff*v(i,Jstr  ,k,nstp)+               &
     &                            Ce *v(i,Jstr+1,k,nout)-               &
     &                            MAX(Cx,0.0_r8)*grad(i  ,Jstr)-        &
     &                            MIN(Cx,0.0_r8)*grad(i+1,Jstr))/       &
     &                           (cff+Ce)

                IF (LBC(isouth,isVvel,ng)%nudging) THEN
                  v(i,Jstr,k,nout)=v(i,Jstr,k,nout)+                    &
     &                             tau*(BOUNDARY(ng)%v_south(i,k)-      &
     &                                  v(i,Jstr,k,nstp))
                END IF
# ifdef MASKING
                v(i,Jstr,k,nout)=v(i,Jstr,k,nout)*                      &
     &                           GRID(ng)%vmask(i,Jstr)
# endif
# ifdef WET_DRY
                v(i,Jstr,k,nout)=v(i,Jstr,k,nout)*                      &
     &                           GRID(ng)%vmask_wet(i,Jstr)
# endif
              END IF
            END DO
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (LBC(isouth,isVvel,ng)%clamped) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                v(i,Jstr,k,nout)=BOUNDARY(ng)%v_south(i,k)
# ifdef MASKING
                v(i,Jstr,k,nout)=v(i,Jstr,k,nout)*                      &
     &                           GRID(ng)%vmask(i,Jstr)
# endif
# ifdef WET_DRY
                v(i,Jstr,k,nout)=v(i,Jstr,k,nout)*                      &
     &                           GRID(ng)%vmask_wet(i,Jstr)
# endif
              END IF
            END DO
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (LBC(isouth,isVvel,ng)%gradient) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                v(i,Jstr,k,nout)=v(i,Jstr+1,k,nout)
# ifdef MASKING
                v(i,Jstr,k,nout)=v(i,Jstr,k,nout)*                      &
     &                           GRID(ng)%vmask(i,Jstr)
# endif
# ifdef WET_DRY
                v(i,Jstr,k,nout)=v(i,Jstr,k,nout)*                      &
     &                           GRID(ng)%vmask_wet(i,Jstr)
# endif
              END IF
            END DO
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE IF (LBC(isouth,isVvel,ng)%closed) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                v(i,Jstr,k,nout)=0.0_r8
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
        IF (LBC(inorth,isVvel,ng)%radiation) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend+1
              grad(i,Jend  )=v(i  ,Jend  ,k,nstp)-                      &
     &                       v(i-1,Jend  ,k,nstp)
              grad(i,Jend+1)=v(i  ,Jend+1,k,nstp)-                      &
     &                       v(i-1,Jend+1,k,nstp)
            END DO
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                dVdt=v(i,Jend,k,nstp)-v(i,Jend  ,k,nout)
                dVde=v(i,Jend,k,nout)-v(i,Jend-1,k,nout)

                IF (LBC(inorth,isVvel,ng)%nudging) THEN
                  IF (LnudgeM3CLM(ng)) THEN
                    obc_out=0.5_r8*                                     &
     &                      (CLIMA(ng)%M3nudgcof(i,Jend  ,k)+           &
     &                       CLIMA(ng)%M3nudgcof(i,Jend+1,k))
                    obc_in =obcfac(ng)*obc_out
                  ELSE
                    obc_out=M3obc_out(ng,inorth)
                    obc_in =M3obc_in (ng,inorth)
                  END IF
                  IF ((dVdt*dVde).lt.0.0_r8) THEN
                    tau=obc_in
                  ELSE
                    tau=obc_out
                  END IF
                  tau=tau*dt(ng)
                END IF

                IF ((dVdt*dVde).lt.0.0_r8) dVdt=0.0_r8
                IF ((dVdt*(grad(i  ,Jend)+                              &
     &                     grad(i+1,Jend))).gt.0.0_r8) THEN
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
# if defined CELERITY_WRITE && defined FORWARD_WRITE
                BOUNDARY(ng)%v_north_Cx(i,k)=Cx
                BOUNDARY(ng)%v_north_Ce(i,k)=Ce
                BOUNDARY(ng)%v_north_C2(i,k)=cff
# endif
                v(i,Jend+1,k,nout)=(cff*v(i,Jend+1,k,nstp)+             &
     &                              Ce *v(i,Jend  ,k,nout)-             &
     &                              MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-    &
     &                              MIN(Cx,0.0_r8)*grad(i+1,Jend+1))/   &
     &                             (cff+Ce)

                IF (LBC(inorth,isVvel,ng)%nudging) THEN
                  v(i,Jend+1,k,nout)=v(i,Jend+1,k,nout)+                &
     &                               tau*(BOUNDARY(ng)%v_north(i,k)-    &
     &                                    v(i,Jend+1,k,nstp))
                END IF
# ifdef MASKING
                v(i,Jend+1,k,nout)=v(i,Jend+1,k,nout)*                  &
     &                             GRID(ng)%vmask(i,Jend+1)
# endif
# ifdef WET_DRY
                v(i,Jend+1,k,nout)=v(i,Jend+1,k,nout)*                  &
     &                             GRID(ng)%vmask_wet(i,Jend+1)
# endif
              END IF
            END DO
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (LBC(inorth,isVvel,ng)%clamped) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                v(i,Jend+1,k,nout)=BOUNDARY(ng)%v_north(i,k)
# ifdef MASKING
                v(i,Jend+1,k,nout)=v(i,Jend+1,k,nout)*                  &
     &                             GRID(ng)%vmask(i,Jend+1)
# endif
# ifdef WET_DRY
                v(i,Jend+1,k,nout)=v(i,Jend+1,k,nout)*                  &
     &                             GRID(ng)%vmask_wet(i,Jend+1)
# endif
              END IF
            END DO
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (LBC(inorth,isVvel,ng)%gradient) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                v(i,Jend+1,k,nout)=v(i,Jend,k,nout)
# ifdef MASKING
                v(i,Jend+1,k,nout)=v(i,Jend+1,k,nout)*                  &
     &                             GRID(ng)%vmask(i,Jend+1)
# endif
# ifdef WET_DRY
                v(i,Jend+1,k,nout)=v(i,Jend+1,k,nout)*                  &
     &                             GRID(ng)%vmask_wet(i,Jend+1)
# endif
              END IF
            END DO
          END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (LBC(inorth,isVvel,ng)%closed) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                v(i,Jend+1,k,nout)=0.0_r8
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
        IF (LBC(iwest,isVvel,ng)%radiation) THEN
          DO k=1,N(ng)
            DO j=JstrV-1,Jend
              grad(Istr-1,j)=v(Istr-1,j+1,k,nstp)-                      &
     &                       v(Istr-1,j  ,k,nstp)
              grad(Istr  ,j)=v(Istr  ,j+1,k,nstp)-                      &
     &                       v(Istr  ,j  ,k,nstp)
            END DO
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                dVdt=v(Istr,j,k,nstp)-v(Istr  ,j,k,nout)
                dVdx=v(Istr,j,k,nout)-v(Istr+1,j,k,nout)

                IF (LBC(iwest,isVvel,ng)%nudging) THEN
                  IF (LnudgeM3CLM(ng)) THEN
                    obc_out=0.5_r8*                                     &
     &                      (CLIMA(ng)%M3nudgcof(Istr-1,j-1,k)+         &
     &                       CLIMA(ng)%M3nudgcof(Istr-1,j  ,k))
                    obc_in =obcfac(ng)*obc_out
                  ELSE
                    obc_out=M3obc_out(ng,iwest)
                    obc_in =M3obc_in (ng,iwest)
                  END IF
                  IF ((dVdt*dVdx).lt.0.0_r8) THEN
                    tau=obc_in
                  ELSE
                    tau=obc_out
                  END IF
                  tau=tau*dt(ng)
                END IF

                IF ((dVdt*dVdx).lt.0.0_r8) dVdt=0.0_r8
                IF ((dVdt*(grad(Istr,j-1)+                              &
     &                     grad(Istr,j  ))).gt.0.0_r8) THEN
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
# if defined CELERITY_WRITE && defined FORWARD_WRITE
                BOUNDARY(ng)%v_west_Cx(j,k)=Cx
                BOUNDARY(ng)%v_west_Ce(j,k)=Ce
                BOUNDARY(ng)%v_west_C2(j,k)=cff
# endif
                v(Istr-1,j,k,nout)=(cff*v(Istr-1,j,k,nstp)+             &
     &                              Cx *v(Istr  ,j,k,nout)-             &
     &                              MAX(Ce,0.0_r8)*grad(Istr-1,j-1)-    &
     &                              MIN(Ce,0.0_r8)*grad(Istr-1,j  ))/   &
     &                             (cff+Cx)

                IF (LBC(iwest,isVvel,ng)%nudging) THEN
                  v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)+                &
     &                               tau*(BOUNDARY(ng)%v_west(j,k)-     &
     &                                    v(Istr-1,j,k,nstp))
                END IF
# ifdef MASKING
                v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                  &
     &                             GRID(ng)%vmask(Istr-1,j)
# endif
# ifdef WET_DRY
                v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                  &
     &                             GRID(ng)%vmask_wet(Istr-1,j)
# endif
              END IF
            END DO
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (LBC(iwest,isVvel,ng)%clamped) THEN
          DO k=1,N(ng)
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                v(Istr-1,j,k,nout)=BOUNDARY(ng)%v_west(j,k)
# ifdef MASKING
                v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                  &
     &                             GRID(ng)%vmask(Istr-1,j)
# endif
# ifdef WET_DRY
                v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                  &
     &                             GRID(ng)%vmask_wet(Istr-1,j)
#   endif
              END IF
            END DO
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (LBC(iwest,isVvel,ng)%gradient) THEN
          DO k=1,N(ng)
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                v(Istr-1,j,k,nout)=v(Istr,j,k,nout)
# ifdef MASKING
                v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                  &
     &                             GRID(ng)%vmask(Istr-1,j)
# endif
# ifdef WET_DRY
                v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                  &
     &                             GRID(ng)%vmask_wet(Istr-1,j)
# endif
              END IF
            END DO
          END DO
!
!  Western edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
        ELSE IF (LBC(iwest,isVvel,ng)%closed) THEN
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
                v(Istr-1,j,k,nout)=gamma2(ng)*v(Istr,j,k,nout)
# ifdef MASKING
                v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                  &
     &                             GRID(ng)%vmask(Istr-1,j)
# endif
# ifdef WET_DRY
                v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                  &
     &                             GRID(ng)%vmask_wet(Istr-1,j)
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
        IF (LBC(ieast,isVvel,ng)%radiation) THEN
          DO k=1,N(ng)
            DO j=JstrV-1,Jend
              grad(Iend  ,j)=v(Iend  ,j+1,k,nstp)-                      &
     &                       v(Iend  ,j  ,k,nstp)
              grad(Iend+1,j)=v(Iend+1,j+1,k,nstp)-                      &
     &                       v(Iend+1,j  ,k,nstp)
            END DO
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                dVdt=v(Iend,j,k,nstp)-v(Iend  ,j,k,nout)
                dVdx=v(Iend,j,k,nout)-v(Iend-1,j,k,nout)

                IF (LBC(ieast,isVvel,ng)%nudging) THEN
                  IF (LnudgeM3CLM(ng)) THEN
                    obc_out=0.5_r8*                                     &
     &                      (CLIMA(ng)%M3nudgcof(Iend+1,j-1,k)+         &
     &                       CLIMA(ng)%M3nudgcof(Iend+1,j  ,k))
                    obc_in =obcfac(ng)*obc_out
                  ELSE
                    obc_out=M3obc_out(ng,ieast)
                    obc_in =M3obc_in (ng,ieast)
                  END IF
                  IF ((dVdt*dVdx).lt.0.0_r8) THEN
                    tau=obc_in
                  ELSE
                    tau=obc_out
                  END IF
                  tau=tau*dt(ng)
                END IF

                IF ((dVdt*dVdx).lt.0.0_r8) dVdt=0.0_r8
                IF ((dVdt*(grad(Iend,j-1)+                              &
     &                     grad(Iend,j  ))).gt.0.0_r8) THEN
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
# if defined CELERITY_WRITE && defined FORWARD_WRITE
                BOUNDARY(ng)%v_east_Cx(j,k)=Cx
                BOUNDARY(ng)%v_east_Ce(j,k)=Ce
                BOUNDARY(ng)%v_east_C2(j,k)=cff
# endif
                v(Iend+1,j,k,nout)=(cff*v(Iend+1,j,k,nstp)+             &
     &                              Cx *v(Iend  ,j,k,nout)-             &
     &                              MAX(Ce,0.0_r8)*grad(Iend+1,j-1)-    &
     &                              MIN(Ce,0.0_r8)*grad(Iend+1,j  ))/   &
     &                             (cff+Cx)

                IF (LBC(ieast,isVvel,ng)%nudging) THEN
                  v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)+                &
     &                               tau*(BOUNDARY(ng)%v_east(j,k)-     &
     &                                    v(Iend+1,j,k,nstp))
                END IF
# ifdef MASKING
                v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                  &
     &                             GRID(ng)%vmask(Iend+1,j)
# endif
# ifdef WET_DRY
                v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                  &
     &                             GRID(ng)%vmask_wet(Iend+1,j)
# endif
              END IF
            END DO
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (LBC(ieast,isVvel,ng)%clamped) THEN
          DO k=1,N(ng)
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                v(Iend+1,j,k,nout)=BOUNDARY(ng)%v_east(j,k)
# ifdef MASKING
                v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                  &
     &                             GRID(ng)%vmask(Iend+1,j)
# endif
# ifdef WET_DRY
                v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                  &
     &                             GRID(ng)%vmask_wet(Iend+1,j)
# endif
              END IF
            END DO
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (LBC(ieast,isVvel,ng)%gradient) THEN
          DO k=1,N(ng)
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                v(Iend+1,j,k,nout)=v(Iend,j,k,nout)
# ifdef MASKING
                v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                  &
     &                             GRID(ng)%vmask(Iend+1,j)
# endif
# ifdef WET_DRY
                v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                  &
     &                             GRID(ng)%vmask_wet(Iend+1,j)
# endif
              END IF
            END DO
          END DO
!
!  Eastern edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
        ELSE IF (LBC(ieast,isVvel,ng)%closed) THEN
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
                v(Iend+1,j,k,nout)=gamma2(ng)*v(Iend,j,k,nout)
# ifdef MASKING
                v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                  &
     &                             GRID(ng)%vmask(Iend+1,j)
# endif
# ifdef WET_DRY
                v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                  &
     &                             GRID(ng)%vmask_wet(Iend+1,j)
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
              v(Istr-1,Jstr,k,nout)=0.5_r8*(v(Istr  ,Jstr  ,k,nout)+    &
     &                                      v(Istr-1,Jstr+1,k,nout))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr  )) THEN
            DO k=1,N(ng)
              v(Iend+1,Jstr,k,nout)=0.5_r8*(v(Iend  ,Jstr  ,k,nout)+    &
     &                                      v(Iend+1,Jstr+1,k,nout))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO k=1,N(ng)
              v(Istr-1,Jend+1,k,nout)=0.5_r8*(v(Istr-1,Jend  ,k,nout)+  &
     &                                        v(Istr  ,Jend+1,k,nout))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO k=1,N(ng)
              v(Iend+1,Jend+1,k,nout)=0.5_r8*(v(Iend+1,Jend  ,k,nout)+  &
     &                                        v(Iend  ,Jend+1,k,nout))
            END DO
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE v3dbc_tile
#endif
      END MODULE v3dbc_mod
