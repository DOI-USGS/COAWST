#include "cppdefs.h"
      MODULE u3dbc_mod
#ifdef SOLVE3D
!
!svn $Id: u3dbc_ex.F 795 2016-05-11 01:42:43Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for total 3D       !
!  U-velocity.                                                         !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: u3dbc_tile
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE u3dbc (ng, tile, nout)
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
      CALL u3dbc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 nstp(ng), nout,                                  &
     &                 OCEAN(ng) % u)

      RETURN
      END SUBROUTINE u3dbc
!
!***********************************************************************
      SUBROUTINE u3dbc_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, UBk,                   &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nstp, nout,                                &
     &                       u)
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
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,UBk,2)
# endif
!
!  Local variable declarations.
!
      integer :: Imin, Imax
      integer :: i, j, k

      real(r8), parameter :: eps = 1.0E-20_r8

      real(r8) :: Ce, Cx, cff, dUde, dUdt, dUdx
      real(r8) :: obc_in, obc_out, tau

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad

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
        IF (LBC(iwest,isUvel,ng)%radiation) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend+1
              grad(Istr  ,j)=u(Istr  ,j  ,k,nstp)-                      &
     &                       u(Istr  ,j-1,k,nstp)
              grad(Istr+1,j)=u(Istr+1,j  ,k,nstp)-                      &
     &                       u(Istr+1,j-1,k,nstp)
            END DO
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                dUdt=u(Istr+1,j,k,nstp)-u(Istr+1,j,k,nout)
                dUdx=u(Istr+1,j,k,nstp)-u(Istr+2,j,k,nstp)

                IF (LBC(iwest,isUvel,ng)%nudging) THEN
                  IF (LnudgeM3CLM(ng)) THEN
                    obc_out=0.5_r8*                                     &
     &                      (CLIMA(ng)%M3nudgcof(Istr-1,j,k)+           &
     &                       CLIMA(ng)%M3nudgcof(Istr  ,j,k))
                    obc_in =obcfac(ng)*obc_out
                  ELSE
                    obc_out=M3obc_out(ng,iwest)
                    obc_in =M3obc_in (ng,iwest)
                  END IF
                  IF ((dUdt*dUdx).lt.0.0_r8) THEN
                    tau=obc_in
                  ELSE
                    tau=obc_out
                  END IF
                  tau=tau*dt(ng)
                END IF

                IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
                IF ((dUdt*(grad(Istr+1,j  )+                            &
     &                     grad(Istr+1,j+1))).gt.0.0_r8) THEN
                  dUde=grad(Istr+1,j  )
                ELSE
                  dUde=grad(Istr+1,j+1)
                END IF
                cff=dUdt/MAX(dUdx*dUdx+dUde*dUde,eps)
                Cx=MIN(1.0_r8,cff*dUdx)
# ifdef RADIATION_2D
                Ce=MIN(1.0_r8,MAX(cff*dUde,-1.0_r8))
# else
                Ce=0.0_r8
# endif
# if defined CELERITY_WRITE && defined FORWARD_WRITE
                BOUNDARY(ng)%u_west_Cx(j,k)=Cx
                BOUNDARY(ng)%u_west_Ce(j,k)=Ce
                BOUNDARY(ng)%u_west_C2(j,k)=cff
# endif
                u(Istr,j,k,nout)=(1.0_r8-Cx)*u(Istr,j,k,nstp)+          &
     &                           Cx*u(Istr+1,j,k,nstp)-                 &
     &                           MAX(Ce,0.0_r8)*grad(Istr,j  )-         &
     &                           MIN(Ce,0.0_r8)*grad(Istr,j+1)

                IF (LBC(iwest,isUvel,ng)%nudging) THEN
                  u(Istr,j,k,nout)=u(Istr,j,k,nout)+                    &
     &                             tau*(BOUNDARY(ng)%u_west(j,k)-       &
     &                                  u(Istr,j,k,nstp))
                END IF
# ifdef MASKING
                u(Istr,j,k,nout)=u(Istr,j,k,nout)*                      &
     &                           GRID(ng)%umask(Istr,j)
# endif
# ifdef WET_DRY
                u(Istr,j,k,nout)=u(Istr,j,k,nout)*                      &
     &                           GRID(ng)%umask_wet(Istr,j)
# endif
              END IF
            END DO
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (LBC(iwest,isUvel,ng)%clamped) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                u(Istr,j,k,nout)=BOUNDARY(ng)%u_west(j,k)
# ifdef MASKING
                u(Istr,j,k,nout)=u(Istr,j,k,nout)*                      &
     &                           GRID(ng)%umask(Istr,j)
# endif
# ifdef WET_DRY
                u(Istr,j,k,nout)=u(Istr,j,k,nout)*                      &
     &                           GRID(ng)%umask_wet(Istr,j)
# endif
              END IF
            END DO
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (LBC(iwest,isUvel,ng)%gradient) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                u(Istr,j,k,nout)=u(Istr+1,j,k,nout)
# ifdef MASKING
                u(Istr,j,k,nout)=u(Istr,j,k,nout)*                      &
     &                           GRID(ng)%umask(Istr,j)
# endif
# ifdef WET_DRY
                u(Istr,j,k,nout)=u(Istr,j,k,nout)*                      &
     &                           GRID(ng)%umask_wet(Istr,j)
# endif
              END IF
            END DO
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (LBC(iwest,isUvel,ng)%closed) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                u(Istr,j,k,nout)=0.0_r8
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
        IF (LBC(ieast,isUvel,ng)%radiation) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend+1
              grad(Iend  ,j)=u(Iend  ,j  ,k,nstp)-                      &
     &                       u(Iend  ,j-1,k,nstp)
              grad(Iend+1,j)=u(Iend+1,j  ,k,nstp)-                      &
     &                       u(Iend+1,j-1,k,nstp)
            END DO
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                dUdt=u(Iend,j,k,nstp)-u(Iend  ,j,k,nout)
                dUdx=u(Iend,j,k,nstp)-u(Iend-1,j,k,nstp)

                IF (LBC(ieast,isUvel,ng)%nudging) THEN
                  IF (LnudgeM3CLM(ng)) THEN
                    obc_out=0.5_r8*                                     &
     &                      (CLIMA(ng)%M3nudgcof(Iend  ,j,k)+           &
     &                       CLIMA(ng)%M3nudgcof(Iend+1,j,k))
                    obc_in =obcfac(ng)*obc_out
                  ELSE
                    obc_out=M3obc_out(ng,ieast)
                    obc_in =M3obc_in (ng,ieast)
                  END IF
                  IF ((dUdt*dUdx).lt.0.0_r8) THEN
                    tau=obc_in
                  ELSE
                    tau=obc_out
                  END IF
                  tau=tau*dt(ng)
                END IF

                IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
                IF ((dUdt*(grad(Iend,j  )+                              &
     &                     grad(Iend,j+1))).gt.0.0_r8) THEN
                  dUde=grad(Iend,j  )
                ELSE
                  dUde=grad(Iend,j+1)
                END IF
                cff=dUdt/MAX(dUdx*dUdx+dUde*dUde,eps)
                Cx=MIN(1.0_r8,cff*dUdx)
# ifdef RADIATION_2D
                Ce=MIN(1.0_r8,MAX(cff*dUde,-1.0_r8))
# else
                Ce=0.0_r8
# endif
# if defined CELERITY_WRITE && defined FORWARD_WRITE
                BOUNDARY(ng)%u_east_Cx(j,k)=Cx
                BOUNDARY(ng)%u_east_Ce(j,k)=Ce
                BOUNDARY(ng)%u_east_C2(j,k)=cff
# endif
                u(Iend+1,j,k,nout)=(1.0_r8-Cx)*u(Iend+1,j,k,nstp)+      &
     &                             Cx*u(Iend,j,k,nstp)-                 &
     &                             MAX(Ce,0.0_r8)*grad(Iend+1,j  )-     &
     &                             MIN(Ce,0.0_r8)*grad(Iend+1,j+1)

                IF (LBC(ieast,isUvel,ng)%nudging) THEN
                  u(Iend+1,j,k,nout)=u(Iend+1,j,k,nout)+                &
     &                               tau*(BOUNDARY(ng)%u_east(j,k)-     &
     &                                    u(Iend+1,j,k,nstp))
                END IF
# ifdef MASKING
                u(Iend+1,j,k,nout)=u(Iend+1,j,k,nout)*                  &
     &                             GRID(ng)%umask(Iend+1,j)
# endif
# ifdef WET_DRY
                u(Iend+1,j,k,nout)=u(Iend+1,j,k,nout)*                  &
     &                             GRID(ng)%umask_wet(Iend+1,j)
# endif
              END IF
            END DO
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (LBC(ieast,isUvel,ng)%clamped) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                u(Iend+1,j,k,nout)=BOUNDARY(ng)%u_east(j,k)
# ifdef MASKING
                u(Iend+1,j,k,nout)=u(Iend+1,j,k,nout)*                  &
     &                             GRID(ng)%umask(Iend+1,j)
# endif
# ifdef WET_DRY
                u(Iend+1,j,k,nout)=u(Iend+1,j,k,nout)*                  &
     &                             GRID(ng)%umask_wet(Iend+1,j)
# endif
              END IF
            END DO
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (LBC(ieast,isUvel,ng)%gradient) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                u(Iend+1,j,k,nout)=u(Iend,j,k,nout)
# ifdef MASKING
                u(Iend+1,j,k,nout)=u(Iend+1,j,k,nout)*                  &
     &                             GRID(ng)%umask(Iend+1,j)
# endif
# ifdef WET_DRY
                u(Iend+1,j,k,nout)=u(Iend+1,j,k,nout)*                  &
     &                             GRID(ng)%umask_wet(Iend+1,j)
# endif
              END IF
            END DO
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (LBC(ieast,isUvel,ng)%closed) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                u(Iend+1,j,k,nout)=0.0_r8
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
        IF (LBC(isouth,isUvel,ng)%radiation) THEN
          DO k=1,N(ng)
            DO i=IstrU-1,Iend
              grad(i,Jstr-1)=u(i+1,Jstr-1,k,nstp)-                      &
     &                       u(i  ,Jstr-1,k,nstp)
              grad(i,Jstr  )=u(i+1,Jstr  ,k,nstp)-                      &
     &                       u(i  ,Jstr  ,k,nstp)
            END DO
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                dUdt=u(i,Jstr,k,nstp)-u(i,Jstr  ,k,nout)
                dUde=u(i,Jstr,k,nstp)-u(i,Jstr+1,k,nstp)

                IF (LBC(isouth,isUvel,ng)%nudging) THEN
                  IF (LnudgeM3CLM(ng)) THEN
                    obc_out=0.5_r8*                                     &
     &                      (CLIMA(ng)%M3nudgcof(i-1,Jstr-1,k)+         &
     &                       CLIMA(ng)%M3nudgcof(i  ,Jstr-1,k))
                    obc_in =obcfac(ng)*obc_out
                  ELSE
                    obc_out=M3obc_out(ng,isouth)
                    obc_in =M3obc_in (ng,isouth)
                  END IF
                  IF ((dUdt*dUde).lt.0.0_r8) THEN
                    tau=obc_in
                  ELSE
                    tau=obc_out
                  END IF
                  tau=tau*dt(ng)
                END IF

                IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
                IF ((dUdt*(grad(i-1,Jstr)+                              &
     &                     grad(i  ,Jstr))).gt.0.0_r8) THEN
                  dUdx=grad(i-1,Jstr)
                ELSE
                  dUdx=grad(i  ,Jstr)
                END IF
                cff=dUdt/MAX(dUdx*dUdx+dUde*dUde,eps)
# ifdef RADIATION_2D
                Cx=MIN(1.0_r8,MAX(cff*dUdx,-1.0_r8))
# else
                Cx=0.0_r8
# endif
                Ce=MIN(1.0_r8,cff*dUde)
# if defined CELERITY_WRITE && defined FORWARD_WRITE
                BOUNDARY(ng)%u_south_Cx(i,k)=Cx
                BOUNDARY(ng)%u_south_Ce(i,k)=Ce
                BOUNDARY(ng)%u_south_C2(i,k)=cff
# endif
                u(i,Jstr-1,k,nout)=(1.0_r8-Ce)*u(i,Jstr-1,k,nstp)+      &
     &                             Ce*u(i,Jstr,k,nstp)-                 &
     &                             MAX(Cx,0.0_r8)*grad(i-1,Jstr-1)-     &
     &                             MIN(Cx,0.0_r8)*grad(i  ,Jstr-1)

                IF (LBC(isouth,isUvel,ng)%nudging) THEN
                  u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)+                &
     &                               tau*(BOUNDARY(ng)%u_south(i,k)-    &
     &                                    u(i,Jstr-1,k,nstp))
                END IF
# ifdef MASKING
                u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)*                  &
     &                             GRID(ng)%umask(i,Jstr-1)
# endif
# ifdef WET_DRY
                u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)*                  &
     &                             GRID(ng)%umask_wet(i,Jstr-1)
# endif
              END IF
            END DO
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (LBC(isouth,isUvel,ng)%clamped) THEN
          DO k=1,N(ng)
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                u(i,Jstr-1,k,nout)=BOUNDARY(ng)%u_south(i,k)
# ifdef MASKING
                u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)*                  &
     &                             GRID(ng)%umask(i,Jstr-1)
# endif
# ifdef WET_DRY
                u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)*                  &
     &                             GRID(ng)%umask_wet(i,Jstr-1)
# endif
              END IF
            END DO
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (LBC(isouth,isUvel,ng)%gradient) THEN
          DO k=1,N(ng)
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                u(i,Jstr-1,k,nout)=u(i,Jstr,k,nout)
# ifdef MASKING
                u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)*                  &
     &                             GRID(ng)%umask(i,Jstr-1)
# endif
# ifdef WET_MASK
                u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)*                  &
     &                             GRID(ng)%umask_wet(i,Jstr-1)
# endif
              END IF
            END DO
          END DO
!
!  Southern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
        ELSE IF (LBC(isouth,isUvel,ng)%closed) THEN
          IF (EWperiodic(ng)) THEN
            Imin=IstrU
            Imax=Iend
          ELSE
            Imin=Istr
            Imax=IendR
          END IF
          DO k=1,N(ng)
            DO i=Imin,Imax
              IF (LBC_apply(ng)%south(i)) THEN
                u(i,Jstr-1,k,nout)=gamma2(ng)*u(i,Jstr,k,nout)
# ifdef MASKING
                u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)*                  &
     &                             GRID(ng)%umask(i,Jstr-1)
# endif
# ifdef WET_DRY
                u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)*                  &
     &                             GRID(ng)%umask_wet(i,Jstr-1)
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
        IF (LBC(inorth,isUvel,ng)%radiation) THEN
          DO k=1,N(ng)
            DO i=IstrU-1,Iend
              grad(i,Jend  )=u(i+1,Jend  ,k,nstp)-                      &
     &                       u(i  ,Jend  ,k,nstp)
              grad(i,Jend+1)=u(i+1,Jend+1,k,nstp)-                      &
     &                       u(i  ,Jend+1,k,nstp)
            END DO
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                dUdt=u(i,Jend,k,nstp)-u(i,Jend  ,k,nout)
                dUde=u(i,Jend,k,nstp)-u(i,Jend-1,k,nstp)

                IF (LBC(inorth,isUvel,ng)%nudging) THEN
                  IF (LnudgeM3CLM(ng)) THEN
                    obc_out=0.5_r8*                                     &
     &                      (CLIMA(ng)%M3nudgcof(i-1,Jend+1,k)+         &
     &                       CLIMA(ng)%M3nudgcof(i  ,Jend+1,k))
                    obc_in =obcfac(ng)*obc_out
                  ELSE
                    obc_out=M3obc_out(ng,inorth)
                    obc_in =M3obc_in (ng,inorth)
                  END IF
                  IF ((dUdt*dUde).lt.0.0_r8) THEN
                    tau=obc_in
                  ELSE
                    tau=obc_out
                  END IF
                  tau=tau*dt(ng)
                END IF

                IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
                IF ((dUdt*(grad(i-1,Jend)+                              &
     &                     grad(i  ,Jend))).gt.0.0_r8) THEN
                  dUdx=grad(i-1,Jend)
                ELSE
                  dUdx=grad(i  ,Jend)
                END IF
                cff=dUdt/MAX(dUdx*dUdx+dUde*dUde,eps)
# ifdef RADIATION_2D
                Cx=MIN(1.0_r8,MAX(cff*dUdx,-1.0_r8))
# else
                Cx=0.0_r8
# endif
                Ce=MIN(1.0_r8,cff*dUde)
# if defined CELERITY_WRITE && defined FORWARD_WRITE
                BOUNDARY(ng)%u_north_Cx(i,k)=Cx
                BOUNDARY(ng)%u_north_Ce(i,k)=Ce
                BOUNDARY(ng)%u_north_C2(i,k)=cff
# endif
                u(i,Jend+1,k,nout)=(1.0_r8-Ce)*u(i,Jend+1,k,nstp)+      &
     &                             Ce*u(i,Jend,k,nstp)-                 &
     &                             MAX(Cx,0.0_r8)*grad(i-1,Jend+1)-     &
     &                             MIN(Cx,0.0_r8)*grad(i  ,Jend+1)

                IF (LBC(inorth,isUvel,ng)%nudging) THEN
                  u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)+                &
     &                               tau*(BOUNDARY(ng)%u_north(i,k)-    &
     &                                    u(i,Jend+1,k,nstp))
                END IF
# ifdef MASKING
                u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)*                  &
     &                             GRID(ng)%umask(i,Jend+1)
# endif
# ifdef WET_DRY
                u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)*                  &
     &                             GRID(ng)%umask_wet(i,Jend+1)
# endif
              END IF
            END DO
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (LBC(inorth,isUvel,ng)%clamped) THEN
          DO k=1,N(ng)
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                u(i,Jend+1,k,nout)=BOUNDARY(ng)%u_north(i,k)
# ifdef MASKING
                u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)*                  &
     &                             GRID(ng)%umask(i,Jend+1)
# endif
# ifdef WET_DRY
                u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)*                  &
     &                             GRID(ng)%umask_wet(i,Jend+1)
# endif
              END IF
            END DO
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (LBC(inorth,isUvel,ng)%gradient) THEN
          DO k=1,N(ng)
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                u(i,Jend+1,k,nout)=u(i,Jend,k,nout)
# ifdef MASKING
                u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)*                  &
     &                             GRID(ng)%umask(i,Jend+1)
# endif
# ifdef WET_DRY
                u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)*                  &
     &                             GRID(ng)%umask_wet(i,Jend+1)
# endif
              END IF
            END DO
          END DO
!
!  Northern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
        ELSE IF (LBC(inorth,isUvel,ng)%closed) THEN
          IF (EWperiodic(ng)) THEN
            Imin=IstrU
            Imax=Iend
          ELSE
            Imin=Istr
            Imax=IendR
          END IF
          DO k=1,N(ng)
            DO i=Imin,Imax
              IF (LBC_apply(ng)%north(i)) THEN
                u(i,Jend+1,k,nout)=gamma2(ng)*u(i,Jend,k,nout)
# ifdef MASKING
                u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)*                  &
     &                             GRID(ng)%umask(i,Jend+1)
# endif
# ifdef WET_DRY
                u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)*                  &
     &                             GRID(ng)%umask_wet(i,Jend+1)
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
          IF (LBC_apply(ng)%south(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
            DO k=1,N(ng)
              u(Istr,Jstr-1,k,nout)=0.5_r8*(u(Istr+1,Jstr-1,k,nout)+    &
     &                                      u(Istr  ,Jstr  ,k,nout))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            DO k=1,N(ng)
              u(Iend+1,Jstr-1,k,nout)=0.5_r8*(u(Iend  ,Jstr-1,k,nout)+  &
     &                                        u(Iend+1,Jstr  ,k,nout))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO k=1,N(ng)
              u(Istr,Jend+1,k,nout)=0.5_r8*(u(Istr  ,Jend  ,k,nout)+    &
     &                                      u(Istr+1,Jend+1,k,nout))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO k=1,N(ng)
              u(Iend+1,Jend+1,k,nout)=0.5_r8*(u(Iend+1,Jend  ,k,nout)+  &
     &                                        u(Iend  ,Jend+1,k,nout))
            END DO
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE u3dbc_tile
#endif
      END MODULE u3dbc_mod
