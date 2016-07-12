#include "cppdefs.h"
      MODULE rp_u3dbc_mod
#if defined TL_IOMS && defined SOLVE3D
!
!svn $Id: rp_u3dbc_im.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine sets  representers tangent linear lateral boundary  !
!  conditions for total 3D U-velocity. It updates the specified "nout" !
!  time index.                                                         !
!                                                                      !
!  BASIC STATE variables needed: u                                     !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: rp_u3dbc, rp_u3dbc_tile
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE rp_u3dbc (ng, tile, nout)
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
      CALL rp_u3dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp(ng), nout,                               &
     &                    OCEAN(ng) % tl_u)

      RETURN
      END SUBROUTINE rp_u3dbc
!
!***********************************************************************
      SUBROUTINE rp_u3dbc_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, UBk,                &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nstp, nout,                             &
     &                          tl_u)
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
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,UBk,2)
# endif
!
!  Local variable declarations.
!
      integer :: Imin, Imax
      integer :: i, j, k

      real(r8) :: Ce, Cx, cff
      real(r8) :: obc_in, obc_out, tau

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_grad

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
        IF (tl_LBC(iwest,isUvel,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend+1
!>              grad(Istr,j)=u(Istr,j  ,k,nstp)-                        &
!>   &                       u(Istr,j-1,k,nstp)
!>
                tl_grad(Istr,j)=0.0_r8
              END DO
              DO j=Jstr,Jend
                IF (LBC_apply(ng)%west(j)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                  IF (tl_LBC(iwest,isUvel,ng)%nudging) THEN
                    IF (LnudgeM3CLM(ng)) THEN
                      obc_out=0.5_r8*                                   &
     &                        (CLIMA(ng)%M3nudgcof(Istr-1,j,k)+         &
     &                         CLIMA(ng)%M3nudgcof(Istr  ,j,k))
                      obc_in =obcfac(ng)*obc_out
                    ELSE
                      obc_out=M3obc_out(ng,iwest)
                      obc_in =M3obc_in (ng,iwest)
                    END IF
                    IF (BOUNDARY(ng)%u_west_Cx(j,k).lt.0.0_r8) THEN
                      tau=obc_in
                    ELSE
                      tau=obc_out
                    END IF
                    tau=tau*dt(ng)
                  END IF
                  Cx=BOUNDARY(ng)%u_west_Cx(j,k)
#  ifdef RADIATION_2D
                  Ce=BOUNDARY(ng)%u_west_Ce(j,k)
#  else
                  Ce=0.0_r8
#  endif
                  cff=BOUNDARY(ng)%u_west_C2(j,k)
# endif
!>                u(Istr,j,k,nout)=(cff*u(Istr  ,j,k,nstp)+             &
!>   &                              Cx *u(Istr+1,j,k,nout)-             &
!>   &                              MAX(Ce,0.0_r8)*grad(Istr,j  )-      &
!>   &                              MIN(Ce,0.0_r8)*grad(Istr,j+1))/     &
!>   &                             (cff+Cx)
!>
                  tl_u(Istr,j,k,nout)=(cff*tl_u(Istr  ,j,k,nstp)+       &
     &                                 Cx *tl_u(Istr+1,j,k,nout)-       &
     &                                 MAX(Ce,0.0_r8)*                  &
     &                                    tl_grad(Istr,j  )-            &
     &                                 MIN(Ce,0.0_r8)*                  &
     &                                    tl_grad(Istr,j+1))/           &
     &                                (cff+Cx)

                  IF (tl_LBC(iwest,isUvel,ng)%nudging) THEN
!>                  u(Istr,j,k,nout)=u(Istr,j,k,nout)+                  &
!>   &                               tau*(BOUNDARY(ng)%u_west(j,k)-     &
!>   &                                    u(Istr,j,k,nstp))
!>
                    tl_u(Istr,j,k,nout)=tl_u(Istr,j,k,nout)-            &
     &                                  tau*tl_u(Istr,j,k,nstp)
                  END IF
# ifdef MASKING
!>                u(Istr,j,k,nout)=u(Istr,j,k,nout)*                    &
!>   &                             GRID(ng)%umask(Istr,j)
!>
                  tl_u(Istr,j,k,nout)=tl_u(Istr,j,k,nout)*              &
     &                                GRID(ng)%umask(Istr,j)
# endif
                END IF
              END DO
            END DO
          END IF
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (tl_LBC(iwest,isUvel,ng)%clamped) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
!>              u(Istr,j,k,nout)=BOUNDARY(ng)%u_west(j,k)
!>
                tl_u(Istr,j,k,nout)=BOUNDARY(ng)%tl_u_west(j,k)
# ifdef MASKING
!>              u(Istr,j,k,nout)=u(Istr,j,k,nout)*                      &
!>   &                           GRID(ng)%umask(Istr,j)
!>
                tl_u(Istr,j,k,nout)=tl_u(Istr,j,k,nout)*                &
     &                              GRID(ng)%umask(Istr,j)
# endif
              END IF
            END DO
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (tl_LBC(iwest,isUvel,ng)%gradient) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
!>              u(Istr,j,k,nout)=u(Istr+1,j,k,nout)
!>
                tl_u(Istr,j,k,nout)=tl_u(Istr+1,j,k,nout)
# ifdef MASKING
!>              u(Istr,j,k,nout)=u(Istr,j,k,nout)*                      &
!>   &                           GRID(ng)%umask(Istr,j)
!>
                tl_u(Istr,j,k,nout)=tl_u(Istr,j,k,nout)*                &
     &                              GRID(ng)%umask(Istr,j)
# endif
              END IF
            END DO
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (tl_LBC(iwest,isUvel,ng)%closed) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
!>              u(Istr,j,k,nout)=0.0_r8
!>
                tl_u(Istr,j,k,nout)=0.0_r8
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
        IF (tl_LBC(ieast,isUvel,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend+1
!>              grad(Iend+1,j)=u(Iend+1,j  ,k,nstp)-                    &
!>   &                         u(Iend+1,j-1,k,nstp)
!>
                tl_grad(Iend+1,j)=0.0_r8
              END DO
              DO j=Jstr,Jend
                IF (LBC_apply(ng)%east(j)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                  IF (tl_LBC(ieast,isUvel,ng)%nudging) THEN
                    IF (LnudgeM3CLM(ng)) THEN
                      obc_out=0.5_r8*                                   &
     &                        (CLIMA(ng)%M3nudgcof(Iend  ,j,k)+         &
     &                         CLIMA(ng)%M3nudgcof(Iend+1,j,k))
                      obc_in =obcfac(ng)*obc_out
                    ELSE
                      obc_out=M3obc_out(ng,ieast)
                      obc_in =M3obc_in (ng,ieast)
                    END IF
                    IF (BOUNDARY(ng)%u_east_Cx(j,k).lt.0.0_r8) THEN
                      tau=obc_in
                    ELSE
                      tau=obc_out
                    END IF
                    tau=tau*dt(ng)
                  END IF
                  Cx=BOUNDARY(ng)%u_east_Cx(j,k)
#  ifdef RADIATION_2D
                  Ce=BOUNDARY(ng)%u_east_Ce(j,k)
#  else
                  Ce=0.0_r8
#  endif
                  cff=BOUNDARY(ng)%u_east_C2(j,k)
# endif
!>                u(Iend+1,j,k,nout)=(cff*u(Iend+1,j,k,nstp)+           &
!>   &                                Cx *u(Iend  ,j,k,nout)-           &
!>   &                                MAX(Ce,0.0_r8)*grad(Iend+1,j  )-  &
!>   &                                MIN(Ce,0.0_r8)*grad(Iend+1,j+1))/ &
!>   &                               (cff+Cx)
!>
                  tl_u(Iend+1,j,k,nout)=(cff*tl_u(Iend+1,j,k,nstp)+     &
     &                                   Cx *tl_u(Iend  ,j,k,nout)-     &
     &                                   MAX(Ce,0.0_r8)*                &
     &                                      tl_grad(Iend+1,j  )-        &
     &                                   MIN(Ce,0.0_r8)*                &
     &                                      tl_grad(Iend+1,j+1))/       &
     &                                  (cff+Cx)

                  IF (tl_LBC(ieast,isUvel,ng)%nudging) THEN
!>                  u(Iend+1,j,k,nout)=u(Iend+1,j,k,nout)+              &
!>   &                                 tau*(BOUNDARY(ng)%u_east(j,k)-   &
!>   &                                      u(Iend+1,j,k,nstp))
!>
                    tl_u(Iend+1,j,k,nout)=tl_u(Iend+1,j,k,nout)-        &
     &                                    tau*tl_u(Iend+1,j,k,nstp)
                  END IF
# ifdef MASKING
!>                u(Iend+1,j,k,nout)=u(Iend+1,j,k,nout)*                &
!>   &                               GRID(ng)%umask(Iend+1,j)
!>
                  tl_u(Iend+1,j,k,nout)=tl_u(Iend+1,j,k,nout)*          &
     &                                  GRID(ng)%umask(Iend+1,j)
# endif
                END IF
              END DO
            END DO
          END IF
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (tl_LBC(ieast,isUvel,ng)%clamped) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
!>              u(Iend+1,j,k,nout)=BOUNDARY(ng)%u_east(j,k)
!>
                tl_u(Iend+1,j,k,nout)=BOUNDARY(ng)%tl_u_east(j,k)
# ifdef MASKING
!>              u(Iend+1,j,k,nout)=u(Iend+1,j,k,nout)*                  &
!>   &                             GRID(ng)%umask(Iend+1,j)
!>
                tl_u(Iend+1,j,k,nout)=tl_u(Iend+1,j,k,nout)*            &
     &                                GRID(ng)%umask(Iend+1,j)
# endif
              END IF
            END DO
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (tl_LBC(ieast,isUvel,ng)%gradient) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
!>              u(Iend+1,j,k,nout)=u(Iend,j,k,nout)
!>
                tl_u(Iend+1,j,k,nout)=tl_u(Iend,j,k,nout)
# ifdef MASKING
!>              u(Iend+1,j,k,nout)=u(Iend+1,j,k,nout)*                  &
!>   &                             GRID(ng)%umask(Iend+1,j)
!>
                tl_u(Iend+1,j,k,nout)=tl_u(Iend+1,j,k,nout)*            &
     &                                GRID(ng)%umask(Iend+1,j)
# endif
              END IF
            END DO
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (tl_LBC(ieast,isUvel,ng)%closed) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
!>              u(Iend+1,j,k,nout)=0.0_r8
!>
                tl_u(Iend+1,j,k,nout)=0.0_r8
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
        IF (tl_LBC(isouth,isUvel,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
            DO k=1,N(ng)
              DO i=IstrU-1,Iend
!>              grad(i,Jstr-1)=u(i+1,Jstr-1,k,nstp)-                    &
!>   &                         u(i  ,Jstr-1,k,nstp)
!>
                tl_grad(i,Jstr-1)=0.0_r8
              END DO
              DO i=IstrU,Iend
                IF (LBC_apply(ng)%south(i)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                  IF (tl_LBC(isouth,isUvel,ng)%nudging) THEN
                    IF (LnudgeM3CLM(ng)) THEN
                      obc_out=0.5_r8*                                   &
     &                        (CLIMA(ng)%M3nudgcof(i-1,Jstr-1,k)+       &
     &                         CLIMA(ng)%M3nudgcof(i  ,Jstr-1,k))
                      obc_in =obcfac(ng)*obc_out
                    ELSE
                      obc_out=M3obc_out(ng,isouth)
                      obc_in =M3obc_in (ng,isouth)
                    END IF
                    IF (BOUNDARY(ng)%u_south_Ce(i,k).lt.0.0_r8) THEN
                      tau=obc_in
                    ELSE
                      tau=obc_out
                    END IF
                    tau=tau*dt(ng)
                  END IF
#  ifdef RADIATION_2D
                  Cx=BOUNDARY(ng)%u_south_Cx(i,k)
#  else
                  Cx=0.0_r8
#  endif
                  Ce=BOUNDARY(ng)%u_south_Ce(i,k)
                  cff=BOUNDARY(ng)%u_south_C2(i,k)
# endif
!>                u(i,Jstr-1,k,nout)=(cff*u(i,Jstr-1,k,nstp)+           &
!>   &                                Ce *u(i,Jstr  ,k,nout)-           &
!>   &                                MAX(Cx,0.0_r8)*grad(i-1,Jstr-1)-  &
!>   &                                MIN(Cx,0.0_r8)*grad(i  ,Jstr-1))/ &
!>   &                               (cff+Ce)
!>
                  tl_u(i,Jstr-1,k,nout)=(cff*tl_u(i,Jstr-1,k,nstp)+     &
     &                                   Ce *tl_u(i,Jstr  ,k,nout)-     &
     &                                   MAX(Cx,0.0_r8)*                &
     &                                      tl_grad(i-1,Jstr-1)-        &
     &                                   MIN(Cx,0.0_r8)*                &
     &                                      tl_grad(i  ,Jstr-1))/       &
     &                                  (cff+Ce)

                  IF (tl_LBC(isouth,isUvel,ng)%nudging) THEN
!>                  u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)+              &
!>   &                                 tau*(BOUNDARY(ng)%u_south(i,k)-  &
!>   &                                      u(i,Jstr-1,k,nstp))
!>
                    tl_u(i,Jstr-1,k,nout)=tl_u(i,Jstr-1,k,nout)-        &
     &                                    tau*tl_u(i,Jstr-1,k,nstp)
                  END IF
# ifdef MASKING
!>                u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)*                &
!>   &                               GRID(ng)%umask(i,Jstr-1)
!>
                  tl_u(i,Jstr-1,k,nout)=tl_u(i,Jstr-1,k,nout)*          &
     &                                  GRID(ng)%umask(i,Jstr-1)
# endif
                END IF
              END DO
            END DO
          END IF
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (tl_LBC(isouth,isUvel,ng)%clamped) THEN
          DO k=1,N(ng)
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%south(i)) THEN
!>              u(i,Jstr-1,k,nout)=BOUNDARY(ng)%u_south(i,k)
!>
                tl_u(i,Jstr-1,k,nout)=BOUNDARY(ng)%tl_u_south(i,k)
# ifdef MASKING
!>              u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)*                  &
!>   &                             GRID(ng)%umask(i,Jstr-1)
!>
                tl_u(i,Jstr-1,k,nout)=tl_u(i,Jstr-1,k,nout)*            &
     &                                GRID(ng)%umask(i,Jstr-1)
# endif
              END IF
            END DO
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (tl_LBC(isouth,isUvel,ng)%gradient) THEN
          DO k=1,N(ng)
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%south(i)) THEN
!>              u(i,Jstr-1,k,nout)=u(i,Jstr,k,nout)
!>
                tl_u(i,Jstr-1,k,nout)=tl_u(i,Jstr,k,nout)
# ifdef MASKING
!>              u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)*                  &
!>   &                             GRID(ng)%umask(i,Jstr-1)
!>
                tl_u(i,Jstr-1,k,nout)=tl_u(i,Jstr-1,k,nout)*            &
     &                                GRID(ng)%umask(i,Jstr-1)
# endif
              END IF
            END DO
          END DO
!
!  Southern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
        ELSE IF (tl_LBC(isouth,isUvel,ng)%closed) THEN
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
!>              u(i,Jstr-1,k,nout)=gamma2(ng)*u(i,Jstr,k,nout)
!>
                tl_u(i,Jstr-1,k,nout)=gamma2(ng)*tl_u(i,Jstr,k,nout)
# ifdef MASKING
!>              u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)*                  &
!>   &                             GRID(ng)%umask(i,Jstr-1)
!>
                tl_u(i,Jstr-1,k,nout)=tl_u(i,Jstr-1,k,nout)*            &
     &                                GRID(ng)%umask(i,Jstr-1)
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
        IF (tl_LBC(inorth,isUvel,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
            DO k=1,N(ng)
              DO i=IstrU-1,Iend
!>              grad(i,Jend+1)=u(i+1,Jend+1,k,nstp)-                    &
!>   &                         u(i  ,Jend+1,k,nstp)
!>
                tl_grad(i,Jend+1)=0.0_r8
              END DO
              DO i=IstrU,Iend
                IF (LBC_apply(ng)%north(i)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                  IF (tl_LBC(inorth,isUvel,ng)%nudging) THEN
                    IF (LnudgeM3CLM(ng)) THEN
                      obc_out=0.5_r8*                                   &
     &                        (CLIMA(ng)%M3nudgcof(i-1,Jend+1,k)+       &
     &                         CLIMA(ng)%M3nudgcof(i  ,Jend+1,k))
                      obc_in =obcfac(ng)*obc_out
                    ELSE
                      obc_out=M3obc_out(ng,inorth)
                      obc_in =M3obc_in (ng,inorth)
                    END IF
                    IF (BOUNDARY(ng)%u_north_Ce(i,k).lt.0.0_r8) THEN
                      tau=obc_in
                    ELSE
                      tau=obc_out
                    END IF
                    tau=tau*dt(ng)
                  END IF
#  ifdef RADIATION_2D
                  Cx=BOUNDARY(ng)%u_north_Cx(i,k)
#  else
                  Cx=0.0_r8
#  endif
                  Ce=BOUNDARY(ng)%u_north_Ce(i,k)
                  cff=BOUNDARY(ng)%u_north_C2(i,k)
# endif
!>                u(i,Jend+1,k,nout)=(cff*u(i,Jend+1,k,nstp)+           &
!>   &                                Ce *u(i,Jend  ,k,nout)-           &
!>   &                                MAX(Cx,0.0_r8)*grad(i-1,Jend+1)-  &
!>   &                                MIN(Cx,0.0_r8)*grad(i  ,Jend+1))/ &
!>   &                               (cff+Ce)
!>
                  tl_u(i,Jend+1,k,nout)=(cff*tl_u(i,Jend+1,k,nstp)+     &
     &                                   Ce *tl_u(i,Jend  ,k,nout)-     &
     &                                   MAX(Cx,0.0_r8)*                &
     &                                      tl_grad(i-1,Jend+1)-        &
     &                                   MIN(Cx,0.0_r8)*                &
     &                                      tl_grad(i  ,Jend+1))/       &
     &                                  (cff+Ce)

                  IF (tl_LBC(inorth,isUvel,ng)%nudging) THEN
!>                  u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)+              &
!>   &                                 tau*(BOUNDARY(ng)%u_north(i,k)-  &
!>   &                                      u(i,Jend+1,k,nstp))
!>
                    tl_u(i,Jend+1,k,nout)=tl_u(i,Jend+1,k,nout)-        &
     &                                    tau*tl_u(i,Jend+1,k,nstp)
                  END IF
# ifdef MASKING
!>                u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)*                &
!>   &                               GRID(ng)%umask(i,Jend+1)
!>
                  tl_u(i,Jend+1,k,nout)=tl_u(i,Jend+1,k,nout)*          &
     &                                  GRID(ng)%umask(i,Jend+1)
# endif
                END IF
              END DO
            END DO
          END IF
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (tl_LBC(inorth,isUvel,ng)%clamped) THEN
          DO k=1,N(ng)
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%north(i)) THEN
!>              u(i,Jend+1,k,nout)=BOUNDARY(ng)%u_north(i,k)
!>
                tl_u(i,Jend+1,k,nout)=BOUNDARY(ng)%tl_u_north(i,k)
# ifdef MASKING
!>              u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)*                  &
!>   &                             GRID(ng)%umask(i,Jend+1)
!>
                tl_u(i,Jend+1,k,nout)=tl_u(i,Jend+1,k,nout)*            &
     &                                GRID(ng)%umask(i,Jend+1)
# endif
              END IF
            END DO
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (tl_LBC(inorth,isUvel,ng)%gradient) THEN
          DO k=1,N(ng)
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%north(i)) THEN
!>              u(i,Jend+1,k,nout)=u(i,Jend,k,nout)
!>
                tl_u(i,Jend+1,k,nout)=tl_u(i,Jend,k,nout)
# ifdef MASKING
!>              u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)*                  &
!>   &                             GRID(ng)%umask(i,Jend+1)
!>
                tl_u(i,Jend+1,k,nout)=tl_u(i,Jend+1,k,nout)*            &
     &                                GRID(ng)%umask(i,Jend+1)
# endif
              END IF
            END DO
          END DO
!
!  Northern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
        ELSE IF (tl_LBC(inorth,isUvel,ng)%closed) THEN
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
!>              u(i,Jend+1,k,nout)=gamma2(ng)*u(i,Jend,k,nout)
!>
                tl_u(i,Jend+1,k,nout)=gamma2(ng)*tl_u(i,Jend,k,nout)
# ifdef MASKING
!>              u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)*                  &
!>   &                             GRID(ng)%umask(i,Jend+1)
!>
                tl_u(i,Jend+1,k,nout)=tl_u(i,Jend+1,k,nout)*            &
     &                                GRID(ng)%umask(i,Jend+1)
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
!>            u(Istr,Jstr-1,k,nout)=0.5_r8*(u(Istr+1,Jstr-1,k,nout)+    &
!>   &                                      u(Istr  ,Jstr  ,k,nout))
!>
              tl_u(Istr,Jstr-1,k,nout)=0.5_r8*                          &
     &                                 (tl_u(Istr+1,Jstr-1,k,nout)+     &
     &                                  tl_u(Istr  ,Jstr  ,k,nout))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            DO k=1,N(ng)
!>            u(Iend+1,Jstr-1,k,nout)=0.5_r8*(u(Iend  ,Jstr-1,k,nout)+  &
!>   &                                        u(Iend+1,Jstr  ,k,nout))
!>
              tl_u(Iend+1,Jstr-1,k,nout)=0.5_r8*                        &
     &                                   (tl_u(Iend  ,Jstr-1,k,nout)+   &
     &                                    tl_u(Iend+1,Jstr  ,k,nout))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO k=1,N(ng)
!>            u(Istr,Jend+1,k,nout)=0.5_r8*(u(Istr  ,Jend  ,k,nout)+    &
!>   &                                      u(Istr+1,Jend+1,k,nout))
!>
              tl_u(Istr,Jend+1,k,nout)=0.5_r8*                          &
     &                                 (tl_u(Istr  ,Jend  ,k,nout)+     &
     &                                  tl_u(Istr+1,Jend+1,k,nout))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO k=1,N(ng)
!>            u(Iend+1,Jend+1,k,nout)=0.5_r8*(u(Iend+1,Jend  ,k,nout)+  &
!>   &                                        u(Iend  ,Jend+1,k,nout))
!>
              tl_u(Iend+1,Jend+1,k,nout)=0.5_r8*                        &
     &                                   (tl_u(Iend+1,Jend  ,k,nout)+   &
     &                                    tl_u(Iend  ,Jend+1,k,nout))
            END DO
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE rp_u3dbc_tile
#endif
      END MODULE rp_u3dbc_mod
