#include "cppdefs.h"
      MODULE tl_v3dbc_mod
#if defined TANGENT && defined SOLVE3D
!
!svn $Id: tl_v3dbc_im.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine sets tangent linear lateral boundary conditions for !
!  total 3D V-velocity. It updates the specified "nout" time index.    !
!                                                                      !
!  BASIC STATE variables needed: v                                     !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: tl_v3dbc, tl_v3dbc_tile
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_v3dbc (ng, tile, nout)
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
      CALL tl_v3dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp(ng), nout,                               &
     &                    OCEAN(ng) % tl_v)

      RETURN
      END SUBROUTINE tl_v3dbc
!
!***********************************************************************
      SUBROUTINE tl_v3dbc_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, UBk,                &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nstp, nout,                             &
     &                          tl_v)
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
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,UBk,2)
# endif
!
!  Local variable declarations.
!
      integer :: Jmin, Jmax
      integer :: i, j, k

      real(r8) :: Ce, Cx, cff
      real(r8) :: obc_in, obc_out, tau

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_grad

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
        IF (tl_LBC(isouth,isVvel,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend+1
!>              grad(i,Jstr)=v(i  ,Jstr,k,nstp)-                        &
!>   &                       v(i-1,Jstr,k,nstp)
!>
                tl_grad(i,Jstr)=0.0_r8
              END DO
              DO i=Istr,Iend
                IF (LBC_apply(ng)%south(i)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                  IF (tl_LBC(isouth,isVvel,ng)%nudging) THEN
                    IF (LnudgeM3CLM(ng)) THEN
                      obc_out=0.5_r8*                                   &
     &                        (CLIMA(ng)%M3nudgcof(i,Jstr-1,k)+         &
     &                         CLIMA(ng)%M3nudgcof(i,Jstr  ,k))
                      obc_in =obcfac(ng)*obc_out
                    ELSE
                      obc_out=M3obc_out(ng,isouth)
                      obc_in =M3obc_in (ng,isouth)
                    END IF
                    IF (BOUNDARY(ng)%v_south_Ce(i,k).lt.0.0_r8) THEN
                      tau=obc_in
                    ELSE
                      tau=obc_out
                    END IF
                    tau=tau*dt(ng)
                  END IF
#  ifdef RADIATION_2D
                  Cx=BOUNDARY(ng)%v_south_Cx(i,k)
#  else
                  Cx=0.0_r8
#  endif
                  Ce=BOUNDARY(ng)%v_south_Ce(i,k)
                  cff=BOUNDARY(ng)%v_south_C2(i,k)
# endif
!>                v(i,Jstr,k,nout)=(cff*v(i,Jstr  ,k,nstp)+             &
!>   &                              Ce *v(i,Jstr+1,k,nout)-             &
!>   &                              MAX(Cx,0.0_r8)*grad(i  ,Jstr)-      &
!>   &                              MIN(Cx,0.0_r8)*grad(i+1,Jstr))/     &
!>   &                             (cff+Ce)
!>
                  tl_v(i,Jstr,k,nout)=(cff*tl_v(i,Jstr  ,k,nstp)+       &
     &                                 Ce *tl_v(i,Jstr+1,k,nout)-       &
     &                                 MAX(Cx,0.0_r8)*                  &
     &                                    tl_grad(i  ,Jstr)-            &
     &                                 MIN(Cx,0.0_r8)*                  &
     &                                    tl_grad(i+1,Jstr))/           &
     &                                (cff+Ce)

                  IF (tl_LBC(isouth,isVvel,ng)%nudging) THEN
!>                  v(i,Jstr,k,nout)=v(i,Jstr,k,nout)+                  &
!>   &                               tau*(BOUNDARY(ng)%v_south(i,k)-    &
!>   &                                    v(i,Jstr,k,nstp))
!>
                    tl_v(i,Jstr,k,nout)=tl_v(i,Jstr,k,nout)-            &
     &                                  tau*tl_v(i,Jstr,k,nstp)
                  END IF
# ifdef MASKING
!>                v(i,Jstr,k,nout)=v(i,Jstr,k,nout)*                    &
!>   &                             GRID(ng)%vmask(i,Jstr)
!>
                  tl_v(i,Jstr,k,nout)=tl_v(i,Jstr,k,nout)*              &
     &                                GRID(ng)%vmask(i,Jstr)
# endif
                END IF
              END DO
            END DO
          END IF
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (tl_LBC(isouth,isVvel,ng)%clamped) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
!>              v(i,Jstr,k,nout)=BOUNDARY(ng)%v_south(i,k)
!>
# ifdef ADJUST_BOUNDARY
                IF (Lobc(isouth,isVvel,ng)) THEN
                  tl_v(i,Jstr,k,nout)=BOUNDARY(ng)%tl_v_south(i,k)
                ELSE
                  tl_v(i,Jstr,k,nout)=0.0_r8
                END IF
# else
                tl_v(i,Jstr,k,nout)=0.0_r8
# endif
# ifdef MASKING
!>              v(i,Jstr,k,nout)=v(i,Jstr,k,nout)*                      &
!>   &                           GRID(ng)%vmask(i,Jstr)
!>
                tl_v(i,Jstr,k,nout)=tl_v(i,Jstr,k,nout)*                &
     &                              GRID(ng)%vmask(i,Jstr)
# endif
              END IF
            END DO
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (tl_LBC(isouth,isVvel,ng)%gradient) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
!>              v(i,Jstr,k,nout)=v(i,Jstr+1,k,nout)
!>
                tl_v(i,Jstr,k,nout)=tl_v(i,Jstr+1,k,nout)
# ifdef MASKING
!>              v(i,Jstr,k,nout)=v(i,Jstr,k,nout)*                      &
!>   &                           GRID(ng)%vmask(i,Jstr)
!>
                tl_v(i,Jstr,k,nout)=tl_v(i,Jstr,k,nout)*                &
     &                              GRID(ng)%vmask(i,Jstr)
# endif
              END IF
            END DO
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE IF (tl_LBC(isouth,isVvel,ng)%closed) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
!>              v(i,Jstr,k,nout)=0.0_r8
!>
                tl_v(i,Jstr,k,nout)=0.0_r8
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
        IF (tl_LBC(inorth,isVvel,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend+1
!>              grad(i,Jend+1)=v(i  ,Jend+1,k,nstp)-                    &
!>   &                         v(i-1,Jend+1,k,nstp)
!>
                tl_grad(i,Jend+1)=0.0_r8
              END DO
              DO i=Istr,Iend
                IF (LBC_apply(ng)%north(i)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                  IF (tl_LBC(inorth,isVvel,ng)%nudging) THEN
                    IF (LnudgeM3CLM(ng)) THEN
                      obc_out=0.5_r8*                                   &
     &                        (CLIMA(ng)%M3nudgcof(i,Jend  ,k)+         &
     &                         CLIMA(ng)%M3nudgcof(i,Jend+1,k))
                      obc_in =obcfac(ng)*obc_out
                    ELSE
                      obc_out=M3obc_out(ng,inorth)
                      obc_in =M3obc_in (ng,inorth)
                    END IF
                    IF (BOUNDARY(ng)%v_south_Ce(i,k).lt.0.0_r8) THEN
                      tau=obc_in
                    ELSE
                      tau=obc_out
                    END IF
                    tau=tau*dt(ng)
                  END IF
#  ifdef RADIATION_2D
                  Cx=BOUNDARY(ng)%v_south_Cx(i,k)
#  else
                  Cx=0.0_r8
#  endif
                  Ce=BOUNDARY(ng)%v_south_Ce(i,k)
                  cff=BOUNDARY(ng)%v_south_C2(i,k)
# endif
!>                v(i,Jend+1,k,nout)=(cff*v(i,Jend+1,k,nstp)+           &
!>   &                                Ce *v(i,Jend  ,k,nout)-           &
!>   &                                MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-  &
!>   &                                MIN(Cx,0.0_r8)*grad(i+1,Jend+1))/ &
!>   &                               (cff+Ce)
!>
                  tl_v(i,Jend+1,k,nout)=(cff*tl_v(i,Jend+1,k,nstp)+     &
     &                                   Ce *tl_v(i,Jend  ,k,nout)-     &
     &                                   MAX(Cx,0.0_r8)*                &
     &                                      tl_grad(i  ,Jend+1)-        &
     &                                   MIN(Cx,0.0_r8)*                &
     &                                      tl_grad(i+1,Jend+1))/       &
     &                                  (cff+Ce)

                  IF (tl_LBC(inorth,isVvel,ng)%nudging) THEN
!>                  v(i,Jend+1,k,nout)=v(i,Jend+1,k,nout)+              &
!>   &                                 tau*(BOUNDARY(ng)%v_north(i,k)-  &
!>   &                                      v(i,Jend+1,k,nstp))
!>
                    tl_v(i,Jend+1,k,nout)=tl_v(i,Jend+1,k,nout)-        &
     &                                    tau*tl_v(i,Jend+1,k,nstp)
                  END IF
# ifdef MASKING
!>                v(i,Jend+1,k,nout)=v(i,Jend+1,k,nout)*                &
!>   &                               GRID(ng)%vmask(i,Jend+1)
!>
                  tl_v(i,Jend+1,k,nout)=tl_v(i,Jend+1,k,nout)*          &
     &                                  GRID(ng)%vmask(i,Jend+1)
# endif
                END IF
              END DO
            END DO
          END IF
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (tl_LBC(inorth,isVvel,ng)%clamped) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
!>              v(i,Jend+1,k,nout)=BOUNDARY(ng)%v_north(i,k)
!>
# ifdef ADJUST_BOUNDARY
                IF (Lobc(inorth,isVvel,ng)) THEN
                  tl_v(i,Jend+1,k,nout)=BOUNDARY(ng)%tl_v_north(i,k)
                ELSE
                  tl_v(i,Jend+1,k,nout)=0.0_r8
                END IF
# else
                tl_v(i,Jend+1,k,nout)=0.0_r8
# endif
# ifdef MASKING
!>              v(i,Jend+1,k,nout)=v(i,Jend+1,k,nout)*                  &
!>   &                             GRID(ng)%vmask(i,Jend+1)
!>
                tl_v(i,Jend+1,k,nout)=tl_v(i,Jend+1,k,nout)*            &
     &                                GRID(ng)%vmask(i,Jend+1)
# endif
              END IF
            END DO
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (tl_LBC(inorth,isVvel,ng)%gradient) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
!>              v(i,Jend+1,k,nout)=v(i,Jend,k,nout)
!>
                tl_v(i,Jend+1,k,nout)=tl_v(i,Jend,k,nout)
# ifdef MASKING
!>              v(i,Jend+1,k,nout)=v(i,Jend+1,k,nout)*                  &
!>   &                             GRID(ng)%vmask(i,Jend+1)
!>
                tl_v(i,Jend+1,k,nout)=tl_v(i,Jend+1,k,nout)*            &
     &                                GRID(ng)%vmask(i,Jend+1)
# endif
              END IF
            END DO
          END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (tl_LBC(inorth,isVvel,ng)%closed) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
!>              v(i,Jend+1,k,nout)=0.0_r8
!>
                tl_v(i,Jend+1,k,nout)=0.0_r8
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
        IF (tl_LBC(iwest,isVvel,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
            DO k=1,N(ng)
              DO j=JstrV-1,Jend
!>              grad(Istr-1,j)=v(Istr-1,j+1,k,nstp)-                    &
!>   &                         v(Istr-1,j  ,k,nstp)
!>
                tl_grad(Istr-1,j)=0.0_r8
              END DO
              DO j=JstrV,Jend
                IF (LBC_apply(ng)%west(j)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                  IF (tl_LBC(iwest,isVvel,ng)%nudging) THEN
                    IF (LnudgeM3CLM(ng)) THEN
                      obc_out=0.5_r8*                                   &
     &                        (CLIMA(ng)%M3nudgcof(Istr-1,j-1,k)+       &
     &                         CLIMA(ng)%M3nudgcof(Istr-1,j  ,k))
                      obc_in =obcfac(ng)*obc_out
                    ELSE
                      obc_out=M3obc_out(ng,iwest)
                      obc_in =M3obc_in (ng,iwest)
                    END IF
                    IF (BOUNDARY(ng)%v_west_Cx(j,k).lt.0.0_r8) THEN
                      tau=obc_in
                    ELSE
                      tau=obc_out
                    END IF
                    tau=tau*dt(ng)
                  END IF
                  Cx=BOUNDARY(ng)%v_west_Cx(j,k)
#  ifdef RADIATION_2D
                  Ce=BOUNDARY(ng)%v_west_Ce(j,k)
#  else
                  Ce=0.0_r8
#  endif
                  cff=BOUNDARY(ng)%v_west_C2(j,k)
# endif
!>                v(Istr-1,j,k,nout)=(cff*v(Istr-1,j,k,nstp)+           &
!>   &                                Cx *v(Istr  ,j,k,nout)-           &
!>   &                                MAX(Ce,0.0_r8)*grad(Istr-1,j-1)-  &
!>   &                                MIN(Ce,0.0_r8)*grad(Istr-1,j  ))/ &
!>   &                               (cff+Cx)
!>
                  tl_v(Istr-1,j,k,nout)=(cff*tl_v(Istr-1,j,k,nstp)+     &
     &                                   Cx *tl_v(Istr  ,j,k,nout)-     &
     &                                   MAX(Ce,0.0_r8)*                &
     &                                      tl_grad(Istr-1,j-1)-        &
     &                                   MIN(Ce,0.0_r8)*                &
     &                                      tl_grad(Istr-1,j  ))/       &
     &                                  (cff+Cx)

                  IF (tl_LBC(iwest,isVvel,ng)%nudging) THEN
!>                  v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)+              &
!>   &                                 tau*(BOUNDARY(ng)%v_west(j,k)-   &
!>   &                                      v(Istr-1,j,k,nstp))
!>
                    tl_v(Istr-1,j,k,nout)=tl_v(Istr-1,j,k,nout)-        &
     &                                    tau*tl_v(Istr-1,j,k,nstp)
                  END IF
# ifdef MASKING
!>                v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                &
!>   &                               GRID(ng)%vmask(Istr-1,j)
!>
                  tl_v(Istr-1,j,k,nout)=tl_v(Istr-1,j,k,nout)*          &
     &                                  GRID(ng)%vmask(Istr-1,j)
# endif
                END IF
              END DO
            END DO
          END IF
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (tl_LBC(iwest,isVvel,ng)%clamped) THEN
          DO k=1,N(ng)
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%west(j)) THEN
!>              v(Istr-1,j,k,nout)=BOUNDARY(ng)%v_west(j,k)
!>
# ifdef ADJUST_BOUNDARY
                IF (Lobc(iwest,isVvel,ng)) THEN
                  tl_v(Istr-1,j,k,nout)=BOUNDARY(ng)%tl_v_west(j,k)
                ELSE
                  tl_v(Istr-1,j,k,nout)=0.0_r8
                END IF
# else
                tl_v(Istr-1,j,k,nout)=0.0_r8
# endif
# ifdef MASKING
!>              v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                  &
!>   &                             GRID(ng)%vmask(Istr-1,j)
!>
                tl_v(Istr-1,j,k,nout)=tl_v(Istr-1,j,k,nout)*            &
     &                                GRID(ng)%vmask(Istr-1,j)
# endif
              END IF
            END DO
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (tl_LBC(iwest,isVvel,ng)%gradient) THEN
          DO k=1,N(ng)
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%west(j)) THEN
!>              v(Istr-1,j,k,nout)=v(Istr,j,k,nout)
!>
                tl_v(Istr-1,j,k,nout)=tl_v(Istr,j,k,nout)
# ifdef MASKING
!>              v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                  &
!>   &                             GRID(ng)%vmask(Istr-1,j)
!>
                tl_v(Istr-1,j,k,nout)=tl_v(Istr-1,j,k,nout)*            &
     &                                GRID(ng)%vmask(Istr-1,j)
# endif
              END IF
            END DO
          END DO
!
!  Western edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
        ELSE IF (tl_LBC(iwest,isVvel,ng)%closed) THEN
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
!>              v(Istr-1,j,k,nout)=gamma2(ng)*v(Istr,j,k,nout)
!>
                tl_v(Istr-1,j,k,nout)=gamma2(ng)*tl_v(Istr,j,k,nout)
# ifdef MASKING
!>              v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                  &
!>   &                             GRID(ng)%vmask(Istr-1,j)
!>
                tl_v(Istr-1,j,k,nout)=tl_v(Istr-1,j,k,nout)*            &
     &                                GRID(ng)%vmask(Istr-1,j)
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
        IF (tl_LBC(ieast,isVvel,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
            DO k=1,N(ng)
              DO j=JstrV-1,Jend
!>              grad(Iend+1,j)=v(Iend+1,j+1,k,nstp)-                    &
!>   &                         v(Iend+1,j  ,k,nstp)
!>
                tl_grad(Iend+1,j)=0.0_r8
              END DO
              DO j=JstrV,Jend
                IF (LBC_apply(ng)%east(j)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                  IF (tl_LBC(ieast,isVvel,ng)%nudging) THEN
                    IF (LnudgeM3CLM(ng)) THEN
                      obc_out=0.5_r8*                                   &
     &                        (CLIMA(ng)%M3nudgcof(Iend+1,j-1,k)+       &
     &                         CLIMA(ng)%M3nudgcof(Iend+1,j  ,k))
                      obc_in =obcfac(ng)*obc_out
                    ELSE
                      obc_out=M3obc_out(ng,ieast)
                      obc_in =M3obc_in (ng,ieast)
                    END IF
                    IF (BOUNDARY(ng)%v_east_Cx(j,k).lt.0.0_r8) THEN
                      tau=obc_in
                    ELSE
                      tau=obc_out
                    END IF
                    tau=tau*dt(ng)
                  END IF
                  Cx=BOUNDARY(ng)%v_east_Cx(j,k)
#  ifdef RADIATION_2D
                  Ce=BOUNDARY(ng)%v_east_Ce(j,k)
#  else
                  Ce=0.0_r8
#  endif
                  cff=BOUNDARY(ng)%v_east_C2(j,k)
# endif
!>                v(Iend+1,j,k,nout)=(cff*v(Iend+1,j,k,nstp)+           &
!>   &                                Cx *v(Iend  ,j,k,nout)-           &
!>   &                                MAX(Ce,0.0_r8)*grad(Iend+1,j-1)-  &
!>   &                                MIN(Ce,0.0_r8)*grad(Iend+1,j  ))/ &
!>   &                               (cff+Cx)
!>
                  tl_v(Iend+1,j,k,nout)=(cff*tl_v(Iend+1,j,k,nstp)+     &
     &                                   Cx *tl_v(Iend  ,j,k,nout)-     &
     &                                   MAX(Ce,0.0_r8)*                &
     &                                      tl_grad(Iend+1,j-1)-        &
     &                                   MIN(Ce,0.0_r8)*                &
     &                                      tl_grad(Iend+1,j  ))/       &
     &                                  (cff+Cx)

                  IF (tl_LBC(ieast,isVvel,ng)%nudging) THEN
!>                  v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)+              &
!>   &                                 tau*(BOUNDARY(ng)%v_east(j,k)-   &
!>   &                                      v(Iend+1,j,k,nstp))
!>
                    tl_v(Iend+1,j,k,nout)=tl_v(Iend+1,j,k,nout)-        &
     &                                  tau*tl_v(Iend+1,j,k,nstp)
                  END IF
# ifdef MASKING
!>                v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                &
!>   &                               GRID(ng)%vmask(Iend+1,j)
!>
                  tl_v(Iend+1,j,k,nout)=tl_v(Iend+1,j,k,nout)*          &
     &                                  GRID(ng)%vmask(Iend+1,j)
# endif
                END IF
              END DO
            END DO
          END IF
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (tl_LBC(ieast,isVvel,ng)%clamped) THEN
          DO k=1,N(ng)
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%east(j)) THEN
!>              v(Iend+1,j,k,nout)=BOUNDARY(ng)%v_east(j,k)
!>
# ifdef ADJUST_BOUNDARY
                IF (Lobc(ieast,isVvel,ng)) THEN
                  tl_v(Iend+1,j,k,nout)=BOUNDARY(ng)%tl_v_east(j,k)
                ELSE
                  tl_v(Iend+1,j,k,nout)=0.0_r8
                END IF
# else
                tl_v(Iend+1,j,k,nout)=0.0_r8
# endif
# ifdef MASKING
!>              v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                  &
!>   &                             GRID(ng)%vmask(Iend+1,j)
!>
                tl_v(Iend+1,j,k,nout)=tl_v(Iend+1,j,k,nout)*            &
     &                                GRID(ng)%vmask(Iend+1,j)
# endif
              END IF
            END DO
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (tl_LBC(ieast,isVvel,ng)%gradient) THEN
          DO k=1,N(ng)
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%east(j)) THEN
!>              v(Iend+1,j,k,nout)=v(Iend,j,k,nout)
!>
                tl_v(Iend+1,j,k,nout)=tl_v(Iend,j,k,nout)
# ifdef MASKING
!>              v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                  &
!>   &                             GRID(ng)%vmask(Iend+1,j)
!>
                tl_v(Iend+1,j,k,nout)=tl_v(Iend+1,j,k,nout)*            &
     &                                GRID(ng)%vmask(Iend+1,j)
# endif
              END IF
            END DO
          END DO
!
!  Eastern edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
        ELSE IF (tl_LBC(ieast,isVvel,ng)%closed) THEN
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
!>              v(Iend+1,j,k,nout)=gamma2(ng)*v(Iend,j,k,nout)
!>
                tl_v(Iend+1,j,k,nout)=gamma2(ng)*tl_v(Iend,j,k,nout)
# ifdef MASKING
!>              v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                  &
!>   &                             GRID(ng)%vmask(Iend+1,j)
!>
                tl_v(Iend+1,j,k,nout)=tl_v(Iend+1,j,k,nout)*            &
     &                                GRID(ng)%vmask(Iend+1,j)
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
!>            v(Istr-1,Jstr,k,nout)=0.5_r8*(v(Istr  ,Jstr  ,k,nout)+    &
!>   &                                      v(Istr-1,Jstr+1,k,nout))
!>
              tl_v(Istr-1,Jstr,k,nout)=0.5_r8*                          &
     &                                 (tl_v(Istr  ,Jstr  ,k,nout)+     &
     &                                  tl_v(Istr-1,Jstr+1,k,nout))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr  )) THEN
            DO k=1,N(ng)
!>            v(Iend+1,Jstr,k,nout)=0.5_r8*(v(Iend  ,Jstr  ,k,nout)+    &
!>   &                                      v(Iend+1,Jstr+1,k,nout))
!>
              tl_v(Iend+1,Jstr,k,nout)=0.5_r8*                          &
     &                                 (tl_v(Iend  ,Jstr  ,k,nout)+     &
     &                                  tl_v(Iend+1,Jstr+1,k,nout))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO k=1,N(ng)
!>            v(Istr-1,Jend+1,k,nout)=0.5_r8*(v(Istr-1,Jend  ,k,nout)+  &
!>   &                                        v(Istr  ,Jend+1,k,nout))
!>
              tl_v(Istr-1,Jend+1,k,nout)=0.5_r8*                        &
     &                                   (tl_v(Istr-1,Jend  ,k,nout)+   &
     &                                    tl_v(Istr  ,Jend+1,k,nout))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO k=1,N(ng)
!>            v(Iend+1,Jend+1,k,nout)=0.5_r8*(v(Iend+1,Jend  ,k,nout)+  &
!>   &                                        v(Iend  ,Jend+1,k,nout))
!>
              tl_v(Iend+1,Jend+1,k,nout)=0.5_r8*                        &
     &                                   (tl_v(Iend+1,Jend  ,k,nout)+   &
     &                                    tl_v(Iend  ,Jend+1,k,nout))
            END DO
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE tl_v3dbc_tile
#endif
      END MODULE tl_v3dbc_mod
