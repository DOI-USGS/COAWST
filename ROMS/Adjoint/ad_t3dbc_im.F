#include "cppdefs.h"
      MODULE ad_t3dbc_mod
#if defined ADJOINT && defined SOLVE3D
!
!svn $Id: ad_t3dbc_im.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine sets adjoint lateral boundary conditions for the    !
!  ITRC-th tracer field. It updates the specified "nout" time index.   !
!                                                                      !
!  BASIC STATE variables needed: t                                     !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ad_t3dbc, ad_t3dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_t3dbc (ng, tile, nout, itrc, ic)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, nout, itrc, ic
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL ad_t3dbc_tile (ng, tile, itrc, ic,                           &
     &                    LBi, UBi, LBj, UBj, N(ng), NT(ng),            &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp(ng), nout,                               &
     &                    OCEAN(ng)% ad_t)
      RETURN
      END SUBROUTINE ad_t3dbc

!
!***********************************************************************
      SUBROUTINE ad_t3dbc_tile (ng, tile, itrc, ic,                     &
     &                          LBi, UBi, LBj, UBj, UBk, UBt,           &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nstp, nout,                             &
     &                          ad_t)
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
      integer, intent(in) :: ng, tile, itrc, ic
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nout
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
# else
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,UBk,3,UBt)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8) :: Ce, Cx, cff
      real(r8) :: obc_in, obc_out, tau
      real(r8) :: adfac

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_grad

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_grad(LBi:UBi,LBj:UBj)=0.0_r8
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO k=1,N(ng)
!>            tl_t(Iend+1,Jend+1,k,nout,itrc)=0.5_r8*                   &
!>   &                                        (tl_t(Iend+1,Jend  ,k,    &
!>   &                                              nout,itrc)+         &
!>   &                                         tl_t(Iend  ,Jend+1,k,    &
!>   &                                              nout,itrc))
!>
              adfac=0.5_r8*ad_t(Iend+1,Jend+1,k,nout,itrc)
              ad_t(Iend+1,Jend  ,k,nout,itrc)=ad_t(Iend+1,Jend  ,k,     &
     &                                             nout,itrc)+          &
     &                                        adfac
              ad_t(Iend  ,Jend+1,k,nout,itrc)=ad_t(Iend  ,Jend+1,k,     &
     &                                             nout,itrc)+          &
     &                                        adfac
              ad_t(Iend+1,Jend+1,k,nout,itrc)=0.0_r8
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO k=1,N(ng)
!>            tl_t(Istr-1,Jend+1,k,nout,itrc)=0.5_r8*                   &
!>   &                                        (tl_t(Istr-1,Jend  ,k,    &
!>   &                                              nout,itrc)+         &
!>   &                                         tl_t(Istr  ,Jend+1,k,    &
!>   &                                              nout,itrc))
!>
              adfac=0.5_r8*ad_t(Istr-1,Jend+1,k,nout,itrc)
              ad_t(Istr-1,Jend  ,k,nout,itrc)=ad_t(Istr-1,Jend  ,k,     &
     &                                             nout,itrc)+          &
     &                                        adfac
              ad_t(Istr  ,Jend+1,k,nout,itrc)=ad_t(Istr  ,Jend+1,k,     &
     &                                             nout,itrc)+          &
     &                                        adfac
              ad_t(Istr-1,Jend+1,k,nout,itrc)=0.0_r8
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            DO k=1,N(ng)
!>            tl_t(Iend+1,Jstr-1,k,nout,itrc)=0.5_r8*                   &
!>   &                                        (tl_t(Iend  ,Jstr-1,k,    &
!>   &                                              nout,itrc)+         &
!>   &                                         tl_t(Iend+1,Jstr  ,k,    &
!>   &                                              nout,itrc))
!>
              adfac=0.5_r8*ad_t(Iend+1,Jstr-1,k,nout,itrc)
              ad_t(Iend  ,Jstr-1,k,nout,itrc)=ad_t(Iend  ,Jstr-1,k,     &
     &                                             nout,itrc)+          &
     &                                        adfac
              ad_t(Iend+1,Jstr  ,k,nout,itrc)=ad_t(Iend+1,Jstr  ,k,     &
     &                                             nout,itrc)+          &
     &                                        adfac
              ad_t(Iend+1,Jstr-1,k,nout,itrc)=0.0_r8
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
            DO k=1,N(ng)
!>            tl_t(Istr-1,Jstr-1,k,nout,itrc)=0.5_r8*                   &
!>   &                                        (tl_t(Istr  ,Jstr-1,k,    &
!>   &                                              nout,itrc)+         &
!>   &                                         tl_t(Istr-1,Jstr  ,k,    &
!>   &                                              nout,itrc))
!>
              adfac=0.5_r8*ad_t(Istr-1,Jstr-1,k,nout,itrc)
              ad_t(Istr  ,Jstr-1,k,nout,itrc)=ad_t(Istr  ,Jstr-1,k,     &
     &                                             nout,itrc)+          &
     &                                        adfac
              ad_t(Istr-1,Jstr  ,k,nout,itrc)=ad_t(Istr-1,Jstr  ,k,     &
     &                                             nout,itrc)+          &
     &                                        adfac
              ad_t(Istr-1,Jstr-1,k,nout,itrc)=0.0_r8
            END DO
          END IF
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
        IF (ad_LBC(inorth,isTvar(itrc),ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                IF (LBC_apply(ng)%north(i)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                  IF (ad_LBC(inorth,isTvar(itrc),ng)%nudging) THEN
                    IF (LnudgeTCLM(itrc,ng)) THEN
                      obc_out=CLIMA(ng)%Tnudgcof(i,Jend+1,k,ic)
                      obc_in =obcfac(ng)*obc_out
                    ELSE
                      obc_out=Tobc_out(itrc,ng,inorth)
                      obc_in =Tobc_in (itrc,ng,inorth)
                    END IF
                    IF (BOUNDARY(ng)%t_north_Ce(i,k,itrc).lt.           &
     &                  0.0_r8) THEN
                      tau=obc_in
                    ELSE
                      tau=obc_out
                    END IF
                    tau=tau*dt(ng)
                  END IF
#  ifdef RADIATION_2D
                  Cx=BOUNDARY(ng)%t_north_Cx(i,k,itrc)
#  else
                  Cx=0.0_r8
#  endif
                  Ce=BOUNDARY(ng)%t_north_Ce(i,k,itrc)
                  cff=BOUNDARY(ng)%t_north_C2(i,k,itrc)
# endif
# ifdef MASKING
!>                tl_t(i,Jend+1,k,nout,itrc)=tl_t(i,Jend+1,k,nout,      &
!>   &                                            itrc)*                &
!>   &                                       GRID(ng)%rmask(i,Jend+1)
!>
                  ad_t(i,Jend+1,k,nout,itrc)=ad_t(i,Jend+1,k,nout,      &
     &                                            itrc)*                &
     &                                       GRID(ng)%rmask(i,Jend+1)
# endif
                  IF (ad_LBC(inorth,isTvar(itrc),ng)%nudging) THEN
!>                  tl_t(i,Jend+1,k,nout,itrc)=tl_t(i,Jend+1,k,nout,    &
!>   &                                              itrc)-              &
!>   &                                         tau*                     &
!>   &                                         tl_t(i,Jend+1,k,nstp,    &
!>   &                                              itrc)
!>
                    ad_t(i,Jend+1,k,nstp,itrc)=ad_t(i,Jend+1,k,nstp,    &
     &                                              itrc)-              &
     &                                         tau*                     &
     &                                         ad_t(i,Jend+1,k,nout,    &
     &                                              itrc)
                  END IF
!>                tl_t(i,Jend+1,k,nout,itrc)=(cff*                      &
!>   &                                        tl_t(i,Jend+1,k,nstp,     &
!>   &                                             itrc)+               &
!>   &                                        Ce*                       &
!>   &                                        tl_t(i,Jend  ,k,nout,     &
!>   &                                             itrc)-               &
!>   &                                        MAX(Cx,0.0_r8)*           &
!>   &                                           tl_grad(i  ,Jend+1)-   &
!>   &                                        MIN(Cx,0.0_r8)*           &
!>   &                                           tl_grad(i+1,Jend+1))/  &
!>   &                                       (cff+Ce)
!>
                  adfac=ad_t(i,Jend+1,k,nout,itrc)/(cff+Ce)
                  ad_grad(i  ,Jend+1)=ad_grad(i  ,Jend+1)-              &
     &                                MAX(Cx,0.0_r8)*adfac
                  ad_grad(i+1,Jend+1)=ad_grad(i+1,Jend+1)-              &
     &                                MIN(Cx,0.0_r8)*adfac
                  ad_t(i,Jend  ,k,nout,itrc)=ad_t(i,Jend  ,k,nout,      &
     &                                            itrc)+                &
     &                                       Ce *adfac
                  ad_t(i,Jend+1,k,nstp,itrc)=ad_t(i,Jend+1,k,nstp,      &
     &                                            itrc)+                &
     &                                       cff*adfac
                  ad_t(i,Jend+1,k,nout,itrc)=0.0_r8
                END IF
              END DO
            END DO
          END IF
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (ad_LBC(inorth,isTvar(itrc),ng)%clamped) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
# ifdef MASKING
!>              tl_t(i,Jend+1,k,nout,itrc)=tl_t(i,Jend+1,k,nout,itrc)*  &
!>   &                                     GRID(ng)%rmask(i,Jend+1)
!>
                ad_t(i,Jend+1,k,nout,itrc)=ad_t(i,Jend+1,k,nout,itrc)*  &
     &                                     GRID(ng)%rmask(i,Jend+1)
# endif
# ifdef ADJUST_BOUNDARY
                IF (Lobc(inorth,isTvar(itrc),ng)) THEN
!>                tl_t(i,Jend+1,k,nout,itrc)=BOUNDARY(ng)%              &
!>                                                 tl_t_north(i,k,itrc)
!>
                  BOUNDARY(ng)%ad_t_north(i,k,itrc)=                    &
     &                              BOUNDARY(ng)%ad_t_north(i,k,itrc)+  &
     &                                              ad_t(i,Jend+1,k,    &
     &                                                   nout,itrc)
                  ad_t(i,Jend+1,k,nout,itrc)=0.0_r8
                ELSE
!>                tl_t(i,Jend+1,k,nout,itrc)=0.0_r8
!>
                  ad_t(i,Jend+1,k,nout,itrc)=0.0_r8
                END IF
# else
!>              tl_t(i,Jend+1,k,nout,itrc)=0.0_r8
!>
                ad_t(i,Jend+1,k,nout,itrc)=0.0_r8
# endif
              END IF
            END DO
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (ad_LBC(inorth,isTvar(itrc),ng)%gradient) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
# ifdef MASKING
!>              tl_t(i,Jend+1,k,nout,itrc)=tl_t(i,Jend+1,k,nout,itrc)*  &
!>   &                                     GRID(ng)%rmask(i,Jend+1)
!>
                ad_t(i,Jend+1,k,nout,itrc)=ad_t(i,Jend+1,k,nout,itrc)*  &
     &                                     GRID(ng)%rmask(i,Jend+1)
# endif
!>              tl_t(i,Jend+1,k,nout,itrc)=tl_t(i,Jend,k,nout,itrc)
!>
                ad_t(i,Jend  ,k,nout,itrc)=ad_t(i,Jend  ,k,nout,itrc)+  &
     &                                     ad_t(i,Jend+1,k,nout,itrc)
                ad_t(i,Jend+1,k,nout,itrc)=0.0_r8
              END IF
            END DO
          END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (ad_LBC(inorth,isTvar(itrc),ng)%closed) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
# ifdef MASKING
!>              tl_t(i,Jend+1,k,nout,itrc)=tl_t(i,Jend+1,k,nout,itrc)*  &
!>   &                                     GRID(ng)%rmask(i,Jend+1)
!>
                ad_t(i,Jend+1,k,nout,itrc)=ad_t(i,Jend+1,k,nout,itrc)*  &
     &                                     GRID(ng)%rmask(i,Jend+1)
# endif
!>              tl_t(i,Jend+1,k,nout,itrc)=tl_t(i,Jend,k,nout,itrc)
!>
                ad_t(i,Jend  ,k,nout,itrc)=ad_t(i,Jend  ,k,nout,itrc)+  &
     &                                     ad_t(i,Jend+1,k,nout,itrc)
                ad_t(i,Jend+1,k,nout,itrc)=0.0_r8
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
        IF (ad_LBC(isouth,isTvar(itrc),ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                IF (LBC_apply(ng)%south(i)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                  IF (ad_LBC(isouth,isTvar(itrc),ng)%nudging) THEN
                    IF (LnudgeTCLM(itrc,ng)) THEN
                      obc_out=CLIMA(ng)%Tnudgcof(i,Jstr-1,k,ic)
                      obc_in =obcfac(ng)*obc_out
                    ELSE
                      obc_out=Tobc_out(itrc,ng,isouth)
                      obc_in =Tobc_in (itrc,ng,isouth)
                    END IF
                    IF (BOUNDARY(ng)%t_south_Ce(i,k,itrc).lt.           &
     &                  0.0_r8) THEN
                      tau=obc_in
                    ELSE
                      tau=obc_out
                    END IF
                    tau=tau*dt(ng)
                  END IF
#  ifdef RADIATION_2D
                  Cx=BOUNDARY(ng)%t_south_Cx(i,k,itrc)
#  else
                  Cx=0.0_r8
#  endif
                  Ce=BOUNDARY(ng)%t_south_Ce(i,k,itrc)
                  cff=BOUNDARY(ng)%t_south_C2(i,k,itrc)
# endif
# ifdef MASKING
!>                tl_t(i,Jstr-1,k,nout,itrc)=tl_t(i,Jstr-1,k,nout,      &
!>   &                                            itrc)*                &
!>   &                                       GRID(ng)%rmask(i,Jstr-1)
!>
                  ad_t(i,Jstr-1,k,nout,itrc)=ad_t(i,Jstr-1,k,nout,      &
     &                                            itrc)*                &
     &                                       GRID(ng)%rmask(i,Jstr-1)
# endif
                  IF (ad_LBC(isouth,isTvar(itrc),ng)%nudging) THEN
!>                  tl_t(i,Jstr-1,k,nout,itrc)=tl_t(i,Jstr-1,k,nout,    &
!>   &                                              itrc)-              &
!>   &                                         tau*                     &
!>   &                                         tl_t(i,Jstr-1,k,nstp,    &
!>   &                                              itrc)
!>
                    ad_t(i,Jstr-1,k,nstp,itrc)=ad_t(i,Jstr-1,k,nstp,    &
     &                                              itrc)-              &
     &                                         tau*                     &
     &                                         ad_t(i,Jstr-1,k,nout,    &
     &                                              itrc)
                  END IF
!>                tl_t(i,Jstr-1,k,nout,itrc)=(cff*                      &
!>   &                                        tl_t(i,Jstr-1,k,nstp,     &
!>   &                                             itrc)+               &
!>   &                                        Ce*                       &
!>   &                                        tl_t(i,Jstr  ,k,nout,     &
!>   &                                             itrc)-               &
!>   &                                        MAX(Cx,0.0_r8)*           &
!>   &                                           tl_grad(i  ,Jstr-1)-   &
!>   &                                        MIN(Cx,0.0_r8)*           &
!>   &                                           tl_grad(i+1,Jstr-1))/  &
!>   &                                       (cff+Ce)
!>
                  adfac=ad_t(i,Jstr-1,k,nout,itrc)/(cff+Ce)
                  ad_grad(i  ,Jstr-1)=ad_grad(i  ,Jstr-1)-              &
     &                                MAX(Cx,0.0_r8)*adfac
                  ad_grad(i+1,Jstr-1)=ad_grad(i+1,Jstr-1)-              &
     &                                MIN(Cx,0.0_r8)*adfac
                  ad_t(i,Jstr-1,k,nstp,itrc)=ad_t(i,Jstr-1,k,nstp,      &
     &                                            itrc)+                &
     &                                       cff*adfac
                  ad_t(i,Jstr  ,k,nout,itrc)=ad_t(i,Jstr  ,k,nout,      &
     &                                            itrc)+                &
     &                                       Ce *adfac
                  ad_t(i,Jstr-1,k,nout,itrc)=0.0_r8
                END IF
              END DO
            END DO
          END IF
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (ad_LBC(isouth,isTvar(itrc),ng)%clamped) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
# ifdef MASKING
!>              tl_t(i,Jstr-1,k,nout,itrc)=tl_t(i,Jstr-1,k,nout,itrc)*  &
!>   &                                     GRID(ng)%rmask(i,Jstr-1)
!>
                ad_t(i,Jstr-1,k,nout,itrc)=ad_t(i,Jstr-1,k,nout,itrc)*  &
     &                                     GRID(ng)%rmask(i,Jstr-1)
# endif
# ifdef ADJUST_BOUNDARY
                IF (Lobc(isouth,isTvar(itrc),ng)) THEN
!>                tl_t(i,Jstr-1,k,nout,itrc)=BOUNDARY(ng)%              &
!>   &                                             tl_t_south(i,k,itrc)
!>
                  BOUNDARY(ng)%ad_t_south(i,k,itrc)=                    &
     &                              BOUNDARY(ng)%ad_t_south(i,k,itrc)+  &
     &                                              ad_t(i,Jstr-1,k,    &
     &                                                   nout,itrc)
                  ad_t(i,Jstr-1,k,nout,itrc)=0.0_r8
                ELSE
!>                tl_t(i,Jstr-1,k,nout,itrc)=0.0_r8
!>
                  ad_t(i,Jstr-1,k,nout,itrc)=0.0_r8
                END IF
# else
!>              tl_t(i,Jstr-1,k,nout,itrc)=0.0_r8
!>
                ad_t(i,Jstr-1,k,nout,itrc)=0.0_r8
# endif
              END IF
            END DO
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (ad_LBC(isouth,isTvar(itrc),ng)%gradient) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
# ifdef MASKING
!>              tl_t(i,Jstr-1,k,nout,itrc)=tl_t(i,Jstr-1,k,nout,itrc)*  &
!>   &                                     GRID(ng)%rmask(i,Jstr-1)
!>
                ad_t(i,Jstr-1,k,nout,itrc)=ad_t(i,Jstr-1,k,nout,itrc)*  &
     &                                     GRID(ng)%rmask(i,Jstr-1)
# endif
!>              tl_t(i,Jstr-1,k,nout,itrc)=tl_t(i,Jstr,k,nout,itrc)
!>
                ad_t(i,Jstr  ,k,nout,itrc)=ad_t(i,Jstr  ,k,nout,itrc)+  &
     &                                     ad_t(i,Jstr-1,k,nout,itrc)
                ad_t(i,Jstr-1,k,nout,itrc)=0.0_r8
              END IF
            END DO
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE IF (ad_LBC(isouth,isTvar(itrc),ng)%closed) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
# ifdef MASKING
!>              tl_t(i,Jstr-1,k,nout,itrc)=tl_t(i,Jstr-1,k,nout,itrc)*  &
!>   &                                     GRID(ng)%rmask(i,Jstr-1)
!>
                ad_t(i,Jstr-1,k,nout,itrc)=ad_t(i,Jstr-1,k,nout,itrc)*  &
     &                                     GRID(ng)%rmask(i,Jstr-1)
# endif
!>              tl_t(i,Jstr-1,k,nout,itrc)=tl_t(i,Jstr,k,nout,itrc)
!>
                ad_t(i,Jstr  ,k,nout,itrc)=ad_t(i,Jstr  ,k,nout,itrc)+  &
     &                                     ad_t(i,Jstr-1,k,nout,itrc)
                ad_t(i,Jstr-1,k,nout,itrc)=0.0_r8
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
        IF (ad_LBC(ieast,isTvar(itrc),ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                IF (LBC_apply(ng)%east(j)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                  IF (ad_LBC(ieast,isTvar(itrc),ng)%nudging) THEN
                    IF (LnudgeTCLM(itrc,ng)) THEN
                      obc_out=CLIMA(ng)%Tnudgcof(Iend+1,j,k,ic)
                      obc_in =obcfac(ng)*obc_out
                    ELSE
                      obc_out=Tobc_out(itrc,ng,ieast)
                      obc_in =Tobc_in (itrc,ng,ieast)
                    END IF
                    IF (BOUNDARY(ng)%t_east_Cx(j,k,itrc).lt.            &
     &                  0.0_r8) THEN
                      tau=obc_in
                    ELSE
                      tau=obc_out
                    END IF
                    tau=tau*dt(ng)
                  END IF
                  Cx=BOUNDARY(ng)%t_east_Cx(j,k,itrc)
#  ifdef RADIATION_2D
                  Ce=BOUNDARY(ng)%t_east_Ce(j,k,itrc)
#  else
                  Ce=0.0_r8
#  endif
                  cff=BOUNDARY(ng)%t_east_C2(j,k,itrc)
# endif
# ifdef MASKING
!>                tl_t(Iend+1,j,k,nout,itrc)=tl_t(Iend+1,j,k,nout,      &
!>   &                                            itrc)*                &
!>   &                                       GRID(ng)%rmask(Iend+1,j)
!>
                  ad_t(Iend+1,j,k,nout,itrc)=ad_t(Iend+1,j,k,nout,      &
     &                                            itrc)*                &
     &                                       GRID(ng)%rmask(Iend+1,j)
# endif
                  IF (ad_LBC(ieast,isTvar(itrc),ng)%nudging) THEN
!>                  tl_t(Iend+1,j,k,nout,itrc)=tl_t(Iend+1,j,k,nout,    &
!>   &                                              itrc)-              &
!>   &                                         tau*                     &
!>   &                                         tl_t(Iend+1,j,k,nstp,    &
!>   &                                              itrc)
!>
                    ad_t(Iend+1,j,k,nstp,itrc)=ad_t(Iend+1,j,k,nstp,    &
     &                                              itrc)-              &
     &                                         tau*                     &
     &                                         ad_t(Iend+1,j,k,nout,    &
     &                                              itrc)
                  END IF
!>                tl_t(Iend+1,j,k,nout,itrc)=(cff*                      &
!>   &                                        tl_t(Iend+1,j,k,nstp,     &
!>   &                                             itrc)+               &
!>   &                                        Cx*                       &
!>   &                                        tl_t(Iend  ,j,k,nout,     &
!>   &                                             itrc)-               &
!>   &                                        MAX(Ce,0.0_r8)*           &
!>   &                                           tl_grad(Iend+1,j  )-   &
!>   &                                        MIN(Ce,0.0_r8)*           &
!>   &                                           tl_grad(Iend+1,j+1))/  &
!>   &                                       (cff+Cx)
!>
                  adfac=ad_t(Iend+1,j,k,nout,itrc)/(cff+Cx)
                  ad_grad(Iend+1,j  )=ad_grad(Iend+1,j  )-              &
     &                                MAX(Ce,0.0_r8)*adfac
                  ad_grad(Iend+1,j+1)=ad_grad(Iend+1,j+1)-              &
     &                                MIN(Ce,0.0_r8)*adfac
                  ad_t(Iend  ,j,k,nout,itrc)=ad_t(Iend  ,j,k,nout,      &
     &                                            itrc)+                &
     &                                       Cx *adfac
                  ad_t(Iend+1,j,k,nstp,itrc)=ad_t(Iend+1,j,k,nstp,      &
     &                                            itrc)+                &
     &                                       cff*adfac
                  ad_t(Iend+1,j,k,nout,itrc)=0.0_r8
                END IF
              END DO
            END DO
          END IF
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (ad_LBC(ieast,isTvar(itrc),ng)%clamped) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
# ifdef MASKING
!>              tl_t(Iend+1,j,k,nout,itrc)=tl_t(Iend+1,j,k,nout,itrc)*  &
!>   &                                     GRID(ng)%rmask(Iend+1,j)
!>
                ad_t(Iend+1,j,k,nout,itrc)=ad_t(Iend+1,j,k,nout,itrc)*  &
     &                                     GRID(ng)%rmask(Iend+1,j)
# endif
# ifdef ADJUST_BOUNDARY
                IF (Lobc(ieast,isTvar(itrc),ng)) THEN
!>                tl_t(Iend+1,j,k,nout,itrc)=BOUNDARY(ng)%              &
!>   &                                             tl_t_east(j,k,itrc)
!>
                  BOUNDARY(ng)%ad_t_east(j,k,itrc)=                     &
     &                              BOUNDARY(ng)%ad_t_east(j,k,itrc)+   &
     &                                             ad_t(Iend+1,j,k,     &
     &                                                  nout,itrc)
                  ad_t(Iend+1,j,k,nout,itrc)=0.0_r8
                ELSE
!>                tl_t(Iend+1,j,k,nout,itrc)=0.0_r8
!>
                  ad_t(Iend+1,j,k,nout,itrc)=0.0_r8
                END IF
# else
!>              tl_t(Iend+1,j,k,nout,itrc)=0.0_r8
!>
                ad_t(Iend+1,j,k,nout,itrc)=0.0_r8
# endif
              END IF
            END DO
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (ad_LBC(ieast,isTvar(itrc),ng)%gradient) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
# ifdef MASKING
!>              tl_t(Iend+1,j,k,nout,itrc)=tl_t(Iend+1,j,k,nout,itrc)*  &
!>   &                                     GRID(ng)%rmask(Iend+1,j)
!>
                ad_t(Iend+1,j,k,nout,itrc)=ad_t(Iend+1,j,k,nout,itrc)*  &
     &                                     GRID(ng)%rmask(Iend+1,j)
# endif
!>              tl_t(Iend+1,j,k,nout,itrc)=tl_t(Iend,j,k,nout,itrc)
!>
                ad_t(Iend  ,j,k,nout,itrc)=ad_t(Iend  ,j,k,nout,itrc)+  &
     &                                     ad_t(Iend+1,j,k,nout,itrc)
                ad_t(Iend+1,j,k,nout,itrc)=0.0_r8
              END IF
            END DO
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (ad_LBC(ieast,isTvar(itrc),ng)%closed) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
# ifdef MASKING
!>              tl_t(Iend+1,j,k,nout,itrc)=tl_t(Iend+1,j,k,nout,itrc)*  &
!>   &                                     GRID(ng)%rmask(Iend+1,j)
!>
                ad_t(Iend+1,j,k,nout,itrc)=ad_t(Iend+1,j,k,nout,itrc)*  &
     &                                     GRID(ng)%rmask(Iend+1,j)
# endif
!>              tl_t(Iend+1,j,k,nout,itrc)=tl_t(Iend,j,k,nout,itrc)
!>
                ad_t(Iend  ,j,k,nout,itrc)=ad_t(Iend  ,j,k,nout,itrc)+  &
     &                                     ad_t(Iend+1,j,k,nout,itrc)
                ad_t(Iend+1,j,k,nout,itrc)=0.0_r8
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
        IF (ad_LBC(iwest,isTvar(itrc),ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                IF (LBC_apply(ng)%west(j)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                  IF (ad_LBC(iwest,isTvar(itrc),ng)%nudging) THEN
                    IF (LnudgeTCLM(itrc,ng)) THEN
                      obc_out=CLIMA(ng)%Tnudgcof(Istr-1,j,k,ic)
                      obc_in =obcfac(ng)*obc_out
                    ELSE
                      obc_out=Tobc_out(itrc,ng,iwest)
                      obc_in =Tobc_in (itrc,ng,iwest)
                    END IF
                    IF (BOUNDARY(ng)%t_west_Cx(j,k,itrc).lt.            &
     &                  0.0_r8) THEN
                      tau=obc_in
                    ELSE
                      tau=obc_out
                    END IF
                    tau=tau*dt(ng)
                  END IF
                  Cx=BOUNDARY(ng)%t_west_Cx(j,k,itrc)
#  ifdef RADIATION_2D
                  Ce=BOUNDARY(ng)%t_west_Ce(j,k,itrc)
#  else
                  Ce=0.0_r8
#  endif
                  cff=BOUNDARY(ng)%t_west_C2(j,k,itrc)
# endif
# ifdef MASKING
!>                tl_t(Istr-1,j,k,nout,itrc)=tl_t(Istr-1,j,k,nout,      &
!>   &                                            itrc)*                &
!>   &                                       GRID(ng)%rmask(Istr-1,j)
!>
                  ad_t(Istr-1,j,k,nout,itrc)=ad_t(Istr-1,j,k,nout,      &
     &                                            itrc)*                &
     &                                       GRID(ng)%rmask(Istr-1,j)
# endif
                  IF (ad_LBC(iwest,isTvar(itrc),ng)%nudging) THEN
!>                  tl_t(Istr-1,j,k,nout,itrc)=tl_t(Istr-1,j,k,nout,    &
!>                                                  itrc)-              &
!>   &                                         tau*                     &
!>   &                                         tl_t(Istr-1,j,k,nstp,    &
!>   &                                              itrc)
!>
                    ad_t(Istr-1,j,k,nstp,itrc)=ad_t(Istr-1,j,k,nstp,    &
     &                                              itrc)-              &
     &                                         tau*                     &
     &                                         ad_t(Istr-1,j,k,nout,    &
     &                                              itrc)
                  END IF
!>                tl_t(Istr-1,j,k,nout,itrc)=(cff*                      &
!>   &                                        tl_t(Istr-1,j,k,nstp,     &
!>   &                                             itrc)+               &
!>   &                                        Cx*                       &
!>   &                                        tl_t(Istr  ,j,k,nout,     &
!>   &                                             itrc)-               &
!>   &                                        MAX(Ce,0.0_r8)*           &
!>   &                                           tl_grad(Istr-1,j  )-   &
!>   &                                        MIN(Ce,0.0_r8)*           &
!>   &                                           tl_grad(Istr-1,j+1))/  &
!>   &                                       (cff+Cx)
!>
                  adfac=ad_t(Istr-1,j,k,nout,itrc)/(cff+Cx)
                  ad_grad(Istr-1,j  )=ad_grad(Istr-1,j  )-              &
     &                                MAX(Ce,0.0_r8)*adfac
                  ad_grad(Istr-1,j+1)=ad_grad(Istr-1,j+1)-              &
     &                                MIN(Ce,0.0_r8)*adfac
                  ad_t(Istr-1,j,k,nstp,itrc)=ad_t(Istr-1,j,k,nstp,      &
     &                                            itrc)+                &
     &                                       cff*adfac
                  ad_t(Istr  ,j,k,nout,itrc)=ad_t(Istr  ,j,k,nout,      &
     &                                            itrc)+                &
     &                                       Cx *adfac
                  ad_t(Istr-1,j,k,nout,itrc)=0.0_r8
                END IF
              END DO
            END DO
          END IF
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (ad_LBC(iwest,isTvar(itrc),ng)%clamped) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
# ifdef MASKING
!>              tl_t(Istr-1,j,k,nout,itrc)=tl_t(Istr-1,j,k,nout,itrc)*  &
!>   &                                     GRID(ng)%rmask(Istr-1,j)
!>
                ad_t(Istr-1,j,k,nout,itrc)=ad_t(Istr-1,j,k,nout,itrc)*  &
     &                                     GRID(ng)%rmask(Istr-1,j)
# endif
# ifdef ADJUST_BOUNDARY
                IF (Lobc(iwest,isTvar(itrc),ng)) THEN
!>                tl_t(Istr-1,j,k,nout,itrc)=BOUNDARY(ng)%              &
!>   &                                             tl_t_west(j,k,itrc)
!>
                  BOUNDARY(ng)%ad_t_west(j,k,itrc)=                     &
     &                              BOUNDARY(ng)%ad_t_west(j,k,itrc)+   &
     &                                             ad_t(Istr-1,j,k,     &
     &                                                  nout,itrc)
                  ad_t(Istr-1,j,k,nout,itrc)=0.0_r8
                ELSE
!>                tl_t(Istr-1,j,k,nout,itrc)=0.0_r8
!>
                  ad_t(Istr-1,j,k,nout,itrc)=0.0_r8
                END IF
# else
!>              tl_t(Istr-1,j,k,nout,itrc)=0.0_r8
!>
                ad_t(Istr-1,j,k,nout,itrc)=0.0_r8
# endif
              END IF
            END DO
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (ad_LBC(iwest,isTvar(itrc),ng)%gradient) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
# ifdef MASKING
!>              tl_t(Istr-1,j,k,nout,itrc)=tl_t(Istr-1,j,k,nout,itrc)*  &
!>   &                                     GRID(ng)%rmask(Istr-1,j)
!>
                ad_t(Istr-1,j,k,nout,itrc)=ad_t(Istr-1,j,k,nout,itrc)*  &
     &                                     GRID(ng)%rmask(Istr-1,j)
# endif
!>              tl_t(Istr-1,j,k,nout,itrc)=tl_t(Istr,j,k,nout,itrc)
!>
                ad_t(Istr  ,j,k,nout,itrc)=ad_t(Istr  ,j,k,nout,itrc)+  &
     &                                     ad_t(Istr-1,j,k,nout,itrc)
                ad_t(Istr-1,j,k,nout,itrc)=0.0_r8
              END IF
            END DO
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (ad_LBC(iwest,isTvar(itrc),ng)%closed) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
# ifdef MASKING
!>              tl_t(Istr-1,j,k,nout,itrc)=tl_t(Istr-1,j,k,nout,itrc)*  &
!>   &                                     GRID(ng)%rmask(Istr-1,j)
!>
                ad_t(Istr-1,j,k,nout,itrc)=ad_t(Istr-1,j,k,nout,itrc)*  &
     &                                     GRID(ng)%rmask(Istr-1,j)
# endif
!>              tl_t(Istr-1,j,k,nout,itrc)=tl_t(Istr,j,k,nout,itrc)
!>
                ad_t(Istr  ,j,k,nout,itrc)=ad_t(Istr  ,j,k,nout,itrc)+  &
     &                                     ad_t(Istr-1,j,k,nout,itrc)
                ad_t(Istr-1,j,k,nout,itrc)=0.0_r8
              END IF
            END DO
          END DO
        END IF
      END IF

      RETURN
      END SUBROUTINE ad_t3dbc_tile
#endif
      END MODULE ad_t3dbc_mod
