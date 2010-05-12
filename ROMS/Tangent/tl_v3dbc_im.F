#include "cppdefs.h"
      MODULE tl_v3dbc_mod
#if defined TANGENT && defined SOLVE3D
!
!svn $Id: tl_v3dbc_im.F 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group       Andrew M. Moore   !
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

      PRIVATE
      PUBLIC  :: tl_v3dbc, tl_v3dbc_tile

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
      integer :: i, j, k

      real(r8) :: Ce, Cx, cff, tau

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_grad

# include "set_bounds.h"

# ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE) THEN

#  if defined SOUTH_M3RADIATION_NOT_YET

        IF (iic(ng).ne.0) THEN
!
!  Southern edge, implicit upstream radiation condition.
!
          DO k=1,N(ng)
            DO i=Istr,Iend+1
!>            grad(i,Jstr)=v(i  ,Jstr,k,nstp)-                          &
!>   &                     v(i-1,Jstr,k,nstp)
!>
              tl_grad(i,Jstr)=0.0_r8
            END DO
            DO i=Istr,Iend
#   ifdef SOUTH_M3NUDGING
              IF (BOUNDARY(ng)%v_south_Ce(i,k).eq.0.0_r8) THEN
                tau=M3obc_in(ng,isouth)
              ELSE
                tau=M3obc_out(ng,isouth)
              END IF
              tau=tau*dt(ng)
#   endif
#   ifdef RADIATION_2D
              Cx=BOUNDARY(ng)%v_south_Cx(i,k)
#   else
              Cx=0.0_r8
#   endif
              Ce=BOUNDARY(ng)%v_south_Ce(i,k)
              cff=BOUNDARY(ng)%v_south_C2(i,k)
!>            v(i,Jstr,k,nout)=(cff*v(i,Jstr  ,k,nstp)+                 &
!>   &                          Ce *v(i,Jstr+1,k,nout)-                 &
!>   &                          MAX(Cx,0.0_r8)*grad(i  ,Jstr)-          &
!>   &                          MIN(Cx,0.0_r8)*grad(i+1,Jstr))/         &
!>   &                         (cff+Ce)
!>
              tl_v(i,Jstr,k,nout)=(cff*tl_v(i,Jstr  ,k,nstp)+           &
     &                             Ce *tl_v(i,Jstr+1,k,nout)-           &
     &                             MAX(Cx,0.0_r8)*                      &
     &                                tl_grad(i  ,Jstr)-                &
     &                             MIN(Cx,0.0_r8)*                      &
     &                                tl_grad(i+1,Jstr))/               &
     &                            (cff+Ce)
#   ifdef SOUTH_M3NUDGING
!>            v(i,Jstr,k,nout)=v(i,Jstr,k,nout)+                        &
!>   &                         tau*(BOUNDARY(ng)%v_south(i,k)-          &
!>   &                              v(i,Jstr,k,nstp))
!>
              tl_v(i,Jstr,k,nout)=tl_v(i,Jstr,k,nout)-                  &
     &                            tau*tl_v(i,Jstr,k,nstp)
#   endif
#   ifdef MASKING
!>            v(i,Jstr,k,nout)=v(i,Jstr,k,nout)*                        &
!>   &                         GRID(ng)%vmask(i,Jstr)
!>
              tl_v(i,Jstr,k,nout)=tl_v(i,Jstr,k,nout)*                  &
     &                            GRID(ng)%vmask(i,Jstr)
#   endif
            END DO
          END DO
        END IF

#  elif defined SOUTH_M3CLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
!>          v(i,Jstr,k,nout)=BOUNDARY(ng)%v_south(i,k)
!>
#   ifdef ADJUST_BOUNDARY
            IF (Lobc(isouth,isVvel,ng)) THEN
              tl_v(i,Jstr,k,nout)=BOUNDARY(ng)%tl_v_south(i,k)
            ELSE
              tl_v(i,Jstr,k,nout)=0.0_r8
            END IF
#   else
            tl_v(i,Jstr,k,nout)=0.0_r8
#   endif
#   ifdef MASKING
!>          v(i,Jstr,k,nout)=v(i,Jstr,k,nout)*                          &
!>   &                       GRID(ng)%vmask(i,Jstr)
!>
            tl_v(i,Jstr,k,nout)=tl_v(i,Jstr,k,nout)*                    &
     &                          GRID(ng)%vmask(i,Jstr)
#   endif
          END DO
        END DO

#  elif defined SOUTH_M3GRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
!>          v(i,Jstr,k,nout)=v(i,Jstr+1,k,nout)
!>
            tl_v(i,Jstr,k,nout)=tl_v(i,Jstr+1,k,nout)
#   ifdef MASKING
!>          v(i,Jstr,k,nout)=v(i,Jstr,k,nout)*                          &
!>   &                       GRID(ng)%vmask(i,Jstr)
!>
            tl_v(i,Jstr,k,nout)=tl_v(i,Jstr,k,nout)*                    &
     &                          GRID(ng)%vmask(i,Jstr)
#   endif
          END DO
        END DO

#  else
!
!  Southern edge, closed boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
!>          v(i,Jstr,k,nout)=0.0_r8
!>
            tl_v(i,Jstr,k,nout)=0.0_r8
          END DO
        END DO
#  endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
       IF (NORTHERN_EDGE) THEN

#  if defined NORTH_M3RADIATION_NOT_YET

         IF (iic(ng).ne.0) THEN
!
!  Northern edge, implicit upstream radiation condition.
!
          DO k=1,N(ng)
            DO i=Istr,Iend+1
!>            grad(i,Jend+1)=v(i  ,Jend+1,k,nstp)-                      &
!>   &                       v(i-1,Jend+1,k,nstp)
!>
              tl_grad(i,Jend+1)=0.0_r8
            END DO
            DO i=Istr,Iend
#   ifdef NORTH_M3NUDGING
              IF (BOUNDARY(ng)%v_south_Ce(i,k).eq.0.0_r8) THEN
                tau=M3obc_in(ng,inorth)
              ELSE
                tau=M3obc_out(ng,inorth)
              END IF
              tau=tau*dt(ng)
#   endif
#   ifdef RADIATION_2D
              Cx=BOUNDARY(ng)%v_south_Cx(i,k)
#   else
              Cx=0.0_r8
#   endif
              Ce=BOUNDARY(ng)%v_south_Ce(i,k)
              cff=BOUNDARY(ng)%v_south_C2(i,k)
!>            v(i,Jend+1,k,nout)=(cff*v(i,Jend+1,k,nstp)+               &
!>   &                            Ce *v(i,Jend  ,k,nout)-               &
!>   &                            MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-      &
!>   &                            MIN(Cx,0.0_r8)*grad(i+1,Jend+1))/     &
!>   &                           (cff+Ce)
!>
              tl_v(i,Jend+1,k,nout)=(cff*tl_v(i,Jend+1,k,nstp)+         &
     &                               Ce *tl_v(i,Jend  ,k,nout)-         &
     &                               MAX(Cx,0.0_r8)*                    &
     &                                  tl_grad(i  ,Jend+1)-            &
     &                               MIN(Cx,0.0_r8)*                    &
     &                                  tl_grad(i+1,Jend+1))/           &
     &                              (cff+Ce)
#   ifdef NORTH_M3NUDGING
!>            v(i,Jend+1,k,nout)=v(i,Jend+1,k,nout)+                    &
!>   &                           tau*(BOUNDARY(ng)%v_north(i,k)-        &
!>   &                                v(i,Jend+1,k,nstp))
!>
              tl_v(i,Jend+1,k,nout)=tl_v(i,Jend+1,k,nout)-              &
     &                              tau*tl_v(i,Jend+1,k,nstp)
#   endif
#   ifdef MASKING
!>            v(i,Jend+1,k,nout)=v(i,Jend+1,k,nout)*                    &
!>   &                           GRID(ng)%vmask(i,Jend+1)
!>
              tl_v(i,Jend+1,k,nout)=tl_v(i,Jend+1,k,nout)*              &
     &                              GRID(ng)%vmask(i,Jend+1)
#   endif
            END DO
          END DO
        END IF

#  elif defined NORTH_M3CLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
!>          v(i,Jend+1,k,nout)=BOUNDARY(ng)%v_north(i,k)
!>
#   ifdef ADJUST_BOUNDARY
            IF (Lobc(inorth,isVvel,ng)) THEN
              tl_v(i,Jend+1,k,nout)=BOUNDARY(ng)%tl_v_north(i,k)
            ELSE
              tl_v(i,Jend+1,k,nout)=0.0_r8
            END IF
#   else
            tl_v(i,Jend+1,k,nout)=0.0_r8
#   endif
#   ifdef MASKING
!>          v(i,Jend+1,k,nout)=v(i,Jend+1,k,nout)*                      &
!>   &                         GRID(ng)%vmask(i,Jend+1)
!>
            tl_v(i,Jend+1,k,nout)=tl_v(i,Jend+1,k,nout)*                &
     &                            GRID(ng)%vmask(i,Jend+1)
#   endif
          END DO
        END DO

#  elif defined NORTH_M3GRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
!>          v(i,Jend+1,k,nout)=v(i,Jend,k,nout)
!>
            tl_v(i,Jend+1,k,nout)=tl_v(i,Jend,k,nout)
#   ifdef MASKING
!>          v(i,Jend+1,k,nout)=v(i,Jend+1,k,nout)*                      &
!>   &                         GRID(ng)%vmask(i,Jend+1)
!>
            tl_v(i,Jend+1,k,nout)=tl_v(i,Jend+1,k,nout)*                &
     &                            GRID(ng)%vmask(i,Jend+1)
#   endif
          END DO
        END DO
#  else
!
!  Northern edge, closed boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
!>          v(i,Jend+1,k,nout)=0.0_r8
!>
            tl_v(i,Jend+1,k,nout)=0.0_r8
          END DO
        END DO
#  endif
      END IF
# endif

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN

#  if defined WEST_M3RADIATION_NOT_YET

        IF (iic(ng).ne.0) THEN
!
!  Western edge, implicit upstream radiation condition.
!
          DO k=1,N(ng)
            DO j=JstrV-1,Jend
!>            grad(Istr-1,j)=v(Istr-1,j+1,k,nstp)-                      &
!>   &                       v(Istr-1,j  ,k,nstp)
!>
              tl_grad(Istr-1,j)=0.0_r8
            END DO
            DO j=JstrV,Jend
#   ifdef WEST_M3NUDGING
              IF (BOUNDARY(ng)%v_west_Cx(j,k).eq.0.0_r8) THEN
                tau=M3obc_in(ng,iwest)
              ELSE
                tau=M3obc_out(ng,iwest)
              END IF
              tau=tau*dt(ng)
#   endif
              Cx=BOUNDARY(ng)%v_west_Cx(j,k)
#   ifdef RADIATION_2D
              Ce=BOUNDARY(ng)%v_west_Ce(j,k)
#   else
              Ce=0.0_r8
#   endif
              cff=BOUNDARY(ng)%v_west_C2(j,k)
!>            v(Istr-1,j,k,nout)=(cff*v(Istr-1,j,k,nstp)+               &
!>   &                            Cx *v(Istr  ,j,k,nout)-               &
!>   &                            MAX(Ce,0.0_r8)*grad(Istr-1,j-1)-      &
!>   &                            MIN(Ce,0.0_r8)*grad(Istr-1,j  ))/     &
!>   &                           (cff+Cx)
!>
              tl_v(Istr-1,j,k,nout)=(cff*tl_v(Istr-1,j,k,nstp)+         &
     &                               Cx *tl_v(Istr  ,j,k,nout)-         &
     &                               MAX(Ce,0.0_r8)*                    &
     &                                  tl_grad(Istr-1,j-1)-            &
     &                               MIN(Ce,0.0_r8)*                    &
     &                                  tl_grad(Istr-1,j  ))/           &
     &                              (cff+Cx)
#   ifdef WEST_M3NUDGING
!>            v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)+                    &
!>   &                           tau*(BOUNDARY(ng)%v_west(j,k)-         &
!>   &                                v(Istr-1,j,k,nstp))
!>
              tl_v(Istr-1,j,k,nout)=tl_v(Istr-1,j,k,nout)-              &
     &                              tau*tl_v(Istr-1,j,k,nstp)
#    endif
#   ifdef MASKING
!>            v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                    &
!>   &                           GRID(ng)%vmask(Istr-1,j)
!>
              tl_v(Istr-1,j,k,nout)=tl_v(Istr-1,j,k,nout)*              &
     &                              GRID(ng)%vmask(Istr-1,j)
#   endif
            END DO
          END DO
        END IF

#  elif defined WEST_M3CLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO j=JstrV,Jend
!>          v(Istr-1,j,k,nout)=BOUNDARY(ng)%v_west(j,k)
!>
#   ifdef ADJUST_BOUNDARY
            IF (Lobc(iwest,isVvel,ng)) THEN
              tl_v(Istr-1,j,k,nout)=BOUNDARY(ng)%tl_v_west(j,k)
            ELSE
              tl_v(Istr-1,j,k,nout)=0.0_r8
            END IF
#   else
            tl_v(Istr-1,j,k,nout)=0.0_r8
#   endif
#   ifdef MASKING
!>          v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                      &
!>   &                         GRID(ng)%vmask(Istr-1,j)
!>
            tl_v(Istr-1,j,k,nout)=tl_v(Istr-1,j,k,nout)*                &
     &                            GRID(ng)%vmask(Istr-1,j)
#   endif
          END DO
        END DO

#  elif defined WEST_M3GRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO j=JstrV,Jend
!>          v(Istr-1,j,k,nout)=v(Istr,j,k,nout)
!>
            tl_v(Istr-1,j,k,nout)=tl_v(Istr,j,k,nout)
#   ifdef MASKING
!>          v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                      &
!>   &                         GRID(ng)%vmask(Istr-1,j)
!>
            tl_v(Istr-1,j,k,nout)=tl_v(Istr-1,j,k,nout)*                &
     &                            GRID(ng)%vmask(Istr-1,j)
#   endif
          END DO
        END DO

#  else
!
!  Western edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
#   ifdef NS_PERIODIC
#    define J_RANGE JstrV,Jend
#   else
#    define J_RANGE Jstr,JendR
#   endif
        DO k=1,N(ng)
          DO j=J_RANGE
!>          v(Istr-1,j,k,nout)=gamma2(ng)*v(Istr,j,k,nout)
!>
            tl_v(Istr-1,j,k,nout)=gamma2(ng)*tl_v(Istr,j,k,nout)
#   ifdef MASKING
!>          v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                      &
!>   &                         GRID(ng)%vmask(Istr-1,j)
!>
            tl_v(Istr-1,j,k,nout)=tl_v(Istr-1,j,k,nout)*                &
     &                            GRID(ng)%vmask(Istr-1,j)
#   endif
          END DO
        END DO
#   undef J_RANGE
#  endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN

#  if defined EAST_M3RADIATION_NOT_YET

        IF (iic(ng).ne.0) THEN
!
!  Eastern edge, implicit upstream radiation condition.
!
          DO k=1,N(ng)
            DO j=JstrV-1,Jend
!>            grad(Iend+1,j)=v(Iend+1,j+1,k,nstp)-                      &
!>   &                       v(Iend+1,j  ,k,nstp)
!>
              tl_grad(Iend+1,j)=0.0_r8
            END DO
            DO j=JstrV,Jend
#   ifdef EAST_M3NUDGING
              IF (BOUNDARY(ng)%v_east_Cx(j,k).eq.0.0_r8) THEN
                tau=M3obc_in(ng,ieast)
              ELSE
                tau=M3obc_out(ng,ieast)
              END IF
              tau=tau*dt(ng)
#   endif
              Cx=BOUNDARY(ng)%v_east_Cx(j,k)
#   ifdef RADIATION_2D
              Ce=BOUNDARY(ng)%v_east_Ce(j,k)
#   else
              Ce=0.0_r8
#   endif
              cff=BOUNDARY(ng)%v_east_C2(j,k)
!>            v(Iend+1,j,k,nout)=(cff*v(Iend+1,j,k,nstp)+               &
!>   &                            Cx *v(Iend  ,j,k,nout)-               &
!>   &                            MAX(Ce,0.0_r8)*grad(Iend+1,j-1)-      &
!>   &                            MIN(Ce,0.0_r8)*grad(Iend+1,j  ))/     &
!>   &                           (cff+Cx)
!>
              tl_v(Iend+1,j,k,nout)=(cff*tl_v(Iend+1,j,k,nstp)+         &
     &                               Cx *tl_v(Iend  ,j,k,nout)-         &
     &                               MAX(Ce,0.0_r8)*                    &
     &                                  tl_grad(Iend+1,j-1)-            &
     &                               MIN(Ce,0.0_r8)*                    &
     &                                  tl_grad(Iend+1,j  ))/           &
     &                              (cff+Cx)
#   ifdef EAST_M3NUDGING
!>            v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)+                    &
!>   &                           tau*(BOUNDARY(ng)%v_east(j,k)-         &
!>   &                                v(Iend+1,j,k,nstp))
!>
              tl_v(Iend+1,j,k,nout)=tl_v(Iend+1,j,k,nout)-              &
     &                              tau*tl_v(Iend+1,j,k,nstp)
#   endif
#   ifdef MASKING
!>            v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                    &
!>   &                           GRID(ng)%vmask(Iend+1,j)
!>
              tl_v(Iend+1,j,k,nout)=tl_v(Iend+1,j,k,nout)*              &
     &                              GRID(ng)%vmask(Iend+1,j)
#   endif
            END DO
          END DO
        END IF

#  elif defined EAST_M3CLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO j=JstrV,Jend
!>          v(Iend+1,j,k,nout)=BOUNDARY(ng)%v_east(j,k)
!>
#   ifdef ADJUST_BOUNDARY
            IF (Lobc(ieast,isVvel,ng)) THEN
              tl_v(Iend+1,j,k,nout)=BOUNDARY(ng)%tl_v_east(j,k)
            ELSE
              tl_v(Iend+1,j,k,nout)=0.0_r8
            END IF
#   else
            tl_v(Iend+1,j,k,nout)=0.0_r8
#   endif
#   ifdef MASKING
!>          v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                      &
!>   &                         GRID(ng)%vmask(Iend+1,j)
!>
            tl_v(Iend+1,j,k,nout)=tl_v(Iend+1,j,k,nout)*                &
     &                            GRID(ng)%vmask(Iend+1,j)
#   endif
          END DO
        END DO

#  elif defined EAST_M3GRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO j=JstrV,Jend
!>          v(Iend+1,j,k,nout)=v(Iend,j,k,nout)
!>
            tl_v(Iend+1,j,k,nout)=tl_v(Iend,j,k,nout)
#   ifdef MASKING
!>          v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                      &
!>   &                         GRID(ng)%vmask(Iend+1,j)
!>
            tl_v(Iend+1,j,k,nout)=tl_v(Iend+1,j,k,nout)*                &
     &                            GRID(ng)%vmask(Iend+1,j)
#   endif
          END DO
        END DO

#  else
!
!  Eastern edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
#   ifdef NS_PERIODIC
#    define J_RANGE JstrV,Jend
#   else
#    define J_RANGE Jstr,JendR
#   endif
        DO k=1,N(ng)
          DO j=J_RANGE
!>          v(Iend+1,j,k,nout)=gamma2(ng)*v(Iend,j,k,nout)
!>
            tl_v(Iend+1,j,k,nout)=gamma2(ng)*tl_v(Iend,j,k,nout)
#   ifdef MASKING
!>          v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                      &
!>   &                         GRID(ng)%vmask(Iend+1,j)
!>
            tl_v(Iend+1,j,k,nout)=tl_v(Iend+1,j,k,nout)*                &
     &                            GRID(ng)%vmask(Iend+1,j)
#   endif
          END DO
        END DO
#   undef J_RANGE
#  endif
      END IF
# endif

# if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=1,N(ng)
!>        v(Istr-1,Jstr,k,nout)=0.5_r8*(v(Istr  ,Jstr  ,k,nout)+        &
!>   &                                  v(Istr-1,Jstr+1,k,nout))
!>
          tl_v(Istr-1,Jstr,k,nout)=0.5_r8*                              &
     &                             (tl_v(Istr  ,Jstr  ,k,nout)+         &
     &                              tl_v(Istr-1,Jstr+1,k,nout))
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,N(ng)
!>        v(Iend+1,Jstr,k,nout)=0.5_r8*(v(Iend  ,Jstr  ,k,nout)+        &
!>   &                                  v(Iend+1,Jstr+1,k,nout))
!>
          tl_v(Iend+1,Jstr,k,nout)=0.5_r8*                              &
     &                             (tl_v(Iend  ,Jstr  ,k,nout)+         &
     &                              tl_v(Iend+1,Jstr+1,k,nout))
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=1,N(ng)
!>        v(Istr-1,Jend+1,k,nout)=0.5_r8*(v(Istr-1,Jend  ,k,nout)+      &
!>   &                                    v(Istr  ,Jend+1,k,nout))
!>
          tl_v(Istr-1,Jend+1,k,nout)=0.5_r8*                            &
     &                               (tl_v(Istr-1,Jend  ,k,nout)+       &
     &                                tl_v(Istr  ,Jend+1,k,nout))
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,N(ng)
!>        v(Iend+1,Jend+1,k,nout)=0.5_r8*(v(Iend+1,Jend  ,k,nout)+      &
!>   &                                    v(Iend  ,Jend+1,k,nout))
!>
          tl_v(Iend+1,Jend+1,k,nout)=0.5_r8*                            &
     &                               (tl_v(Iend+1,Jend  ,k,nout)+       &
     &                                tl_v(Iend  ,Jend+1,k,nout))
        END DO
      END IF
# endif
      RETURN
      END SUBROUTINE tl_v3dbc_tile
#endif
      END MODULE tl_v3dbc_mod
