#include "cppdefs.h"
      MODULE ad_u3dbc_mod
#if defined ADJOINT && defined SOLVE3D
!
!svn $Id: ad_u3dbc_im.F 694 2008-08-08 18:33:05Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine sets adjoint lateral boundary conditions for total  !
!  3D U-velocity. It updates the specified "nout" time index.          !
!                                                                      !
!  BASIC STATE variables needed: u                                     !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ad_u3dbc, ad_u3dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_u3dbc (ng, tile, nout)
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
      CALL ad_u3dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp(ng), nout,                               &
     &                    OCEAN(ng) % ad_u)
      RETURN
      END SUBROUTINE ad_u3dbc

!
!***********************************************************************
      SUBROUTINE ad_u3dbc_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, UBk,                &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nstp, nout,                             &
     &                          ad_u)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
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
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,UBk,2)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8) :: Ce, Cx, cff, tau
      real(r8) :: adfac

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_grad

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_grad(LBi:UBi,LBj:UBj)=0.0_r8

# if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,N(ng)
!>        tl_u(Iend+1,Jend+1,k,nout)=0.5_r8*                            &
!>   &                               (tl_u(Iend+1,Jend  ,k,nout)+       &
!>   &                                tl_u(Iend  ,Jend+1,k,nout))
!>
          adfac=0.5_r8*ad_u(Iend+1,Jend+1,k,nout)
          ad_u(Iend+1,Jend  ,k,nout)=ad_u(Iend+1,Jend  ,k,nout)+adfac
          ad_u(Iend  ,Jend+1,k,nout)=ad_u(Iend  ,Jend+1,k,nout)+adfac
          ad_u(Iend+1,Jend+1,k,nout)=0.0_r8
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=1,N(ng)
!>        tl_u(Istr,Jend+1,k,nout)=0.5_r8*                              &
!>   &                             (tl_u(Istr  ,Jend  ,k,nout)+         &
!>   &                              tl_u(Istr+1,Jend+1,k,nout))
!>
          adfac=0.5_r8*ad_u(Istr,Jend+1,k,nout)
          ad_u(Istr  ,Jend  ,k,nout)=ad_u(Istr  ,Jend  ,k,nout)+adfac
          ad_u(Istr+1,Jend+1,k,nout)=ad_u(Istr+1,Jend+1,k,nout)+adfac
          ad_u(Istr  ,Jend+1,k,nout)=0.0_r8
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,N(ng)
!>        tl_u(Iend+1,Jstr-1,k,nout)=0.5_r8*                            &
!>   &                               (tl_u(Iend  ,Jstr-1,k,nout)+       &
!>   &                                tl_u(Iend+1,Jstr  ,k,nout))
!>
          adfac=0.5_r8*ad_u(Iend+1,Jstr-1,k,nout)
          ad_u(Iend  ,Jstr-1,k,nout)=ad_u(Iend  ,Jstr-1,k,nout)+adfac
          ad_u(Iend+1,Jstr  ,k,nout)=ad_u(Iend+1,Jstr  ,k,nout)+adfac
          ad_u(Iend+1,Jstr-1,k,nout)=0.0_r8
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=1,N(ng)
!>        tl_u(Istr,Jstr-1,k,nout)=0.5_r8*                              &
!>   &                             (tl_u(Istr+1,Jstr-1,k,nout)+         &
!>   &                              tl_u(Istr  ,Jstr  ,k,nout))
!>
          adfac=0.5_r8*ad_u(Istr,Jstr-1,k,nout)
          ad_u(Istr+1,Jstr-1,k,nout)=ad_u(Istr+1,Jstr-1,k,nout)+adfac
          ad_u(Istr  ,Jstr  ,k,nout)=ad_u(Istr  ,Jstr  ,k,nout)+adfac
          ad_u(Istr  ,Jstr-1,k,nout)=0.0_r8
        END DO
      END IF
# endif

# ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (NORTHERN_EDGE) THEN

#  if defined NORTH_M3RADIATION

        IF (iic(ng).ne.0) THEN
!
!  Northern edge, implicit upstream radiation condition.
!
          DO k=1,N(ng)
            DO i=IstrU,Iend
#   ifdef NORTH_M3NUDGING
              IF (BOUNDARY(ng)%u_north_Ce(i,k).eq.0.0_r8) THEN
                tau=M3obc_in(ng,inorth)
              ELSE
                tau=M3obc_out(ng,inorth)
              END IF
              tau=tau*dt(ng)
#   endif
#   ifdef RADIATION_2D
              Cx=BOUNDARY(ng)%u_north_Cx(i,k)
#   else
              Cx=0.0_r8
#   endif
              Ce=BOUNDARY(ng)%u_north_Ce(i,k)
              cff=BOUNDARY(ng)%u_north_C2(i,k)
#   ifdef MASKING
!>            tl_u(i,Jend+1,k,nout)=tl_u(i,Jend+1,k,nout)*              &
!>   &                              GRID(ng)%umask(i,Jend+1)
!>
              ad_u(i,Jend+1,k,nout)=ad_u(i,Jend+1,k,nout)*              &
     &                              GRID(ng)%umask(i,Jend+1)
#   endif
#   ifdef NORTH_M3NUDGING
!>            tl_u(i,Jend+1,k,nout)=tl_u(i,Jend+1,k,nout)-              &
!>   &                              tau*tl_u(i,Jend+1,k,nstp)
!>
              ad_u(i,Jend+1,k,nstp)=ad_u(i,Jend+1,k,nstp)-              &
     &                              tau*ad_u(i,Jend+1,k,nout)
#    endif
!>            tl_u(i,Jend+1,k,nout)=(cff*tl_u(i,Jend+1,k,nstp)+         &
!>   &                               Ce *tl_u(i,Jend  ,k,nout)-         &
!>   &                               MAX(Cx,0.0_r8)*                    &
!>   &                                  tl_grad(i-1,Jend+1)-            &
!>   &                               MIN(Cx,0.0_r8)*                    &
!>   &                                  tl_grad(i  ,Jend+1))/           &
!>   &                              (cff+Ce) 
!>
              adfac=ad_u(i,Jend+1,k,nout)/(cff+Ce)
              ad_grad(i-1,Jend+1)=ad_grad(i-1,Jend+1)-                  &
     &                            MAX(Cx,0.0_r8)*adfac
              ad_grad(i  ,Jend+1)=ad_grad(i  ,Jend+1)-                  &
     &                            MIN(Cx,0.0_r8)*adfac
              ad_u(i,Jend  ,k,nout)=ad_u(i,Jend  ,k,nout)+Ce *adfac
              ad_u(i,Jend+1,k,nstp)=ad_u(i,Jend+1,k,nstp)+cff*adfac
              ad_u(i,Jend+1,k,nout)=0.0_r8
            END DO
          END DO
        END IF

#  elif defined NORTH_M3CLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
!>          tl_u(i,Jend+1,k,nout)=0.0_r8
!>
            ad_u(i,Jend+1,k,nout)=0.0_r8
          END DO
        END DO

#  elif defined NORTH_M3GRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
#   ifdef MASKING
!>          tl_u(i,Jend+1,k,nout)=tl_u(i,Jend+1,k,nout)*                &
!>   &                            GRID(ng)%umask(i,Jend+1)
!>
            ad_u(i,Jend+1,k,nout)=ad_u(i,Jend+1,k,nout)*                &
     &                            GRID(ng)%umask(i,Jend+1)
#   endif
!>          tl_u(i,Jend+1,k,nout)=tl_u(i,Jend,k,nout)
!>
            ad_u(i,Jend  ,k,nout)=ad_u(i,Jend  ,k,nout)+                &
     &                            ad_u(i,Jend+1,k,nout)
            ad_u(i,Jend+1,k,nout) = 0.0_r8
          END DO
        END DO

#  else
!
!  Northern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
#   ifdef EW_PERIODIC
#    define I_RANGE IstrU,Iend
#   else
#    define I_RANGE Istr,IendR
#   endif
        DO k=1,N(ng)
          DO i=I_RANGE
#   ifdef MASKING
!>          tl_u(i,Jend+1,k,nout)=tl_u(i,Jend+1,k,nout)*                &
!>   &                            GRID(ng)%umask(i,Jend+1)
!>
            ad_u(i,Jend+1,k,nout)=ad_u(i,Jend+1,k,nout)*                &
     &                            GRID(ng)%umask(i,Jend+1)
#   endif
!>          tl_u(i,Jend+1,k,nout)=gamma2(ng)*tl_u(i,Jend,k,nout)
!>
            ad_u(i,Jend  ,k,nout)=ad_u(i,Jend  ,k,nout)+                &
     &                            gamma2(ng)*ad_u(i,Jend+1,k,nout)
            ad_u(i,Jend+1,k,nout)=0.0_r8
          END DO
        END DO
#   undef I_RANGE
#  endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE) THEN

#  if defined SOUTH_M3RADIATION

        IF (iic(ng).ne.0) THEN
!
!  Southern edge, implicit upstream radiation condition.
!             
          DO k=1,N(ng)
            DO i=IstrU,Iend
#   ifdef SOUTH_M3NUDGING
              IF (BOUNDARY(ng)%u_south_Ce(i,k).eq.0.0_r8) THEN
                tau=M3obc_in(ng,isouth)
              ELSE
                tau=M3obc_out(ng,isouth)
              END IF
              tau=tau*dt(ng)
#   endif
#   ifdef RADIATION_2D
              Cx=BOUNDARY(ng)%u_south_Cx(i,k)
#   else
              Cx=0.0_r8
#   endif
              Ce=BOUNDARY(ng)%u_south_Ce(i,k)
              cff=BOUNDARY(ng)%u_south_C2(i,k)
#   ifdef MASKING
!>            tl_u(i,Jstr-1,k,nout)=tl_u(i,Jstr-1,k,nout)*              &
!>   &                              GRID(ng)%umask(i,Jstr-1)
!>
              ad_u(i,Jstr-1,k,nout)=ad_u(i,Jstr-1,k,nout)*              &
     &                              GRID(ng)%umask(i,Jstr-1)
#   endif
#   ifdef SOUTH_M3NUDGING
!>            tl_u(i,Jstr-1,k,nout)=tl_u(i,Jstr-1,k,nout)-              &
!>   &                              tau*tl_u(i,Jstr-1,k,nstp)
!>
              ad_u(i,Jstr-1,k,nstp)=ad_u(i,Jstr-1,k,nstp)-              &
     &                              tau*ad_u(i,Jstr-1,k,nout)
#   endif
!>            tl_u(i,Jstr-1,k,nout)=(cff*tl_u(i,Jstr-1,k,nstp)+         &
!>   &                               Ce *tl_u(i,Jstr  ,k,nout)-         &
!>   &                               MAX(Cx,0.0_r8)*                    &
!>   &                                  tl_grad(i-1,Jstr-1)-            &
!>   &                               MIN(Cx,0.0_r8)*                    &
!>   &                                  tl_grad(i  ,Jstr-1))/           &
!>   &                              (cff+Ce) 
!>
              adfac=ad_u(i,Jstr-1,k,nout)/(cff+Ce)
              ad_grad(i-1,Jstr-1)=ad_grad(i-1,Jstr-1)-                  &
     &                            MAX(Cx,0.0_r8)*adfac
              ad_grad(i  ,Jstr-1)=ad_grad(i  ,Jstr-1)-                  &
     &                            MIN(Cx,0.0_r8)*adfac
              ad_u(i,Jstr-1,k,nstp)=ad_u(i,Jstr-1,k,nstp)+cff*adfac
              ad_u(i,Jstr  ,k,nout)=ad_u(i,Jstr  ,k,nout)+Ce *adfac
              ad_u(i,Jstr-1,k,nout)=0.0_r8
            END DO
          END DO
        END IF

#  elif defined SOUTH_M3CLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
!>          tl_u(i,Jstr-1,k,nout)=0.0_r8
!>
            ad_u(i,Jstr-1,k,nout)=0.0_r8
          END DO
        END DO

#  elif defined SOUTH_M3GRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
#   ifdef MASKING
!>          tl_u(i,Jstr-1,k,nout)=tl_u(i,Jstr-1,k,nout)*                &
!>   &                            GRID(ng)%umask(i,Jstr-1)
!>
            ad_u(i,Jstr-1,k,nout)=ad_u(i,Jstr-1,k,nout)*                &
     &                            GRID(ng)%umask(i,Jstr-1)
#   endif
!>          tl_u(i,Jstr-1,k,nout)=tl_u(i,Jstr,k,nout)
!>
            ad_u(i,Jstr  ,k,nout)=ad_u(i,Jstr  ,k,nout)+                &
     &                            ad_u(i,Jstr-1,k,nout)
            ad_u(i,Jstr-1,k,nout)=0.0_r8
          END DO
        END DO

#  else
!
!  Southern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
#   ifdef EW_PERIODIC
#    define I_RANGE IstrU,Iend
#   else
#    define I_RANGE Istr,IendR
#   endif
        DO k=1,N(ng)
          DO i=I_RANGE
#   ifdef MASKING
!>          tl_u(i,Jstr-1,k,nout)=tl_u(i,Jstr-1,k,nout)*                &
!>   &                            GRID(ng)%umask(i,Jstr-1)
!>
            ad_u(i,Jstr-1,k,nout)=ad_u(i,Jstr-1,k,nout)*                &
     &                            GRID(ng)%umask(i,Jstr-1)
#   endif
!>          tl_u(i,Jstr-1,k,nout)=gamma2(ng)*tl_u(i,Jstr,k,nout)
!>
            ad_u(i,Jstr  ,k,nout)=ad_u(i,Jstr  ,k,nout)+                &
     &                            gamma2(ng)*ad_u(i,Jstr-1,k,nout)
            ad_u(i,Jstr-1,k,nout)=0.0_r8
          END DO
        END DO
#   undef I_RANGE
#  endif
      END IF
# endif
# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN

#  if defined EAST_M3RADIATION

        IF (iic(ng).ne.0) THEN
!
!  Eastern edge, implicit upstream radiation condition.
! 
          DO k=1,N(ng)
            DO j=Jstr,Jend
#   ifdef EAST_M3NUDGING
              IF (BOUNDARY(ng)%u_east_Cx(j,k).eq.0.0_r8) THEN
                tau=M3obc_in(ng,ieast)
              ELSE
                tau=M3obc_out(ng,ieast)
              END IF
              tau=tau*dt(ng)
#   endif
              Cx=BOUNDARY(ng)%u_east_Cx(j,k)
#   ifdef RADIATION_2D
              Ce=BOUNDARY(ng)%u_east_Ce(j,k)
#   else
              Ce=0.0_r8
#   endif
              cff=BOUNDARY(ng)%u_east_C2(j,k)
#   ifdef MASKING
!>            tl_u(Iend+1,j,k,nout)=tl_u(Iend+1,j,k,nout)*              &
!>   &                              GRID(ng)%umask(Iend+1,j)
!>
              ad_u(Iend+1,j,k,nout)=ad_u(Iend+1,j,k,nout)*              &
     &                              GRID(ng)%umask(Iend+1,j)
#   endif
#   ifdef EAST_M3NUDGING
!>            tl_u(Iend+1,j,k,nout)=tl_u(Iend+1,j,k,nout)-              &
!>   &                              tau*tl_u(Iend+1,j,k,nstp)
!>
              ad_u(Iend+1,j,k,nstp)=ad_u(Iend+1,j,k,nstp)-              &
     &                              tau*ad_u(Iend+1,j,k,nout)
#   endif
!>            tl_u(Iend+1,j,k,nout)=(cff*tl_u(Iend+1,j,k,nstp)+         &
!>   &                               Cx *tl_u(Iend  ,j,k,nout)-         &
!>   &                               MAX(Ce,0.0_r8)*                    &
!>   &                                  tl_grad(Iend+1,j  )-            &
!>   &                               MIN(Ce,0.0_r8)*                    &
!>   &                                  tl_grad(Iend+1,j+1))/           &
!>   &                              (cff+Cx) 
!>
              adfac=ad_u(Iend+1,j,k,nout)/(cff+Cx)
              ad_grad(Iend+1,j  )=ad_grad(Iend+1,j  )-                  &
     &                            MAX(Ce,0.0_r8)*adfac
              ad_grad(Iend+1,j+1)=ad_grad(Iend+1,j+1)-                  &
     &                            MIN(Ce,0.0_r8)*adfac
              ad_u(Iend  ,j,k,nout)=ad_u(Iend  ,j,k,nout)+Cx *adfac
              ad_u(Iend+1,j,k,nstp)=ad_u(Iend+1,j,k,nstp)+cff*adfac
              ad_u(Iend+1,j,k,nout)=0.0_r8
            END DO
          END DO
        END IF

#  elif defined EAST_M3CLAMPED
!
!  Eastern edge, clamped boundary condition.
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
!>        tl_u(Iend+1,j,k,nout)=0.0_r8
!>
          ad_u(Iend+1,j,k,nout)=0.0_r8
        END DO
      END DO

#  elif defined EAST_M3GRADIENT
!
!  Eastern edge, gradient boundary condition.
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
#   ifdef MASKING
!>        tl_u(Iend+1,j,k,nout)=tl_u(Iend+1,j,k,nout)*                  &
!>   &                          GRID(ng)%umask(Iend+1,j)
!>
          ad_u(Iend+1,j,k,nout)=ad_u(Iend+1,j,k,nout)*                  &
     &                          GRID(ng)%umask(Iend+1,j)
#   endif
!>        tl_u(Iend+1,j,k,nout)=tl_u(Iend,j,k,nout)
!>
          ad_u(Iend  ,j,k,nout)=ad_u(Iend  ,j,k,nout)+                  &
     &                          ad_u(Iend+1,j,k,nout)
          ad_u(Iend+1,j,k,nout)=0.0_r8
        END DO
      END DO

#  else
!
!  Eastern edge, closed boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
!>          tl_u(Iend+1,j,k,nout)=0.0_r8
!>
            ad_u(Iend+1,j,k,nout)=0.0_r8
          END DO
        END DO
#  endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN

#  if defined WEST_M3RADIATION

        IF (iic(ng).ne.0) THEN
!
!  Western edge, implicit upstream radiation condition.
!
          DO k=1,N(ng)
            DO j=Jstr,Jend
#   ifdef WEST_M2NUDGING
              IF (BOUNDARY(ng)%u_west_Cx(j,k).eq.0.0_r8) THEN
     &          tau=M3obc_in(ng,iwest)
              ELSE
                tau=M3obc_out(ng,iwest)
              END IF
              tau=tau*dt(ng)
#   endif
              Cx=BOUNDARY(ng)%u_west_Cx(j,k)
#   ifdef RADIATION_2D
              Ce=BOUNDARY(ng)%u_west_Ce(j,k)
#   else
              Ce=0.0_r8
#   endif
              cff=BOUNDARY(ng)%u_west_C2(j,k)
#   ifdef MASKING
!>            tl_u(Istr,j,k,nout)=tl_u(Istr,j,k,nout)*                  &
!>   &                            GRID(ng)%umask(Istr,j)
!>
              ad_u(Istr,j,k,nout)=ad_u(Istr,j,k,nout)*                  &
     &                            GRID(ng)%umask(Istr,j)
#   endif
#   ifdef WEST_M3NUDGING
!>            tl_u(Istr,j,k,nout)=tl_u(Istr,j,k,nout)-                  &
!>   &                            tau*tl_u(Istr,j,k,nstp)
!>
              ad_u(Istr,j,k,nstp)=ad_u(Istr,j,k,nstp)-                  &
     &                            tau*ad_u(Istr,j,k,nout)
#   endif
!>            tl_u(Istr,j,k,nout)=(cff*tl_u(Istr  ,j,k,nstp)+           &
!>   &                             Cx *tl_u(Istr+1,j,k,nout)-           &
!>   &                             MAX(Ce,0.0_r8)*                      &
!>   &                                tl_grad(Istr,j  )-                &
!>   &                             MIN(Ce,0.0_r8)*                      &
!>   &                                tl_grad(Istr,j+1))/               &
!>   &                            (cff+Cx) 
!>
              adfac=ad_u(Istr,j,k,nout)/(cff+Cx)
              ad_grad(Istr,j  )=ad_grad(Istr,j  )-                      &
     &                          MAX(Ce,0.0_r8)*adfac
              ad_grad(Istr,j+1)=ad_grad(Istr,j+1)-                      &
     &                          MIN(Ce,0.0_r8)*adfac
              ad_u(Istr  ,j,k,nstp)=ad_u(Istr  ,j,k,nstp)+cff*adfac
              ad_u(Istr+1,j,k,nout)=ad_u(Istr+1,j,k,nout)+Cx *adfac
              ad_u(Istr  ,j,k,nout)=0.0_r8
            END DO
          END DO
        END IF

#  elif defined WEST_M3CLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
!>          tl_u(Istr,j,k,nout)=0.0_r8
!>
            ad_u(Istr,j,k,nout)=0.0_r8
          END DO
        END DO

#  elif defined WEST_M3GRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
#   ifdef MASKING
!>          tl_u(Istr,j,k,nout)=tl_u(Istr,j,k,nout)*                    &
!>   &                          GRID(ng)%umask(Istr,j)
!>
            ad_u(Istr  ,j,k,nout)=ad_u(Istr  ,j,k,nout)*                &
     &                            GRID(ng)%umask(Istr,j)
#   endif
!>          tl_u(Istr,j,k,nout)=tl_u(Istr+1,j,k,nout)
!>
            ad_u(Istr+1,j,k,nout)=ad_u(Istr+1,j,k,nout)+                &
     &                            ad_u(Istr  ,j,k,nout)
            ad_u(Istr  ,j,k,nout)=0.0_r8
          END DO
        END DO

#  else
!
!  Western edge, closed boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
!>          tl_u(Istr,j,k,nout)=0.0_r8
!>
            ad_u(Istr,j,k,nout)=0.0_r8
          END DO
        END DO
#  endif
      END IF 
# endif
      RETURN
      END SUBROUTINE ad_u3dbc_tile
#endif
      END MODULE ad_u3dbc_mod
