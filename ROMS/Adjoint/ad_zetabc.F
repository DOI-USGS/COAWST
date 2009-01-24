#include "cppdefs.h"
      MODULE ad_zetabc_mod
#ifdef ADJOINT
!
!svn $Id: ad_zetabc.F 694 2008-08-08 18:33:05Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine sets adjoint lateral boundary conditions for        !
!  free-surface. It updates the specified "kout" index.                !
!                                                                      !
!  BASIC STATE variables needed: zeta                                  !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ad_zetabc, ad_zetabc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_zetabc (ng, tile, kout)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, kout
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL ad_zetabc_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     krhs(ng), kstp(ng), kout,                    &
     &                     OCEAN(ng) % zeta,                            &
     &                     OCEAN(ng) % ad_zeta)
      RETURN
      END SUBROUTINE ad_zetabc
!
!***********************************************************************
      SUBROUTINE ad_zetabc_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           krhs, kstp, kout,                      &
     &                           zeta, ad_zeta)
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
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: krhs, kstp, kout
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)

      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)

      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, know

      real(r8) :: Ce, Cx
      real(r8) :: cff, cff1, cff2, dt2d, tau

      real(r8) :: ad_Ce, ad_Cx
      real(r8) :: ad_cff1, ad_cff2
      real(r8) :: adfac

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_grad

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_Ce=0.0_r8
      ad_Cx=0.0_r8
      ad_cff1=0.0_r8
      ad_cff2=0.0_r8

      ad_grad(LBi:UBi,LBj:UBj)=0.0_r8
!
!-----------------------------------------------------------------------
!  Set time-indices
!-----------------------------------------------------------------------
!
      IF (FIRST_2D_STEP) THEN
        know=krhs
        dt2d=dtfast(ng)
      ELSE IF (PREDICTOR_2D_STEP(ng)) THEN
        know=krhs
        dt2d=2.0_r8*dtfast(ng)
      ELSE
        know=kstp
        dt2d=dtfast(ng)
      END IF

# if defined WET_DRY
!
!-----------------------------------------------------------------------
! Ensure that water level on boundary cells is above bed elevation.
!-----------------------------------------------------------------------
!
      cff=Dcrit(ng)-eps

#  if !defined EW_PERIODIC && !defined NS_PERIODIC
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        IF (zeta(Iend+1,Jend+1,kout).le.                                &
     &      (Dcrit(ng)-GRID(ng)%h(Iend+1,Jend+1))) THEN
!>        tl_zeta(Iend+1,Jend+1,kout)=-GRID(ng)%tl_h(Iend+1,Jend+1)
!>
          GRID(ng)%ad_h(Iend+1,Jend+1)=GRID(ng)%ad_h(Iend+1,Jend+1)-    &
     &                                 ad_zeta(Iend+1,Jend+1,kout)
          ad_zeta(Iend+1,Jend+1,kout)=0.0_r8
        END IF
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        IF (zeta(Istr-1,Jend+1,kout).le.                                &
     &      (Dcrit(ng)-GRID(ng)%h(Istr-1,Jend+1))) THEN
!>        tl_zeta(Istr-1,Jend+1,kout)=-GRID(ng)%tl_h(Istr-1,Jend+1)
!>
          GRID(ng)%ad_h(Istr-1,Jend+1)=GRID(ng)%ad_h(Istr-1,Jend+1)-    &
     &                                 ad_zeta(Istr-1,Jend+1,kout)
          ad_zeta(Istr-1,Jend+1,kout)=0.0_r8
        END IF
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        IF (zeta(Iend+1,Jstr-1,kout).le.                                &
     &      (Dcrit(ng)-GRID(ng)%h(Iend+1,Jstr-1))) THEN
!>        tl_zeta(Iend+1,Jstr-1,kout)=-GRID(ng)%tl_h(Iend+1,Jstr-1)
!>
          GRID(ng)%ad_h(Iend+1,Jstr-1)=GRID(ng)%ad_h(Iend+1,Jstr-1)-    &
     &                                 ad_zeta(Iend+1,Jstr-1,kout)
          tl_zeta(Iend+1,Jstr-1,kout)=0.0_r8
        END IF
      END IF
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        IF (zeta(Istr-1,Jstr-1,kout).le.                                &
     &      (Dcrit(ng)-GRID(ng)%h(Istr-1,Jstr-1))) THEN
!>        tl_zeta(Istr-1,Jstr-1,kout)=-GRID(ng)%tl_h(Istr-1,Jstr-1)
!>
          GRID(ng)%ad_h(Istr-1,Jstr-1)=GRID(ng)%ad_h(Istr-1,Jstr-1)-    &
     &                                 ad_zeta(Istr-1,Jstr-1,kout)
          ad_zeta(Istr-1,Jstr-1,kout)=0.0_r8
        END IF
      END IF
#  endif
#  ifndef NS_PERIODIC
      IF (NORTHERN_EDGE) THEN
        DO i=Istr,Iend
          IF (zeta(i,Jend+1,kout).le.                                   &
     &        (Dcrit(ng)-GRID(ng)%h(i,Jend+1))) THEN
!>          tl_zeta(i,Jend+1,kout)=-GRID(ng)%tl_h(i,Jend+1)
!>
            GRID(ng)%ad_h(i,Jend+1)=GRID(ng)%ad_h(i,Jend+1)-            &
     &                              ad_zeta(i,Jend+1,kout)
            ad_zeta(i,Jend+1,kout)=0.0_r8
          END IF
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO i=Istr,Iend
          IF (zeta(i,Jstr-1,kout).le.                                   &
     &        (Dcrit(ng)-GRID(ng)%h(i,Jstr-1))) THEN
!>          tl_zeta(i,Jstr-1,kout)=-GRID(ng)%tl_h(i,Jstr-1)
!>
            GRID(ng)%ad_h(i,Jstr-1)=GRID(ng)%ad_h(i,Jstr-1)-            &
     &                              ad_zeta(i,Jstr-1,kout)
            ad_zeta(i,Jstr-1,kout)=0.0_r8
          END IF
        END DO
      END IF
#  endif
#  ifndef EW_PERIODIC
      IF (EASTERN_EDGE) THEN
        DO j=Jstr,Jend
          IF (zeta(Iend+1,j,kout).le.                                   &
     &        (Dcrit(ng)-GRID(ng)%h(Iend+1,j))) THEN
!>          tl_zeta(Iend+1,j,kout)=-GRID(ng)%tl_h(Iend+1,j)
!>
            GRID(ng)%ad_h(Iend+1,j)=GRID(ng)%ad_h(Iend+1,j)-            &
     &                              ad_zeta(Iend+1,j,kout)
            ad_zeta(Iend+1,j,kout)=0.0_r8
          END IF
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO j=Jstr,Jend
          IF (zeta(Istr-1,j,kout).le.                                   &
     &        (Dcrit(ng)-GRID(ng)%h(Istr-1,j))) THEN
!>          tl_zeta(Istr-1,j,kout)=-GRID(ng)%tl_h(Istr-1,j)
!>
            GRID(ng)%ad_h(Istr-1,j)=GRID(ng)%ad_h(Istr-1,j)-            &
     &                              ad_zeta(Istr-1,j,kout)  
            ad_zeta(Istr-1,j,kout)=0.0_r8
          END IF
        END DO
      END IF
#  endif
# endif

# if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
!>      tl_zeta(Iend+1,Jend+1,kout)=0.5_r8*                             &
!>   &                              (tl_zeta(Iend+1,Jend  ,kout)+       &
!>   &                               tl_zeta(Iend  ,Jend+1,kout))
!>
        adfac=0.5_r8*ad_zeta(Iend+1,Jend+1,kout)
        ad_zeta(Iend  ,Jend+1,kout)=ad_zeta(Iend  ,Jend+1,kout)+adfac
        ad_zeta(Iend+1,Jend  ,kout)=ad_zeta(Iend+1,Jend  ,kout)+adfac
        ad_zeta(Iend+1,Jend+1,kout)=0.0_r8
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
!>      tl_zeta(Istr-1,Jend+1,kout)=0.5_r8*                             &
!>   &                              (tl_zeta(Istr-1,Jend  ,kout)+       &
!>   &                               tl_zeta(Istr  ,Jend+1,kout))
!>
        adfac=0.5_r8*ad_zeta(Istr-1,Jend+1,kout)
        ad_zeta(Istr-1,Jend  ,kout)=ad_zeta(Istr-1,Jend  ,kout)+adfac
        ad_zeta(Istr  ,Jend+1,kout)=ad_zeta(Istr  ,Jend+1,kout)+adfac
        ad_zeta(Istr-1,Jend+1,kout)=0.0_r8
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
!>      tl_zeta(Iend+1,Jstr-1,kout)=0.5_r8*                             &
!>   &                              (tl_zeta(Iend  ,Jstr-1,kout)+       &
!>   &                               tl_zeta(Iend+1,Jstr  ,kout))
!>
        adfac=0.5_r8*ad_zeta(Iend+1,Jstr-1,kout)
        ad_zeta(Iend  ,Jstr-1,kout)=ad_zeta(Iend  ,Jstr-1,kout)+adfac
        ad_zeta(Iend+1,Jstr  ,kout)=ad_zeta(Iend+1,Jstr  ,kout)+adfac
        ad_zeta(Iend+1,Jstr-1,kout)=0.0_r8
      END IF
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
!>      tl_zeta(Istr-1,Jstr-1,kout)=0.5_r8*                             &
!>   &                              (tl_zeta(Istr  ,Jstr-1,kout)+       &
!>   &                               tl_zeta(Istr-1,Jstr  ,kout))
!>
        adfac=0.5_r8*ad_zeta(Istr-1,Jstr-1,kout)
        ad_zeta(Istr  ,Jstr-1,kout)=ad_zeta(Istr  ,Jstr-1,kout)+adfac
        ad_zeta(Istr-1,Jstr  ,kout)=ad_zeta(Istr-1,Jstr  ,kout)+adfac
        ad_zeta(Istr-1,Jstr-1,kout)=0.0_r8
      END IF
# endif

# ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (NORTHERN_EDGE) THEN

#  if defined NORTH_FSRADIATION

        IF (iic(ng).ne.0) THEN
!
!  Northern edge, implicit upstream radiation condition.
!
          DO i=Istr,Iend
#   ifdef NORTH_FSNUDGING
            IF (BOUNDARY(ng)%zeta_north_Ce(i).eq.0.0_r8) THEN
              tau=FSobc_in(ng,inorth)
            ELSE
              tau=FSobc_out(ng,inorth)
            END IF
            tau=tau*dt2d
#   endif
#   ifdef RADIATION_2D
            Cx=BOUNDARY(ng)%zeta_north_Cx(i)
#   else
            Cx=0.0_r8
#   endif
            Ce=BOUNDARY(ng)%zeta_north_Ce(i)
            cff=BOUNDARY(ng)%zeta_north_C2(i)
#   ifdef MASKING
!>          tl_zeta(i,Jend+1,kout)=tl_zeta(i,Jend+1,kout)*              &
!>   &                             GRID(ng)%rmask(i,Jend+1)
!>
            ad_zeta(i,Jend+1,kout)=ad_zeta(i,Jend+1,kout)*              &
     &                             GRID(ng)%rmask(i,Jend+1)
#   endif
#   ifdef NORTH_FSNUDGING
!>          tl_zeta(i,Jend+1,kout)=tl_zeta(i,Jend+1,kout)-              &
!>   &                             tau*tl_zeta(i,Jend+1,know)
!>
            ad_zeta(i,Jend+1,know)=ad_zeta(i,Jend+1,know)-              &
     &                             tau*ad_zeta(i,Jend+1,kout)
#   endif
!>          tl_zeta(i,Jend+1,kout)=(cff*tl_zeta(i,Jend+1,know)+         &
!>   &                              Ce *tl_zeta(i,Jend  ,kout)-         &
!>   &                              MAX(Cx,0.0_r8)*tl_grad(i  ,Jend+1)- &
!>   &                              MIN(Cx,0.0_r8)*tl_grad(i+1,Jend+1))/&
!>   &                             (cff+Ce)
!>
            adfac=ad_zeta(i,Jend+1,kout)/(cff+Ce)
            ad_grad(i  ,Jend+1)=ad_grad(i  ,Jend+1)-MAX(Cx,0.0_r8)*adfac
            ad_grad(i+1,Jend+1)=ad_grad(i+1,Jend+1)-MIN(Cx,0.0_r8)*adfac
            ad_zeta(i,Jend  ,kout)=ad_zeta(i,Jend  ,kout)+Ce *adfac
            ad_zeta(i,Jend+1,know)=ad_zeta(i,Jend+1,know)+cff*adfac
            ad_zeta(i,Jend+1,kout)=0.0_r8
          END DO
        END IF

#  elif defined NORTH_FSCHAPMAN
!
!  Northern edge, Chapman boundary condition.
!
        DO i=Istr,Iend
          cff=dt2d*GRID(ng)%pn(i,Jend)
          cff1=SQRT(g*(GRID(ng)%h(i,Jend)+                              &
     &                 zeta(i,Jend,know)))
          Ce=cff*cff1
          cff2=1.0_r8/(1.0_r8+Ce)
#   ifdef MASKING
!>        tl_zeta(i,Jend+1,kout)=tl_zeta(i,Jend+1,kout)*                &
!>   &                           GRID(ng)%rmask(i,Jend+1)
!>
          ad_zeta(i,Jend+1,kout)=ad_zeta(i,Jend+1,kout)*                &
     &                           GRID(ng)%rmask(i,Jend+1)
#   endif
!>        tl_zeta(i,Jend+1,kout)=tl_cff2*(zeta(i,Jend+1,know)+          &
!>   &                                    Ce*zeta(i,Jend,kout))+        &
!>   &                           cff2*(tl_zeta(i,Jend+1,know)+          &
!>   &                                 tl_Ce*zeta(i,Jend,kout)+         &
!>   &                                 Ce*tl_zeta(i,Jend,kout))
!>
          adfac=cff2*ad_zeta(i,Jend+1,kout)
          ad_zeta(i,Jend  ,kout)=ad_zeta(i,Jend  ,kout)+Ce*adfac
          ad_zeta(i,Jend+1,know)=ad_zeta(i,Jend+1,know)+adfac
          ad_Ce=ad_Ce+zeta(i,Jend,kout)*adfac
          ad_cff2=ad_cff2+                                              &
     &            (zeta(i,Jend+1,know)+                                 &
     &             Ce*zeta(i,Jend,kout))*ad_zeta(i,Jend+1,kout)
          ad_zeta(i,Jend+1,kout)=0.0_r8
!>        tl_cff2=-cff2*cff2*tl_Ce
!>
          ad_Ce=ad_Ce-cff2*cff2*ad_cff2
          ad_cff2=0.0_r8
!>        tl_Ce=cff*tl_cff1
!>
          ad_cff1=ad_cff1+cff*ad_Ce
          ad_Ce=0.0_r8
!>        tl_cff1=0.5_r8*g*(GRID(ng)%tl_h(i,Jend)+                      &
!>   &                      tl_zeta(i,Jend,know))/cff1
!>
          adfac=0.5_r8*g*ad_cff1/cff1
          GRID(ng)%ad_h(i,Jend)=GRID(ng)%ad_h(i,Jend)+adfac
          ad_zeta(i,Jend,know)=ad_zeta(i,Jend,know)+adfac
          ad_cff1=0.0_r8
        END DO

#  elif defined NORTH_FSCLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=Istr,Iend
!>        tl_zeta(i,Jend+1,kout)=0.0_r8
!>
          ad_zeta(i,Jend+1,kout)=0.0_r8
        END DO

#  elif defined NORTH_FSGRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=Istr,Iend
#   ifdef MASKING
!>        tl_zeta(i,Jend+1,kout)=tl_zeta(i,Jend+1,kout)*                &
!>   &                           GRID(ng)%rmask(i,Jend+1)
!>
          ad_zeta(i,Jend+1,kout)=ad_zeta(i,Jend+1,kout)*                &
     &                           GRID(ng)%rmask(i,Jend+1)
#   endif
!>        tl_zeta(i,Jend+1,kout)=tl_zeta(i,Jend,kout)
!>
          ad_zeta(i,Jend  ,kout)=ad_zeta(i,Jend  ,kout)+                &
     &                           ad_zeta(i,Jend+1,kout)
          ad_zeta(i,Jend+1,kout)=0.0_r8
        END DO

#  else
!
!  Northern edge, closed boundary condition.
!
        DO i=Istr,Iend
#   ifdef MASKING
!>        tl_zeta(i,Jend+1,kout)=tl_zeta(i,Jend+1,kout)*                &
!>   &                           GRID(ng)%rmask(i,Jend+1)
!>
          ad_zeta(i,Jend+1,kout)=ad_zeta(i,Jend+1,kout)*                &
     &                           GRID(ng)%rmask(i,Jend+1)
#   endif
!>        tl_zeta(i,Jend+1,kout)=tl_zeta(i,Jend,kout)
!>
          ad_zeta(i,Jend  ,kout)=ad_zeta(i,Jend  ,kout)+                &
     &                           ad_zeta(i,Jend+1,kout)
          ad_zeta(i,Jend+1,kout)=0.0_r8
        END DO
#  endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE) THEN

#  if defined SOUTH_FSRADIATION

        IF (iic(ng).ne.0) THEN
!
!  Southern edge, implicit upstream radiation condition.
!
          DO i=Istr,Iend
#   ifdef SOUTH_FSNUDGING
            IF (BOUNDARY(ng)%zeta_south_Ce(i).eq.0.0_r8) THEN
              tau=FSobc_in(ng,isouth)
            ELSE
              tau=FSobc_out(ng,isouth)
            END IF
            tau=tau*dt2d
#   endif
#   ifdef RADIATION_2D
            Cx=BOUNDARY(ng)%zeta_south_Cx(i)
#   else
            Cx=0.0_r8
#   endif
            Ce=BOUNDARY(ng)%zeta_south_Ce(i)
            cff=BOUNDARY(ng)%zeta_south_C2(i)
#   ifdef MASKING
!>          tl_zeta(i,Jstr-1,kout)=tl_zeta(i,Jstr-1,kout)*              &
!>   &                             GRID(ng)%rmask(i,Jstr-1)
!>
            ad_zeta(i,Jstr-1,kout)=ad_zeta(i,Jstr-1,kout)*              &
     &                             GRID(ng)%rmask(i,Jstr-1)
#   endif
#   ifdef SOUTH_FSNUDGING
!>          tl_zeta(i,Jstr-1,kout)=tl_zeta(i,Jstr-1,kout)-              &
!>   &                             tau*tl_zeta(i,Jstr-1,know)
!>
            ad_zeta(i,Jstr-1,know)=ad_zeta(i,Jstr-1,know)-              &
     &                             tau*ad_zeta(i,Jstr-1,kout)
#   endif
!>          tl_zeta(i,Jstr-1,kout)=(cff*tl_zeta(i,Jstr-1,know)+         &
!>   &                              Ce *tl_zeta(i,Jstr  ,kout)-         &
!>   &                              MAX(Cx,0.0_r8)*tl_grad(i  ,Jstr-1)- &
!>   &                              MIN(Cx,0.0_r8)*tl_grad(i+1,Jstr-1))/&
!>   &                             (cff+Ce)
!>
            adfac=ad_zeta(i,Jstr-1,kout)/(cff+Ce)
            ad_grad(i  ,Jstr-1)=ad_grad(i  ,Jstr-1)-MAX(Cx,0.0_r8)*adfac
            ad_grad(i+1,Jstr-1)=ad_grad(i+1,Jstr-1)-MIN(Cx,0.0_r8)*adfac
            ad_zeta(i,Jstr-1,know)=ad_zeta(i,Jstr-1,know)+cff*adfac
            ad_zeta(i,Jstr  ,kout)=ad_zeta(i,Jstr  ,kout)+Ce *adfac
            ad_zeta(i,Jstr-1,kout)=0.0_r8
          END DO
        END IF

#  elif defined SOUTH_FSCHAPMAN
!
!  Southern edge, Chapman boundary condition.
!
        DO i=Istr,Iend
          cff=dt2d*GRID(ng)%pn(i,Jstr)
          cff1=SQRT(g*(GRID(ng)%h(i,Jstr)+                              &
     &                 zeta(i,Jstr,know)))
          Ce=cff*cff1
          cff2=1.0_r8/(1.0_r8+Ce)
#  ifdef MASKING
!>        tl_zeta(i,Jstr-1,kout)=tl_zeta(i,Jstr-1,kout)*                &
!>   &                           GRID(ng)%rmask(i,Jstr-1)
!>
          ad_zeta(i,Jstr-1,kout)=ad_zeta(i,Jstr-1,kout)*                &
     &                           GRID(ng)%rmask(i,Jstr-1)
#  endif
!>        tl_zeta(i,Jstr-1,kout)=tl_cff2*(zeta(i,Jstr-1,know)+          &
!>   &                                    Ce*zeta(i,Jstr,kout))+        &
!>   &                           cff2*(tl_zeta(i,Jstr-1,know)+          &
!>   &                                 tl_Ce*zeta(i,Jstr,kout)+         &
!>   &                                 Ce*tl_zeta(i,Jstr,kout))
!>
          adfac=cff2*ad_zeta(i,Jstr-1,kout)
          ad_zeta(i,Jstr-1,know)=ad_zeta(i,Jstr-1,know)+adfac
          ad_zeta(i,Jstr  ,kout)=ad_zeta(i,Jstr  ,kout)+Ce*adfac
          ad_Ce=ad_Ce+zeta(i,Jstr,kout)*adfac
          ad_cff2=ad_cff2+                                              &
     &            (zeta(i,Jstr-1,know)+                                 &
     &             Ce*zeta(i,Jstr,kout))*ad_zeta(i,Jstr-1,kout)
          ad_zeta(i,Jstr-1,kout)=0.0_r8
!>        tl_cff2=-cff2*cff2*tl_Ce
!>
          ad_Ce=ad_Ce-cff2*cff2*ad_cff2
          ad_cff2=0.0_r8
!>        tl_Ce=cff*tl_cff1
!>
          ad_cff1=ad_cff1+cff*ad_Ce
          ad_Ce=0.0_r8
!>        tl_cff1=0.5_r8*g*(GRID(ng)%tl_h(i,Jstr)+                      &
!>   &                      tl_zeta(i,Jstr,know))/cff1
!>
          adfac=0.5_r8*g*ad_cff1/cff1
          GRID(ng)%ad_h(i,Jstr)=GRID(ng)%ad_h(i,Jstr)+adfac
          ad_zeta(i,Jstr,know)=ad_zeta(i,Jstr,know)+adfac
          ad_cff1=0.0_r8
        END DO

#  elif defined SOUTH_FSCLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=Istr,Iend
!>        tl_zeta(i,Jstr-1,kout)=0.0_r8
!>
          ad_zeta(i,Jstr-1,kout)=0.0_r8
        END DO

#  elif defined SOUTH_FSGRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=Istr,Iend
#   ifdef MASKING
!>        tl_zeta(i,Jstr-1,kout)=tl_zeta(i,Jstr-1,kout)*                &
!>   &                           GRID(ng)%rmask(i,Jstr-1)
!>
          ad_zeta(i,Jstr-1,kout)=ad_zeta(i,Jstr-1,kout)*                &
     &                           GRID(ng)%rmask(i,Jstr-1)
#   endif
!>        tl_zeta(i,Jstr-1,kout)=tl_zeta(i,Jstr,kout)
!>
          ad_zeta(i,Jstr  ,kout)=ad_zeta(i,Jstr  ,kout)+                &
     &                           ad_zeta(i,Jstr-1,kout)
          ad_zeta(i,Jstr-1,kout)=0.0_r8
        END DO

#  else
!
!  Southern edge, closed boundary condition.
!
        DO i=Istr,Iend
#   ifdef MASKING
!>        tl_zeta(i,Jstr-1,kout)=tl_zeta(i,Jstr-1,kout)*                &
!>   &                           GRID(ng)%rmask(i,Jstr-1)
!>
          ad_zeta(i,Jstr-1,kout)=ad_zeta(i,Jstr-1,kout)*                &
     &                           GRID(ng)%rmask(i,Jstr-1)
#   endif
!>        tl_zeta(i,Jstr-1,kout)=tl_zeta(i,Jstr,kout)
!>
          ad_zeta(i,Jstr  ,kout)=ad_zeta(i,Jstr  ,kout)+                &
     &                           ad_zeta(i,Jstr-1,kout)
          ad_zeta(i,Jstr-1,kout)=0.0_r8
        END DO
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

#  if defined EAST_FSRADIATION

        IF (iic(ng).ne.0) THEN
!
!  Eastern edge, implicit upstream radiation condition.
!
          DO j=Jstr,Jend
#   ifdef EAST_FSNUDGING
            IF (BOUNDARY(ng)%zeta_east_Cx(j).eq.0.0_r8) THEN
              tau=FSobc_in(ieast)
            ELSE
              tau=FSobc_out(ieast)
            END IF
            tau=tau*dt2d
#   endif
            Cx=BOUNDARY(ng)%zeta_east_Cx(j)
#   ifdef RADIATION_2D
            Ce=BOUNDARY(ng)%zeta_east_Ce(j)
#   else
            Ce=0.0_r8
#   endif
            cff=BOUNDARY(ng)%zeta_east_C2(j)
#   ifdef MASKING
!>           tl_zeta(Iend+1,j,kout)=tl_zeta(Iend+1,j,kout)*             &
!>   &                              GRID(ng)%rmask(Iend+1
!>
             ad_zeta(Iend+1,j,kout)=ad_zeta(Iend+1,j,kout)*             &
     &                              GRID(ng)%rmask(Iend+1
#   endif
#   ifdef EAST_FSNUDGING
!>          tl_zeta(Iend+1,j,kout)=tl_zeta(Iend+1,j,kout)-              &
!>   &                             tau*tl_zeta(Iend+1,j,know)
!>
            ad_zeta(Iend+1,j,know)=ad_zeta(Iend+1,j,know)-              &
     &                              tau*ad_zeta(Iend+1,j,kout)
#   endif
!>          tl_zeta(Iend+1,j,kout)=(cff*tl_zeta(Iend+1,j,know)+         &
!>   &                              Cx *tl_zeta(Iend  ,j,kout)-         &
!>   &                              MAX(Ce,0.0_r8)*tl_grad(Iend+1,j  )- &
!>   &                              MIN(Ce,0.0_r8)*tl_grad(Iend+1,j+1))/&
!>   &                             (cff+Cx)
!>
            adfac=ad_zeta(Iend+1,j,kout)/(cff+Cx)
            ad_grad(Iend+1,j  )=ad_grad(Iend+1,j  )-MAX(Ce,0.0_r8)*adfac
            ad_grad(Iend+1,j+1)=ad_grad(Iend+1,j+1)-MIN(Ce,0.0_r8)*adfac
            ad_zeta(Iend  ,j,kout)=ad_zeta(Iend  ,j,kout)+Cx *adfac
            ad_zeta(Iend+1,j,know)=ad_zeta(Iend+1,j,know)+cff*adfac
            ad_zeta(Iend+1,j,kout)=0.0_r8
          END DO
        END IF

#  elif defined EAST_FSCHAPMAN
!
!  Eastern edge, Chapman boundary condition.
!
        DO j=Jstr,Jend
          cff=dt2d*GRID(ng)%pm(Iend,j)
          cff1=SQRT(g*(GRID(ng)%h(Iend,j)+                              &
     &                 zeta(Iend,j,know)))
          Cx=cff*cff1
          cff2=1.0_r8/(1.0_r8+Cx)
#   ifdef MASKING
!>        tl_zeta(Iend+1,j,kout)=tl_zeta(Iend+1,j,kout)*                &
!>   &                           GRID(ng)%rmask(Iend+1,j)
!>
          ad_zeta(Iend+1,j,kout)=ad_zeta(Iend+1,j,kout)*                &
     &                           GRID(ng)%rmask(Iend+1,j)
#   endif
!>        tl_zeta(Iend+1,j,kout)=tl_cff2*(zeta(Iend+1,j,know)+          &
!>   &                                    Cx*zeta(Iend,j,kout))+        &
!>   &                           cff2*(tl_zeta(Iend+1,j,know)+          &
!>   &                                 tl_Cx*zeta(Iend,j,kout)+         &
!>   &                                 Cx*tl_zeta(Iend,j,kout))
!>
          adfac=cff2*ad_zeta(Iend+1,j,kout)
          ad_zeta(Iend  ,j,kout)=ad_zeta(Iend  ,j,kout)+Cx*adfac
          ad_zeta(Iend+1,j,know)=ad_zeta(Iend+1,j,know)+adfac
          ad_Cx=ad_Cx+zeta(Iend,j,kout)*adfac
          ad_cff2=ad_cff2+                                              &
     &            (zeta(Iend+1,j,know)+                                 &
     &             Cx*zeta(Iend,j,kout))*ad_zeta(Iend+1,j,kout)
          ad_zeta(Iend+1,j,kout)=0.0_r8
!>        tl_cff2=-cff2*cff2*tl_Cx
!>
          ad_Cx=ad_Cx-cff2*cff2*ad_cff2
          ad_cff2=0.0_r8
!>        tl_Cx=cff*tl_cff1
!>
          ad_cff1=ad_cff1+cff*ad_Cx
          ad_Cx=0.0_r8
!>        tl_cff1=0.5_r8*g*(GRID(ng)%tl_h(Iend,j)+                      &
!>   &                      tl_zeta(Iend,j,know))/cff1
!>
          adfac=0.5_r8*g*ad_cff1/cff1
          GRID(ng)%ad_h(Iend,j)=GRID(ng)%ad_h(Iend,j)+adfac
          ad_zeta(Iend,j,know)=ad_zeta(Iend,j,know)+adfac
          ad_cff1=0.0_r8
        END DO

#  elif defined EAST_FSCLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=Jstr,Jend
!>        tl_zeta(Iend+1,j,kout)=0.0_r8
!>
          ad_zeta(Iend+1,j,kout)=0.0_r8
        END DO

#  elif defined EAST_FSGRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=Jstr,Jend
#  ifdef MASKING
!>        tl_zeta(Iend+1,j,kout)=tl_zeta(Iend+1,j,kout)*                &
!>   &                           GRID(ng)%rmask(Iend+1,j)
!>
          ad_zeta(Iend+1,j,kout)=ad_zeta(Iend+1,j,kout)*                &
     &                           GRID(ng)%rmask(Iend+1,j)
#  endif
!>        tl_zeta(Iend+1,j,kout)=tl_zeta(Iend,j,kout)
!>
          ad_zeta(Iend  ,j,kout)=ad_zeta(Iend  ,j,kout)+                &
     &                           ad_zeta(Iend+1,j,kout)
          ad_zeta(Iend+1,j,kout)=0.0_r8
        END DO

#  else
!
!  Eastern edge, closed boundary condition.
!
        DO j=Jstr,Jend
#  ifdef MASKING
!>        tl_zeta(Iend+1,j,kout)=tl_zeta(Iend+1,j,kout)*                &
!>   &                           GRID(ng)%rmask(Iend+1,j)
!>
          ad_zeta(Iend+1,j,kout)=ad_zeta(Iend+1,j,kout)*                &
     &                           GRID(ng)%rmask(Iend+1,j)
#  endif
!>        tl_zeta(Iend+1,j,kout)=tl_zeta(Iend,j,kout)
!>
          ad_zeta(Iend  ,j,kout)=ad_zeta(Iend  ,j,kout)+                &
     &                           ad_zeta(Iend+1,j,kout)
          ad_zeta(Iend+1,j,kout)=0.0_r8
        END DO
#  endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN

#  if defined WEST_FSRADIATION

        IF (iic(ng).ne.0) THEN
!
!  Western edge, implicit upstream radiation condition.
!
          DO j=Jstr,Jend
#   ifdef WEST_FSNUDGING
            IF (BOUNDARY(ng)%zeta_west_Cx(j).eq.0.0_r8) THEN
     &        tau=FSobc_in(ng,iwest)
            ELSE
              tau=FSobc_out(ng,iwest)
            END IF
            tau=tau*dt2d
#   endif
            Cx=BOUNDARY(ng)%zeta_west_Cx(j)
#   ifdef RADIATION_2D
            Ce=BOUNDARY(ng)%zeta_west_Ce(j)
#   else
            Ce=0.0_r8
#   endif
            cff=BOUNDARY(ng)%zeta_west_C2(j)
#   ifdef MASKING
!>          tl_zeta(Istr-1,j,kout)=tl_zeta(Istr-1,j,kout)*              &
!>   &                             GRID(ng)%rmask(Istr-1,j)
!>
            ad_zeta(Istr-1,j,kout)=ad_zeta(Istr-1,j,kout)*              &
     &                             GRID(ng)%rmask(Istr-1,j)
#   endif
#   ifdef WEST_FSNUDGING
!>          tl_zeta(Istr-1,j,kout)=tl_zeta(Istr-1,j,kout)-              &
!>   &                             tau*tl_zeta(Istr-1,j,know)
!>
            ad_zeta(Istr-1,j,know)=ad_zeta(Istr-1,j,know)-              &
     &                            tau*ad_zeta(Istr-1,j,kout)

#   endif
!>          tl_zeta(Istr-1,j,kout)=(cff*tl_zeta(Istr-1,j,know)+         &
!>   &                              Cx *tl_zeta(Istr  ,j,kout)-         &
!>   &                              MAX(Ce,0.0_r8)*tl_grad(Istr-1,j  )- &
!>   &                              MIN(Ce,0.0_r8)*tl_grad(Istr-1,j+1))/&
!>   &                             (cff+Cx)
!>
            adfac=ad_zeta(Istr-1,j,kout)/(cff+Cx)
            ad_grad(Istr-1,j  )=ad_grad(Istr-1,j  )-MAX(Ce,0.0_r8)*adfac
            ad_grad(Istr-1,j+1)=ad_grad(Istr-1,j+1)-MIN(Ce,0.0_r8)*adfac
            ad_zeta(Istr-1,j,know)=ad_zeta(Istr-1,j,know)+cff*adfac
            ad_zeta(Istr  ,j,kout)=ad_zeta(Istr  ,j,kout)+Cx *adfac
            ad_zeta(Istr-1,j,kout)=0.0_r8
          END DO
        END IF

#  elif defined WEST_FSCHAPMAN
!
!  Western edge, Chapman boundary condition.
!
        DO j=Jstr,Jend
          cff=dt2d*GRID(ng)%pm(Istr,j)
          cff1=SQRT(g*(GRID(ng)%h(Istr,j)+                              &
     &                 zeta(Istr,j,know)))
          Cx=cff*cff1
          cff2=1.0_r8/(1.0_r8+Cx)
#  ifdef MASKING
!>        tl_zeta(Istr-1,j,kout)=tl_zeta(Istr-1,j,kout)*                &
!>   &                           GRID(ng)%rmask(Istr-1,j)
!>
          ad_zeta(Istr-1,j,kout)=ad_zeta(Istr-1,j,kout)*                &
     &                           GRID(ng)%rmask(Istr-1,j)
#  endif
!>        tl_zeta(Istr-1,j,kout)=tl_cff2*(zeta(Istr-1,j,know)+          &
!>   &                                    Cx*zeta(Istr,j,kout))+        &
!>   &                           cff2*(tl_zeta(Istr-1,j,know)+          &
!>   &                                 tl_Cx*zeta(Istr,j,kout)+         &
!>   &                                 Cx*tl_zeta(Istr,j,kout))
!>
          adfac=cff2*ad_zeta(Istr-1,j,kout)
          ad_zeta(Istr-1,j,know)=ad_zeta(Istr-1,j,know)+adfac
          ad_zeta(Istr  ,j,kout)=ad_zeta(Istr  ,j,kout)+Cx*adfac
          ad_Cx=ad_Cx+zeta(Istr,j,kout)*adfac
          ad_cff2=ad_cff2+                                              &
     &            (zeta(Istr-1,j,know)+                                 &
     &             Cx*zeta(Istr,j,kout))*ad_zeta(Istr-1,j,kout)
          ad_zeta(Istr-1,j,kout)=0.0_r8
!>        tl_cff2=-cff2*cff2*tl_Cx
!>
          ad_Cx=ad_Cx-cff2*cff2*ad_cff2
          ad_cff2=0.0_r8
!>        tl_Cx=cff*tl_cff1
!>
          ad_cff1=ad_cff1+cff*ad_Cx
          ad_Cx=0.0_r8
!>        tl_cff1=0.5_r8*g*(GRID(ng)%tl_h(Istr,j)+                      &
!>   &                      tl_zeta(Istr,j,know))/cff1
!>
          adfac=0.5_r8*g*ad_cff1/cff1
          GRID(ng)%ad_h(Istr,j)=GRID(ng)%ad_h(Istr,j)+adfac
          ad_zeta(Istr,j,know)=ad_zeta(Istr,j,know)+adfac
          ad_cff1=0.0_r8
        END DO

#  elif defined WEST_FSCLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=Jstr,Jend
!>        tl_zeta(Istr-1,j,kout)=0.0_r8
!>
          ad_zeta(Istr-1,j,kout)=0.0_r8
        END DO

#  elif defined WEST_FSGRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=Jstr,Jend
#  ifdef MASKING
!>        tl_zeta(Istr-1,j,kout)=tl_zeta(Istr-1,j,kout)*                &
!>   &                           GRID(ng)%rmask(Istr-1,j)
!>
          ad_zeta(Istr-1,j,kout)=ad_zeta(Istr-1,j,kout)*                &
     &                           GRID(ng)%rmask(Istr-1,j)
#  endif
!>        tl_zeta(Istr-1,j,kout)=tl_zeta(Istr,j,kout)
!>
          ad_zeta(Istr  ,j,kout)=ad_zeta(Istr  ,j,kout)+                &
     &                           ad_zeta(Istr-1,j,kout)
          ad_zeta(Istr-1,j,kout)=0.0_r8
        END DO

#  else
!
!  Western edge, closed boundary condition.
!
        DO j=Jstr,Jend
#  ifdef MASKING
!>        tl_zeta(Istr-1,j,kout)=tl_zeta(Istr-1,j,kout)*                &
!>   &                           GRID(ng)%rmask(Istr-1,j)
!>
          ad_zeta(Istr-1,j,kout)=ad_zeta(Istr-1,j,kout)*                &
     &                           GRID(ng)%rmask(Istr-1,j)
#  endif
!>        tl_zeta(Istr-1,j,kout)=tl_zeta(Istr,j,kout)
!>
          ad_zeta(Istr  ,j,kout)=ad_zeta(Istr  ,j,kout)+                &
     &                           ad_zeta(Istr-1,j,kout)
          ad_zeta(Istr-1,j,kout)=0.0_r8
        END DO
#  endif
      END IF
# endif
      RETURN
      END SUBROUTINE ad_zetabc_tile
#endif
      END MODULE ad_zetabc_mod
