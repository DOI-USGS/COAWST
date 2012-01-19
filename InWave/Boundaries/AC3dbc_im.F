#include "cppdefs.h"
      MODULE AC3dbc_mod
#ifdef SOLVE3D
!
!svn $Id: t3dbc_im.F 732 2008-09-07 01:55:51Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for the ITRC-th    !
!  tracer field.                                                       !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: AC3dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE AC3dbc (ng, tile, nout)
!***********************************************************************
!
      USE mod_param
      USE mod_inwave_vars
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
      CALL AC3dbc_tile (ng, tile,                                       &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 nstp(ng),nout,                                   &
     &                 WAVEP(ng)% AC)
      RETURN
      END SUBROUTINE AC3dbc

!
!***********************************************************************
      SUBROUTINE AC3dbc_tile (ng, tile,                                 &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nstp,nout,                                 &
     &                       AC)
!***********************************************************************

      USE mod_param
      USE mod_inwave_params
      USE mod_inwave_bound
      USE mod_inwave_vars
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
# ifdef REFINED_GRID
      USE mod_stepping
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nout
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: AC(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: AC(LBi:UBi,LBj:UBj,ND(ng),3)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k, d_bnd

      real(r8), parameter :: eps =1.0E-20_r8

      real(r8) :: Ce, Cx, cff, dTde, dTdt, dTdx, tau

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: grad
#ifdef REFINED_GRID
      integer :: t1, t2, ngp
      real(r8) :: fac, fac1, fac2, my_val
#endif

# include "set_bounds.h"

#ifdef REFINED_GRID
!set this for Parent Grid ONLY
      IF (ng.eq.1) THEN
#endif

# if !defined EW_AC_PERIODIC && !defined COMPOSED_GRID
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN

#  if defined WEST_AC_RADIATION
!
!  Western edge, implicit upstream radiation condition.
!
        DO k=1,ND(ng)
          DO j=Jstr,Jend+1
            grad(Istr-1,j)=AC(Istr-1,j  ,k,nstp)-                       &
     &                     AC(Istr-1,j-1,k,nstp)
#   ifdef MASKING
            grad(Istr-1,j)=grad(Istr-1,j)*                              &
     &                     GRID(ng)%vmask(Istr-1,j)
#   endif
            grad(Istr  ,j)=AC(Istr  ,j  ,k,nstp)-                       &
     &                     AC(Istr  ,j-1,k,nstp)
#   ifdef MASKING
            grad(Istr  ,j)=grad(Istr  ,j)*                              &
     &                     GRID(ng)%vmask(Istr  ,j)
#   endif
          END DO
          DO j=Jstr,Jend
            dTdt=AC(Istr,j,k,nstp)-AC(Istr  ,j,k,nout)
            dTdx=AC(Istr,j,k,nout)-AC(Istr+1,j,k,nout)

            IF ((dTdt*dTdx).lt.0.0_r8) dTdt=0.0_r8
            IF ((dTdt*(grad(Istr,j)+grad(Istr,j+1))).gt.0.0_r8) THEN
              dTde=grad(Istr,j  )
            ELSE
              dTde=grad(Istr,j+1)
            END IF
            cff=MAX(dTdx*dTdx+dTde*dTde,eps)
            Cx=dTdt*dTdx
            Ce=MIN(cff,MAX(dTdt*dTde,-cff))

            AC(Istr-1,j,k,nout)=(cff*AC(Istr-1,j,k,nstp)+               &
     &                               Cx *AC(Istr  ,j,k,nout)-           &
     &                               MAX(Ce,0.0_r8)*                    &
     &                                  grad(Istr-1,j  )-               &
     &                               MIN(Ce,0.0_r8)*                    &
     &                                  grad(Istr-1,j+1))/              &
     &                              (cff+Cx)

#   ifdef MASKING
            AC(Istr-1,j,k,nout)=AC(Istr-1,j,k,nout)*                    &
     &                          GRID(ng)%rmask(Istr-1,j)
#   endif
          END DO
        END DO

#  elif defined WEST_AC_CLAMPED
!
!  Western edge, clamped boundary condition.
!
      DO k=1,ND(ng)
       DO d_bnd=1,WAVEB(ng)%ND_BND
        IF(WAVEB(ng)%WD_BND(d_bnd).EQ.WAVEG(ng)%wd(k))then
            DO j=Jstr-1,Jend+1
              AC(Istr-1,j,k,nout)=WAVEB(ng)%AC_west(j,d_bnd)
#   ifdef MASKING
              AC(Istr-1,j,k,nout)=AC(Istr-1,j,k,nout)*                  &
     &                          GRID(ng)%rmask(Istr-1,j)
#   endif
            END DO
          ENDIF
        END DO
      ENDDO

#  elif defined WEST_AC_GRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO k=1,ND(ng)
          DO j=Jstr-1,Jend+1
            AC(Istr-1,j,k,nout)=AC(Istr,j,k,nout)
#   ifdef MASKING
            AC(Istr-1,j,k,nout)=AC(Istr-1,j,k,nout)*                    &
     &                          GRID(ng)%rmask(Istr-1,j)
#   endif
          END DO
        END DO

#  else
!
!  Western edge, closed boundary condition.
!
        DO k=1,ND(ng)
          DO j=Jstr-1,Jend+1
            AC(Istr-1,j,k,nout)=AC(Istr,j,k,nout)
#   ifdef MASKING
            AC(Istr-1,j,k,nout)=AC(Istr-1,j,k,nout)*                    &
     &                          GRID(ng)%rmask(Istr-1,j)
#   endif
          END DO
        END DO
#  endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN

#  if defined EAST_AC_RADIATION
!
!  Eastern edge, implicit upstream radiation condition.
!
        DO k=1,ND(ng)
          DO j=Jstr,Jend+1
           grad(Iend  ,j)=AC(Iend  ,j  ,k,nstp)-                        &
     &                    AC(Iend  ,j-1,k,nstp)
#   ifdef MASKING
           grad(Iend  ,j)=grad(Iend  ,j)*                               &
     &                    GRID(ng)%vmask(Iend  ,j)
#   endif
           grad(Iend+1,j)=AC(Iend+1,j  ,k,nstp)-                        &
     &                    AC(Iend+1,j-1,k,nstp)
#   ifdef MASKING
           grad(Iend+1,j)=grad(Iend+1,j)*                               &
     &                    GRID(ng)%vmask(Iend+1,j)
#   endif
          END DO
          DO j=Jstr,Jend
            dTdt=AC(Iend,j,k,nstp)-AC(Iend  ,j,k,nout)
            dTdx=AC(Iend,j,k,nout)-AC(Iend-1,j,k,nout)
            IF ((dTdt*dTdx).lt.0.0_r8) dTdt=0.0_r8
            IF ((dTdt*(grad(Iend,j)+grad(Iend,j+1))).gt.0.0_r8) THEN
              dTde=grad(Iend,j  )
            ELSE
              dTde=grad(Iend,j+1)
            END IF
            cff=MAX(dTdx*dTdx+dTde*dTde,eps)
            Cx=dTdt*dTdx
            Ce=MIN(cff,MAX(dTdt*dTde,-cff))
            AC(Iend+1,j,k,nout)=(cff*AC(Iend+1,j,k,nstp)+               &
     &                               Cx *AC(Iend  ,j,k,nout)-           &
     &                               MAX(Ce,0.0_r8)*                    &
     &                                  grad(Iend+1,j  )-               &
     &                               MIN(Ce,0.0_r8)*                    &
     &                                  grad(Iend+1,j+1))/              &
     &                              (cff+Cx)

#   ifdef MASKING
            AC(Iend+1,j,k,nout)=AC(Iend+1,j,k,nout)*                    &
     &                          GRID(ng)%rmask(Iend+1,j)
#   endif
          END DO
        END DO

#  elif defined EAST_AC_CLAMPED
!
!  Eastern edge, clamped boundary condition.
!
      DO d_bnd=1,WAVEB(ng)%ND_BND
        DO k=1,ND(ng)
        IF(WAVEB(ng)%WD_BND(d_bnd).EQ.WAVEG(ng)%wd(k))then
          DO j=Jstr-1,Jend+1
            AC(Iend+1,j,k,nout)=WAVEB(ng)%AC_east(j,d_bnd)
#   ifdef MASKING
            AC(Iend+1,j,k,nout)=AC(Iend+1,j,k,nout)*                    &
     &                          GRID(ng)%rmask(Iend+1,j)
#   endif
          END DO
        ENDIF
        END DO
      ENDDO


#  elif defined EAST_AC_GRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO k=1,ND(ng)
          DO j=Jstr-1,Jend+1
            AC(Iend+1,j,k,nout)=AC(Iend,j,k,nout)
#   ifdef MASKING
            AC(Iend+1,j,k,nout)=AC(Iend+1,j,k,nout)*                    &
     &                          GRID(ng)%rmask(Iend+1,j)
#   endif
          END DO
        END DO

#  else
!
!  Eastern edge, closed boundary condition.
!
        DO k=1,ND(ng)
          DO j=Jstr-1,Jend+1
            AC(Iend+1,j,k,nout)=AC(Iend,j,k,nout)
#   ifdef MASKING
            AC(Iend+1,j,k,nout)=AC(Iend+1,j,k,nout)*                    &
     &                          GRID(ng)%rmask(Iend+1,j)
#   endif
          END DO
        END DO
#  endif
      END IF
# endif

# if !defined NS_AC_PERIODIC && !defined COMPOSED_GRID
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE) THEN

#  if defined SOUTH_AC_RADIATION
!
!  Southern edge, implicit upstream radiation condition.
!
        DO k=1,ND(ng)
          DO i=Istr,Iend+1
            grad(i,Jstr  )=AC(i  ,Jstr  ,k,nstp)-                       &
     &                     AC(i-1,Jstr  ,k,nstp)
#   ifdef MASKING
            grad(i,Jstr  )=grad(i,Jstr  )*                              &
     &                     GRID(ng)%umask(i,Jstr  )
#   endif
            grad(i,Jstr-1)=AC(i  ,Jstr-1,k,nstp)-                       &
     &                     AC(i-1,Jstr-1,k,nstp)
#   ifdef MASKING
            grad(i,Jstr-1)=grad(i,Jstr-1)*                              &
     &                     GRID(ng)%umask(i,Jstr-1)
#   endif
          END DO
          DO i=Istr,Iend
            dTdt=AC(i,Jstr,k,nstp)-AC(i,Jstr  ,k,nout)
            dTde=AC(i,Jstr,k,nout)-AC(i,Jstr+1,k,nout)
            IF ((dTdt*dTde).lt.0.0_r8) dTdt=0.0_r8
            IF ((dTdt*(grad(i,Jstr)+grad(i+1,Jstr))).gt.0.0_r8) THEN
              dTdx=grad(i  ,Jstr)
            ELSE
              dTdx=grad(i+1,Jstr)
            END IF
            cff=MAX(dTdx*dTdx+dTde*dTde,eps)
            Cx=MIN(cff,MAX(dTdt*dTdx,-cff))
            Ce=dTdt*dTde
            AC(i,Jstr-1,k,nout)=(cff*AC(i,Jstr-1,k,nstp)+               &
     &                               Ce *AC(i,Jstr  ,k,nout)-           &
     &                               MAX(Cx,0.0_r8)*                    &
     &                                  grad(i  ,Jstr-1)-               &
     &                               MIN(Cx,0.0_r8)*                    &
     &                                  grad(i+1,Jstr-1))/              &
     &                              (cff+Ce)
#   ifdef MASKING
            AC(i,Jstr-1,k,nout)=AC(i,Jstr-1,k,nout)*                    &
     &                          GRID(ng)%rmask(i,Jstr-1)
#   endif
          END DO
        END DO

#  elif defined SOUTH_AC_CLAMPED
!
!  Southern edge, clamped boundary condition.
!
      DO d_bnd=1,WAVEB(ng)%ND_BND
        DO k=1,ND(ng)
        IF(WAVEB(ng)%WD_BND(d_bnd).EQ.WAVEG(ng)%wd(k))then
          DO i=Istr-1,Iend+1
            AC(i,Jstr-1,k,nout)=WAVEB(ng)%AC_south(i,d_bnd)
#   ifdef MASKING
            AC(i,Jstr-1,k,nout)=AC(i,Jstr-1,k,nout)*                    &
     &                          GRID(ng)%rmask(i,Jstr-1)
#   endif
            END DO
          ENDIF
          END DO
        END DO

#  elif defined SOUTH_AC_GRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO k=1,ND(ng)
          DO i=Istr-1,Iend+1
            AC(i,Jstr-1,k,nout)=AC(i,Jstr,k,nout)
#   ifdef MASKING
            AC(i,Jstr-1,k,nout)=AC(i,Jstr-1,k,nout)*                    &
     &                          GRID(ng)%rmask(i,Jstr-1)
#   endif
          END DO
        END DO

#  else
!
!  Southern edge, closed boundary condition.
!
        DO k=1,ND(ng)
          DO i=Istr-1,Iend+1
            AC(i,Jstr-1,k,nout)=AC(i,Jstr,k,nout)
#   ifdef MASKING
            AC(i,Jstr-1,k,nout)=AC(i,Jstr-1,k,nout)*                    &
     &                          GRID(ng)%rmask(i,Jstr-1)
#   endif
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

#  if defined NORTH_AC_RADIATION
!
!  Northern edge, implicit upstream radiation condition.
!
        DO k=1,ND(ng)
          DO i=Istr,Iend+1
            grad(i,Jend  )=AC(i  ,Jend  ,k,nstp)-                       &
     &                     AC(i-1,Jend  ,k,nstp)
#   ifdef MASKING
            grad(i,Jend  )=grad(i,Jend  )*                              &
     &                     GRID(ng)%umask(i,Jend  )
#   endif
            grad(i,Jend+1)=AC(i  ,Jend+1,k,nstp)-                       &
     &                     AC(i-1,Jend+1,k,nstp)
#   ifdef MASKING
            grad(i,Jend+1)=grad(i,Jend+1)*                              &
     &                     GRID(ng)%umask(i,Jend+1)
#   endif
          END DO
          DO i=Istr,Iend
            dTdt=t(i,Jend,k,nstp,itrc)-t(i,Jend  ,k,nout,itrc)
            dTde=t(i,Jend,k,nout,itrc)-t(i,Jend-1,k,nout,itrc)
            IF ((dTdt*dTde).lt.0.0_r8) dTdt=0.0_r8
            IF ((dTdt*(grad(i,Jend)+grad(i+1,Jend))).gt.0.0_r8) THEN
              dTdx=grad(i  ,Jend)
            ELSE
              dTdx=grad(i+1,Jend)
            END IF
            cff=MAX(dTdx*dTdx+dTde*dTde,eps)
            Cx=MIN(cff,MAX(dTdt*dTdx,-cff))
            Ce=dTdt*dTde
            AC(i,Jend+1,k,nout)=(cff*AC(i,Jend+1,k,nstp)+               &
     &                               Ce *AC(i,Jend  ,k,nout)-           &
     &                               MAX(Cx,0.0_r8)*                    &
     &                                  grad(i  ,Jend+1)-               &
     &                               MIN(Cx,0.0_r8)*                    &
     &                                  grad(i+1,Jend+1))/              &
     &                              (cff+Ce)
#   ifdef MASKING
            AC(i,Jend+1,k,nout)=AC(i,Jend+1,k,nout)*                    &
     &                          GRID(ng)%rmask(i,Jend+1)
#   endif
          END DO
        END DO

#  elif defined NORTH_AC_CLAMPED
!
!  Northern edge, clamped boundary condition.
!
      DO d_bnd=1,WAVEB(ng)%ND_BND
        DO k=1,ND(ng)
        IF(WAVEB(ng)%WD_BND(d_bnd).EQ.WAVEG(ng)%wd(k))then
          DO i=Istr-1,Iend+1
            AC(i,Jend+1,k,nout)=WAVEB(ng)%AC_north(i,d_bnd)
#   ifdef MASKING
            AC(i,Jend+1,k,nout)=AC(i,Jend+1,k,nout)*                    &
     &                          GRID(ng)%rmask(i,Jend+1)
#   endif
          END DO
         ENDIF
        END DO
      ENDDO

#  elif defined NORTH_AC_GRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO k=1,ND(ng)
          DO i=Istr-1,Iend+1
            AC(i,Jend+1,k,nout)=AC(i,Jend,k,nout)
#   ifdef MASKING
            AC(i,Jend+1,k,nout)=AC(i,Jend+1,k,nout)*                    &
     &                          GRID(ng)%rmask(i,Jend+1)
#   endif
          END DO
        END DO

#  else
!
!  Northern edge, closed boundary condition.
!
        DO k=1,ND(ng)
          DO i=Istr-1,Iend+1
            AC(i,Jend+1,k,nout)=AC(i,Jend,k,nout)
#   ifdef MASKING
            AC(i,Jend+1,k,nout)=AC(i,Jend+1,k,nout)*                    &
     &                          GRID(ng)%rmask(i,Jend+1)
#   endif
          END DO
        END DO
#  endif
      END IF
# endif

# if !defined EW_AC_PERIODIC && !defined NS_AC_PERIODIC && !defined COMPOSED_GRID
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=1,ND(ng)
          AC(Istr-1,Jstr-1,k,nout)=0.5_r8*                              &
     &                                 (AC(Istr  ,Jstr-1,k,nout)+       &
     &                                  AC(Istr-1,Jstr  ,k,nout))
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,ND(ng)
          AC(Iend+1,Jstr-1,k,nout)=0.5_r8*                              &
     &                                 (AC(Iend  ,Jstr-1,k,nout)+       &
     &                                  AC(Iend+1,Jstr  ,k,nout))
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=1,ND(ng)
          AC(Istr-1,Jend+1,k,nout)=0.5_r8*                              &
     &                                 (AC(Istr-1,Jend  ,k,nout)+       &
     &                                  AC(Istr  ,Jend+1,k,nout))
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,ND(ng)
          AC(Iend+1,Jend+1,k,nout)=0.5_r8*                              &
     &                                 (AC(Iend+1,Jend  ,k,nout)+       &
     &                                  AC(Iend  ,Jend+1,k,nout))
        END DO
      END IF
# endif

# ifdef REFINED_GRID
!set this for Parent Grid ONLY
      END IF
# endif

      RETURN
      END SUBROUTINE AC3dbc_tile
#endif
      END MODULE AC3dbc_mod
