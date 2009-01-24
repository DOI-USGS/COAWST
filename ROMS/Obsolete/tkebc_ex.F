#include "cppdefs.h"
      MODULE tkebc
#if defined SOLVE3D && (defined MY25_MIXING || defined GLS_MIXING)
!
!svn $Id: tkebc_ex.F 503 2008-01-10 00:11:51Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!======================================================================= 
!                                                                      !
!  This subroutine sets lateral boundary conditions for turbulent      !
!  kinetic energy and turbulent length scale variables associated      !
!  with the Mellor and Yamada or GOTM closures.                        !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: tkebc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE tkebc (ng, tile, nout)
!***********************************************************************
!
      USE mod_param
      USE mod_mixing

      integer, intent(in) :: ng, tile, nout

# include "tile.h"
!
      CALL tkebc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 nout, nstp(ng),                                  &
     &                 MIXING(ng)% gls,                                 &
     &                 MIXING(ng)% tke)
      RETURN
      END SUBROUTINE tkebc
!
!***********************************************************************
      SUBROUTINE tkebc_tile (ng, Istr, Iend, Jstr, Jend,                &
     &                       LBi, UBi, LBj, UBj, UBk,                   &
     &                       nout, nstp,                                &
     &                       gls, tke)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk
      integer, intent(in) :: nout, nstp

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: gls(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: tke(LBi:,LBj:,0:,:)
# else
      real(r8), intent(inout) :: gls(LBi:UBi,LBj:UBj,0:UBk,3)
      real(r8), intent(inout) :: tke(LBi:UBi,LBj:UBj,0:UBk,3)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

      real(r8), parameter :: eps = 1.0e-20_r8

      real(r8) :: Ce, Cx, cff, dKde, dKdt, dKdx

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: grad
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: gradL

# include "set_bounds.h"

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN
!
#  if defined WEST_KRADIATION
!
!  Western edge, implicit upstream radiation condition.
!
        DO k=0,N(ng)
          DO j=Jstr,Jend+1
            grad(Istr-1,j)=tke(Istr-1,j  ,k,nstp)-                      &
     &                     tke(Istr-1,j-1,k,nstp)
#   ifdef MASKING
            grad(Istr-1,j)=grad (Istr-1,j)*GRID(ng)%vmask(Istr-1,j)
#   endif
            grad(Istr  ,j)=tke(Istr  ,j  ,k,nstp)-                      &
     &                     tke(Istr  ,j-1,k,nstp)
#   ifdef MASKING
            grad(Istr  ,j)=grad(Istr  ,j)*GRID(ng)%vmask(Istr  ,j)
#   endif
            gradL(Istr-1,j)=gls(Istr-1,j  ,k,nstp)-                     &
     &                      gls(Istr-1,j-1,k,nstp)
#   ifdef MASKING
            gradL(Istr-1,j)=gradL(Istr-1,j)*GRID(ng)%vmask(Istr-1,j)
#   endif
            gradL(Istr  ,j)=gls(Istr  ,j  ,k,nstp)-                     &
     &                      gls(Istr  ,j-1,k,nstp)
#   ifdef MASKING
            gradL(Istr,j)=gradL(Istr,j)*GRID(ng)%vmask(Istr,j)
#   endif
          END DO
          DO j=Jstr,Jend
            dKdt=tke(Istr,j,k,nstp)-tke(Istr  ,j,k,nout)
            dKdx=tke(Istr,j,k,nstp)-tke(istr+1,j,k,nstp)
            IF ((dKdt*dKdx).lt.0.0_r8) dKdt=0.0_r8
            IF ((dKdt*(grad(Istr,j)+grad(Istr,j+1))).gt.0.0_r8) THEN
              dKde=grad(Istr,j  )
            ELSE
              dKde=grad(Istr,j+1)
            END IF
            cff=dKdt/MAX(dKdx*dKdx+dKde*dKde,eps)
            Cx=MIN(1.0_r8,cff*dKdx)
#   ifdef RADIATION_2D
            Ce=MIN(1.0_r8,MAX(cff*dKde,-1.0_r8))
#   else
            Ce=0.0_r8
#   endif
            tke(Istr-1,j,k,nout)=(1.0_r8-Cx)*tke(Istr-1,j,k,nstp)+      &
     &                           Cx*tke(Istr,j,k,nstp)-                 &
     &                           MAX(Ce,0.0_r8)*grad(Istr-1,j  )-       &
     &                           MIN(Ce,0.0_r8)*grad(Istr-1,j+1)
#   ifdef MASKING
            tke(Istr-1,j,k,nout)=tke(Istr-1,j,k,nout)*                  &
     &                           GRID(ng)%rmask(Istr-1,j)
#   endif
            dKdt=gls(Istr,j,k,nstp)-gls(Istr  ,j,k,nout)
            dKdx=gls(Istr,j,k,nstp)-gls(Istr+1,j,k,nstp)
            IF ((dKdt*dKdx).lt.0.0_r8) dKdt=0.0_r8
            IF ((dKdt*(gradL(Istr,j)+gradL(Istr,j+1))).gt.0.0_r8) THEN
              dKde=gradL(Istr,j  )
            ELSE
              dKde=gradL(Istr,j+1)
            END IF
            cff=dKdt/MAX(dKdx*dKdx+dKde*dKde,eps)
            Cx=MIN(1.0_r8,cff*dKdx)
#   ifdef RADIATION_2D
            Ce=MIN(1.0_r8,MAX(cff*dKde,-1.0_r8))
#   else
            Ce=0.0_r8
#   endif
            gls(Istr-1,j,k,nout)=(1.0_r8-Cx)*gls(Istr-1,j,k,nstp)+      &
     &                           Cx*gls(Istr,j,k,nstp)-                 &
     &                           MAX(Ce,0.0_r8)*gradL(Istr-1,j  )-      &
     &                           MIN(Ce,0.0_r8)*gradL(Istr-1,j+1)
#   ifdef MASKING
            gls(Istr-1,j,k,nout)=gls(Istr-1,j,k,nout)*                  &
     &                           GRID(ng)%rmask(Istr-1,j)
#   endif
          END DO
        END DO

#  elif defined WEST_KGRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO k=0,N(ng)
          DO j=Jstr,Jend
            tke(Istr-1,j,k,nout)=tke(Istr,j,k,nout)
#   ifdef MASKING
            tke(Istr-1,j,k,nout)=tke(Istr-1,j,k,nout)*                  &
     &                           GRID(ng)%rmask(Istr-1,j)
#   endif
            gls(Istr-1,j,k,nout)=gls(Istr,j,k,nout)
#   ifdef MASKING
            gls(Istr-1,j,k,nout)=gls(Istr-1,j,k,nout)*                  &
     &                           GRID(ng)%rmask(Istr-1,j)
#   endif
          END DO
        END DO

#  else
!
!  Western edge, closed boundary condition.
!
        DO k=0,N(ng)
          DO j=Jstr,Jend
            tke(Istr-1,j,k,nout)=tke(Istr,j,k,nout)
#   ifdef MASKING
            tke(Istr-1,j,k,nout)=tke(Istr-1,j,k,nout)*                  &
     &                           GRID(ng)%rmask(Istr-1,j)
#   endif
            gls(Istr-1,j,k,nout)=gls(Istr,j,k,nout)
#   ifdef MASKING
            gls(Istr-1,j,k,nout)=gls(Istr-1,j,k,nout)*                  &
     &                           GRID(ng)%rmask(Istr-1,j)
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
!
#  if defined EAST_KRADIATION
!
!  Eastern edge, implicit upstream radiation condition.
!
        DO k=0,N(ng)
          DO j=Jstr,Jend+1
           grad(Iend  ,j)=tke(Iend  ,j  ,k,nstp)-                       &
     &                    tke(Iend  ,j-1,k,nstp)
#   ifdef MASKING
           grad(Iend,j)=grad(Iend,j)*GRID(ng)%vmask(Iend,j)
#   endif
           grad(Iend+1,j)=tke(Iend+1,j  ,k,nstp)-                       &
     &                    tke(Iend+1,j-1,k,nstp)
#   ifdef MASKING
           grad(Iend+1,j)=grad(Iend+1,j)*                               &
     &                    GRID(ng)%vmask(Iend+1,j)
#   endif
           gradL(Iend  ,j)=gls(Iend  ,j  ,k,nstp)-                      &
     &                     gls(Iend  ,j-1,k,nstp)
#   ifdef MASKING
           gradL(Iend  ,j)=gradL(Iend  ,j)*GRID(ng)%vmask(Iend  ,j)
#   endif
           gradL(Iend+1,j)=gls(Iend+1,j  ,k,nstp)-                      &
     &                     gls(Iend+1,j-1,k,nstp)
#   ifdef MASKING
           gradL(Iend+1,j)=gradL(Iend+1,j)*GRID(ng)%vmask(Iend+1,j)
#   endif
          END DO
          DO j=Jstr,Jend
            dKdt=tke(Iend,j,k,nstp)-tke(Iend  ,j,k,nout)
            dKdx=tke(Iend,j,k,nstp)-tke(Iend-1,j,k,nstp)
            IF ((dKdt*dKdx).lt.0.0_r8) dKdt=0.0_r8
            IF ((dKdt*(grad(Iend,j)+grad(Iend,j+1))).gt.0.0_r8) THEN
              dKde=grad(Iend,j  )
            ELSE
              dKde=grad(Iend,j+1)
            END IF
            cff=dKdt/MAX(dKdx*dKdx+dKde*dKde,eps)
            Cx=MIN(1.0_r8,cff*dKdx)
#   ifdef RADIATION_2D
            Ce=MIN(1.0_r8,MAX(cff*dKde,-1.0_r8))
#   else
            Ce=0.0_r8
#   endif
            tke(Iend+1,j,k,nout)=(1.0_r8-Cx)*tke(Iend+1,j,k,nstp)+      &
     &                           Cx*tke(Iend,j,k,nstp)-                 &
     &                           MAX(Ce,0.0_r8)*grad(Iend+1,j  )-       &
     &                           MIN(Ce,0.0_r8)*grad(Iend+1,j+1)
#   ifdef MASKING
            tke(Iend+1,j,k,nout)=tke(Iend+1,j,k,nout)*                  &
     &                           GRID(ng)%rmask(Iend+1,j)
#   endif
            dKdt=gls(Iend,j,k,nstp)-gls(Iend  ,j,k,nout)
            dKdx=gls(Iend,j,k,nstp)-gls(Iend-1,j,k,nstp)
            IF ((dKdt*dKdx).lt.0.0_r8) dKdt=0.0_r8
            IF ((dKdt*(gradL(Iend,j)+gradL(Iend,j+1))).gt.0.0_r8) THEN
              dKde=gradL(Iend,j  )
            ELSE
              dKde=gradL(Iend,j+1)
            END IF
            cff=dKdt/MAX(dKdx*dKdx+dKde*dKde,eps)
            Cx=MIN(1.0_r8,cff*dKdx)
#   ifdef RADIATION_2D
            Ce=MIN(1.0_r8,MAX(cff*dKde,-1.0_r8))
#   else
            Ce=0.0_r8
#   endif
            gls(Iend+1,j,k,nout)=(1.0_r8-Cx)*gls(Iend+1,j,k,nstp)+      &
     &                           Cx*gls(Iend,j,k,nstp)-                 &
     &                           MAX(Ce,0.0_r8)*gradL(Iend+1,j  )-      &
     &                           MIN(Ce,0.0_r8)*gradL(Iend+1,j+1)
#   ifdef MASKING
            gls(Iend+1,j,k,nout)=gls(Iend+1,j,k,nout)*                  &
     &                           GRID(ng)%rmask(Iend+1,j)
#   endif
          END DO
        END DO

#  elif defined EAST_KGRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO k=0,N(ng)
          DO j=Jstr,Jend
            tke(Iend+1,j,k,nout)=tke(Iend,j,k,nout)
#   ifdef MASKING
            tke(Iend+1,j,k,nout)=tke(Iend+1,j,k,nout)*                  &
     &                           GRID(ng)%rmask(Iend+1,j)
#   endif
            gls(Iend+1,j,k,nout)=gls(Iend,j,k,nout)
#   ifdef MASKING
            gls(Iend+1,j,k,nout)=gls(Iend+1,j,k,nout)*                  &
     &                           GRID(ng)%rmask(Iend+1,j)
#   endif
          END DO
        END DO

#  else
!
!  Eastern edge, closed boundary condition.
!
        DO k=0,N(ng)
          DO j=Jstr,Jend
            tke(Iend+1,j,k,nout)=tke(Iend,j,k,nout)
#   ifdef MASKING
            tke(Iend+1,j,k,nout)=tke(Iend+1,j,k,nout)*                  &
     &                           GRID(ng)%rmask(Iend+1,j)
#   endif
            gls(Iend+1,j,k,nout)=gls(Iend,j,k,nout)
#   ifdef MASKING
            gls(Iend+1,j,k,nout)=gls(Iend+1,j,k,nout)*                  &
     &                          GRID(ng)%rmask(Iend+1,j)
#   endif
          END DO
        END DO
#  endif
      END IF
# endif

# ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE) THEN
!
#  if defined SOUTH_KRADIATION
!
!  Southern edge, implicit upstream radiation condition.
!
        DO k=0,N(ng)
          DO i=Istr,Iend+1
            grad(i,Jstr  )=tke(i  ,Jstr  ,k,nstp)-                      &
     &                     tke(i-1,Jstr  ,k,nstp)
#   ifdef MASKING
            grad(i,Jstr  )=grad(i,Jstr  )*GRID(ng)%umask(i,Jstr  )
#   endif
            grad(i,Jstr-1)=tke(i  ,Jstr-1,k,nstp)-                      &
     &                     tke(i-1,Jstr-1,k,nstp)
#   ifdef MASKING
            grad(i,Jstr-1)=grad(i,Jstr-1)*GRID(ng)%umask(i,Jstr-1)
#   endif
            gradL(i,Jstr  )=gls(i  ,Jstr  ,k,nstp)-                     &
     &                      gls(i-1,Jstr  ,k,nstp)
#   ifdef MASKING
            gradL(i,Jstr  )=gradL(i,Jstr  )*GRID(ng)%umask(i,Jstr  )
#   endif
            gradL(i,Jstr-1)=gls(i  ,Jstr-1,k,nstp)-                     &
     &                      gls(i-1,Jstr-1,k,nstp)
#   ifdef MASKING
            gradL(i,Jstr-1)=gradL(i,Jstr-1)*GRID(ng)%umask(i,Jstr-1)
#   endif
          END DO
          DO i=Istr,Iend
            dKdt=tke(i,Jstr,k,nstp)-tke(i,Jstr  ,k,nout)
            dKde=tke(i,Jstr,k,nstp)-tke(i,Jstr+1,k,nstp)
            IF ((dKdt*dKde).lt.0.0_r8) dKdt=0.0_r8
            IF ((dKdt*(grad(i,Jstr)+grad(i+1,Jstr))).gt.0.0_r8) THEN
              dKdx=grad(i  ,Jstr)
            ELSE
              dKdx=grad(i+1,Jstr)
            END IF
            cff=dKdt/MAX(dKdx*dKdx+dKde*dKde,eps)
#   ifdef RADIATION_2D
            Cx=MIN(1.0_r8,MAX(cff*dKdx,-1.0_r8))
#   else
            Cx=0.0_r8
#   endif
            Ce=MIN(1.0_r8,cff*dKde)
            tke(i,Jstr-1,k,nout)=(1.0_r8-Ce)*tke(i,Jstr-1,k,nstp)+      &
     &                           Ce*tke(i,Jstr,k,nstp)-                 &
     &                           MAX(Cx,0.0_r8)*grad(i  ,Jstr-1)-       &
     &                           MIN(Cx,0.0_r8)*grad(i+1,Jstr-1)
#   ifdef MASKING
            tke(i,Jstr-1,k,nout)=tke(i,Jstr-1,k,nout)*                  &
     &                           GRID(ng)%rmask(i,Jstr-1)
#   endif
            dKdt=gls(i,Jstr,k,nstp)-gls(i,Jstr  ,k,nout)
            dKde=gls(i,Jstr,k,nstp)-gls(i,Jstr+1,k,nstp)
            IF ((dKdt*dKde).lt.0.0_r8) dKdt=0.0_r8
            IF ((dKdt*(gradL(i,Jstr)+gradL(i+1,Jstr))).gt.0.0_r8) THEN
              dKdx=gradL(i  ,Jstr)
            ELSE
              dKdx=gradL(i+1,Jstr)
            END IF
            cff=dKdt/MAX(dKdx*dKdx+dKde*dKde,eps)
#   ifdef RADIATION_2D
            Cx=MIN(1.0_r8,MAX(cff*dKdx,-1.0_r8))
#   else
            Cx=0.0_r8
#   endif
            Ce=MIN(1.0_r8,cff*dKde)
            gls(i,Jstr-1,k,nout)=(1.0_r8-Ce)*gls(i,Jstr-1,k,nstp)+      &
     &                           Ce*gls(i,Jstr,k,nstp)-                 &
     &                           MAX(Cx,0.0_r8)*gradL(i  ,Jstr-1)-      &
     &                           MIN(Cx,0.0_r8)*gradL(i+1,Jstr-1)
#   ifdef MASKING
            gls(i,Jstr-1,k,nout)=gls(i,Jstr-1,k,nout)*                  &
     &                           GRID(ng)%rmask(i,Jstr-1)
#   endif
          END DO
        END DO

#  elif defined SOUTH_KGRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO k=0,N(ng)
          DO i=Istr,Iend
            tke(i,Jstr-1,k,nout)=tke(i,Jstr,k,nout)
#   ifdef MASKING
            tke(i,Jstr-1,k,nout)=tke(i,Jstr-1,k,nout)*                  &
     &                           GRID(ng)%rmask(i,Jstr-1)
#   endif
            gls(i,Jstr-1,k,nout)=gls(i,Jstr,k,nout)
#   ifdef MASKING
            gls(i,Jstr-1,k,nout)=gls(i,Jstr-1,k,nout)*                  &
     &                           GRID(ng)%rmask(i,Jstr-1)
#   endif
          END DO
        END DO

#  else
!
!  Southern edge, closed boundary condition.
!
        DO k=0,N(ng)
          DO i=Istr,Iend
            tke(i,Jstr-1,k,nout)=tke(i,Jstr,k,nout)
#   ifdef MASKING
            tke(i,Jstr-1,k,nout)=tke(i,Jstr-1,k,nout)*                  &
     &                           GRID(ng)%rmask(i,Jstr-1)
#   endif
            gls(i,Jstr-1,k,nout)=gls(i,Jstr,k,nout)
#   ifdef MASKING
            gls(i,Jstr-1,k,nout)=gls(i,Jstr-1,k,nout)*                  &
     &                           GRID(ng)%rmask(i,Jstr-1)
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
!
#  if defined NORTH_KRADIATION
!
!  Northern edge, implicit upstream radiation condition.
!
        DO k=0,N(ng)
          DO i=Istr,Iend+1
            grad(i,Jend  )=tke(i  ,Jend  ,k,nstp)-                      &
     &                     tke(i-1,Jend  ,k,nstp)
#   ifdef MASKING
            grad(i,Jend  )=grad(i,Jend  )*GRID(ng)%umask(i,Jend  )
#   endif
            grad(i,Jend+1)=tke(i  ,Jend+1,k,nstp)-                      &
     &                     tke(i-1,Jend+1,k,nstp)
#   ifdef MASKING
            grad(i,Jend+1)=grad(i,Jend+1)*GRID(ng)%umask(i,Jend+1)
#   endif
            gradL(i,Jend  )=gls(i  ,Jend  ,k,nstp)-                     &
     &                      gls(i-1,Jend  ,k,nstp)
#   ifdef MASKING
            gradL(i,Jend  )=gradL(i,Jend  )*GRID(ng)%umask(i,Jend  )
#   endif
            gradL(i,Jend+1)=gls(i  ,Jend+1,k,nstp)-                     &
     &                      gls(i-1,Jend+1,k,nstp)
#   ifdef MASKING
            gradL(i,Jend+1)=gradL(i,Jend+1)*GRID(ng)%umask(i,Jend+1)
#   endif
          END DO
          DO i=Istr,Iend
            dKdt=tke(i,Jend,k,nstp)-tke(i,Jend  ,k,nout)
            dKde=tke(i,Jend,k,nstp)-tke(i,Jend-1,k,nstp)
            IF ((dKdt*dKde).lt.0.0_r8) dKdt=0.0_r8
            IF ((dKdt*(grad(i,Jend)+grad(i+1,Jend))).gt.0.0_r8) THEN
              dKdx=grad(i  ,Jend)
            ELSE
              dKdx=grad(i+1,Jend)
            END IF
            cff=dKdt/MAX(dKdx*dKdx+dKde*dKde,eps)
#   ifdef RADIATION_2D
            Cx=MIN(1.0_r8,MAX(cff*dKdx,-1.0_r8))
#   else
            Cx=0.0_r8
#   endif
            Ce=MIN(1.0_r8,cff*dKde)
            tke(i,Jend+1,k,nout)=(1.0_r8-Ce)*tke(i,Jend+1,k,nstp)+      &
     &                           Ce*tke(i,Jend,k,nstp)-                 &
     &                           MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-       &
     &                           MIN(Cx,0.0_r8)*grad(i+1,Jend+1)
#   ifdef MASKING
            tke(i,Jend+1,k,nout)=tke(i,Jend+1,k,nout)*                  &
     &                           GRID(ng)%rmask(i,Jend+1)
#   endif
            dKdt=gls(i,Jend,k,nstp)-gls(i,Jend  ,k,nout)
            dKde=gls(i,Jend,k,nstp)-gls(i,Jend-1,k,nstp)
            IF ((dKdt*dKde).lt.0.0_r8) dKdt=0.0_r8
            IF ((dKdt*(gradL(i  ,Jend)+gradL(i+1,Jend))).gt.0.0_r8) THEN
              dKdx=gradL(i  ,Jend)
            ELSE
              dKdx=gradL(i+1,Jend)
            END IF
            cff=dKdt/MAX(dKdx*dKdx+dKde*dKde,eps)
#   ifdef RADIATION_2D
            Cx=MIN(1.0_r8,MAX(cff*dKdx,-1.0_r8))
#   else
            Cx=0.0_r8
#   endif
            Ce=MIN(1.0_r8,cff*dKde)
            gls(i,Jend+1,k,nout)=(1.0_r8-Ce)*gls(i,Jend+1,k,nstp)+      &
     &                           Ce*gls(i,Jend,k,nstp)-                 &
     &                           MAX(Cx,0.0_r8)*gradL(i  ,Jend+1)-      &
     &                           MIN(Cx,0.0_r8)*gradL(i+1,Jend+1)
#   ifdef MASKING
            gls(i,Jend+1,k,nout)=gls(i,Jend+1,k,nout)*                  &
     &                           GRID(ng)%rmask(i,Jend+1)
#   endif
          END DO
        END DO

#  elif defined NORTH_KGRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO k=0,N(ng)
          DO i=Istr,Iend
            tke(i,Jend+1,k,nout)=tke(i,Jend,k,nout)
#   ifdef MASKING
            tke(i,Jend+1,k,nout)=tke(i,Jend+1,k,nout)*                  &
     &                           GRID(ng)%rmask(i,Jend+1)
#   endif
            gls(i,Jend+1,k,nout)=gls(i,Jend,k,nout)
#   ifdef MASKING
            gls(i,Jend+1,k,nout)=gls(i,Jend+1,k,nout)*                  &
     &                           GRID(ng)%rmask(i,Jend+1)
#   endif
          END DO
        END DO

#  else
!
!  Northern edge, closed boundary condition.
!
        DO k=0,N(ng)
          DO i=Istr,Iend
            tke(i,Jend+1,k,nout)=tke(i,Jend,k,nout)
#   ifdef MASKING
            tke(i,Jend+1,k,nout)=tke(i,Jend+1,k,nout)*                  &
     &                           GRID(ng)%rmask(i,Jend+1)
#   endif
            gls(i,Jend+1,k,nout)=gls(i,Jend,k,nout)
#   ifdef MASKING
            gls(i,Jend+1,k,nout)=gls(i,Jend+1,k,nout)*                  &
     &                           GRID(ng)%rmask(i,Jend+1)
#   endif
          END DO
        END DO
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
        DO k=0,N(ng)
          tke(Istr-1,Jstr-1,k,nout)=0.5_r8*(tke(Istr  ,Jstr-1,k,nout)+  &
     &                                      tke(Istr-1,Jstr  ,k,nout))
          gls(Istr-1,Jstr-1,k,nout)=0.5_r8*(gls(Istr  ,Jstr-1,k,nout)+  &
     &                                      gls(Istr-1,Jstr  ,k,nout))
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=0,N(ng)
          tke(Iend+1,Jstr-1,k,nout)=0.5_r8*(tke(Iend  ,Jstr-1,k,nout)+  &
     &                                      tke(Iend+1,Jstr  ,k,nout))
          gls(Iend+1,Jstr-1,k,nout)=0.5_r8*(gls(Iend  ,Jstr-1,k,nout)+  &
     &                                      gls(Iend+1,Jstr  ,k,nout))
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=0,N(ng)
          tke(Istr-1,Jend+1,k,nout)=0.5_r8*(tke(Istr  ,Jend+1,k,nout)+  &
     &                                      tke(Istr-1,Jend  ,k,nout))
          gls(Istr-1,Jend+1,k,nout)=0.5_r8*(gls(Istr  ,Jend+1,k,nout)+  &
     &                                      gls(Istr-1,Jend  ,k,nout))
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=0,N(ng)
          tke(Iend+1,Jend+1,k,nout)=0.5_r8*(tke(Iend  ,Jend+1,k,nout)+  &
     &                                      tke(Iend+1,Jend  ,k,nout))
          gls(Iend+1,Jend+1,k,nout)=0.5_r8*(gls(Iend  ,Jend+1,k,nout)+  &
     &                                      gls(Iend+1,Jend  ,k,nout))
        END DO
      END IF
# endif
      RETURN
      END SUBROUTINE tkebc
#endif
      END MODULE tkebc_mod
