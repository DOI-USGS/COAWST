#include "cppdefs.h"
      MODULE extract_sta_mod
!
!svn $Id: extract_sta.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine extracts field at the requested (Xpos,Ypos,Zpos)       !
!  positions.  The extraction is done using linear interpolation.      !
!  The (Xpos,Ypos) positions are in fractional grid coordinates.       !
!  Zpos is in fractional grid coordinates (Zpos >= 0) or actual        !
!  depths (Zpos < 0), if applicable.                                   !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Cgrid      Switch to interpolate at native C-grid (TRUE) or to   !
!                  interpolate at RHO-points (FALSE).                  !
!     ifield     Field ID.                                             !
!     gtype      Grid type.                                            !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     LBk        K-dimension Lower bound, if any. Otherwise, a value   !
!                  of one is expected.                                 !
!     LBk        K-dimension Upper bound, if any. Otherwise, a value   !
!                  of one is expected.                                 !
!     UBk        K-dimension Upper bound.                              !
!     Ascl       Factor to scale field after extraction.               !
!     A          Tile array (2D or 3D) to process.                     !
!     Npos       Number of values to extract.                          !
!     Xpos       X-extraction positions (grid coordinates).            !
!     Ypos       Y-extraction positions (grid coordinates).            !
!     Zpos       Z-extraction positions (grid coordinates or depth).   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Apos       Extracted values.                                     !
!                                                                      !
!  Note:                                                               !
!                                                                      !
!  Starting F95 zero values can be signed (-0 or +0) following the     !
!  IEEE 754 floating point standard. This can be advantageous in       !
!  some computations but not here when "Ascl" is negative and "Apos"   !
!  is zero.  This will produce different output files in serial        !
!  and distributed memory applications. Since comparing serial and     !
!  parallel output is essential for tracking parallel partition        !
!  bugs, "positive zero" is enforced.                                  !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE extract_sta2d (ng, model, Cgrid, ifield, gtype,        &
     &                          LBi, UBi, LBj, UBj, Ascl, A,            &
     &                          Npos, Xpos, Ypos, Apos)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_collect
# endif
!
!  Imported variable declarations.
!
      logical, intent(in) :: Cgrid

      integer, intent(in) :: ng, model, ifield, gtype, Npos
      integer, intent(in) :: LBi, UBi, LBj, UBj
      real(r8), intent(in) :: Ascl

#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: A(LBi:,LBj:)
      real(r8), intent(in) :: Xpos(:), Ypos(:)
      real(r8), intent(out) :: Apos(Npos)
#else
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Xpos(Npos), Ypos(Npos)
      real(r8), intent(out) :: Apos(Npos)
#endif
!
!  Local variable declarations.
!
      integer :: i1, i2, j1, j2, np

      real(r8), parameter :: Aspv = 0.0_r8

      real(r8) :: Xmin, Xmax, Ymin, Ymax
      real(r8) :: Xgrd, Xoff, Ygrd, Yoff
      real(r8) :: p1, p2, q1, q2, r1, r2, wsum
      real(r8) :: w111, w211, w121, w221

      real(r8), dimension(Npos) :: bounded
!
!-----------------------------------------------------------------------
!  Interpolate from 2D field at RHO-points.
!-----------------------------------------------------------------------
!
      IF (gtype.eq.r2dvar) THEN
        Xmin=rXmin(ng)
        Xmax=rXmax(ng)
        Ymin=rYmin(ng)
        Ymax=rYmax(ng)
        DO np=1,Npos
          Xgrd=Xpos(np)
          Ygrd=Ypos(np)
          bounded(np)=0.0_r8
          IF (((Xmin.le.Xgrd).and.(Xgrd.lt.Xmax)).and.                  &
     &        ((Ymin.le.Ygrd).and.(Ygrd.lt.Ymax))) THEN
            i1=INT(Xgrd)
            j1=INT(Ygrd)
            i2=i1+1
            j2=j1+1
            IF (i2.gt.Lm(ng)+1) THEN
              i2=i1                   ! station at the eastern boundary
            END IF
            IF (j2.gt.Mm(ng)+1) THEN
              j2=j1                   ! station at the northern boundary
            END IF
            bounded(np)=1.0_r8
            p2=REAL(i2-i1,r8)*(Xgrd-REAL(i1,r8))
            q2=REAL(j2-j1,r8)*(Ygrd-REAL(j1,r8))
            p1=1.0_r8-p2
            q1=1.0_r8-q2
            w111=p1*q1
            w211=p2*q1
            w121=p1*q2
            w221=p2*q2
#ifdef MASKING
            w111=w111*GRID(ng)%rmask(i1,j1)
            w211=w211*GRID(ng)%rmask(i2,j1)
            w121=w121*GRID(ng)%rmask(i1,j2)
            w221=w221*GRID(ng)%rmask(i2,j2)
            wsum=w111+w211+w121+w221
            IF (wsum.gt.0.0_r8) THEN
              wsum=1.0_r8/wsum
              w111=w111*wsum
              w211=w211*wsum
              w121=w121*wsum
              w221=w221*wsum
            ELSE
              bounded(np)=0.0_r8
            ENDIF
#endif
            Apos(np)=Ascl*(w111*A(i1,j1)+                               &
     &                     w211*A(i2,j1)+                               &
     &                     w121*A(i1,j2)+                               &
     &                     w221*A(i2,j2))
            IF (ABS(Apos(np)).eq.0.0_r8) Apos(np)=0.0_r8      ! positive
          ELSE                                                ! zero
            Apos(np)=Aspv
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Interpolate from 2D field at U-points.
!-----------------------------------------------------------------------
!
      ELSE IF (gtype.eq.u2dvar) THEN
        IF (Cgrid) THEN
          Xmin=uXmin(ng)+0.5_r8
          Xmax=uXmax(ng)+0.5_r8
          Ymin=uYmin(ng)
          Ymax=uYmax(ng)
          Xoff=0.0_r8
        ELSE
          Xmin=rXmin(ng)
          Xmax=rXmax(ng)
          Ymin=rYmin(ng)
          Ymax=rYmax(ng)
          Xoff=0.5_r8
        END IF
        DO np=1,Npos
          Xgrd=Xpos(np)+Xoff
          Ygrd=Ypos(np)
          bounded(np)=0.0_r8
          IF (((Xmin.le.Xgrd).and.(Xgrd.lt.Xmax)).and.                  &
     &        ((Ymin.le.Ygrd).and.(Ygrd.lt.Ymax))) THEN
            i1=INT(Xgrd)
            j1=INT(Ygrd)
            i2=i1+1
            j2=j1+1
            IF (i2.gt.Lm(ng)+1) THEN
              i2=i1                   ! station at the eastern boundary
            END IF
            IF (j2.gt.Mm(ng)+1) THEN
              j2=j1                   ! station at the northern boundary
            END IF
            bounded(np)=1.0_r8
            p2=REAL(i2-i1,r8)*(Xgrd-REAL(i1,r8))
            q2=REAL(j2-j1,r8)*(Ygrd-REAL(j1,r8))
            p1=1.0_r8-p2
            q1=1.0_r8-q2
            w111=p1*q1
            w211=p2*q1
            w121=p1*q2
            w221=p2*q2
#ifdef MASKING
            w111=w111*GRID(ng)%umask(i1,j1)
            w211=w211*GRID(ng)%umask(i2,j1)
            w121=w121*GRID(ng)%umask(i1,j2)
            w221=w221*GRID(ng)%umask(i2,j2)
            wsum=w111+w211+w121+w221
            IF (wsum.gt.0.0_r8) THEN
              wsum=1.0_r8/wsum
              w111=w111*wsum
              w211=w211*wsum
              w121=w121*wsum
              w221=w221*wsum
            ELSE
              bounded(np)=0.0_r8
            END IF
#endif
            Apos(np)=Ascl*(w111*A(i1,j1)+                               &
     &                     w211*A(i2,j1)+                               &
     &                     w121*A(i1,j2)+                               &
     &                     w221*A(i2,j2))
            IF (ABS(Apos(np)).eq.0.0_r8) Apos(np)=0.0_r8      ! positive
          ELSE                                                ! zero
            Apos(np)=Aspv
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Interpolate from 2D field at V-points.
!-----------------------------------------------------------------------
!
      ELSE IF (gtype.eq.v2dvar) THEN
        IF (Cgrid) THEN
          Xmin=vXmin(ng)
          Xmax=vXmax(ng)
          Ymin=vYmin(ng)+0.5_r8
          Ymax=vYmax(ng)+0.5_r8
          Yoff=0.0_r8
        ELSE
          Xmin=rXmin(ng)
          Xmax=rXmax(ng)
          Ymin=rYmin(ng)
          Ymax=rYmax(ng)
          Yoff=0.5_r8
        END IF
        DO np=1,Npos
          Xgrd=Xpos(np)
          Ygrd=Ypos(np)+Yoff
          bounded(np)=0.0_r8
          IF (((Xmin.le.Xgrd).and.(Xgrd.lt.Xmax)).and.                  &
     &        ((Ymin.le.Ygrd).and.(Ygrd.lt.Ymax))) THEN
            i1=INT(Xgrd)
            j1=INT(Ygrd)
            i2=i1+1
            j2=j1+1
            IF (i2.gt.Lm(ng)+1) THEN
              i2=i1                   ! station at the eastern boundary
            END IF
            IF (j2.gt.Mm(ng)+1) THEN
              j2=j1                   ! station at the northern boundary
            END IF
            bounded(np)=1.0_r8
            p2=REAL(i2-i1,r8)*(Xgrd-REAL(i1,r8))
            q2=REAL(j2-j1,r8)*(Ygrd-REAL(j1,r8))
            p1=1.0_r8-p2
            q1=1.0_r8-q2
            w111=p1*q1
            w211=p2*q1
            w121=p1*q2
            w221=p2*q2
#ifdef MASKING
            w111=w111*GRID(ng)%vmask(i1,j1)
            w211=w211*GRID(ng)%vmask(i2,j1)
            w121=w121*GRID(ng)%vmask(i1,j2)
            w221=w221*GRID(ng)%vmask(i2,j2)
            wsum=w111+w211+w121+w221
            IF (wsum.gt.0.0_r8) THEN
              wsum=1.0_r8/wsum
              w111=w111*wsum
              w211=w211*wsum
              w121=w121*wsum
              w221=w221*wsum
            ELSE
              bounded(np)=0.0_r8
            END IF
#endif
            Apos(np)=Ascl*(w111*A(i1,j1)+                               &
     &                     w211*A(i2,j1)+                               &
     &                     w121*A(i1,j2)+                               &
     &                     w221*A(i2,j2))
            IF (ABS(Apos(np)).eq.0.0_r8) Apos(np)=0.0_r8      ! positive
          ELSE                                                ! zero
            Apos(np)=Aspv
          END IF
        END DO
      END IF
#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Collect all extracted data.
!-----------------------------------------------------------------------
!
      CALL mp_collect (ng, model, Npos, Aspv, Apos)
      CALL mp_collect (ng, model, Npos, 0.0_r8, bounded)
#endif
!
!-----------------------------------------------------------------------
!  Set unbounded data to special value.
!-----------------------------------------------------------------------
!
      DO np=1,Npos
        IF (bounded(np).lt.1.0_r8) THEN
          Apos(np)=spval
        END IF
      END DO
      RETURN
      END SUBROUTINE extract_sta2d

#ifdef SOLVE3D
!
!***********************************************************************
      SUBROUTINE extract_sta3d (ng, model, Cgrid, ifield, gtype,        &
     &                          LBi, UBi, LBj, UBj, LBk, UBk, Ascl, A,  &
     &                          Npos, Xpos, Ypos, Zpos, Apos)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_collect
# endif
!
!  Imported variable declarations.
!
      logical, intent(in) :: Cgrid

      integer, intent(in) :: ng, model, ifield, gtype, Npos
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      real(r8), intent(in) :: Ascl
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: A(LBi:,LBj:,LBk:)
      real(r8), intent(in) :: Xpos(:), Ypos(:), Zpos(:)
      real(r8), intent(out) :: Apos(:)
# else
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
      real(r8), intent(in) :: Xpos(Npos), Ypos(Npos), Zpos(Npos)
      real(r8), intent(out) :: Apos(Npos)
# endif
!
!  Local variable declarations.
!
      integer :: i1, i2, j1, j2, k, k1, k2, np

      real(r8), parameter :: Aspv = 0.0_r8

      real(r8) :: Xmin, Xmax, Ymin, Ymax
      real(r8) :: Xgrd, Xoff, Ygrd, Yoff, Zgrd, Zbot, Ztop
      real(r8) :: dz, p1, p2, q1, q2, r1, r2, wsum
      real(r8) :: w111, w211, w121, w221, w112, w212, w122, w222

      real(r8), dimension(Npos) :: bounded
!
!-----------------------------------------------------------------------
!  Interpolate from 3D field at RHO-points.
!-----------------------------------------------------------------------
!
      IF (gtype.eq.r3dvar) THEN
        Xmin=rXmin(ng)
        Xmax=rXmax(ng)
        Ymin=rYmin(ng)
        Ymax=rYmax(ng)
        DO np=1,Npos
          Xgrd=Xpos(np)
          Ygrd=Ypos(np)
          Zgrd=Zpos(np)
          bounded(np)=0.0_r8
          IF (((Xmin.le.Xgrd).and.(Xgrd.lt.Xmax)).and.                  &
     &        ((Ymin.le.Ygrd).and.(Ygrd.lt.Ymax))) THEN
            i1=INT(Xgrd)
            j1=INT(Ygrd)
            i2=i1+1
            j2=j1+1
            IF (i2.gt.Lm(ng)+1) THEN
              i2=i1                   ! station at the eastern boundary
            END IF
            IF (j2.gt.Mm(ng)+1) THEN
              j2=j1                   ! station at the northern boundary
            END IF
            bounded(np)=1.0_r8
            p2=REAL(i2-i1,r8)*(Xgrd-REAL(i1,r8))
            q2=REAL(j2-j1,r8)*(Ygrd-REAL(j1,r8))
            p1=1.0_r8-p2
            q1=1.0_r8-q2
            w111=p1*q1
            w211=p2*q1
            w121=p1*q2
            w221=p2*q2
            w112=0.0_r8
            w212=0.0_r8
            w122=0.0_r8
            w222=0.0_r8
            IF (Zgrd.ge.0.0_r8) THEN
              k1=INT(Zgrd)
              k2=INT(Zgrd)
              r1=1.0_r8
              r2=0.0_r8
            ELSE
              Ztop=GRID(ng)%z_r(i1,j1,N(ng))
              Zbot=GRID(ng)%z_r(i1,j1,1)
              IF (Zgrd.ge.Ztop) THEN
                k1=N(ng)                     ! If shallower, assign
                k2=N(ng)                     ! station to surface
                r1=1.0_r8                    ! level
                r2=0.0_r8
              ELSE IF (Zbot.ge.Zgrd) THEN
                k1=1                         ! If deeper, assign
                k2=1                         ! station to bottom
                r1=1.0_r8                    ! level
                r2=0.0_r8
              ELSE
                DO k=N(ng),2,-1
                  Ztop=GRID(ng)%z_r(i1,j1,k)
                  Zbot=GRID(ng)%z_r(i1,j1,k-1)
                  IF ((Ztop.gt.Zgrd).and.(Zgrd.ge.Zbot)) THEN
                    k1=k-1
                    k2=k
                  END IF
                END DO
                dz=GRID(ng)%z_r(i1,j1,k2)-GRID(ng)%z_r(i1,j1,k1)
                r2=(Zgrd-GRID(ng)%z_r(i1,j1,k1))/dz
                r1=1.0_r8-r2
              END IF
            END IF
            w112=w111*r2
            w212=w211*r2
            w122=w121*r2
            w222=w221*r2
            w111=w111*r1
            w211=w211*r1
            w121=w121*r1
            w221=w221*r1
# ifdef MASKING
            w111=w111*GRID(ng)%rmask(i1,j1)
            w211=w211*GRID(ng)%rmask(i2,j1)
            w121=w121*GRID(ng)%rmask(i1,j2)
            w221=w221*GRID(ng)%rmask(i2,j2)
            w112=w112*GRID(ng)%rmask(i1,j1)
            w212=w212*GRID(ng)%rmask(i2,j1)
            w122=w122*GRID(ng)%rmask(i1,j2)
            w222=w222*GRID(ng)%rmask(i2,j2)
            wsum=w111+w211+w121+w221+w112+w212+w122+w222
            IF (wsum.gt.0.0_r8) THEN
              wsum=1.0_r8/wsum
              w111=w111*wsum
              w211=w211*wsum
              w121=w121*wsum
              w221=w221*wsum
              w112=w112*wsum
              w212=w212*wsum
              w122=w122*wsum
              w222=w222*wsum
            ELSE
              bounded(np)=0.0_r8
            END IF
# endif
            Apos(np)=Ascl*(w111*A(i1,j1,k1)+                            &
     &                     w211*A(i2,j1,k1)+                            &
     &                     w121*A(i1,j2,k1)+                            &
     &                     w221*A(i2,j2,k1)+                            &
     &                     w112*A(i1,j1,k2)+                            &
     &                     w212*A(i2,j1,k2)+                            &
     &                     w122*A(i1,j2,k2)+                            &
     &                     w222*A(i2,j2,k2))
            IF (ABS(Apos(np)).eq.0.0_r8) Apos(np)=0.0_r8      ! positive
          ELSE                                                ! zero
            Apos(np)=Aspv
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Interpolate from 3D field at U-points.
!-----------------------------------------------------------------------
!
      ELSE IF (gtype.eq.u3dvar) THEN
        IF (Cgrid) THEN
          Xmin=uXmin(ng)+0.5_r8
          Xmax=uXmax(ng)+0.5_r8
          Ymin=uYmin(ng)
          Ymax=uYmax(ng)
          Xoff=0.0_r8
        ELSE
          Xmin=rXmin(ng)
          Xmax=rXmax(ng)
          Ymin=rYmin(ng)
          Ymax=rYmax(ng)
          Xoff=0.5_r8
        END IF
        DO np=1,Npos
          Xgrd=Xpos(np)+Xoff
          Ygrd=Ypos(np)
          Zgrd=Zpos(np)
          bounded(np)=0.0_r8
          IF (((Xmin.le.Xgrd).and.(Xgrd.lt.Xmax)).and.                  &
     &        ((Ymin.le.Ygrd).and.(Ygrd.lt.Ymax))) THEN
            i1=INT(Xgrd)
            j1=INT(Ygrd)
            i2=i1+1
            j2=j1+1
            IF (i2.gt.Lm(ng)+1) THEN
              i2=i1                   ! station at the eastern boundary
            END IF
            IF (j2.gt.Mm(ng)+1) THEN
              j2=j1                   ! station at the northern boundary
            END IF
            bounded(np)=1.0_r8
            p2=REAL(i2-i1,r8)*(Xgrd-REAL(i1,r8))
            q2=REAL(j2-j1,r8)*(Ygrd-REAL(j1,r8))
            p1=1.0_r8-p2
            q1=1.0_r8-q2
            w111=p1*q1
            w211=p2*q1
            w121=p1*q2
            w221=p2*q2
            w112=0.0_r8
            w212=0.0_r8
            w122=0.0_r8
            w222=0.0_r8
            IF (Zgrd.ge.0.0_r8) THEN
              k1=INT(Zgrd)
              k2=INT(Zgrd)
              r1=1.0_r8
              r2=0.0_r8
            ELSE
              Ztop=0.5_r8*(GRID(ng)%z_r(i1-1,j1,N(ng))+                 &
     &                     GRID(ng)%z_r(i1  ,j1,N(ng)))
              Zbot=0.5_r8*(GRID(ng)%z_r(i1-1,j1,1)+                     &
     &                     GRID(ng)%z_r(i1  ,j1,1))
              IF (Zgrd.ge.Ztop) THEN
                k1=N(ng)                     ! If shallower, assign
                k2=N(ng)                     ! station to surface
                r1=1.0_r8                    ! level
                r2=0.0_r8
              ELSE IF (Zbot.ge.Zgrd) THEN
                k1=1                         ! If deeper, assign
                k2=1                         ! station to bottom
                r1=1.0_r8                    ! level
                r2=0.0_r8
              ELSE
                DO k=N(ng),2,-1
                  Ztop=0.5_r8*(GRID(ng)%z_r(i1-1,j1,k)+                 &
     &                         GRID(ng)%z_r(i1  ,j1,k))
                  Zbot=0.5_r8*(GRID(ng)%z_r(i1-1,j1,k-1)+               &
     &                         GRID(ng)%z_r(i1  ,j1,k-1))
                  IF ((Ztop.gt.Zgrd).and.(Zgrd.ge.Zbot)) THEN
                    k1=k-1
                    k2=k
                  END IF
                END DO
                dz=0.5_r8*((GRID(ng)%z_r(i1-1,j1,k2)+                   &
     &                      GRID(ng)%z_r(i1  ,j1,k2))-                  &
     &                     (GRID(ng)%z_r(i1-1,j1,k1)+                   &
     &                      GRID(ng)%z_r(i1  ,j1,k1)))
                r2=(Zgrd-0.5_r8*(GRID(ng)%z_r(i1-1,j1,k1)+              &
     &                           GRID(ng)%z_r(i1  ,j1,k1)))/dz
                r1=1.0_r8-r2
              END IF
            END IF
            w112=w111*r2
            w212=w211*r2
            w122=w121*r2
            w222=w221*r2
            w111=w111*r1
            w211=w211*r1
            w121=w121*r1
            w221=w221*r1
# ifdef MASKING
            w111=w111*GRID(ng)%umask(i1,j1)
            w211=w211*GRID(ng)%umask(i2,j1)
            w121=w121*GRID(ng)%umask(i1,j2)
            w221=w221*GRID(ng)%umask(i2,j2)
            w112=w112*GRID(ng)%umask(i1,j1)
            w212=w212*GRID(ng)%umask(i2,j1)
            w122=w122*GRID(ng)%umask(i1,j2)
            w222=w222*GRID(ng)%umask(i2,j2)
            wsum=w111+w211+w121+w221+w112+w212+w122+w222
            IF (wsum.gt.0.0_r8) THEN
              wsum=1.0_r8/wsum
              w111=w111*wsum
              w211=w211*wsum
              w121=w121*wsum
              w221=w221*wsum
              w112=w112*wsum
              w212=w212*wsum
              w122=w122*wsum
              w222=w222*wsum
            ELSE
              bounded(np)=0.0_r8
            END IF
# endif
            Apos(np)=Ascl*(w111*A(i1,j1,k1)+                            &
     &                     w211*A(i2,j1,k1)+                            &
     &                     w121*A(i1,j2,k1)+                            &
     &                     w221*A(i2,j2,k1)+                            &
     &                     w112*A(i1,j1,k2)+                            &
     &                     w212*A(i2,j1,k2)+                            &
     &                     w122*A(i1,j2,k2)+                            &
     &                     w222*A(i2,j2,k2))
            IF (ABS(Apos(np)).eq.0.0_r8) Apos(np)=0.0_r8      ! positive
          ELSE                                                ! zero
            Apos(np)=Aspv
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Interpolate from 3D field at V-points.
!-----------------------------------------------------------------------
!
      ELSE IF (gtype.eq.v3dvar) THEN
        IF (Cgrid) THEN
          Xmin=vXmin(ng)
          Xmax=vXmax(ng)
          Ymin=vYmin(ng)+0.5_r8
          Ymax=vYmax(ng)+0.5_r8
          Yoff=0.0_r8
        ELSE
          Xmin=rXmin(ng)
          Xmax=rXmax(ng)
          Ymin=rYmin(ng)
          Ymax=rYmax(ng)
          Yoff=0.5_r8
        END IF
        DO np=1,Npos
          Xgrd=Xpos(np)
          Ygrd=Ypos(np)+Yoff
          Zgrd=Zpos(np)
          bounded(np)=0.0_r8
          IF (((Xmin.le.Xgrd).and.(Xgrd.lt.Xmax)).and.                  &
     &        ((Ymin.le.Ygrd).and.(Ygrd.lt.Ymax))) THEN
            i1=INT(Xgrd)
            j1=INT(Ygrd)
            i2=i1+1
            j2=j1+1
            IF (i2.gt.Lm(ng)+1) THEN
              i2=i1                   ! station at the eastern boundary
            END IF
            IF (j2.gt.Mm(ng)+1) THEN
              j2=j1                   ! station at the northern boundary
            END IF
            bounded(np)=1.0_r8
            p2=REAL(i2-i1,r8)*(Xgrd-REAL(i1,r8))
            q2=REAL(j2-j1,r8)*(Ygrd-REAL(j1,r8))
            p1=1.0_r8-p2
            q1=1.0_r8-q2
            w111=p1*q1
            w211=p2*q1
            w121=p1*q2
            w221=p2*q2
            w112=0.0_r8
            w212=0.0_r8
            w122=0.0_r8
            w222=0.0_r8
            IF (Zgrd.ge.0.0_r8) THEN
              k1=INT(Zgrd)
              k2=INT(Zgrd)
              r1=1.0_r8
              r2=0.0_r8
            ELSE
              Ztop=0.5_r8*(GRID(ng)%z_r(i1,j1-1,N(ng))+                 &
     &                     GRID(ng)%z_r(i1,j1,  N(ng)))
              Zbot=0.5_r8*(GRID(ng)%z_r(i1,j1-1,1)+                     &
     &                     GRID(ng)%z_r(i1,j1  ,1))
              IF (Zgrd.ge.Ztop) THEN
                k1=N(ng)                     ! If shallower, assign
                k2=N(ng)                     ! station to surface
                r1=1.0_r8                    ! level
                r2=0.0_r8
              ELSE IF (Zbot.ge.Zgrd) THEN
                k1=1                         ! If deeper, assign
                k2=1                         ! station to bottom
                r1=1.0_r8                    ! level
                r2=0.0_r8
              ELSE
                DO k=N(ng),2,-1
                  Ztop=0.5_r8*(GRID(ng)%z_r(i1,j1-1,k)+                 &
     &                         GRID(ng)%z_r(i1,j1  ,k))
                  Zbot=0.5_r8*(GRID(ng)%z_r(i1,j1-1,k-1)+               &
     &                         GRID(ng)%z_r(i1,j1  ,k-1))
                  IF ((Ztop.gt.Zgrd).and.(Zgrd.ge.Zbot)) THEN
                    k1=k-1
                    k2=k
                  END IF
                END DO
                dz=0.5_r8*((GRID(ng)%z_r(i1,j1-1,k2)+                   &
     &                      GRID(ng)%z_r(i1,j1  ,k2))-                  &
     &                     (GRID(ng)%z_r(i1,j1-1,k1)+                   &
     &                      GRID(ng)%z_r(i1,j1  ,k1)))
                r2=(Zgrd-0.5_r8*(GRID(ng)%z_r(i1,j1-1,k1)+              &
     &                           GRID(ng)%z_r(i1,j1  ,k1)))/dz
                r1=1.0_r8-r2
              END IF
            END IF
            w112=w111*r2
            w212=w211*r2
            w122=w121*r2
            w222=w221*r2
            w111=w111*r1
            w211=w211*r1
            w121=w121*r1
            w221=w221*r1
# ifdef MASKING
            w111=w111*GRID(ng)%vmask(i1,j1)
            w211=w211*GRID(ng)%vmask(i2,j1)
            w121=w121*GRID(ng)%vmask(i1,j2)
            w221=w221*GRID(ng)%vmask(i2,j2)
            w112=w112*GRID(ng)%vmask(i1,j1)
            w212=w212*GRID(ng)%vmask(i2,j1)
            w122=w122*GRID(ng)%vmask(i1,j2)
            w222=w222*GRID(ng)%vmask(i2,j2)
            wsum=w111+w211+w121+w221+w112+w212+w122+w222
            IF (wsum.gt.0.0_r8) THEN
              wsum=1.0_r8/wsum
              w111=w111*wsum
              w211=w211*wsum
              w121=w121*wsum
              w221=w221*wsum
              w112=w112*wsum
              w212=w212*wsum
              w122=w122*wsum
              w222=w222*wsum
            ELSE
              bounded(np)=0.0_r8
            END IF
# endif
            Apos(np)=Ascl*(w111*A(i1,j1,k1)+                            &
     &                     w211*A(i2,j1,k1)+                            &
     &                     w121*A(i1,j2,k1)+                            &
     &                     w221*A(i2,j2,k1)+                            &
     &                     w112*A(i1,j1,k2)+                            &
     &                     w212*A(i2,j1,k2)+                            &
     &                     w122*A(i1,j2,k2)+                            &
     &                     w222*A(i2,j2,k2))
            IF (ABS(Apos(np)).eq.0.0_r8) Apos(np)=0.0_r8      ! positive
          ELSE                                                ! zero
            Apos(np)=Aspv
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Interpolate from 3D field at W-points.
!-----------------------------------------------------------------------
!
      ELSE IF (gtype.eq.w3dvar) THEN
        Xmin=rXmin(ng)
        Xmax=rXmax(ng)
        Ymin=rYmin(ng)
        Ymax=rYmax(ng)
        DO np=1,Npos
          Xgrd=Xpos(np)
          Ygrd=Ypos(np)
          Zgrd=Zpos(np)
          bounded(np)=0.0_r8
          IF (((Xmin.le.Xgrd).and.(Xgrd.lt.Xmax)).and.                  &
     &        ((Ymin.le.Ygrd).and.(Ygrd.lt.Ymax))) THEN
            i1=INT(Xgrd)
            j1=INT(Ygrd)
            i2=i1+1
            j2=j1+1
            IF (i2.gt.Lm(ng)+1) THEN
              i2=i1                   ! station at the eastern boundary
            END IF
            IF (j2.gt.Mm(ng)+1) THEN
              j2=j1                   ! station at the northern boundary
            END IF
            bounded(np)=1.0_r8
            p2=REAL(i2-i1,r8)*(Xgrd-REAL(i1,r8))
            q2=REAL(j2-j1,r8)*(Ygrd-REAL(j1,r8))
            p1=1.0_r8-p2
            q1=1.0_r8-q2
            w111=p1*q1
            w211=p2*q1
            w121=p1*q2
            w221=p2*q2
            w112=0.0_r8
            w212=0.0_r8
            w122=0.0_r8
            w222=0.0_r8
            IF (Zgrd.ge.0.0_r8) THEN
              k1=INT(Zgrd)
              k2=INT(Zgrd)
              r1=1.0_r8
              r2=0.0_r8
            ELSE
              Ztop=GRID(ng)%z_w(i1,j1,N(ng))
              Zbot=GRID(ng)%z_w(i1,j1,0)
              IF (Zgrd.ge.Ztop) THEN
                k1=N(ng)                     ! If shallower, assign
                k2=N(ng)                     ! station to surface
                r1=1.0_r8                    ! level
                r2=0.0_r8
              ELSE IF (Zbot.ge.Zgrd) THEN
                k1=0                         ! If deeper, assign
                k2=0                         ! station to bottom
                r1=1.0_r8                    ! level
                r2=0.0_r8
              ELSE
                DO k=N(ng),2,-1
                  Ztop=GRID(ng)%z_w(i1,j1,k)
                  Zbot=GRID(ng)%z_w(i1,j1,k-1)
                  IF ((Ztop.gt.Zgrd).and.(Zgrd.ge.Zbot)) THEN
                    k1=k-1
                    k2=k
                  END IF
                END DO
                dz=GRID(ng)%z_w(i1,j1,k2)-GRID(ng)%z_w(i1,j1,k1)
                r2=(Zgrd-GRID(ng)%z_w(i1,j1,k1))/dz
                r1=1.0_r8-r2
              END IF
            END IF
            w112=w111*r2
            w212=w211*r2
            w122=w121*r2
            w222=w221*r2
            w111=w111*r1
            w211=w211*r1
            w121=w121*r1
            w221=w221*r1
# ifdef MASKING
            w111=w111*GRID(ng)%rmask(i1,j1)
            w211=w211*GRID(ng)%rmask(i2,j1)
            w121=w121*GRID(ng)%rmask(i1,j2)
            w221=w221*GRID(ng)%rmask(i2,j2)
            w112=w112*GRID(ng)%rmask(i1,j1)
            w212=w212*GRID(ng)%rmask(i2,j1)
            w122=w122*GRID(ng)%rmask(i1,j2)
            w222=w222*GRID(ng)%rmask(i2,j2)
            wsum=w111+w211+w121+w221+w112+w212+w122+w222
            IF (wsum.gt.0.0_r8) THEN
              wsum=1.0_r8/wsum
              w111=w111*wsum
              w211=w211*wsum
              w121=w121*wsum
              w221=w221*wsum
              w112=w112*wsum
              w212=w212*wsum
              w122=w122*wsum
              w222=w222*wsum
            ELSE
              bounded(np)=0.0_r8
            END IF
# endif
            Apos(np)=Ascl*(w111*A(i1,j1,k1)+                            &
     &                     w211*A(i2,j1,k1)+                            &
     &                     w121*A(i1,j2,k1)+                            &
     &                     w221*A(i2,j2,k1)+                            &
     &                     w112*A(i1,j1,k2)+                            &
     &                     w212*A(i2,j1,k2)+                            &
     &                     w122*A(i1,j2,k2)+                            &
     &                     w222*A(i2,j2,k2))
            IF (ABS(Apos(np)).eq.0.0_r8) Apos(np)=0.0_r8      ! positive
          ELSE                                                ! zero
            Apos(np)=Aspv
          END IF
        END DO
      END IF
# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Collect all extracted data.
!-----------------------------------------------------------------------
!
      CALL mp_collect (ng, model, Npos, Aspv, Apos)
      CALL mp_collect (ng, model, Npos, 0.0_r8, bounded)
# endif
!
!-----------------------------------------------------------------------
!  Set unbounded data to special value.
!-----------------------------------------------------------------------
!
      DO np=1,Npos
        IF (bounded(np).lt.1.0_r8) THEN
          Apos(np)=spval
        END IF
      END DO
      RETURN
      END SUBROUTINE extract_sta3d
#endif
      END MODULE extract_sta_mod

