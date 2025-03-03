#include "cppdefs.h"
#if defined FLOATS || defined STATIONS
       SUBROUTINE grid_coords (ng, model)
!
!git $Id$
!svn $Id: grid_coords.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine converts initial locations to fractional grid (I,J)    !
!  coordinates, if appropriate.                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef FLOATS
      USE mod_floats
# endif
      USE mod_grid
      USE mod_scalars
      USE roms_interpolate_mod
!
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_collect
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer :: IstrR, Iend, JstrR, Jend
      integer :: LBi, UBi, LBj, UBj
      integer :: i, j, k, l, mc

      real(r8), parameter :: spv = 0.0_r8

# ifdef FLOATS
      real(r8) :: Xstr, Xend, Ystr, Yend, zfloat
      logical, dimension(Nfloats(ng)) :: my_thread
      real(r8), dimension(Nfloats(ng)) :: Iflt, Jflt
#  ifdef SOLVE3D
      real(r8), dimension(Nfloats(ng)) :: Kflt
#  endif
# endif

# ifdef STATIONS
      real(r8), dimension(Nstation(ng)) :: Slon, Slat
      real(r8), dimension(Nstation(ng)) :: Ista, Jsta
# endif
!
!-----------------------------------------------------------------------
!  Determine searching model grid box and arrays bounds.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
      IstrR=BOUNDS(ng)%IstrR(MyRank)
      Iend =BOUNDS(ng)%Iend (MyRank)
      JstrR=BOUNDS(ng)%JstrR(MyRank)
      Jend =BOUNDS(ng)%Jend (MyRank)
# else
      IstrR=0
      Iend =Lm(ng)
      JstrR=0
      Jend =Mm(ng)
# endif
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)

# ifdef FLOATS
!
      Xstr=REAL(BOUNDS(ng)%Istr(MyRank),r8)-0.5_r8
      Xend=REAL(BOUNDS(ng)%Iend(MyRank),r8)+0.5_r8
      Ystr=REAL(BOUNDS(ng)%Jstr(MyRank),r8)-0.5_r8
      Yend=REAL(BOUNDS(ng)%Jend(MyRank),r8)+0.5_r8
!
!-----------------------------------------------------------------------
!  If applicable, convert initial floats locations (Flon,Flat) to
!  fractional grid coordinates.
!-----------------------------------------------------------------------
!
      IF (spherical) THEN
        IF (Lfloats(ng)) THEN
          mc=DRIFTER(ng)%Findex(0)
          IF (DRIFTER(ng)%Findex(0).gt.0) THEN
            CALL hindices (ng, LBi, UBi, LBj, UBj,                      &
     &                     IstrR, Iend+1, JstrR, Jend+1,                &
     &                     GRID(ng)%angler,                             &
     &                     GRID(ng)%lonr,                               &
     &                     GRID(ng)%latr,                               &
     &                     1, mc, 1, 1,                                 &
     &                     1, mc, 1, 1,                                 &
     &                     DRIFTER(ng)%Flon,                            &
     &                     DRIFTER(ng)%Flat,                            &
     &                     Iflt, Jflt, spv, .FALSE.)
#  ifdef DISTRIBUTE
            CALL mp_collect (ng, model, mc, spv, Iflt)
            CALL mp_collect (ng, model, mc, spv, Jflt)
#  endif
            DO i=1,mc
              l=DRIFTER(ng)%Findex(i)
              DRIFTER(ng)%Tinfo(ixgrd,l)=Iflt(i)
              DRIFTER(ng)%Tinfo(iygrd,l)=Jflt(i)
            END DO
          END IF
        END IF
      END IF

#  ifdef SOLVE3D
!
!  Determine which node bounds the initial float location.
!
#   ifdef DISTRIBUTE
      IF (Lfloats(ng)) THEN
        DO l=1,Nfloats(ng)
          IF ((Xstr.le.DRIFTER(ng)%Tinfo(ixgrd,l)).and.                 &
     &        (DRIFTER(ng)%Tinfo(ixgrd,l).lt.Xend).and.                 &
     &        (Ystr.le.DRIFTER(ng)%Tinfo(iygrd,l)).and.                 &
     &        (DRIFTER(ng)%Tinfo(iygrd,l).lt.Yend)) THEN
            my_thread(l)=.TRUE.
          ELSE
            my_thread(l)=.FALSE.
          END IF
        END DO
      END IF
#   else
      DO l=1,Nfloats(ng)
        my_thread(l)=.TRUE.
      END DO
#   endif
#  endif
!
!-----------------------------------------------------------------------
!  Set float initial vertical level position, if inside application
!  grid.  If the initial float depth (in meters) is not found, release
!  float at the surface model level.
!-----------------------------------------------------------------------
!
      DO l=1,Nfloats(ng)
        IF (Lfloats(ng)) THEN
#  ifdef SOLVE3D
          DRIFTER(ng)%Fz0(l)=spv
          IF (my_thread(l).and.                                         &
     &        ((DRIFTER(ng)%Tinfo(ixgrd,l).ge.0.5_r8).and.              &
     &         (DRIFTER(ng)%Tinfo(iygrd,l).ge.0.5_r8).and.              &
     &         (DRIFTER(ng)%Tinfo(ixgrd,l).le.                          &
     &          REAL(Lm(ng),r8)+0.5_r8).and.                            &
     &         (DRIFTER(ng)%Tinfo(iygrd,l).le.                          &
     &          REAL(Mm(ng),r8)+0.5_r8))) THEN
            zfloat=DRIFTER(ng)%Tinfo(izgrd,l)
            DRIFTER(ng)%Fz0(l)=zfloat           ! Save original value
            Kflt(l)=zfloat
            IF (zfloat.le.0.0_r8) THEN
              i=INT(DRIFTER(ng)%Tinfo(ixgrd,l)) ! Fractional positions
              j=INT(DRIFTER(ng)%Tinfo(iygrd,l)) ! are still in this cell
              IF (zfloat.lt.GRID(ng)%z_w(i,j,0)) THEN
                zfloat=GRID(ng)%z_w(i,j,0)+5.0_r8
                DRIFTER(ng)%Fz0(l)=zfloat
              END IF
              DRIFTER(ng)%Tinfo(izgrd,l)=REAL(N(ng),r8)
              DO k=N(ng),1,-1
                IF ((GRID(ng)%z_w(i,j,k)-zfloat)*                       &
     &              (zfloat-GRID(ng)%z_w(i,j,k-1)).ge.0.0_r8) THEN
                  Kflt(l)=REAL(k-1,r8)+                                 &
     &                    (zfloat-GRID(ng)%z_w(i,j,k-1))/               &
     &                    GRID(ng)%Hz(i,j,k)
                END IF
              END DO
            END IF
          ELSE
            Kflt(l)=spv
          END IF
#  else
          DRIFTER(ng)%Tinfo(izgrd,l)=0.0_r8
#  endif
        END IF
      END DO
#  ifdef SOLVE3D
      IF (Lfloats(ng)) THEN
#   ifdef DISTRIBUTE
        CALL mp_collect (ng, model, Nfloats(ng), spv, DRIFTER(ng)%Fz0)
        CALL mp_collect (ng, model, Nfloats(ng), spv, Kflt)
#   endif
        DO l=1,Nfloats(ng)
          DRIFTER(ng)%Tinfo(izgrd,l)=Kflt(l)
        END DO
      END IF
#  endif
# endif
# ifdef STATIONS
!
!-----------------------------------------------------------------------
!  If applicable, convert station locations (SposX,SposY) to fractional
!  grid coordinates.
!-----------------------------------------------------------------------
!
      IF (spherical) THEN
        mc=0
        DO l=1,Nstation(ng)
          IF (SCALARS(ng)%Sflag(l).gt.0) THEN
            mc=mc+1
            Slon(mc)=SCALARS(ng)%SposX(l)
            Slat(mc)=SCALARS(ng)%SposY(l)
          END IF
        END DO
        IF (mc.gt.0) THEN
          CALL hindices (ng, LBi, UBi, LBj, UBj,                        &
     &                   IstrR, Iend+1, JstrR, Jend+1,                  &
     &                   GRID(ng)%angler,                               &
     &                   GRID(ng)%lonr,                                 &
     &                   GRID(ng)%latr,                                 &
     &                   1, mc, 1, 1,                                   &
     &                   1, mc, 1, 1,                                   &
     &                   Slon, Slat,                                    &
     &                   Ista, Jsta,                                    &
     &                   spv, .FALSE.)
#  ifdef DISTRIBUTE
          CALL mp_collect (ng, model, mc, spv, Ista)
          CALL mp_collect (ng, model, mc, spv, Jsta)
#  endif
          mc=0
          DO l=1,Nstation(ng)
            IF (SCALARS(ng)%Sflag(l).gt.0) THEN
              mc=mc+1
              SCALARS(ng)%SposX(l)=Ista(mc)
              SCALARS(ng)%SposY(l)=Jsta(mc)
            END IF
          END DO
        END IF
      END IF
# endif
      RETURN
      END SUBROUTINE grid_coords
#else
      SUBROUTINE grid_coords
      RETURN
      END SUBROUTINE grid_coords
#endif
