      SUBROUTINE ana_mask (ng, tile, model)
!
!! svn $Id: ana_mask.h 737 2008-09-07 02:06:44Z jcwarner $
!!======================================================================
!! Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine sets analytical Land/Sea masking.                   !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_mask_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % pmask,                             &
     &                    GRID(ng) % rmask,                             &
     &                    GRID(ng) % umask,                             &
     &                    GRID(ng) % vmask)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        ANANAME(15)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_mask
!
!***********************************************************************
      SUBROUTINE ana_mask_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmask, rmask, umask, vmask)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
#if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
#endif
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: pmask(LBi:,LBj:)
      real(r8), intent(out) :: rmask(LBi:,LBj:)
      real(r8), intent(out) :: umask(LBi:,LBj:)
      real(r8), intent(out) :: vmask(LBi:,LBj:)
#else
      real(r8), intent(out) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: vmask(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
#ifdef DISTRIBUTE
# ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
# else
      logical :: EWperiodic=.FALSE.
# endif
# ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
# else
      logical :: NSperiodic=.FALSE.
# endif
#endif
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j
      real(r8) :: mask(PRIVATE_2D_SCRATCH_ARRAY)

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set Land/Sea mask of RHO-points: Land=0, Sea=1.
!-----------------------------------------------------------------------
!
!  Notice that private scratch array "mask" is used to allow
!  computation within a parallel loop.
!
#ifdef DOUBLE_GYRE
      Imin=-2+(Lm(ng)+1)/2
      Imax=Imin+2
      Jmin=-2+(Mm(ng)+1)/2
      Jmax=Jmin+2
      DO j=Jstr-2,Jend+2
        DO i=Istr-2,Iend+2
          mask(i,j)=1.0_r8
          IF (((Imin.le.i).and.(i.le.Imax)).and.                        &
     &        ((Jmin.le.j).and.(j.le.Jmax))) THEN
            mask(i,j)=0.0_r8
          END IF
        END DO
      END DO
#elif defined FLT_TEST
      DO j=Jstr-2,Jend+2
        DO i=Istr-2,Iend+2
          mask(i,j)=1.0_r8
          IF (j.eq.1 ) mask(i,j)=0.0_r8
          IF (j.eq.Mm(ng)) mask(i,j)=0.0_r8
          IF ((i.ge.((Lm(ng)+1)/2)).and.                                &
     &        (i.le.((Lm(ng)+1)/2+1)).and.                              &
     &        (j.ge.((Mm(ng)+1)/2)).and.                                &
     &        (j.le.((Mm(ng)+1)/2+1))) mask(i,j)=0.0_r8
        END DO
      END DO
#elif defined LAKE_SIGNELL
      DO j=Jstr-2,Jend+2
        DO i=Istr-2,Iend+2
          mask(i,j)=1.0_r8
        END DO
      END DO
      IF (WESTERN_EDGE) THEN
        DO j=Jstr-1,Jend+1
          mask(Istr-1,j)=0.0_r8
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO j=Jstr-1,Jend+1
          mask(Iend+1,j)=0.0_r8
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO i=Istr-1,Iend+1
          mask(i,Jstr-1)=0.0_r8
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO i=Istr-1,Iend+1
          mask(i,Jend+1)=0.0_r8
        END DO
      END IF
#elif defined RIVERPLUME1
      DO j=Jstr-2,Jend+2
        DO i=Istr-2,Iend+2
          mask(i,j)=1.0_r8
        END DO
      END DO
      DO i=Istr-2,MIN(5,Iend+2)
        DO j=Jstr-2,MIN(Mm(ng)-18,Jend+2)
          mask(i,j)=0.0_r8
        END DO
        DO j=MAX(Jstr-2,Mm(ng)-16),Jend+2
          mask(i,j)=0.0_r8
        END DO
      END DO
#elif defined RIVERPLUME2
      DO j=Jstr-2,Jend+2
        DO i=Istr-2,Iend+2
          mask(i,j)=1.0_r8
        END DO
      END DO
      DO i=Istr-2,MIN(5,Iend+2)
        DO j=Jstr-2,MIN(Mm(ng)-11,Jend+2)
          mask(i,j)=0.0_r8
        END DO
        DO j=MAX(Jstr-2,Mm(ng)-9),Jend+2
          mask(i,j)=0.0_r8
        END DO
      END DO
#elif defined SHOREFACE
      DO j=Jstr-2,Jend+2
        DO i=Istr-2,Iend+2
          mask(i,j)=1.0_r8
        END DO
      END DO
#else
      ana_mask.h: no values provided for mask.
#endif
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          rmask(i,j)=mask(i,j)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute Land/Sea mask of U- and V-points.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          umask(i,j)=mask(i-1,j)*mask(i,j)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          vmask(i,j)=mask(i,j-1)*mask(i,j)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute Land/Sea mask of PSI-points.
!-----------------------------------------------------------------------
!
      DO j=Jstr,JendR
        DO i=Istr,IendR
          pmask(i,j)=mask(i-1,j-1)*mask(i,j-1)*                         &
     &               mask(i-1,j  )*mask(i,j  )
        END DO
      END DO
#if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        rmask)
      CALL exchange_p2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        pmask)
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        umask)
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        vmask)
#endif
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    rmask, pmask, umask, vmask)
#endif
      RETURN
      END SUBROUTINE ana_mask_tile
