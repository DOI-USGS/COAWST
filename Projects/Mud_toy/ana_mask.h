      SUBROUTINE ana_mask (ng, tile, model)
!
!! svn $Id: ana_mask.h 727 2007-04-27 05:03:38Z jcwarner $
!!======================================================================
!! Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
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
      CALL ana_mask_tile (ng, model, Istr, Iend, Jstr, Jend,            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % pmask,                             &
     &                    GRID(ng) % rmask,                             &
     &                    GRID(ng) % umask,                             &
     &                    GRID(ng) % vmask)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        ANANAME(15)='ROMS/Functionals/ana_mask.h'
      END IF

      RETURN
      END SUBROUTINE ana_mask
!
!***********************************************************************
      SUBROUTINE ana_mask_tile (ng, model, Istr, Iend, Jstr, Jend,      &
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
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
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
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
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
#ifdef MUD_TOY
      DO j=Jstr-2,Jend+2
        DO i=Istr-2,Iend+2
          mask(i,j)=1.0_r8
        END DO
      END DO
      DO j=Jstr-2,Jend+2
        mask(Istr-0,j)=0.0_r8
        mask(Istr-1,j)=0.0_r8
        mask(Istr-2,j)=0.0_r8
        mask(Iend+1,j)=0.0_r8
        mask(Iend+2,j)=0.0_r8
      END DO
      DO i=Istr-2,Iend+2
        mask(i,Jstr-0)=0.0_r8
        mask(i,Jstr-1)=0.0_r8
        mask(i,Jstr-2)=0.0_r8
        mask(i,Jend+1)=0.0_r8
        mask(i,Jend+2)=0.0_r8
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
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        rmask)
      CALL exchange_p2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        pmask)
      CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        umask)
      CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        vmask)
#endif
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, model, 4, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    rmask, pmask, umask, vmask)
#endif
      RETURN
      END SUBROUTINE ana_mask_tile
