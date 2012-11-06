#include "cppdefs.h"
      MODULE set_masks_mod
#ifdef MASKING
!
!svn $Id: set_masks.F 585 2012-01-03 18:44:28Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines set internal Land/Sea masking arrays that are used   !
!  to process fields into output NetCDF files.  The Land grid points   !
!  are replaced by the _FillValue in the output files to  facilitate   !
!  post-processing with generic tools.                                 !
!                                                                      !
# ifdef UV_PSOURCE
!  If point sources, insure that masks at point source locations are   !
!  set to water to avoid masking with _FillValue at those locations.   !
# endif
# ifdef WET_DRY
!
!  If wetting and drying, masks array are time dependent and changed   !
!  at every time-step in routine "wetdry".  Notice that time-average   !
!  masks are needed for the selected time window.
# endif
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: set_masks
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_masks (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
# ifdef UV_PSOURCE
      USE mod_sources
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 2)
# endif
      CALL set_masks_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
# ifdef UV_PSOURCE
     &                     Msrc(ng), Nsrc(ng),                          &
     &                     SOURCES(ng) % Isrc,                          &
     &                     SOURCES(ng) % Jsrc,                          &
     &                     SOURCES(ng) % Dsrc,                          &
# endif
     &                     GRID(ng) % pmask,                            &
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
     &                     GRID(ng) % pmask_io,                         &
     &                     GRID(ng) % rmask_io,                         &
     &                     GRID(ng) % umask_io,                         &
     &                     GRID(ng) % vmask_io)
# ifdef PROFILE
      CALL wclock_off (ng, model, 2)
# endif

      RETURN
      END SUBROUTINE set_masks
!
!***********************************************************************
      SUBROUTINE set_masks_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
# ifdef UV_PSOURCE
     &                           Msrc, Nsrc,                            &
     &                           Isrc, Jsrc, Dsrc,                      &
# endif
     &                           pmask, rmask,                          &
     &                           umask, vmask,                          &
     &                           pmask_io, rmask_io,                    &
     &                           umask_io, vmask_io)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
# ifdef UV_PSOURCE
      integer, intent(in) :: Msrc, Nsrc
# endif
!
# ifdef ASSUMED_SHAPE
#  ifdef UV_PSOURCE
      integer, intent(in) :: Isrc(:)
      integer, intent(in) :: Jsrc(:)

      real(r8), intent(in) :: Dsrc(:)
#  endif
      real(r8), intent(in) :: pmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(inout) :: pmask_io(LBi:,LBj:)
      real(r8), intent(inout) :: rmask_io(LBi:,LBj:)
      real(r8), intent(inout) :: umask_io(LBi:,LBj:)
      real(r8), intent(inout) :: vmask_io(LBi:,LBj:)
# else
#  ifdef UV_PSOURCE
      integer, intent(in) :: Isrc(Msrc)
      integer, intent(in) :: Jsrc(Msrc)

      real(r8), intent(in) :: Dsrc(Msrc)
#  endif
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: pmask_io(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rmask_io(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: umask_io(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vmask_io(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
# ifdef UV_PSOURCE
      integer :: is
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize internal history files Land/Sea masks with its respective
!  application grid mask.
!-----------------------------------------------------------------------
!
      DO j=Jstr,JendR
        DO i=Istr,IendR
          pmask_io(i,j)=pmask(i,j)
        END DO
      END DO
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          rmask_io(i,j)=rmask(i,j)
        END DO
      END DO
      DO j=JstrR,JendR
        DO i=Istr,IendR
          umask_io(i,j)=umask(i,j)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          vmask_io(i,j)=vmask(i,j)
        END DO
      END DO

# ifdef UV_PSOURCE
!
!  Insure that masks at mass point source locations are set to water
!  to avoid masking with _FillValue at those locations.
!
      DO is=1,Nsrc
        i=Isrc(is)
        j=Jsrc(is)
        IF (((IstrR.le.i).and.(i.le.IendR)).and.                        &
     &      ((JstrR.le.j).and.(j.le.JendR))) THEN
          IF (INT(Dsrc(is)).eq.0) THEN
            umask_io(i,j)=1.0_r8
          ELSE
            vmask_io(i,j)=1.0_r8
          END IF
        END IF
      END DO
# endif
!
#  if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_p2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        pmask_io)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        rmask_io)
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        umask_io)
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        vmask_io)
# endif
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    pmask_io, rmask_io, umask_io, vmask_io)
# endif

      RETURN
      END SUBROUTINE set_masks_tile
#endif
      END MODULE set_masks_mod
