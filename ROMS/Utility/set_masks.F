#include "cppdefs.h"
      MODULE set_masks_mod
#ifdef MASKING
!
!svn $Id: set_masks.F 585 2012-01-03 18:44:28Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines set internal Land/Sea masking arrays that are used   !
!  to process fields into output NetCDF files.  The Land grid points   !
!  are replaced by the _FillValue in the output files to  facilitate   !
!  post-processing with generic tools.                                 !
!                                                                      !
!  If point sources, insure that masks at point source locations are   !
!  set to water to avoid masking with _FillValue at those locations.   !
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
# if defined WET_DRY     && \
     defined AVERAGES    || \
    (defined AD_AVERAGES && defined ADJOINT) || \
    (defined RP_AVERAGES && defined TL_IOMS) || \
    (defined TL_AVERAGES && defined TANGENT)
      PUBLIC :: set_avg_masks
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_masks (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
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
     &                     GRID(ng) % pmask,                            &
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
# if defined PERFECT_RESTART && defined WET_DRY
     &                     GRID(ng) % rmask_wet,                        &
     &                     GRID(ng) % umask_wet,                        &
     &                     GRID(ng) % vmask_wet,                        &
# endif
# ifdef OUTFLOW_MASK
     &                     GRID(ng) % mask_outflow,                     &
# endif
# if defined AVERAGES    || \
    (defined AD_AVERAGES && defined ADJOINT) || \
    (defined RP_AVERAGES && defined TL_IOMS) || \
    (defined TL_AVERAGES && defined TANGENT)
     &                     GRID(ng) % pmask_avg,                        &
     &                     GRID(ng) % rmask_avg,                        &
     &                     GRID(ng) % umask_avg,                        &
     &                     GRID(ng) % vmask_avg,                        &
# endif
# if defined AVERAGES2
     &                     GRID(ng) % pmask_avg2,                       &
     &                     GRID(ng) % rmask_avg2,                       &
     &                     GRID(ng) % umask_avg2,                       &
     &                     GRID(ng) % vmask_avg2,                       &
# endif
# ifdef DIAGNOSTICS
     &                     GRID(ng) % pmask_dia,                        &
     &                     GRID(ng) % rmask_dia,                        &
     &                     GRID(ng) % umask_dia,                        &
     &                     GRID(ng) % vmask_dia,                        &
# endif
     &                     GRID(ng) % pmask_full,                       &
     &                     GRID(ng) % rmask_full,                       &
     &                     GRID(ng) % umask_full,                       &
     &                     GRID(ng) % vmask_full)
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
     &                           pmask, rmask,                          &
     &                           umask, vmask,                          &
# if defined PERFECT_RESTART && defined WET_DRY
     &                           rmask_wet,                             &
     &                           umask_wet, vmask_wet,                  &
# endif
# ifdef OUTFLOW_MASK
     &                           mask_outflow,                          &
# endif
# if defined AVERAGES    || \
    (defined AD_AVERAGES && defined ADJOINT) || \
    (defined RP_AVERAGES && defined TL_IOMS) || \
    (defined TL_AVERAGES && defined TANGENT)
     &                           pmask_avg, rmask_avg,                  &
     &                           umask_avg, vmask_avg,                  &
# endif
# if defined AVERAGES2
     &                           pmask_avg2, rmask_avg2,                &
     &                           umask_avg2, vmask_avg2,                &
# endif
# ifdef DIAGNOSTICS
     &                           pmask_dia, rmask_dia,                  &
     &                           umask_dia, vmask_dia,                  &
# endif
     &                           pmask_full, rmask_full,                &
     &                           umask_full, vmask_full)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sources
!
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  if defined PERFECT_RESTART && defined WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
#  endif
#  ifdef OUTFLOW_MASK
      real(r8), intent(out) :: mask_outflow(LBi:,LBj:)
#  endif
#  if defined AVERAGES    || \
     (defined AD_AVERAGES && defined ADJOINT) || \
     (defined RP_AVERAGES && defined TL_IOMS) || \
     (defined TL_AVERAGES && defined TANGENT)
      real(r8), intent(inout) :: pmask_avg(LBi:,LBj:)
      real(r8), intent(inout) :: rmask_avg(LBi:,LBj:)
      real(r8), intent(inout) :: umask_avg(LBi:,LBj:)
      real(r8), intent(inout) :: vmask_avg(LBi:,LBj:)
#  endif
#  if defined AVERAGES2
      real(r8), intent(inout) :: pmask_avg2(LBi:,LBj:)
      real(r8), intent(inout) :: rmask_avg2(LBi:,LBj:)
      real(r8), intent(inout) :: umask_avg2(LBi:,LBj:)
      real(r8), intent(inout) :: vmask_avg2(LBi:,LBj:)
#  endif
#  ifdef DIAGNOSTICS
      real(r8), intent(inout) :: pmask_dia(LBi:,LBj:)
      real(r8), intent(inout) :: rmask_dia(LBi:,LBj:)
      real(r8), intent(inout) :: umask_dia(LBi:,LBj:)
      real(r8), intent(inout) :: vmask_dia(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: pmask_full(LBi:,LBj:)
      real(r8), intent(inout) :: rmask_full(LBi:,LBj:)
      real(r8), intent(inout) :: umask_full(LBi:,LBj:)
      real(r8), intent(inout) :: vmask_full(LBi:,LBj:)
# else
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  if defined PERFECT_RESTART && defined WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
#  endif
#  ifdef OUTFLOW_MASK
      real(r8), intent(out) :: mask_outflow(LBi:UBi,LBj:UBj)
#  endif
#  if defined AVERAGES    || \
     (defined AD_AVERAGES && defined ADJOINT) || \
     (defined RP_AVERAGES && defined TL_IOMS) || \
     (defined TL_AVERAGES && defined TANGENT)
      real(r8), intent(inout) :: pmask_avg(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rmask_avg(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: umask_avg(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vmask_avg(LBi:UBi,LBj:UBj)
#  endif
#  if defined AVERAGES2
      real(r8), intent(inout) :: pmask_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rmask_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: umask_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vmask_avg2(LBi:UBi,LBj:UBj)
#  endif
#  ifdef DIAGNOSTICS
      real(r8), intent(inout) :: pmask_dia(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rmask_dia(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: umask_dia(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vmask_dia(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: pmask_full(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rmask_full(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: umask_full(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vmask_full(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, is, j
# ifdef OUTFLOW_MASK
      real(r8) :: cff
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize internal full Land/Sea masks with its respective
!  application time-indrpendent values.
!-----------------------------------------------------------------------
!
!  The full mask values are updated with time-dependent values in
!  file "wetdry.F" if wetting and drying is activated.
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          pmask_full(i,j)=pmask(i,j)
        END DO
      END DO
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          rmask_full(i,j)=rmask(i,j)
        END DO
      END DO
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          umask_full(i,j)=umask(i,j)
        END DO
      END DO
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          vmask_full(i,j)=vmask(i,j)
        END DO
      END DO
!# endif

# ifdef OUTFLOW_MASK
!
!  Set outflow mask for the ice model.
!  The ice strength, etc. get multiplied by this mask.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          mask_outflow(i,j)=1.0_r8
        END DO
      END DO
#  ifdef CHUKCHI
      cff = 0.95_r8
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=JstrR,JendR
          mask_outflow(0,j) = cff
          mask_outflow(1,j) = cff
        END DO
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=JstrR,JendR
          mask_outflow(Lm(ng),j) = cff
          mask_outflow(Lm(ng)+1,j) = cff
        END DO
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=IstrR,IendR
          mask_outflow(i,0) = cff
          mask_outflow(i,1) = cff
        END DO
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=IstrR,IendR
          mask_outflow(i,Mm(ng)) = cff
          mask_outflow(i,Mm(ng)+1) = cff
        END DO
      END IF
#  elif defined ARCTIC
! Just off Greenland
      cff = 0.5_r8
! Old ARCTIC domain
!      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
!        DO j=MAX(JstrR,Mm(ng)/2),JendR
!          mask_outflow(Lm(ng),j) = cff
!          mask_outflow(Lm(ng)+1,j) = cff
!        END DO
!      END IF
! New ARCTIC domain
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=IstrR,IendR
          mask_outflow(i,Mm(ng)) = cff
          mask_outflow(i,Mm(ng)+1) = cff
        END DO
      END IF
#  elif defined SVALBARD
      cff = 0.5_r8
! We want outflow south of Greenland but not north of it.
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=IstrR,MIN(IendR,65)
          mask_outflow(i,Mm(ng)) = cff
          mask_outflow(i,Mm(ng)+1) = cff
        END DO
      END IF
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=JstrR,JendR
          mask_outflow(0,j) = cff
          mask_outflow(1,j) = cff
        END DO
      END IF
#  else
      Need to set outflow mask in set_masks.F for your application.
#  endif
# endif

!
!  Insure that masks at mass point source locations are set to water
!  to avoid masking with _FillValue at those locations.
!
      IF (LuvSrc(ng)) THEN
        DO is=1,Nsrc(ng)
          i=SOURCES(ng)%Isrc(is)
          j=SOURCES(ng)%Jsrc(is)
          IF (((IstrT.le.i).and.(i.le.IendT)).and.                        &
     &        ((JstrT.le.j).and.(j.le.JendT))) THEN
            IF (INT(SOURCES(ng)%Dsrc(is)).eq.0) THEN
              umask_full(i,j)=1.0_r8
            ELSE
              vmask_full(i,j)=1.0_r8
            END IF
          END IF
        END DO
      END IF
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmask_full)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_full)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          umask_full)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vmask_full)
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pmask_full, rmask_full, umask_full, vmask_full)
# endif

# if defined AVERAGES    || \
    (defined AD_AVERAGES && defined ADJOINT) || \
    (defined RP_AVERAGES && defined TL_IOMS) || \
    (defined TL_AVERAGES && defined TANGENT)
!
!-----------------------------------------------------------------------
!  Initialize average file Land/Sea masks for time-averaged fields.
!-----------------------------------------------------------------------
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
#  ifdef WET_DRY
          pmask_avg(i,j)=0.0_r8
#  else
          pmask_avg(i,j)=pmask_full(i,j)
#  endif
        END DO
      END DO

      DO j=JstrT,JendT
        DO i=IstrT,IendT
#  ifdef WET_DRY
          rmask_avg(i,j)=0.0_r8
#  else
          rmask_avg(i,j)=rmask_full(i,j)
#  endif
        END DO
      END DO

      DO j=JstrT,JendT
        DO i=IstrP,IendT
#  ifdef WET_DRY
          umask_avg(i,j)=0.0_r8
#  else
          umask_avg(i,j)=umask_full(i,j)
#  endif
        END DO
      END DO

      DO j=JstrP,JendT
        DO i=IstrT,IendT
#  ifdef WET_DRY
          vmask_avg(i,j)=0.0_r8
#  else
          vmask_avg(i,j)=vmask_full(i,j)
#  endif
        END DO
      END DO
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmask_avg)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_avg)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          umask_avg)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vmask_avg)
      END IF

#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pmask_avg, rmask_avg, umask_avg, vmask_avg)
#  endif
# endif

# if defined AVERAGES2
!
!-----------------------------------------------------------------------
!  Initialize average2 file Land/Sea masks for time-averaged fields.
!-----------------------------------------------------------------------
!
      DO j=Jstr,JendR
        DO i=Istr,IendR
#  ifdef WET_DRY
          pmask_avg2(i,j)=0.0_r8
#  else
          pmask_avg2(i,j)=pmask_full(i,j)
#  endif
        END DO
      END DO

      DO j=JstrR,JendR
        DO i=IstrR,IendR
#  ifdef WET_DRY
          rmask_avg2(i,j)=0.0_r8
#  else
          rmask_avg2(i,j)=rmask_full(i,j)
#  endif
        END DO
      END DO

      DO j=JstrR,JendR
        DO i=Istr,IendR
#  ifdef WET_DRY
          umask_avg2(i,j)=0.0_r8
#  else
          umask_avg2(i,j)=umask_full(i,j)
#  endif
        END DO
      END DO

      DO j=Jstr,JendR
        DO i=IstrR,IendR
#  ifdef WET_DRY
          vmask_avg2(i,j)=0.0_r8
#  else
          vmask_avg2(i,j)=vmask_full(i,j)
#  endif
        END DO
      END DO
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmask_avg2)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_avg2)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          umask_avg2)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vmask_avg2)
      END IF

#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pmask_avg2, rmask_avg2,                       &
     &                    umask_avg2, vmask_avg2)
#  endif
# endif

# ifdef DIAGNOSTICS
!
!-----------------------------------------------------------------------
!  Initialize diagnostic file Land/Sea masks for time-averaged fields.
!-----------------------------------------------------------------------
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
#  ifdef WET_DRY
          pmask_dia(i,j)=0.0_r8
#  else
          pmask_dia(i,j)=pmask_full(i,j)
#  endif
        END DO
      END DO

      DO j=JstrT,JendT
        DO i=IstrT,IendT
#  ifdef WET_DRY
          rmask_dia(i,j)=0.0_r8
#  else
          rmask_dia(i,j)=rmask_full(i,j)
#  endif
        END DO
      END DO

      DO j=JstrT,JendT
        DO i=IstrP,IendT
#  ifdef WET_DRY
          umask_dia(i,j)=0.0_r8
#  else
          umask_dia(i,j)=umask_full(i,j)
#  endif
        END DO
      END DO

      DO j=JstrP,JendT
        DO i=IstrT,IendT
#  ifdef WET_DRY
          vmask_dia(i,j)=0.0_r8
#  else
          vmask_dia(i,j)=vmask_full(i,j)
#  endif
        END DO
      END DO
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmask_dia)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_dia)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          umask_dia)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vmask_dia)
      END IF

#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pmask_dia, rmask_dia, umask_dia, vmask_dia)
#  endif

# endif

      RETURN
      END SUBROUTINE set_masks_tile

# if defined WET_DRY     && \
     defined AVERAGES    || \
    (defined AD_AVERAGES && defined ADJOINT) || \
    (defined RP_AVERAGES && defined TL_IOMS) || \
    (defined TL_AVERAGES && defined TANGENT)
!
!***********************************************************************
      SUBROUTINE set_avg_masks (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          pmask_avg, rmask_avg,                   &
     &                          umask_avg, vmask_avg)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: pmask_avg(LBi:,LBj:)
      real(r8), intent(inout) :: rmask_avg(LBi:,LBj:)
      real(r8), intent(inout) :: umask_avg(LBi:,LBj:)
      real(r8), intent(inout) :: vmask_avg(LBi:,LBj:)
#  else
      real(r8), intent(inout) :: pmask_avg(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rmask_avg(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: umask_avg(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vmask_avg(LBi:UBi,LBj:UBj)
#  endif
!
!
!  Local variable declarations.
!
      integer :: i, j

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Return if time-averaging window is zero.
!-----------------------------------------------------------------------
!
      IF (nAVG(ng).eq.0) RETURN
!
!-----------------------------------------------------------------------
!  If last time-step of average window, convert time dependent counters
!  for wet points to time-averaged Land/Sea masks (dry=0, wet=1) for
!  the current average window period. Notice that a grid point is wet
!  if the count is greater than zero for the current time average
!  window.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).gt.ntsAVG(ng)).and.                                  &
     &    (MOD(iic(ng)-1,nAVG(ng)).eq.0).and.                           &
     &    ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))) THEN

        DO j=JstrP,JendP
          DO i=IstrP,IendP
            pmask_avg(i,j)=MIN(1.0_r8, pmask_avg(i,j))
          END DO
        END DO
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            rmask_avg(i,j)=MIN(1.0_r8, rmask_avg(i,j))
          END DO
        END DO
        DO j=JstrT,JendT
          DO i=IstrP,IendT
            umask_avg(i,j)=MIN(1.0_r8, umask_avg(i,j))
          END DO
        END DO
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            vmask_avg(i,j)=MIN(1.0_r8, vmask_avg(i,j))
          END DO
        END DO
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_p2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            pmask_avg)
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            rmask_avg)
          CALL exchange_u2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            umask_avg)
          CALL exchange_v2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            vmask_avg)
        END IF

#  ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, model, 4,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      pmask_avg, rmask_avg, umask_avg, vmask_avg)
#  endif

      END IF

      RETURN

      END SUBROUTINE set_avg_masks
# endif
#endif
      END MODULE set_masks_mod
