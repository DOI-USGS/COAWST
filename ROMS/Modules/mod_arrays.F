#include "cppdefs.h"
      SUBROUTINE mod_arrays (allocate_vars)
!
!svn $Id: mod_arrays.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine routine allocates and initializa model state arrays    !
!  for each nested and/or multiple connected grids.                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
!
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
      USE mod_average, ONLY : allocate_average, initialize_average
#endif
#ifdef AVERAGES2
      USE mod_average2, ONLY : allocate_average2, initialize_average2
#endif
#ifdef FILTERED
      USE mod_filter, ONLY : allocate_filter, initialize_filter
#endif
      USE mod_boundary, ONLY : allocate_boundary, initialize_boundary
      USE mod_clima, ONLY : allocate_clima, initialize_clima
#ifdef SOLVE3D
      USE mod_coupling, ONLY : allocate_coupling, initialize_coupling
#endif
#ifdef DIAGNOSTICS
      USE mod_diags, ONLY : allocate_diags, initialize_diags
#endif
      USE mod_forces, ONLY : allocate_forces, initialize_forces
      USE mod_grid, ONLY : allocate_grid, initialize_grid
      USE mod_mixing, ONLY : allocate_mixing, initialize_mixing
#ifdef NESTING
      USE mod_nesting, ONLY : allocate_nesting, initialize_nesting
#endif
#if defined ASSIMILATION || defined NUDGING
      USE mod_obs, ONLY : allocate_obs, initialize_obs
#endif
      USE mod_ocean, ONLY : allocate_ocean, initialize_ocean
#if defined SEDIMENT || defined BBL_MODEL
      USE mod_sedbed, ONLY : allocate_sedbed, initialize_sedbed
#endif
#if defined SEDIMENT && defined SED_FLOCS
      USE mod_sedflocs, ONLY : allocate_sedflocs, initialize_sedflocs
#endif
#if defined VEGETATION 
      USE mod_vegarr,ONLY:allocate_vegarr, initialize_vegarr
#endif 
      USE mod_sources, ONLY : allocate_sources
#if defined TRC_PSOURCE
      USE mod_trc_sources, ONLY : allocate_trc_sources
#endif
#if defined SSH_TIDES || defined UV_TIDES || defined POT_TIDES
      USE mod_tides, ONLY : allocate_tides, initialize_tides
#endif
#ifdef BBL_MODEL
      USE mod_bbl, ONLY : allocate_bbl, initialize_bbl
#endif
#if defined ICE_MODEL || defined CICE_MODEL
      USE mod_ice, ONLY : allocate_ice, initialize_ice
#endif
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: allocate_vars
!
!  Local variable declarations.
!
      logical :: LallocateClima

      integer :: ng, thread, tile
      integer :: LBi, UBi, LBj, UBj, LBij, UBij

      integer, parameter :: model = 0

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on allocation time wall clock.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO thread=THREAD_RANGE
            CALL wclock_on (ng, iNLM, 1)
          END DO
!$OMP BARRIER
        END DO
#endif
!
!-----------------------------------------------------------------------
!  Allocate model structures.
!-----------------------------------------------------------------------
!
      IF (allocate_vars) then
#ifdef DISTRIBUTE
        tile=MyRank
#else
        tile=0
#endif

#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI
        LallocateClima=.TRUE.
#else
        LallocateClima=.FALSE.
#endif

        DO ng=1,Ngrids
!$OMP MASTER
          LBi=BOUNDS(ng)%LBi(tile)
          UBi=BOUNDS(ng)%UBi(tile)
          LBj=BOUNDS(ng)%LBj(tile)
          UBj=BOUNDS(ng)%UBj(tile)
          LBij=BOUNDS(ng)%LBij
          UBij=BOUNDS(ng)%UBij
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
          CALL allocate_average (ng, LBi, UBi, LBj, UBj)
#endif
#ifdef AVERAGES2
          CALL allocate_average2 (ng, LBi, UBi, LBj, UBj)
#endif
#ifdef FILTERED
          CALL allocate_filter (ng, LBi, UBi, LBj, UBj)
#endif
          CALL allocate_boundary (ng)
#ifdef BBL_MODEL
          CALL allocate_bbl (ng, LBi, UBi, LBj, UBj)
#endif
#if defined ICE_MODEL || defined CICE_MODEL
          CALL allocate_ice (ng, LBi, UBi, LBj, UBj)
#endif
          IF (LallocateClima.or.Lclimatology(ng)) THEN
            CALL allocate_clima (ng, LBi, UBi, LBj, UBj)
          END IF
#ifdef SOLVE3D
          CALL allocate_coupling (ng, LBi, UBi, LBj, UBj)
#endif
#ifdef DIAGNOSTICS
          CALL allocate_diags (ng, LBi, UBi, LBj, UBj)
#endif
          CALL allocate_forces (ng, LBi, UBi, LBj, UBj)
          CALL allocate_grid (ng, LBi, UBi, LBj, UBj, LBij, UBij)
          CALL allocate_mixing (ng, LBi, UBi, LBj, UBj)
#if defined ASSIMILATION || defined NUDGING
          CALL allocate_obs (ng, LBi, UBi, LBj, UBj)
#endif
          CALL allocate_ocean (ng, LBi, UBi, LBj, UBj)
#if defined SEDIMENT || defined BBL_MODEL
          CALL allocate_sedbed (ng, LBi, UBi, LBj, UBj)
#endif
#if defined SEDIMENT && defined SED_FLOCS
          CALL allocate_sedflocs(ng, LBi, UBi, LBj, UBj)
#endif
#if defined VEGETATION
          CALL allocate_vegarr(ng, LBi, UBi, LBj, UBj)
#endif
#if defined SSH_TIDES || defined UV_TIDES || defined POT_TIDES
          CALL allocate_tides (ng, LBi, UBi, LBj, UBj)
#endif
          IF (LuvSrc(ng).or.LwSrc(ng).or.ANY(LtracerSrc(:,ng))) THEN
            CALL allocate_sources (ng)
          END IF
#if defined TRC_PSOURCE
          CALL allocate_trc_sources (ng)
#endif
!$OMP END MASTER
!$OMP BARRIER
        END DO
#ifdef NESTING
!
!  Allocate and initialized contact points boundaty structure. It
!  needs to be delayed to the end because we need "LBC_apply" to
!  allocated in "mod_boundary" for all nested grid.
!
!$OMP MASTER
        CALL allocate_nesting
!$OMP END MASTER
!$OMP BARRIER
#endif
      END IF
!
!-----------------------------------------------------------------------
!  Intialize variables within structures for each grid.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
#ifdef NESTING
        IF (ng.eq.1) THEN
          CALL initialize_nesting
        END IF
#endif
        DO tile=first_tile(ng),last_tile(ng),+1
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
          CALL initialize_average (ng, tile)
#endif
#ifdef AVERAGES2
          CALL initialize_average2 (ng, tile)
#endif
#ifdef FILTERED
          CALL initialize_filter (ng, tile)
#endif
#ifdef BBL_MODEL
          CALL initialize_bbl (ng, tile)
#endif
#if defined ICE_MODEL || defined CICE_MODEL
          CALL initialize_ice (ng, tile)
#endif
          CALL initialize_boundary (ng, tile, model)

          IF (LallocateClima.or.Lclimatology(ng)) THEN
            CALL initialize_clima (ng, tile)
          END IF
#ifdef SOLVE3D
          CALL initialize_coupling (ng, tile, model)
#endif
#ifdef DIAGNOSTICS
          CALL initialize_diags (ng, tile)
#endif
          CALL initialize_forces (ng, tile, model)
          CALL initialize_grid (ng, tile, model)
          CALL initialize_mixing (ng, tile, model)
#if defined ASSIMILATION || defined NUDGING
          CALL initialize_obs (ng, tile)
#endif
          CALL initialize_ocean (ng, tile, model)
#if defined SEDIMENT || defined BBL_MODEL
          CALL initialize_sedbed (ng, tile, model)
#endif
#if defined SEDIMENT && defined SED_FLOCS
            CALL initialize_sedflocs (ng, TILE, model)
#endif
#if defined VEGETATION 
            CALL initialize_vegarr (ng, TILE, model)
#endif
#if defined SSH_TIDES || defined UV_TIDES || defined POT_TIDES
          CALL initialize_tides (ng, tile)
#endif
        END DO
!$OMP BARRIER
      END DO

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off allocation time wall clock.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iNLM, 1)
        END DO
!$OMP BARRIER
      END DO
#endif

      RETURN
      END SUBROUTINE mod_arrays
