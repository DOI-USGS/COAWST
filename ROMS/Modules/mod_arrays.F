#include "cppdefs.h"
      SUBROUTINE mod_arrays (allocate_vars)
!
!svn $Id: mod_arrays.F 1589 2008-07-24 04:02:51Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
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
!
#ifdef AVERAGES
      USE mod_average, ONLY : allocate_average, initialize_average
#endif
#ifdef OBC
      USE mod_boundary, ONLY : allocate_boundary, initialize_boundary
#endif
#if defined AD_SENSITIVITY      || defined CLIMATOLOGY      || \
    defined IS4DVAR_SENSITIVITY || defined OPT_OBSERVATIONS || \
    defined SENSITIVITY_4DVAR   || defined SO_SEMI
      USE mod_clima, ONLY : allocate_clima, initialize_clima
#endif
#ifdef SOLVE3D
      USE mod_coupling, ONLY : allocate_coupling, initialize_coupling
#endif
#ifdef DIAGNOSTICS
      USE mod_diags, ONLY : allocate_diags, initialize_diags
#endif
      USE mod_forces, ONLY : allocate_forces, initialize_forces
      USE mod_grid, ONLY : allocate_grid, initialize_grid
      USE mod_mixing, ONLY : allocate_mixing, initialize_mixing
#if defined ASSIMILATION || defined NUDGING
      USE mod_obs, ONLY : allocate_obs, initialize_obs
#endif
      USE mod_ocean, ONLY : allocate_ocean, initialize_ocean
#if defined SEDIMENT || defined BBL_MODEL
      USE mod_sedbed, ONLY : allocate_sedbed, initialize_sedbed
#endif
#if defined UV_PSOURCE || defined TS_PSOURCE || defined Q_PSOURCE
      USE mod_sources, ONLY : allocate_sources
#endif
#if defined SSH_TIDES || defined UV_TIDES
      USE mod_tides, ONLY : allocate_tides, initialize_tides
#endif
#ifdef BBL_MODEL
      USE mod_bbl, ONLY : allocate_bbl, initialize_bbl
#endif
#ifdef REFINED_GRID
      USE mod_refined, ONLY : allocate_refined, initialize_refined
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
      integer :: ng
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
      integer :: tile, subs, thread

      integer, parameter :: model = 0
!
!-----------------------------------------------------------------------
!  Allocate model structures.
!-----------------------------------------------------------------------
!
      IF (allocate_vars) then
        tile=0
        DO ng=1,Ngrids
          LBi=BOUNDS(ng)%LBi(TILE)
          UBi=BOUNDS(ng)%UBi(TILE)
          LBj=BOUNDS(ng)%LBj(TILE)
          UBj=BOUNDS(ng)%UBj(TILE)
          LBij=BOUNDS(ng)%LBij
          UBij=BOUNDS(ng)%UBij
#ifdef AVERAGES
          CALL allocate_average (ng, LBi, UBi, LBj, UBj)
#endif
#ifdef OBC
          CALL allocate_boundary (ng)
#endif
#ifdef BBL_MODEL
          CALL allocate_bbl (ng, LBi, UBi, LBj, UBj)
#endif
#if defined AD_SENSITIVITY      || defined CLIMATOLOGY      || \
    defined IS4DVAR_SENSITIVITY || defined OPT_OBSERVATIONS || \
    defined SENSITIVITY_4DVAR   || defined SO_SEMI
          CALL allocate_clima (ng, LBi, UBi, LBj, UBj)
#endif
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
#if defined SSH_TIDES || defined UV_TIDES
          CALL allocate_tides (ng, LBi, UBi, LBj, UBj)
#endif
#if defined UV_PSOURCE || defined TS_PSOURCE || defined Q_PSOURCE
          CALL allocate_sources (ng)
#endif
#ifdef REFINED_GRID
          CALL allocate_refined (ng, LBi, UBi, LBj, UBj)
#endif
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Allocate and intialize variables within structures for each grid.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile) SHARED(numthreads)
      DO thread=0,numthreads-1
        DO ng=1,Ngrids
#if defined _OPENMP || defined DISTRIBUTE
          subs=NtileX(ng)*NtileE(ng)/numthreads
#else
          subs=1
#endif
          DO tile=subs*thread,subs*(thread+1)-1
#ifdef AVERAGES
            CALL initialize_average (ng, TILE)
#endif
#ifdef BBL_MODEL
            CALL initialize_bbl (ng, TILE)
#endif
#ifdef OBC
            CALL initialize_boundary (ng, TILE, model)
#endif
#if defined AD_SENSITIVITY      || defined CLIMATOLOGY      || \
    defined IS4DVAR_SENSITIVITY || defined OPT_OBSERVATIONS || \
    defined SENSITIVITY_4DVAR   || defined SO_SEMI
            CALL initialize_clima (ng, TILE)
#endif
#ifdef SOLVE3D
            CALL initialize_coupling (ng, TILE, model)
#endif
#ifdef DIAGNOSTICS
            CALL initialize_diags (ng, TILE)
#endif
            CALL initialize_forces (ng, TILE, model)
            CALL initialize_grid (ng, TILE, model)
            CALL initialize_mixing (ng, TILE, model)
#if defined ASSIMILATION || defined NUDGING
            CALL initialize_obs (ng, TILE)
#endif
            CALL initialize_ocean (ng, TILE, model)
#if defined SEDIMENT || defined BBL_MODEL
            CALL initialize_sedbed (ng, TILE, model)
#endif
#if defined SSH_TIDES || defined UV_TIDES
            CALL initialize_tides (ng, TILE)
#endif
#ifdef REFINED_GRID
            CALL initialize_refined (ng, TILE, model)
#endif
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE mod_arrays
