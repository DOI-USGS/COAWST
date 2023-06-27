#include "cppdefs.h"
      MODULE mod_arrays
!
!git $Id$
!svn $Id: mod_arrays.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module is used to allocate, initialize, and deallocate ROMS    !
!  state arrays for each nested and/or multiple connected grids.       !
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
      USE mod_average,  ONLY : allocate_average,                        &
     &                         deallocate_average,                      &
     &                         initialize_average
#endif
      USE mod_boundary, ONLY : allocate_boundary,                       &
     &                         deallocate_boundary,                     &
     &                         initialize_boundary
      USE mod_clima,    ONLY : allocate_clima,                          &
     &                         deallocate_clima,                        &
     &                         initialize_clima
#ifdef SOLVE3D
      USE mod_coupling, ONLY : allocate_coupling,                       &
     &                         deallocate_coupling,                     &
     &                         initialize_coupling
#endif
#ifdef DIAGNOSTICS
      USE mod_diags,    ONLY : allocate_diags,                          &
     &                         deallocate_diags,                        &
     &                         initialize_diags
#endif
      USE mod_forces,   ONLY : allocate_forces,                         &
     &                         deallocate_forces,                       &
     &                         initialize_forces
#if defined FOUR_DVAR || defined VERIFICATION
      USE mod_fourdvar, ONLY : deallocate_fourdvar,                     &
     &                         initialize_fourdvar
#endif
      USE mod_grid,     ONLY : allocate_grid,                           &
     &                         deallocate_grid,                         &
     &                         initialize_grid
      USE mod_iounits,  ONLY : deallocate_iounits
      USE mod_mixing,   ONLY : allocate_mixing,                         &
     &                         deallocate_mixing,                       &
     &                         initialize_mixing
#ifdef NESTING
      USE mod_nesting,  ONLY : allocate_nesting,                        &
     &                         deallocate_nesting,                      &
     &                         initialize_nesting
#endif
      USE mod_ocean,    ONLY : allocate_ocean,                          &
     &                         deallocate_ocean,                        &
     &                         initialize_ocean
#if defined SEDIMENT || defined BBL_MODEL
      USE mod_sedbed,   ONLY : allocate_sedbed,                         &
     &                         deallocate_sedbed,                       &
     &                         initialize_sedbed
#endif
#if defined SEDIMENT && defined SED_FLOCS
      USE mod_sedflocs, ONLY : allocate_sedflocs,                       &
     &                         initialize_sedflocs
!    &                         deallocate_sedflocs,                     &
#endif
#if defined VEGETATION 
      USE mod_vegarr,   ONLY : allocate_vegarr,                         &
     &                         initialize_vegarr
!    &                         deallocate_vegarr,                       &
#endif 
      USE mod_sources,  ONLY : allocate_sources,                        &
     &                         deallocate_sources
#if defined SSH_TIDES || defined UV_TIDES
      USE mod_tides,    ONLY : allocate_tides,                          &
     &                         deallocate_tides,                        &
     &                         initialize_tides
#endif
#ifdef BBL_MODEL
      USE mod_bbl,      ONLY : allocate_bbl,                            &
     &                         deallocate_bbl,                          &
     &                         initialize_bbl
#endif
#if defined PIO_LIB && defined DISTRIBUTE
      USE set_pio_mod,  ONLY : set_iodecomp
#endif
!
      implicit none
!
      PUBLIC :: ROMS_allocate_arrays
      PUBLIC :: ROMS_deallocate_arrays
      PUBLIC :: ROMS_initialize_arrays
!
#if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR      || \
    defined SO_SEMI
      logical :: LallocateClima = .TRUE.
#else
      logical :: LallocateClima = .FALSE.
#endif
!
      CONTAINS
!
      SUBROUTINE ROMS_allocate_arrays (allocate_vars)
!
!=======================================================================
!                                                                      !
!  This routine allocates ROMS state variables.                        !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations
!
      logical, intent(in) :: allocate_vars
!
!  Local variable declarations.
!
      integer :: ng, thread, tile
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_allocate_arrays"

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on allocation time wall clock.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_on (ng, iNLM, 1, __LINE__, MyFile)
        END DO
!$OMP BARRIER
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Allocate model type-derived structures and its variables.
!-----------------------------------------------------------------------
!
      IF (allocate_vars) then
#ifdef DISTRIBUTE
        tile=MyRank
#else
        tile=0
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
          CALL allocate_boundary (ng)
#ifdef BBL_MODEL
          CALL allocate_bbl (ng, LBi, UBi, LBj, UBj)
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
#if defined SSH_TIDES || defined UV_TIDES
          CALL allocate_tides (ng, LBi, UBi, LBj, UBj)
#endif
          IF (LuvSrc(ng).or.LwSrc(ng).or.ANY(LtracerSrc(:,ng))) THEN
            CALL allocate_sources (ng)
          END IF
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

#if defined PIO_LIB && defined DISTRIBUTE
!
!  Allocate and initialize PIO single and double precision data
!  decomposition for mapping between computational and I/O processes.
!  It needs to done after all parameters have been allocated and
!  initialize. In particular, we need tidal constituents and state
!  propagator  parameters.
!
        CALL set_iodecomp
#endif

        LallocatedMemory=.TRUE.

      END IF

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off allocation time wall clock.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iNLM, 1, __LINE__, MyFile)
        END DO
!$OMP BARRIER
      END DO
#endif
!
      RETURN
      END SUBROUTINE ROMS_allocate_arrays
!
      SUBROUTINE ROMS_deallocate_arrays
!
!=======================================================================
!                                                                      !
!  This routine deallocates ROMS state objects, arrays and vectors.    !
!                                                                      !
!=======================================================================
!
      USE mod_param, ONLY : Ngrids
!
!  Local variable declarations.
!
      integer :: ng
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_deallocate_arrays"
!
!-----------------------------------------------------------------------
!  Deallocate all structures.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!$OMP MASTER
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
        CALL deallocate_average (ng)
#endif
        CALL deallocate_boundary (ng)
#ifdef BBL_MODEL
        CALL deallocate_bbl (ng)
#endif
        IF (LallocateClima.or.Lclimatology(ng)) THEN
          CALL deallocate_clima (ng)
        END IF
#ifdef SOLVE3D
        CALL deallocate_coupling (ng)
#endif
#ifdef DIAGNOSTICS
        CALL deallocate_diags (ng)
#endif
        CALL deallocate_forces (ng)
        CALL deallocate_grid (ng)
        CALL deallocate_mixing (ng)
        CALL deallocate_ocean (ng)
#if defined SEDIMENT || defined BBL_MODEL
        CALL deallocate_sedbed (ng)
#endif
#if defined SEDIMENT && defined SED_FLOCS
!       CALL deallocate_sedflocs (ng)
#endif
!#if defined VEGETATION
!       CALL deallocate_vegarr(ng)
!#endif
#if defined SSH_TIDES || defined UV_TIDES
        CALL deallocate_tides (ng)
#endif
        IF (LuvSrc(ng).or.LwSrc(ng).or.ANY(LtracerSrc(:,ng))) THEN
          CALL deallocate_sources (ng)
        END IF
!$OMP END MASTER
!$OMP BARRIER
      END DO

#ifdef NESTING
!
!  Deallocate nesting variables and structures.
!
!$OMP MASTER
      CALL deallocate_nesting
!$OMP END MASTER
!$OMP BARRIER
#endif

#if defined FOUR_DVAR || defined VERIFICATION
!
!  Deallocate observations arrays and object.
!
      CALL deallocate_fourdvar
#endif
!
!  Deallocate I/O derived-type structures.
!
      CALL deallocate_iounits
!
!  Deallocate main configuration dimensions and associated parameters.
!
      CALL deallocate_param
!
      RETURN
      END SUBROUTINE ROMS_deallocate_arrays
!
      SUBROUTINE ROMS_initialize_arrays
!
!=======================================================================
!                                                                      !
!  This routine initialize ROMS state variables. In shared-memory it   !
!  important for the first-touch policy in memory.                     !
!                                                                      !
!=======================================================================
!
!  Local variable declarations.
!
      integer :: ng, thread, tile
!
      integer, parameter :: model = 0
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_initialize_arrays"

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on allocation time wall clock.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_on (ng, iNLM, 1, __LINE__, MyFile)
        END DO
!$OMP BARRIER
      END DO
#endif
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
#ifdef BBL_MODEL
          CALL initialize_bbl (ng, tile)
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
          CALL initialize_ocean (ng, tile, model)
#if defined SEDIMENT || defined BBL_MODEL
          CALL initialize_sedbed (ng, tile, model)
#endif
#if defined SEDIMENT && defined SED_FLOCS
            CALL initialize_sedflocs (ng, tile, model)
#endif
#if defined VEGETATION 
            CALL initialize_vegarr (ng, tile, model)
#endif
#if defined SSH_TIDES || defined UV_TIDES
          CALL initialize_tides (ng, tile)
#endif
        END DO
!$OMP BARRIER
      END DO

#if defined FOUR_DVAR || defined VERIFICATION
!
!  Finish allocating observation arrays.  Then, initialize observation
!  arrays. Notice that some of the module variables are allocated in
!  the call to "allocate_fourdvar" in "read_phypar".
!
      CALL initialize_fourdvar
#endif

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off allocation time wall clock.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iNLM, 1, __LINE__, MyFile)
        END DO
!$OMP BARRIER
      END DO
#endif
!
      RETURN
      END SUBROUTINE ROMS_initialize_arrays
!
      END MODULE mod_arrays
