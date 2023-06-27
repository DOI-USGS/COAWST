#include "cppdefs.h"
      SUBROUTINE set_grid (ng, model)
!
!git $Id$
!svn $Id: set_grid.F 1151 2023-02-09 03:08:53Z arango $
!=======================================================================
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine sets application grid and associated variables and     !
!  parameters. It called only once during the initialization stage.    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#ifdef NESTING
      USE mod_nesting
#endif
      USE mod_scalars
!
      USE analytical_mod
#ifdef DISTRIBUTE
      USE distribute_mod,       ONLY : mp_bcasti
#endif
#ifdef TIDE_GENERATING_FORCES
      USE equilibrium_tide_mod, ONLY : harmonic_constituents
#endif
#ifndef ANA_GRID
      USE get_grid_mod,         ONLY : get_grid
#endif
#ifndef ANA_NUDGCOEF
      USE get_nudgcoef_mod,     ONLY : get_nudgcoef
#endif
      USE metrics_mod,          ONLY : metrics
#ifdef NESTING
      USE nesting_mod,          ONLY : nesting
#endif
      USE strings_mod,          ONLY : FoundError
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer :: tile
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Set horizontal grid, bathymetry, and Land/Sea masking (if any).
!  Use analytical functions or read in from a grid NetCDF.
!-----------------------------------------------------------------------
!
#ifdef ANA_GRID
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL ana_grid (ng, tile, model)
# ifdef MASKING
        CALL ana_mask (ng, tile, model)
# endif
# if defined AD_SENSITIVITY    || defined I4DVAR_ANA_SENSITIVITY || \
     defined OPT_OBSERVATIONS  || defined SO_SEMI                || \
     defined SENSITIVITY_4DVAR
        CALL ana_scope (ng, tile, model)
# endif
      END DO
!$OMP BARRIER
#else
!$OMP MASTER
# ifdef DISTRIBUTE
      CALL get_grid (ng, MyRank, model)
# else
      CALL get_grid (ng, -1, model)
# endif
!$OMP END MASTER
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, model, exit_flag)
# endif
!$OMP BARRIER
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Set vertical terrain-following coordinate transformation function.
!-----------------------------------------------------------------------
!
!$OMP MASTER
      CALL set_scoord (ng)
!$OMP END MASTER
!$OMP BARRIER
#endif

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Set barotropic time-steps average weighting function.
!-----------------------------------------------------------------------
!
!$OMP MASTER
      CALL set_weights (ng)
!$OMP END MASTER
!$OMP BARRIER
#endif

!
!-----------------------------------------------------------------------
!  Compute various metric term combinations.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL metrics (ng, tile, model)
      END DO
!$OMP BARRIER

#ifdef NESTING
!
!-----------------------------------------------------------------------
!  If nesting, initialize grid spacing (on_u and om_v) in REFINED(:)
!  structure.  They are used to impose mass flux at the finer grid
!  physical boundaries.
!-----------------------------------------------------------------------
!
      CALL nesting (ng, model, ndxdy)
#endif

#if defined WTYPE_GRID && defined ANA_WTYPE     && \
   (defined LMD_SKPP   || defined SOLAR_SOURCE)
!
!-----------------------------------------------------------------------
!  Set spatially varying Jerlov water type.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL ana_wtype (ng, tile, iNLM)
      END DO
!$OMP BARRIER
#endif

#ifndef CORRELATION
!
!-----------------------------------------------------------------------
!  If appropriate, set spatially varying nudging coefficients time
!  scales.
!-----------------------------------------------------------------------
!
# ifdef ANA_NUDGCOEF
      IF (Lnudging(ng)) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_nudgcoef (ng, tile, model)
        END DO
!$OMP BARRIER
      END IF
# else
      IF (Lnudging(ng)) THEN
!$OMP MASTER
#  ifdef DISTRIBUTE
        CALL get_nudgcoef (ng, MyRank, model)
#  else
        CALL get_nudgcoef (ng, -1, model)
#  endif
!$OMP END MASTER
#  ifdef DISTRIBUTE
        CALL mp_bcasti (ng, model, exit_flag)
#  endif
!$OMP BARRIER
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif

# ifdef TIDE_GENERATING_FORCES
!
!-----------------------------------------------------------------------
!  If applying tide generating forces, compute harmonic constituent
!  parameters of the equilibrium tide on Greenwich meridian for the
!  reference tide date number.
!-----------------------------------------------------------------------
!
      CALL harmonic_constituents (Lnodal)
# endif
#endif
!
      RETURN
      END SUBROUTINE set_grid
