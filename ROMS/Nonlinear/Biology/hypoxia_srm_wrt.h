/*
** git $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2026 The ROMS Group                             **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.md                                              **
************************************************************************
**                                                                    **
**  Writes Hypoxia Simple Respiration Model input parameters into     **
**  output NetCDF files. It is included in routine "wrt_info.F".      **
**                                                                    **
************************************************************************
*/

!
!  Write out Hypoxia Simple Respiration Model parameters.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'BioIter',               &
     &                      BioIter(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ResRate',               &
     &                      ResRate(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
