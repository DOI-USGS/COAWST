/*
** git $Id$
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2023 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
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
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'BioIter',           &
     &                          BioIter(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'ResRate',           &
     &                          ResRate(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
