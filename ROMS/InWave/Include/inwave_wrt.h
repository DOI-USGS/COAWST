/*
** svn $Id: inwave_wrt.h 429 2009-12-20 17:30:26Z arango $
*************************************************** John C. Warner   ***
** Copyright (c) 2002-2016 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Writes inwave model input parameters into output NetCDF files.    **
**  It is included in routine "wrt_info.F".                           **
**                                                                    **
************************************************************************
*/

!
!  Write out InWave model parameters.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'energy_angle',          &
     &                      WAVEG(ng)%WD*180.0d0/pi, (/1/), (/ND/),     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

!      CALL netcdf_put_fvar (ng, model, ncname, 'energy_angle_c',        &
!     &                      WAVEB(ng)%WD_bnd, (/1/),                    &
!     &                      (/WAVEB(ng)%ND_bnd/),                       &
!     &                      ncid = ncid)
!      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

