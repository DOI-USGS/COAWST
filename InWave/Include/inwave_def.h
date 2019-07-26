/*
** svn $Id: inwave_def.h 429 2009-12-20 17:30:26Z jcwarner $
*************************************************** John C. Warner   ***
** Copyright (c) 2002-2016 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Defines inwave model input parameters in output NetCDF files.     **
**  It is included in routine "def_info.F".                           **
**                                                                    **
************************************************************************
*/

!
!  Define InWave model parameters.
!
      Vinfo( 1)='energy_angle'
      Vinfo( 2)='direction respect to the north of the bin'
      Vinfo( 3)='degrees'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/inwdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

!      Vinfo( 1)='energy_angle_c'
!      Vinfo( 2)='direction respect to the north of the bin'
!      Vinfo( 3)='degrees'
!      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
!     &               1, (/inwdim/), Aval, Vinfo, ncname,                &
!     &               SetParAccess = .FALSE.)
!      IF (exit_flag.ne.NoError) RETURN
