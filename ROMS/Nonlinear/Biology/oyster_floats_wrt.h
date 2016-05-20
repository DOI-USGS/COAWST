/*
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2016 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Writes oyster floats model input parameters in output NetCDF      **
**  files.  It is included in routine "wrt_info.F".                   **
**                                                                    **
************************************************************************
*/

!
!  Write out oyster float model (Dekshenieks et al. 1993, 1996, 1997;
!  Narvaez et al. 2012a,b).
!
      CALL netcdf_put_fvar (ng, model, ncname, 'Larvae_GR0',            &
     &                      Larvae_GR0(ng), (/0/), (/0/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Larvae_size0',          &
     &                      Larvae_size0(ng), (/0/), (/0/),             &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'food_supply',           &
     &                      food_supply(ng), (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'settle_size',           &
     &                      settle_size(ng), (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Turbidity effect parameters on planktonic oyster larvae growth.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'turb_ambi',             &
     &                      turb_ambi(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'turb_axis',             &
     &                      turb_axis(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'turb_base',             &
     &                      turb_base(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'turb_crit',             &
     &                      turb_crit(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'turb_mean',             &
     &                      turb_mean(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'turb_rate',             &
     &                      turb_rate(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'turb_size',             &
     &                      turb_size(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'turb_slop',             &
     &                      turb_slop(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Planktonic oyster larvae vertical migration (swimming) parameters.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'slope_Sdec',            &
     &                      slope_Sdec(ng), (/0/), (/0/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'slope_Sinc',            &
     &                      slope_Sinc(ng), (/0/), (/0/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'swim_Sdec',             &
     &                      swim_Sdec(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'swim_Sinc',             &
     &                      swim_Sinc(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'swim_Tmin',             &
     &                      swim_Tmin(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'swim_Tmax',             &
     &                      swim_Tmax(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Planktonic larvae sinking parameters.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'sink_base',             &
     &                      sink_base(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'sink_rate',             &
     &                      sink_rate(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'sink_size',             &
     &                      sink_size(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
