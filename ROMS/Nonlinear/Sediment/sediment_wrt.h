/*
** svn $Id: sediment_wrt.h 429 2009-12-20 17:30:26Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2016 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Writes sediment model input parameters into output NetCDF files.  **
**  It is included in routine "wrt_info.F".                           **
**                                                                    **
************************************************************************
*/

!
!  Write out sediment model parameters.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'minlayer_thick',        &
     &                      minlayer_thick(ng), (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'newlayer_thick',        &
     &                      newlayer_thick(ng), (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

#ifdef BEDLOAD
      CALL netcdf_put_fvar (ng, model, ncname, 'bedload_coeff',         &
     &                      bedload_coeff(ng), (/0/), (/0/),            &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#endif

!#ifdef ANA_SEDIMENT
      CALL netcdf_put_fvar (ng, model, ncname, 'Sd50',                  &
     &                      Sd50(:,ng), (/1/), (/NST/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Srho',                  &
     &                      Srho(:,ng), (/1/), (/NST/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Csed',                  &
     &                      Csed(:,ng), (/1/), (/NST/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
!#endif

      CALL netcdf_put_fvar (ng, model, ncname, 'Wsed',                  &
     &                      Wsed(:,ng), (/1/), (/NST/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Erate',                 &
     &                      Erate(:,ng), (/1/), (/NST/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'tau_ce',                &
     &                      tau_ce(:,ng), (/1/), (/NST/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'tau_cd',                &
     &                      tau_cd(:,ng), (/1/), (/NST/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'poros',                 &
     &                      poros(:,ng), (/1/), (/NST/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
