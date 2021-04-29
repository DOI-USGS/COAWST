/*
** svn $Id: sediment_wrt.h 1054 2021-03-06 19:47:12Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2021 The ROMS/TOMS Group                        **
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
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'newlayer_thick',        &
     &                      newlayer_thick(ng), (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#ifdef BEDLOAD
      CALL netcdf_put_fvar (ng, model, ncname, 'bedload_coeff',         &
     &                      bedload_coeff(ng), (/0/), (/0/),            &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!# ifdef BEDLOAD_VANDERA
      CALL netcdf_put_fvar (ng, model, ncname, 'sg_zwbl',               &
     &                      sg_zwbl(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'sedslope_crit_wet',     &
     &                      sedslope_crit_wet(ng), (/0/), (/0/),        &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'sedslope_crit_dry',     &
     &                      sedslope_crit_dry(ng), (/0/), (/0/),        &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'slopefac_wet',          &
     &                      slopefac_wet(ng), (/0/), (/0/),             &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'slopefac_dry',          &
     &                      slopefac_dry(ng), (/0/), (/0/),             &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'bedload_vandera_alphaw',&
     &                      bedload_vandera_alphaw(ng), (/0/), (/0/),   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'bedload_vandera_alphac',&
     &                      bedload_vandera_alphac(ng), (/0/), (/0/),   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!# endif
#endif

!#ifdef ANA_SEDIMENT
      CALL netcdf_put_fvar (ng, model, ncname, 'Sd50',                  &
     &                      Sd50(:,ng), (/1/), (/NST/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Srho',                  &
     &                      Srho(:,ng), (/1/), (/NST/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Csed',                  &
     &                      Csed(:,ng), (/1/), (/NST/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!#endif

      CALL netcdf_put_fvar (ng, model, ncname, 'Wsed',                  &
     &                      Wsed(:,ng), (/1/), (/NST/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Erate',                 &
     &                      Erate(:,ng), (/1/), (/NST/),                &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'tau_ce',                &
     &                      tau_ce(:,ng), (/1/), (/NST/),               &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'tau_cd',                &
     &                      tau_cd(:,ng), (/1/), (/NST/),               &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'poros',                 &
     &                      poros(:,ng), (/1/), (/NST/),                &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
