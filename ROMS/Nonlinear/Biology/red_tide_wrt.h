/*
** svn $Id: red_tide_wrt.h 859 2017-08-02 01:45:30Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2017 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Writes red tide (Stock et al., 2005; He et al., 2008) biological  **
**  model input parameters into output NetCDF files.  It is included  **
**  in routine "wrt_info.F".                                          **
**                                                                    **
************************************************************************
*/

!
!  Write out red tide (Stock et al., 2005; He et al., 2008) biological
!  model parameters.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'BioIter',               &
     &                      BioIter(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Gmax',                  &
     &                      Gmax(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Dg',                    &
     &                      Dg(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Kn',                    &
     &                      Kn(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'G_eff',                 &
     &                      G_eff(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'G_r',                   &
     &                      G_r(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'AttW',                  &
     &                      AttW(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'AttS',                  &
     &                      AttS(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'E_light',               &
     &                      E_light(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'E_dark',                &
     &                      E_dark(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Tmin_growth',           &
     &                      Tmin_growth(ng), (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'DIN_Cdepth',            &
     &                      DIN_Cdepth(ng), (/0/), (/0/),               &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'wDino',                 &
     &                      wDino(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Mor_a',                 &
     &                      Mor_a(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Mor_b',                 &
     &                      Mor_b(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Mor_Q10',               &
     &                      Mor_Q10(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Mor_T0',                &
     &                      Mor_T0(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
