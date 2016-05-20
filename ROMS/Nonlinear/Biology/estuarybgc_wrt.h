/*
** svn $Id: estuarybgc_wrt.h 429 2009-12-20 17:30:26Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2010 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Writes Fennel et al. (2006) ecosystem model input parameters into **
**  output NetCDF files. It is included in routine "wrt_info.F".      **
**                                                                    **
************************************************************************
*/

!
!  Write out Fennel et al. (2006) ecosystem model parameters.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'BioIter',               &
     &                      BioIter(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'AttSW',                 &
     &                      AttSW(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'AttChl',                &
     &                      AttChl(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'PARfrac',               &
     &                      PARfrac(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Vp0',                   &
     &                      Vp0(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'I_thNH4',               &
     &                      I_thNH4(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'D_p5NH4',               &
     &                      D_p5NH4(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'NitriR',                &
     &                      NitriR(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'K_NO3',                 &
     &                      K_NO3(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'K_NH4',                 &
     &                      K_NH4(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'K_Phy',                 &
     &                      K_Phy(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Chl2C_m',               &
     &                      Chl2C_m(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ChlMin',                &
     &                      ChlMin(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'PhyCN',                 &
     &                      PhyCN(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'PhyIP',                 &
     &                      PhyIP(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'PhyIS',                 &
     &                      PhyIS(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'PhyMin',                &
     &                      PhyMin(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'PhyMR',                 &
     &                      PhyMR(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ZooAE_N',               &
     &                      ZooAE_N(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ZooBM',                 &
     &                      ZooBM(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ZooCN',                 &
     &                      ZooCN(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ZooER',                 &
     &                      ZooER(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ZooGR',                 &
     &                      ZooGR(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ZooMin',                &
     &                      ZooMin(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ZooMR',                 &
     &                      ZooMR(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'LDeRRN',                &
     &                      LDeRRN(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'LDeRRC',                &
     &                      LDeRRC(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'CoagR',                 &
     &                      CoagR(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'SDeRRN',                &
     &                      SDeRRN(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'SDeRRC',                &
     &                      SDeRRC(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'wPhy',                  &
     &                      wPhy(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'wLDet',                 &
     &                      wLDet(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'wSDet',                 &
     &                      wSDet(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'pCO2air',               &
     &                      pCO2air(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

#ifdef SPECTRAL_LIGHT
      CALL netcdf_put_fvar (ng, model, ncname, 'SIGATRB',               &
     &                      SIGATRB(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'STRB',                  &
     &                      STRB(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'BLTRB',                 &
     &                      BLTRB(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'SIGBTRB',               &
     &                      SIGBTRB(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ETASPEC',               &
     &                      ETASPEC(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'BB2B',                  &
     &                      BB2B(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#endif
#ifdef SAV_BIOMASS 
      CALL netcdf_put_ivar (ng, model, ncname, 'GMODopt',               &
     &                      GMODopt(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'KNSED',                 &
     &                      KNSED(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'KNWC',                  &
     &                      KNWC(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'TOPT',                  &
     &                      TOPT(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'THTA',                  &
     &                      THTA(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'THTA2',                 &
     &                      THTA2(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'SCL',                   &
     &                      SCL(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'SCL2',                  &
     &                      SCL2(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'KI',                    &
     &                      KI(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'SR',                    &
     &                      SR(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'LMBAMX',                &
     &                      LMBAMX(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'KMAG',                  &
     &                      KMAG(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ARSC',                  &
     &                      ARSC(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ARC',                   &
     &                      ARC(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'BSRC',                  &
     &                      BSRC(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'RC',                    &
     &                      RC(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'RtStTL',                &
     &                      RtStTL(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'DOWNt',                 &
     &                      DOWNt(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'TRNS',                  &
     &                      TRNS(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'TCRIT',                 &
     &                      TCRIT(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'KM',                    &
     &                      KM(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#endif 
