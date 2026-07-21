/*
** git $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2026 The ROMS Group                             **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.md                                              **
************************************************************************
**                                                                    **
**  Writes the Estuarine Carbon Biogeochemistry model (Feng et al.,   **
**  2015, and follow-ups) parameters into output NetCDF files using   **
**  the PIO library. It is included in routine "wrt_info.F".          **
**                                                                    **
************************************************************************
*/

!
!  Write out the Estuarine Carbon Biogeochemistry model (Feng et al.,
!  2015, and follow-ups) configuration parameters.
!                 
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'BioIter',           &
     &                          BioIter(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'AttSW',             &
     &                          AttSW(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'AttChl',            &
     &                          AttChl(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'PARfrac',           &
     &                          PARfrac(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Vp0',               &
     &                          Vp0(ng), (/0/), (/0/),                  &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'I_thNH4',           &
     &                          I_thNH4(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'D_p5NH4',           &
     &                          D_p5NH4(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'NitriR',            &
     &                          NitriR(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'K_NO3',             &
     &                          K_NO3(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'K_NH4',             &
     &                          K_NH4(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'K_PO4',             &
     &                          K_PO4(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'K_Phy',             &
     &                          K_Phy(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Chl2C_m',           &
     &                          Chl2C_m(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'ChlMin',            &
     &                          ChlMin(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'PhyCN',             &
     &                          PhyCN(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'R_P2N',             &
     &                          R_P2N(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'PhyIP',             &
     &                          PhyIP(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'PhyIS',             &
     &                          PhyIS(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'PhyMin',            &
     &                          PhyMin(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'PhyMR',             &
     &                          PhyMR(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'ZooAE_N',           &
     &                          ZooAE_N(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'ZooBM',             &
     &                          ZooBM(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'ZooCN',             &
     &                          ZooCN(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'ZooER',             &
     &                          ZooER(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'ZooGR',             &
     &                          ZooGR(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'ZooMin',            &
     &                          ZooMin(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'ZooMR',             &
     &                          ZooMR(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'LDeRRN',            &
     &                          LDeRRN(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'LDeRRC',            &
     &                          LDeRRC(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'CoagR',             &
     &                          CoagR(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'SDeRRN',            &
     &                          SDeRRN(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'SDeRRC',            &
     &                          SDeRRC(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'RDeRRN',            &
     &                          RDeRRN(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'RDeRRC',            &
     &                          RDeRRC(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'wPhy',              &
     &                          wPhy(ng), (/0/), (/0/),                 &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'wLDet',             &
     &                          wLDet(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'wSDet',             &
     &                          wSDet(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'rkd1',              &
     &                          rkd1(ng), (/0/), (/0/),                 &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'rkdTSS1',           &
     &                          rkdTSS1(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'rkdS1',             &
     &                          rkdS1(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'ElDON',             &
     &                          ElDON(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'EsDON',             &
     &                          EsDON(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'deltN',             &
     &                          deltN(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'deltC',             &
     &                          deltC(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gammaC',            &
     &                          gammaC(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'slCexc',            &
     &                          slCexc(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'pCO2air',           &
     &                          pCO2air(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
