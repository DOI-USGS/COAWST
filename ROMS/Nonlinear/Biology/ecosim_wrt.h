/*
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2016 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Writes EcoSim bio-optical model input parameters into output      **
**  NetCDF files. It is included in routine "wrt_info.F".             **
**                                                                    **
************************************************************************
*/

!
!  Write out EcoSim bio-optical model parameters.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'BioIter',               &
     &                      BioIter(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_lvar (ng, model, ncname, 'RtUVR_flag',            &
     &                      RtUVR_flag(ng), (/0/), (/0/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_lvar (ng, model, ncname, 'NFIX_flag',             &
     &                      NFIX_flag(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_lvar (ng, model, ncname, 'Regen_flag',            &
     &                      Regen_flag(ng), (/0/), (/0/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'HsNO3',                 &
     &                      HsNO3(:,ng), (/1/), (/Nphy/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'HsNH4',                 &
     &                      HsNH4(:,ng), (/1/), (/Nphy/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'HsSiO',                 &
     &                      HsSiO(:,ng), (/1/), (/Nphy/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'HsPO4',                 &
     &                      HsPO4(:,ng), (/1/), (/Nphy/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'HsFe',                  &
     &                      HsFe(:,ng), (/1/), (/Nphy/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'GtALG_max',             &
     &                      GtALG_max(:,ng), (/1/), (/Nphy/),           &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'PhyTbase',              &
     &                      PhyTbase(:,ng), (/1/), (/Nphy/),            &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'PhyTfac',               &
     &                      PhyTfac(:,ng), (/1/), (/Nphy/),             &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'BET_',                  &
     &                      BET_(:,ng), (/1/), (/Nphy/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'maxC2nALG',             &
     &                      maxC2nALG(:,ng), (/1/), (/Nphy/),           &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'minC2nALG',             &
     &                      minC2nALG(:,ng), (/1/), (/Nphy/),           &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'C2nALGminABS',          &
     &                      C2nALGminABS(:,ng), (/1/), (/Nphy/),        &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'maxC2SiALG',            &
     &                      maxC2SiALG(:,ng), (/1/), (/Nphy/),          &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'minC2SiALG',            &
     &                      minC2SiALG(:,ng), (/1/), (/Nphy/),          &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'C2SiALGminABS',         &
     &                      C2SiALGminABS(:,ng), (/1/), (/Nphy/),       &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'maxC2pALG',             &
     &                      maxC2pALG(:,ng), (/1/), (/Nphy/),           &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'minC2pALG',             &
     &                      minC2pALG(:,ng), (/1/), (/Nphy/),           &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'C2pALGminABS',          &
     &                      C2pALGminABS(:,ng), (/1/), (/Nphy/),        &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'maxC2FeALG',            &
     &                      maxC2FeALG(:,ng), (/1/), (/Nphy/),          &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'minC2FeALG',            &
     &                      minC2FeALG(:,ng), (/1/), (/Nphy/),          &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'C2FeALGminABS',         &
     &                      C2FeALGminABS(:,ng), (/1/), (/Nphy/),       &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'qu_yld',                &
     &                      qu_yld(:,ng), (/1/), (/Nphy/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'E0_comp',               &
     &                      E0_comp(:,ng), (/1/), (/Nphy/),             &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'E0_inhib',              &
     &                      E0_inhib(:,ng), (/1/), (/Nphy/),            &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'inhib_fac',             &
     &                      inhib_fac(:,ng), (/1/), (/Nphy/),           &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'C2CHL_max',             &
     &                      C2CHL_max(:,ng), (/1/), (/Nphy/),           &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'mxC2Cl',                &
     &                      mxC2Cl(:,ng), (/1/), (/Nphy/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'b_C2Cl',                &
     &                      b_C2Cl(:,ng), (/1/), (/Nphy/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'mxC2Cn',                &
     &                      mxC2Cn(:,ng), (/1/), (/Nphy/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'b_C2Cn',                &
     &                      b_C2Cn(:,ng), (/1/), (/Nphy/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'mxPacEff',              &
     &                      mxPacEff(:,ng), (/1/), (/Nphy/),            &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'b_PacEff',              &
     &                      b_PacEff(:,ng), (/1/), (/Nphy/),            &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'mxChlB',                &
     &                      mxChlB(:,ng), (/1/), (/Nphy/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'b_ChlB',                &
     &                      b_ChlB(:,ng), (/1/), (/Nphy/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'mxChlC',                &
     &                      mxChlC(:,ng), (/1/), (/Nphy/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'b_ChlC',                &
     &                      b_ChlC(:,ng), (/1/), (/Nphy/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'mxPSC',                 &
     &                      mxPSC(:,ng), (/1/), (/Nphy/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'b_PSC',                 &
     &                      b_PSC(:,ng), (/1/), (/Nphy/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'mxPPC',                 &
     &                      mxPPC(:,ng), (/1/), (/Nphy/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'b_PPC',                 &
     &                      b_PPC(:,ng), (/1/), (/Nphy/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'mxLPUb',                &
     &                      mxLPUb(:,ng), (/1/), (/Nphy/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'b_LPUb',                &
     &                      b_LPUb(:,ng), (/1/), (/Nphy/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'mxHPUb',                &
     &                      mxHPUb(:,ng), (/1/), (/Nphy/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'b_HPUb',                &
     &                      b_HPUb(:,ng), (/1/), (/Nphy/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'FecDOC',                &
     &                      FecDOC(:,ng), (/1/), (/Nphy/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'FecPEL',                &
     &                      FecPEL(:,:,ng), (/1,1/), (/Nphy,Nfec/),     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'FecCYC',                &
     &                      FecCYC(:,ng), (/1/), (/Nphy/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ExALG',                 &
     &                      ExALG(:,ng), (/1/), (/Nphy/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'WS',                    &
     &                      WS(:,ng), (/1/), (/Nphy/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'HsGRZ',                 &
     &                      HsGRZ(:,ng), (/1/), (/Nphy/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'MinRefuge',             &
     &                      MinRefuge(:,ng),  (/1/), (/Nphy/),          &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'RefugeDep',             &
     &                      RefugeDep(:,ng), (/1/), (/Nphy/),           &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Norm_Vol',              &
     &                      Norm_Vol(:,ng), (/1/), (/Nphy/),            &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Norm_Surf',             &
     &                      Norm_Surf(:,ng), (/1/), (/Nphy/),           &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'HsDOP',                 &
     &                      HsDOP(:,ng), (/1/), (/Nphy/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'C2pALKPHOS',            &
     &                      C2pALKPHOS(:,ng), (/1/), (/Nphy/),          &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'HsDON',                 &
     &                      HsDON(:,ng), (/1/), (/Nphy/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'C2nNupDON',             &
     &                      C2nNupDON(:,ng), (/1/), (/Nphy/),           &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'HsDOC_ba',              &
     &                      HsDOC_ba(:,ng), (/1/), (/Nbac/),            &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'GtBAC_max',             &
     &                      GtBAC_max(:,ng), (/1/), (/Nbac/),           &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'BacTbase',              &
     &                      BacTbase(:,ng), (/1/), (/Nbac/),            &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'BacTfac',               &
     &                      BacTfac(:,ng), (/1/), (/Nbac/),             &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'C2nBAC',                &
     &                      C2nBAC(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'C2pBAC',                &
     &                      C2pBAC(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'C2FeBAC',               &
     &                      C2FeBAC(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'BacDOC',                &
     &                      BacDOC(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'BacPEL',                &
     &                      BacPEL(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'BacCYC',                &
     &                      BacCYC(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ExBAC_c',               &
     &                      ExBAC_c(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ExBacC2N',              &
     &                      ExBacC2N(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Bac_Ceff',              &
     &                      Bac_Ceff(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'RtNIT',                 &
     &                      RtNIT(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'HsNIT',                 &
     &                      HsNIT(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'cDOCfrac_c',            &
     &                      cDOCfrac_c(:,ng), (/1/), (/Ndom/),          &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'RtUVR_DIC',             &
     &                      RtUVR_DIC(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'RtUVR_DOC',             &
     &                      RtUVR_DIC(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'WF',                    &
     &                      WF(:,ng), (/1/), (/Nfec/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'RegTbase',              &
     &                      RegTbase(:,ng), (/1/), (/Nfec/),            &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'RegTfac',               &
     &                      RegTfac(:,ng), (/1/), (/Nfec/),             &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'RegCR',                 &
     &                      RegCR(:,ng), (/1/), (/Nfec/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'RegNR',                 &
     &                      RegNR(:,ng), (/1/), (/Nfec/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'RegSR',                 &
     &                      RegSR(:,ng), (/1/), (/Nfec/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'RegPR',                 &
     &                      RegPR(:,ng), (/1/), (/Nfec/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'RegFR',                 &
     &                      RegFR(:,ng), (/1/), (/Nfec/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
