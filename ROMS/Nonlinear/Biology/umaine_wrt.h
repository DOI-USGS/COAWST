!
!  Write out UMaine CoSiNE ecosystem model parameters.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'BioIter',               &
     &                      BioIter(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'reg1',                  &
     &                      reg1(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'reg2',                  &
     &                      reg2(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gmaxs1',                &
     &                      gmaxs1(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gmaxs2',                &
     &                      gmaxs2(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gmaxs3',                &
     &                      gmaxs3(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'beta1',                 &
     &                      beta1(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'beta2',                 &
     &                      beta2(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'akz1',                  &
     &                      akz1(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'akz2',                  &
     &                      akz2(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'PARfrac',               &
     &                      PARfrac(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'alphachl_s1',           &
     &                      alphachl_s1(ng), (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'alphachl_s2',           &
     &                      alphachl_s2(ng), (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'alphachl_s3',           &
     &                      alphachl_s3(ng), (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'pis1',                  &
     &                      pis1(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'pis2',                  &
     &                      pis2(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'pis3',                  &
     &                      pis3(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'akno3s1',               &
     &                      akno3s1(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'akno3s2',               &
     &                      akno3s2(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'akno3s3',               &
     &                      akno3s3(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'aknh4s1',               &
     &                      aknh4s1(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'aknh4s2',               &
     &                      aknh4s2(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'aknh4s3',               &
     &                      aknh4s3(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'akpo4s1',               &
     &                      akpo4s1(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'akpo4s2',               &
     &                      akpo4s2(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'akpo4s3',               &
     &                      akpo4s3(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'akco2s1',               &
     &                      akco2s1(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'akco2s2',               &
     &                      akco2s2(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'akco2s3',               &
     &                      akco2s3(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'aksio4s2',              &
     &                      aksio4s2(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ES1',                   &
     &                      ES1(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ES2',                   &
     &                      ES2(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ES3',                   &
     &                      ES3(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ak1',                   &
     &                      ak1(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ak2',                   &
     &                      ak2(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Qmax',                  &
     &                      Qmax(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Qmin',                  &
     &                      Qmin(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN	

      CALL netcdf_put_fvar (ng, model, ncname, 'lambdano3_s1',          &
     &                      lambdano3_s1(ng), (/0/), (/0/),             &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'lambdano3_s2',          &
     &                      lambdano3_s2(ng), (/0/), (/0/),             &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'lambdano3_s3',          &
     &                      lambdano3_s3(ng), (/0/), (/0/),             &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN	

      CALL netcdf_put_fvar (ng, model, ncname, 'thetaNmax_s1',          &
     &                      thetaNmax_s1(ng), (/0/), (/0/),             &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'thetaNmax_s2',          &
     &                      thetaNmax_s2(ng), (/0/), (/0/),             &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'thetaNmax_s3',          &
     &                      thetaNmax_s3(ng), (/0/), (/0/),             &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'bgamma',                &
     &                      bgamma(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'bgamma1',               &
     &                      bgamma1(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'bgamma2',               &
     &                      bgamma2(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'bgamma22',              &
     &                      bgamma22(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'bgamma3',               &
     &                      bgamma3(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'bgamma4',               &
     &                      bgamma4(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'bgamma10',              &
     &                      bgamma10(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'bgamma12',              &
     &                      bgamma12(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'bgamma5',               &
     &                      bgamma5(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'bgamma7',               &
     &                      bgamma7(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      	
      CALL netcdf_put_fvar (ng, model, ncname, 'bgamma11',              &
     &                      bgamma11(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN	

      CALL netcdf_put_fvar (ng, model, ncname, 'bgamma13',              &
     &                      bgamma13(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'mtos1',                &
     &                      mtos1(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN	

      CALL netcdf_put_fvar (ng, model, ncname, 'mtos2',                &
     &                      mtos2(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'mtos3',                &
     &                      mtos3(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'flz1',                 &
     &                      flz1(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'flz2',                 &
     &                      flz2(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'lk1',                  &
     &                      lk1(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'lk2',                  &
     &                      lk2(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'lk3',                  &
     &                      lk3(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ratiol1',              &
     &                      ratiol1(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ratiol2',              &
     &                      ratiol2(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'wsdn',                  &
     &                      wsdn(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'wsdc',                  &
     &                      wsdc(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'wsdsi',                 &
     &                      wsdsi(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'wsdca',                 &
     &                      wsdca(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'wsp1',                  &
     &                      wsp1(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'wsp2',                  &
     &                      wsp2(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'wsp3',                  &
     &                      wsp3(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'pco2a',                 &
     &                      pco2a(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'p2n',                   &
     &                      p2n(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'o2no',                  &
     &                      o2no(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'o2nh',                  &
     &                      o2nh(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'cnb',                   &
     &                      cnb(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'apsilon',               &
     &                      apsilon(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ro5',                   &
     &                      ro5(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ro6',                   &
     &                      ro6(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ro7',                   &
     &                      ro7(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ro10',                  &
     &                      ro10(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'rop',                   &
     &                      rop(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'rob',                   &
     &                      rob(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'kabac',                 &
     &                      kabac(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'klbac',                 &
     &                      klbac(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ksdoc',                 &
     &                      ksdoc(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ksdon',                 &
     &                      ksdon(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ratiob',                &
     &                      ratiob(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ratiobc',               &
     &                      ratiobc(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'RtUVLDOC',              &
     &                      RtUVLDOC(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'RtUVSDOC',              &
     &                      RtUVSDOC(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'RtUVLDIC',              &
     &                      RtUVLDIC(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN 	

      CALL netcdf_put_fvar (ng, model, ncname, 'RtUVSDIC',              &
     &                      RtUVSDIC(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'colorFR1',              &
     &                      colorFR1(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN	

      CALL netcdf_put_fvar (ng, model, ncname, 'colorFR2',              &
     &                      colorFR2(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

#ifdef IRON_LIMIT
      CALL netcdf_put_fvar (ng, model, ncname, 'T_Fe',                  &
     &                      T_Fe(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'A_Fe',                  &
     &                      A_Fe(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'B_Fe',                  &
     &                      B_Fe(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'S1_FeC',                &
     &                      S1_FeC(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'S2_FeC',                &
     &                      S2_FeC(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'S3_FeC',                &
     &                      S3_FeC(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'FeRR',                  &
     &                      FeRR(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#endif
