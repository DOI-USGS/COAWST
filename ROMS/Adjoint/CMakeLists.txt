# git $Id$
# svn $Id: CMakeLists.txt 1151 2023-02-09 03:08:53Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "ROMS/Adjoint"

list( APPEND _files
      ROMS/Adjoint/ad_balance.F
      ROMS/Adjoint/ad_bc_2d.F
      ROMS/Adjoint/ad_bc_3d.F
      ROMS/Adjoint/ad_bc_bry2d.F
      ROMS/Adjoint/ad_bc_bry3d.F
      ROMS/Adjoint/ad_bulk_flux.F
      ROMS/Adjoint/ad_congrad.F
      ROMS/Adjoint/ad_conv_2d.F
      ROMS/Adjoint/ad_conv_3d.F
      ROMS/Adjoint/ad_conv_bry2d.F
      ROMS/Adjoint/ad_conv_bry3d.F
      ROMS/Adjoint/ad_convolution.F
      ROMS/Adjoint/ad_def_his.F
      ROMS/Adjoint/ad_diag.F
      ROMS/Adjoint/ad_exchange_2d.F
      ROMS/Adjoint/ad_exchange_3d.F
      ROMS/Adjoint/ad_extract_obs.F
      ROMS/Adjoint/ad_force_dual.F
      ROMS/Adjoint/ad_forcing.F
      ROMS/Adjoint/ad_frc_adjust.F
      ROMS/Adjoint/ad_get_data.F
      ROMS/Adjoint/ad_get_idata.F
      ROMS/Adjoint/ad_htobs.F
      ROMS/Adjoint/ad_ini_fields.F
      ROMS/Adjoint/ad_initial.F
      ROMS/Adjoint/ad_lmd_swfrac.F
      ROMS/Adjoint/ad_main2d.F
      ROMS/Adjoint/ad_main3d.F
      ROMS/Adjoint/ad_misfit.F
      ROMS/Adjoint/ad_nesting.F
      ROMS/Adjoint/ad_obc_adjust.F
      ROMS/Adjoint/ad_obc_volcons.F
      ROMS/Adjoint/ad_omega.F
      ROMS/Adjoint/ad_output.F
      ROMS/Adjoint/ad_pre_step3d.F
      ROMS/Adjoint/ad_prsgrd.F
      ROMS/Adjoint/ad_rho_eos.F
      ROMS/Adjoint/ad_rhs3d.F
      ROMS/Adjoint/ad_rpcg_lanczos.F
      ROMS/Adjoint/adsen_force.F
      ROMS/Adjoint/adsen_initial.F
      ROMS/Adjoint/ad_set_avg.F
      ROMS/Adjoint/ad_set_data.F
      ROMS/Adjoint/ad_set_depth.F
      ROMS/Adjoint/ad_set_massflux.F
      ROMS/Adjoint/ad_set_vbc.F
      ROMS/Adjoint/ad_set_zeta.F
      ROMS/Adjoint/ad_sqlq.F
      ROMS/Adjoint/ad_step2d.F
      ROMS/Adjoint/ad_step3d_t.F
      ROMS/Adjoint/ad_step3d_uv.F
      ROMS/Adjoint/ad_t3dbc_im.F
      ROMS/Adjoint/ad_t3dmix.F
      ROMS/Adjoint/ad_t3drelax.F
      ROMS/Adjoint/ad_u2dbc_im.F
      ROMS/Adjoint/ad_u3dbc_im.F
      ROMS/Adjoint/ad_uv3dmix.F
      ROMS/Adjoint/ad_uv3drelax.F
      ROMS/Adjoint/ad_v2dbc_im.F
      ROMS/Adjoint/ad_v3dbc_im.F
      ROMS/Adjoint/ad_variability.F
      ROMS/Adjoint/ad_wrt_his.F
      ROMS/Adjoint/ad_wvelocity.F
ROMS/Adjoint/ad_zetabc.F
)

set( ROMS_Adjoint_files
     ${_files}
     PARENT_SCOPE
)
