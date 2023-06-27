# git $Id$
# svn $Id: CMakeLists.txt 1151 2023-02-09 03:08:53Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "ROMS/Tangent"

list( APPEND _files
      ROMS/Tangent/tl_balance.F
      ROMS/Tangent/tl_bulk_flux.F
      ROMS/Tangent/tl_congrad.F
      ROMS/Tangent/tl_conv_2d.F
      ROMS/Tangent/tl_conv_3d.F
      ROMS/Tangent/tl_conv_bry2d.F
      ROMS/Tangent/tl_conv_bry3d.F
      ROMS/Tangent/tl_convolution.F
      ROMS/Tangent/tl_def_his.F
      ROMS/Tangent/tl_def_ini.F
      ROMS/Tangent/tl_diag.F
      ROMS/Tangent/tl_forcing.F
      ROMS/Tangent/tl_frc_adjust.F
      ROMS/Tangent/tl_get_data.F
      ROMS/Tangent/tl_get_idata.F
      ROMS/Tangent/tl_ini_fields.F
      ROMS/Tangent/tl_initial.F
      ROMS/Tangent/tl_lmd_swfrac.F
      ROMS/Tangent/tl_main2d.F
      ROMS/Tangent/tl_main3d.F
      ROMS/Tangent/tl_nesting.F
      ROMS/Tangent/tl_obc_adjust.F
      ROMS/Tangent/tl_obc_volcons.F
      ROMS/Tangent/tl_omega.F
      ROMS/Tangent/tl_output.F
      ROMS/Tangent/tl_pre_step3d.F
      ROMS/Tangent/tl_prsgrd.F
      ROMS/Tangent/tl_rho_eos.F
      ROMS/Tangent/tl_rhs3d.F
      ROMS/Tangent/tl_rpcg_lanczos.F
      ROMS/Tangent/tl_set_avg.F
      ROMS/Tangent/tl_set_data.F
      ROMS/Tangent/tl_set_depth.F
      ROMS/Tangent/tl_set_massflux.F
      ROMS/Tangent/tl_set_vbc.F
      ROMS/Tangent/tl_set_zeta.F
      ROMS/Tangent/tl_sqlq.F
      ROMS/Tangent/tl_step2d.F
      ROMS/Tangent/tl_step3d_t.F
      ROMS/Tangent/tl_step3d_uv.F
      ROMS/Tangent/tl_t3dbc_im.F
      ROMS/Tangent/tl_t3dmix.F
      ROMS/Tangent/tl_t3drelax.F
      ROMS/Tangent/tl_u2dbc_im.F
      ROMS/Tangent/tl_u3dbc_im.F
      ROMS/Tangent/tl_uv3dmix.F
      ROMS/Tangent/tl_uv3drelax.F
      ROMS/Tangent/tl_v2dbc_im.F
      ROMS/Tangent/tl_v3dbc_im.F
      ROMS/Tangent/tl_variability.F
      ROMS/Tangent/tl_wrt_his.F
      ROMS/Tangent/tl_wrt_ini.F
      ROMS/Tangent/tl_zetabc.F
)

set( ROMS_Tangent_files
     ${_files}
     PARENT_SCOPE
)
