# svn $Id: CMakeLists.txt 1054 2021-03-06 19:47:12Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2021 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "ROMS/Nonlinear"

list( APPEND _files
      ROMS/Nonlinear/bbl.F
      ROMS/Nonlinear/bc_2d.F
      ROMS/Nonlinear/bc_3d.F
      ROMS/Nonlinear/bc_4d.F
      ROMS/Nonlinear/bc_bry2d.F
      ROMS/Nonlinear/bc_bry3d.F
      ROMS/Nonlinear/bulk_flux.F
      ROMS/Nonlinear/bvf_mix.F
      ROMS/Nonlinear/conv_2d.F
      ROMS/Nonlinear/conv_3d.F
      ROMS/Nonlinear/conv_bry2d.F
      ROMS/Nonlinear/conv_bry3d.F
      ROMS/Nonlinear/diag.F
      ROMS/Nonlinear/exchange_2d.F
      ROMS/Nonlinear/exchange_3d.F
      ROMS/Nonlinear/exchange_4d.F
      ROMS/Nonlinear/forcing.F
      ROMS/Nonlinear/frc_adjust.F
      ROMS/Nonlinear/get_data.F
      ROMS/Nonlinear/get_idata.F
      ROMS/Nonlinear/gls_corstep.F
      ROMS/Nonlinear/gls_prestep.F
      ROMS/Nonlinear/hmixing.F
      ROMS/Nonlinear/ini_fields.F
      ROMS/Nonlinear/initial.F
      ROMS/Nonlinear/interp_floats.F
      ROMS/Nonlinear/lmd_bkpp.F
      ROMS/Nonlinear/lmd_skpp.F
      ROMS/Nonlinear/lmd_swfrac.F
      ROMS/Nonlinear/lmd_vmix.F
      ROMS/Nonlinear/main2d.F
      ROMS/Nonlinear/main3d.F
      ROMS/Nonlinear/mpdata_adiff.F
      ROMS/Nonlinear/my25_corstep.F
      ROMS/Nonlinear/my25_prestep.F
      ROMS/Nonlinear/nesting.F
      ROMS/Nonlinear/obc_adjust.F
      ROMS/Nonlinear/obc_volcons.F
      ROMS/Nonlinear/omega.F
      ROMS/Nonlinear/output.F
      ROMS/Nonlinear/pre_step3d.F
      ROMS/Nonlinear/prsgrd.F
      ROMS/Nonlinear/radiation_stress.F
      ROMS/Nonlinear/rho_eos.F
      ROMS/Nonlinear/rhs3d.F
      ROMS/Nonlinear/set_avg.F
      ROMS/Nonlinear/set_data.F
      ROMS/Nonlinear/set_depth.F
      ROMS/Nonlinear/set_massflux.F
      ROMS/Nonlinear/set_tides.F
      ROMS/Nonlinear/set_vbc.F
      ROMS/Nonlinear/set_zeta.F
      ROMS/Nonlinear/step2d.F
      ROMS/Nonlinear/step3d_t.F
      ROMS/Nonlinear/step3d_uv.F
      ROMS/Nonlinear/step_floats.F
      ROMS/Nonlinear/t3dbc_im.F
      ROMS/Nonlinear/t3dmix.F
      ROMS/Nonlinear/tkebc_im.F
      ROMS/Nonlinear/u2dbc_im.F
      ROMS/Nonlinear/u3dbc_im.F
      ROMS/Nonlinear/uv3dmix.F
      ROMS/Nonlinear/v2dbc_im.F
      ROMS/Nonlinear/v3dbc_im.F
      ROMS/Nonlinear/vwalk_floats.F
      ROMS/Nonlinear/wetdry.F
      ROMS/Nonlinear/wvelocity.F
      ROMS/Nonlinear/zetabc.F
)

set ( ROMS_Nonlinear_files
      ${_files}
      PARENT_SCOPE
)
