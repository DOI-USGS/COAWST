# git $Id$
# svn $Id: CMakeLists.txt 1151 2023-02-09 03:08:53Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "ROMS/Representer"

list( APPEND _files
      ROMS/Representer/rp_bulk_flux.F
      ROMS/Representer/rp_def_ini.F
      ROMS/Representer/rp_diag.F
      ROMS/Representer/rp_frc_adjust.F
      ROMS/Representer/rp_get_data.F
      ROMS/Representer/rp_get_idata.F
      ROMS/Representer/rp_ini_fields.F
      ROMS/Representer/rp_initial.F
      ROMS/Representer/rp_lmd_swfrac.F
      ROMS/Representer/rp_main2d.F
      ROMS/Representer/rp_main3d.F
      ROMS/Representer/rp_obc_adjust.F
      ROMS/Representer/rp_obc_volcons.F
      ROMS/Representer/rp_omega.F
      ROMS/Representer/rp_output.F
      ROMS/Representer/rp_pre_step3d.F
      ROMS/Representer/rp_prsgrd.F
      ROMS/Representer/rp_rho_eos.F
      ROMS/Representer/rp_rhs3d.F
      ROMS/Representer/rp_set_data.F
      ROMS/Representer/rp_set_depth.F
      ROMS/Representer/rp_set_massflux.F
      ROMS/Representer/rp_set_vbc.F
      ROMS/Representer/rp_set_zeta.F
      ROMS/Representer/rp_step2d.F
      ROMS/Representer/rp_step3d_t.F
      ROMS/Representer/rp_step3d_uv.F
      ROMS/Representer/rp_t3dbc_im.F
      ROMS/Representer/rp_t3dmix.F
      ROMS/Representer/rp_t3drelax.F
      ROMS/Representer/rp_u2dbc_im.F
      ROMS/Representer/rp_u3dbc_im.F
      ROMS/Representer/rp_uv3dmix.F
      ROMS/Representer/rp_uv3drelax.F
      ROMS/Representer/rp_v2dbc_im.F
      ROMS/Representer/rp_v3dbc_im.F
      ROMS/Representer/rp_wrt_ini.F
      ROMS/Representer/rp_zetabc.F
)

set( ROMS_Representer_files
     ${_files}
     PARENT_SCOPE
)
