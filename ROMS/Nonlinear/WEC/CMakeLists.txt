# git $Id$
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2026 The ROMS Group                                :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "ROMS/Nonlinear/WEC"

list( APPEND _files
      ROMS/Nonlinear/Wec/wec_banihashemi
      ROMS/Nonlinear/Wec/wec_dissip.F
      ROMS/Nonlinear/Wec/wec_output.F
      ROMS/Nonlinear/Wec/wec_roller.F
      ROMS/Nonlinear/Wec/wec_stokes.F
      ROMS/Nonlinear/Wec/wec_streaming.F
      ROMS/Nonlinear/Wec/wec_u2dbc_im.F
      ROMS/Nonlinear/Wec/wec_u3dbc_im.F
      ROMS/Nonlinear/Wec/wec_v2dbc_im.F
      ROMS/Nonlinear/Wec/wec_v3dbc_im.F
      ROMS/Nonlinear/Wec/wec_vf.F
      ROMS/Nonlinear/Wec/wec_wave_mix.F
      ROMS/Nonlinear/Wec/wec_wvelocity.F
)

set ( ROMS_Nonlinear_Wec_files
      ${_files}
      PARENT_SCOPE
)
