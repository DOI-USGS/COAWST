# git $Id$
# svn $Id: CMakeLists.txt 1151 2023-02-09 03:08:53Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "ROMS/Nonlinear/Wec"

list( APPEND _files
      ROMS/Nonlinear/Wec/radiation_stress.F
      ROMS/Nonlinear/Wec/us2dbc_im.F
      ROMS/Nonlinear/Wec/us3dbc_im.F
      ROMS/Nonlinear/Wec/vs2dbc_im.F
      ROMS/Nonlinear/Wec/vs3dbc_im.F
      ROMS/Nonlinear/Wec/wec_dissip.F
      ROMS/Nonlinear/Wec/wec_roller.F
      ROMS/Nonlinear/Wec/wec_stokes.F
      ROMS/Nonlinear/Wec/wec_streaming.F
      ROMS/Nonlinear/Wec/wec_vf.F
      ROMS/Nonlinear/Wec/wec_wave_mix.F
      ROMS/Nonlinear/Wec/wstvelocity.F
)

set ( ROMS_Nonlinear_Wec_files
      ${_files}
      PARENT_SCOPE
)
