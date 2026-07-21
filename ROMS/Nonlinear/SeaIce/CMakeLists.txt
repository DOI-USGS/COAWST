# git $Id$
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2026 The ROMS Group                                :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "ROMS/Nonlinear/SeaIce"

list( APPEND _files
      ROMS/Nonlinear/SeaIce/ice_advect.F
      ROMS/Nonlinear/SeaIce/ice_albedo.F
      ROMS/Nonlinear/SeaIce/ice_bc2d.F
      ROMS/Nonlinear/SeaIce/ice_elastic.F
      ROMS/Nonlinear/SeaIce/ice_enthalpy.F
      ROMS/Nonlinear/SeaIce/ice_evp.F
      ROMS/Nonlinear/SeaIce/ice_evp_sig.F
      ROMS/Nonlinear/SeaIce/ice_frazil.F
      ROMS/Nonlinear/SeaIce/ice_limit.F
      ROMS/Nonlinear/SeaIce/ice_output.F
      ROMS/Nonlinear/SeaIce/ice_set_avg.F
      ROMS/Nonlinear/SeaIce/ice_spdiw.F
      ROMS/Nonlinear/SeaIce/ice_thermo.F
      ROMS/Nonlinear/SeaIce/ice_tibc.F
      ROMS/Nonlinear/SeaIce/ice_uibc.F
      ROMS/Nonlinear/SeaIce/ice_vbc.F
      ROMS/Nonlinear/SeaIce/ice_vibc.F
      ROMS/Nonlinear/SeaIce/seaice.F
)

set ( ROMS_Nonlinear_SeaIce_files
      ${_files}
      PARENT_SCOPE
)