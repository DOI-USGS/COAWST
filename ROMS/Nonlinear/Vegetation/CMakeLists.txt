# git $Id$
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2026 The ROMS Group                                :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "ROMS/Nonlinear/Vegetation"

list( APPEND _files
      ROMS/Nonlinear/Vegetation/marsh_dynamics.F
      ROMS/Nonlinear/Vegetation/marsh_sed_erosion.F
      ROMS/Nonlinear/Vegetation/marsh_tidal_range.F
      ROMS/Nonlinear/Vegetation/marsh_vert_growth.F
      ROMS/Nonlinear/Vegetation/marsh_wave_thrust.F
      ROMS/Nonlinear/Vegetation/vegetation_biomass.F
      ROMS/Nonlinear/Vegetation/vegetation_drag.F
      ROMS/Nonlinear/Vegetation/vegetation_output.F
      ROMS/Nonlinear/Vegetation/vegetation_stream.F
      ROMS/Nonlinear/Vegetation/vegetation_turb.F     
)

set ( ROMS_Nonlinear_Vegetation_files
      ${_files}
      PARENT_SCOPE
)
