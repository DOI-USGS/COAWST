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
      ROMS/Nonlinear/Vegetation/marsh_dynamics.F
      ROMS/Nonlinear/Vegetation/marsh_sed_erosion.F
      ROMS/Nonlinear/Vegetation/marsh_tidal_range_calc.F
      ROMS/Nonlinear/Vegetation/marsh_vert_growth.F
      ROMS/Nonlinear/Vegetation/marsh_wave_thrust.F
      ROMS/Nonlinear/Vegetation/vegetation_biomass.F
      ROMS/Nonlinear/Vegetation/vegetation_drag.F
      ROMS/Nonlinear/Vegetation/vegetation_stream.F
      ROMS/Nonlinear/Vegetation/vegetation_turb_cal.F     
)

set ( ROMS_Nonlinear_Vegetation_files
      ${_files}
      PARENT_SCOPE
)
