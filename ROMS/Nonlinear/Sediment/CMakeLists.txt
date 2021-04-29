# svn $Id: CMakeLists.txt 1054 2021-03-06 19:47:12Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2021 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "ROMS/Nonlinear/Sediment"

list( APPEND _files
      ROMS/Nonlinear/Sediment/sed_bed.F
      ROMS/Nonlinear/Sediment/sed_bedload.F
      ROMS/Nonlinear/Sediment/sed_fluxes.F
      ROMS/Nonlinear/Sediment/sediment.F
      ROMS/Nonlinear/Sediment/sed_settling.F
      ROMS/Nonlinear/Sediment/sed_surface.F
)

set ( ROMS_Nonlinear_Sediment_files
      ${_files}
      PARENT_SCOPE
)
