# svn $Id: CMakeLists.txt 1054 2021-03-06 19:47:12Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2021 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "ROMS/Modules"

list ( APPEND _files
       ROMS/Modules/mod_arrays.F
       ROMS/Modules/mod_average.F
       ROMS/Modules/mod_bbl.F
       ROMS/Modules/mod_behavior.F
       ROMS/Modules/mod_biology.F
       ROMS/Modules/mod_boundary.F
       ROMS/Modules/mod_clima.F
       ROMS/Modules/mod_coupler.F
       ROMS/Modules/mod_coupling.F
       ROMS/Modules/mod_diags.F
       ROMS/Modules/mod_eclight.F
       ROMS/Modules/mod_eoscoef.F
       ROMS/Modules/mod_floats.F
       ROMS/Modules/mod_forces.F
       ROMS/Modules/mod_fourdvar.F
       ROMS/Modules/mod_grid.F
       ROMS/Modules/mod_iounits.F
       ROMS/Modules/mod_kinds.F
       ROMS/Modules/mod_mixing.F
       ROMS/Modules/mod_ncparam.F
       ROMS/Modules/mod_nesting.F
       ROMS/Modules/mod_netcdf.F
       ROMS/Modules/mod_ocean.F
       ROMS/Modules/mod_parallel.F
       ROMS/Modules/mod_param.F
       ROMS/Modules/mod_scalars.F
       ROMS/Modules/mod_sedbed.F
       ROMS/Modules/mod_sediment.F
       ROMS/Modules/mod_sources.F
       ROMS/Modules/mod_stepping.F
       ROMS/Modules/mod_storage.F
       ROMS/Modules/mod_strings.F
       ROMS/Modules/mod_tides.F
)

set ( ROMS_Modules_files
      ${_files}
      PARENT_SCOPE
)
