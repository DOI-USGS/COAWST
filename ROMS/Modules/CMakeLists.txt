# git $Id$
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2026 The ROMS Group                                :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
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
       ROMS/Modules/mod_extract.F
       ROMS/Modules/mod_floats.F
       ROMS/Modules/mod_forces.F
       ROMS/Modules/mod_fourdvar.F
       ROMS/Modules/mod_grid.F
       ROMS/Modules/mod_ice.F
       ROMS/Modules/mod_iounits.F
       ROMS/Modules/mod_kinds.F
       ROMS/Modules/mod_mixing.F
       ROMS/Modules/mod_ncparam.F
       ROMS/Modules/mod_nesting.F
       ROMS/Modules/mod_netcdf.F
       ROMS/Modules/mod_ocean.F
       ROMS/Modules/mod_parallel.F
       ROMS/Modules/mod_param.F
       ROMS/Modules/mod_pio_netcdf.F
       ROMS/Modules/mod_scalars.F
       ROMS/Modules/mod_sedbed.F
       ROMS/Modules/mod_sedflocs.F
       ROMS/Modules/mod_sediment.F
       ROMS/Modules/mod_sources.F
       ROMS/Modules/mod_stepping.F
       ROMS/Modules/mod_storage.F
       ROMS/Modules/mod_strings.F
       ROMS/Modules/mod_tides.F
       ROMS/Modules/mod_vegetation.F
)

set ( ROMS_Modules_files
      ${_files}
      PARENT_SCOPE
)
