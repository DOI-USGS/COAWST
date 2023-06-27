# git $Id$
# svn $Id: CMakeLists.txt 1151 2023-02-09 03:08:53Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "ROMS/Nonlinear/Sediment"

list( APPEND _files
      ROMS/Nonlinear/Sediment/mod_vandera_funcs.F
      ROMS/Nonlinear/Sediment/sed_bed.F
      ROMS/Nonlinear/Sediment/sed_bed2.F
      ROMS/Nonlinear/Sediment/sed_bed_cohesive.F
      ROMS/Nonlinear/Sediment/sed_bedload.F
      ROMS/Nonlinear/Sediment/sed_bedload_vandera.F
      ROMS/Nonlinear/Sediment/sed_biodiff.F
      ROMS/Nonlinear/Sediment/sed_flocs.F
      ROMS/Nonlinear/Sediment/sed_fluxes.F
      ROMS/Nonlinear/Sediment/sed_settling.F
      ROMS/Nonlinear/Sediment/sed_surface.F
      ROMS/Nonlinear/Sediment/sediment.F
      ROMS/Nonlinear/Sediment/sedtr_decay.F
      ROMS/Nonlinear/Sediment/sedtr_reactions_pom.F
      ROMS/Nonlinear/Sediment/sedtr_reactions_sed_decay.F
)

set ( ROMS_Nonlinear_Sediment_files
      ${_files}
      PARENT_SCOPE
)
