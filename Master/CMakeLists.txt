# svn $Id: CMakeLists.txt 1054 2021-03-06 19:47:12Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2021 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "Master"

list( APPEND _files
      Master/coupler.F
      Master/esmf_roms.F
      Master/ocean_control.F
      Master/propagator.F
      Master/roms_export.F
      Master/roms_import.F
)

set( Master_files
     ${_files}
     PARENT_SCOPE
)
