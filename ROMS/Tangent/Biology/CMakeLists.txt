# git $Id$
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2026 The ROMS Group                                :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "ROMS/Tangent/Biology"

list( APPEND _files
      ROMS/Tangent/Biology/tl_biology.F
)

set( ROMS_Tangent_Biology_files
     ${_files}
     PARENT_SCOPE
)
