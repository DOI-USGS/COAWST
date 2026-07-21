# git $Id$
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2026 The ROMS Group                                :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "ROMS/Adjoint/Biology"

list( APPEND _files
      ROMS/Adjoint/Biology/ad_biology.F
)

set ( ROMS_Adjoint_Biology_files
      ${_files}
      PARENT_SCOPE
)
