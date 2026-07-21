# git $Id$
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2026 The ROMS Group                                :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "ROMS/Functionals"

list( APPEND _files
      ROMS/Functionals/analytical.F
)

set( ROMS_Functionals_files
     ${_files}
     PARENT_SCOPE
)
