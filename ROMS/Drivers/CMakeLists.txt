# git $Id$
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2026 The ROMS Group                                :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "ROMS/Drivers"

list( APPEND _files
      ROMS/Drivers/i4dvar.F
      ROMS/Drivers/r4dvar.F
      ROMS/Drivers/rbl4dvar.F
)

set ( ROMS_Drivers_files
      ${_files}
      PARENT_SCOPE
)
