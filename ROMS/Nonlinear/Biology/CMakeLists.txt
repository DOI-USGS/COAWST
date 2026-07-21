# git $Id$
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2026 The ROMS Group                                :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "ROMS/Nonlinear/Biology"

list( APPEND _files
      ROMS/Nonlinear/Biology/biology.F
      ROMS/Nonlinear/Biology/biology_floats.F
)

set ( ROMS_Nonlinear_Biology_files
      ${_files}
      PARENT_SCOPE
)
