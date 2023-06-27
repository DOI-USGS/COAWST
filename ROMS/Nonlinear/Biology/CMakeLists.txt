# git $Id$
# svn $Id: CMakeLists.txt 1151 2023-02-09 03:08:53Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
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
