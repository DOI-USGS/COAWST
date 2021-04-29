# svn $Id: CMakeLists.txt 1054 2021-03-06 19:47:12Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2021 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
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
