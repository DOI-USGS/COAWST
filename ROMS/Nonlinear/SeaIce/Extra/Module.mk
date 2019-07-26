# svn $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
# Copyright (c) 2002-2016 The ROMS/TOMS Group             Kate Hedstrom :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

local_sub  := SeaIce/Extra

local_lib  := libExtra.a
local_src  := $(wildcard $(local_sub)/*.F)

$(eval $(call make-library,$(local_lib),$(local_src)))

$(eval $(compile-rules))
