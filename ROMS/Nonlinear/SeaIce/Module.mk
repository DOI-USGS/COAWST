# svn $Id: Module.mk 889 2018-02-10 03:32:52Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
# Copyright (c) 2002-2019 The ROMS/TOMS Group             Kate Hedstrom :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

local_sub  := ROMS/SeaIce

local_lib  := libICE.a
local_src  := $(wildcard $(local_sub)/*.F)

$(eval $(call make-library,$(local_lib),$(local_src)))

$(eval $(compile-rules))
