# svn $Id: Module.mk 1311 2008-01-10 04:13:52Z jcwarner $

local_sub  := InWave/Boundaries

local_lib  := libInWave_Boundaries.a
local_src  := $(wildcard $(local_sub)/*.F)

$(eval $(call make-library,$(local_lib),$(local_src)))

$(eval $(compile-rules))
