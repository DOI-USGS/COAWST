# svn $Id: Module.mk 503 2008-01-10 00:11:51Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
# Copyright (c) 2002-2016 The ROMS/TOMS Group             Kate Hedstrom :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

local_sub  := Master

local_src  := $(wildcard $(local_sub)/*.F)
#local_objs := $(call source-to-object,$(local_src))
local_objs := $(subst .F,.o,$(local_src))
local_objs := $(addprefix $(SCRATCH_DIR)/, $(notdir $(local_objs)))

sources    += $(local_src)
ifeq "$(ROMS_APPLICATION)" "CIRCLE"
local_objs += $(addprefix $(SCRATCH_DIR)/, bessi.o)
endif

ifdef LD_WINDOWS
$(BIN):	$(libraries) $(local_objs)
	$(LD) $(FFLAGS) $(local_objs) -o $@ $(libraries) $(LIBS_WIN32) $(LDFLAGS)
else
$(BIN):	$(libraries) $(local_objs)
	$(LD) $(FFLAGS) $(LDFLAGS) $(local_objs) -o $@ $(libraries) $(LIBS)
endif

$(eval $(compile-rules))
