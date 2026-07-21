# git $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
# Copyright (c) 2002-2026 The ROMS Group                  Kate Hedstrom :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

local_sub  := Master

ifdef EXEC
  local_src  := $(wildcard $(local_sub)/*.F)
  local_objs := $(subst .F,.o,$(local_src))
  local_objs := $(addprefix $(BUILD_DIR)/, $(notdir $(local_objs)))

  sources    += $(local_src)

  ifdef LD_WINDOWS
  $(BIN):	$(libraries) $(local_objs)
		$(LD) $(FFLAGS) $(local_objs) -o $@ $(ROMS_LIB) $(LIBS_WIN32) $(LDFLAGS)
  else
  $(BIN):	$(libraries) $(local_objs)
		$(LD) $(FFLAGS) $(LDFLAGS) $(local_objs) -o $@ $(ROMS_LIB) $(LIBS)
  endif
else
  local_src  := $(local_sub)/roms_kernel.F
  local_objs := $(subst .F,.o,$(local_src))
  local_objs := $(addprefix $(BUILD_DIR)/, $(notdir $(local_objs)))

  sources    += $(local_src)
endif

$(eval $(compile-rules))
