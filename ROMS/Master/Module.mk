# git $Id$
# svn $Id: Module.mk 1151 2023-02-09 03:08:53Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
# Copyright (c) 2002-2023 The ROMS/TOMS Group             Kate Hedstrom :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

local_sub  := Master

ifdef EXEC
  local_src  := $(wildcard $(local_sub)/*.F)
  local_objs := $(subst .F,.o,$(local_src))
  local_objs := $(addprefix $(SCRATCH_DIR)/, $(notdir $(local_objs)))

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
  local_objs := $(addprefix $(SCRATCH_DIR)/, $(notdir $(local_objs)))

  sources    += $(local_src)
endif

$(eval $(compile-rules))
