# ======================================================================
# Noah-MP unified Makefile.in (final version)
# Automatic Fortran dependency tracking + correct build order
# ======================================================================

# --- Root paths ---
NOAHMP_HOME     ?= ../..
NOAHMP_OBJ_DIR  ?= $(NOAHMP_HOME)/drivers/erf/obj
NOAHMP_LIB_DIR  ?= $(NOAHMP_HOME)/drivers/erf/lib
NOAHMP_INC_DIR  ?= $(NOAHMP_HOME)/drivers/erf/include
NOAHMP_LIB      ?= $(NOAHMP_LIB_DIR)/libnoahmp.a

# --- Compilers ---
MPIFF    ?= mpif90
CXX      ?= mpicxx
AR       ?= ar

# --- Flags ---
FFLAGS_NOAHMP   += -O3 -cpp -fPIC -ffree-line-length-none -J$(abspath $(NOAHMP_INC_DIR)) -I$(abspath $(NOAHMP_INC_DIR))
CXXFLAGS_NOAHMP += -O3 -fPIC -I$(abspath $(NOAHMP_INC_DIR)) -I$(abspath $(NOAHMP_HOME)/drivers/erf)

# --- NetCDF detection ---
ifneq (, $(shell which nf-config 2>/dev/null))
    NETCDF_INCFLAGS := $(shell nf-config --fflags)
    NETCDF_LIBS     := $(shell nf-config --flibs)
else
    NETCDF_INCFLAGS ?= $(shell pkg-config --cflags netcdf-fortran 2>/dev/null)
    NETCDF_LIBS     ?= $(shell pkg-config --libs netcdf-fortran 2>/dev/null)
endif

FFLAGS_NOAHMP += $(NETCDF_INCFLAGS)

# --- Ensure directories exist ---
$(shell mkdir -p $(NOAHMP_OBJ_DIR) $(NOAHMP_LIB_DIR) $(NOAHMP_INC_DIR))

# ======================================================================
# Source discovery
# ======================================================================
SRC_DIRS := $(NOAHMP_HOME)/src $(NOAHMP_HOME)/utility $(NOAHMP_HOME)/drivers/erf

# Map lowercase module names to actual source files (case-insensitive match)
define find_f90_source
$(firstword $(wildcard $(foreach d,$(SRC_DIRS),$(d)/$(1).F90 $(d)/$(1).f90 $(d)/$(shell echo $(1) | tr a-z A-Z).F90)))
endef

F90_SRCS := $(foreach d,$(SRC_DIRS),$(wildcard $(d)/*.F90))
CPP_SRCS := $(foreach d,$(SRC_DIRS),$(wildcard $(d)/*.cpp))

OBJS := $(addprefix $(NOAHMP_OBJ_DIR)/, $(notdir $(F90_SRCS:.F90=.o))) \
        $(addprefix $(NOAHMP_OBJ_DIR)/, $(notdir $(CPP_SRCS:.cpp=.o)))

DEPS := $(OBJS:.o=.d)

# ======================================================================
# Primary targets
# ======================================================================
.PHONY: all_noahmp clean_noahmp noahmp_includes

all_noahmp: $(NOAHMP_LIB) noahmp_includes

$(NOAHMP_LIB): $(OBJS)
	@mkdir -p $(NOAHMP_LIB_DIR)
	@echo "Archiving Noah-MP library: $@"
	$(AR) rcs $@ $^

# ======================================================================
# Automatic Fortran dependency generation (correct + path-safe)
# ======================================================================
# Each .d file will contain lines like:
#   obj/foo.o: obj/bar.o obj/baz.o
define make-dep
@echo "Generating dependencies for $<"
@outfile=$(NOAHMP_OBJ_DIR)/$*.d; \
mkdir -p $(dir $$outfile); \
mods=$$( $(MPIFF) $(FFLAGS) $(FFLAGS_NOAHMP) -cpp -E $< 2>/dev/null | \
sed 's/!.*//' | \
grep -E '(^|[^A-Za-z0-9_])use[[:space:]]+[A-Za-z0-9_]+' | \
grep -vE 'netcdf|mpi|iso_fortran_env|omp_lib|iso_c_binding|intrinsic' | \
sed -E 's|.*use[[:space:]]+([A-Za-z0-9_]+).*|\1|g' | sort -u ); \
printf "$(NOAHMP_OBJ_DIR)/$*.o:" > $$outfile; \
for m in $$mods; do printf " $(NOAHMP_OBJ_DIR)/$${m}.o"; done >> $$outfile; \
echo "" >> $$outfile; \
test -s $$outfile || echo "$(NOAHMP_OBJ_DIR)/$*.o:" > $$outfile
endef


# Generate .d files from any Fortran source
$(NOAHMP_OBJ_DIR)/%.d: $(NOAHMP_HOME)/src/%.F90 ; $(make-dep)
$(NOAHMP_OBJ_DIR)/%.d: $(NOAHMP_HOME)/utility/%.F90 ; $(make-dep)
$(NOAHMP_OBJ_DIR)/%.d: $(NOAHMP_HOME)/drivers/erf/%.F90 ; $(make-dep)

# Include generated dependency files
ifeq ($(MAKECMDGOALS),all_noahmp)
  -include $(DEPS)
endif

# ======================================================================
# Compilation rules
# ======================================================================

# --- Fortran src modules ---
$(NOAHMP_OBJ_DIR)/%.o: $(NOAHMP_HOME)/src/%.F90 $(NOAHMP_OBJ_DIR)/%.d
	@echo "Compiling $(NOAHMP_HOME)/src: $<"
	$(MPIFF) $(FFLAGS) $(FFLAGS_NOAHMP) -c $< -o $@

# --- Fortran utility modules ---
$(NOAHMP_OBJ_DIR)/%.o: $(NOAHMP_HOME)/utility/%.F90 $(NOAHMP_OBJ_DIR)/%.d
	@echo "Compiling $(NOAHMP_HOME)/util: $<"
	$(MPIFF) $(FFLAGS) $(FFLAGS_NOAHMP) -c $< -o $@

# --- Fortran driver modules ---
$(NOAHMP_OBJ_DIR)/%.o: $(NOAHMP_HOME)/drivers/erf/%.F90 $(NOAHMP_OBJ_DIR)/%.d
	@echo "Compiling $(NOAHMP_HOME)/driver: $<"
	$(MPIFF) $(FFLAGS) $(FFLAGS_NOAHMP) -c $< -o $@

# --- C++ sources ---
$(NOAHMP_OBJ_DIR)/%.o: $(NOAHMP_HOME)/drivers/erf/%.cpp
	@echo "Compiling Noah-MP C++: $<"
	$(CXX) $(CXXFLAGS) $(CXXFLAGS_NOAHMP) -c $< -o $@

$(NOAHMP_OBJ_DIR)/%.o: $(NOAHMP_HOME)/src/%.cpp
	@echo "Compiling Noah-MP C++: $<"
	$(CXX) $(CXXFLAGS) $(CXXFLAGS_NOAHMP) -c $< -o $@

# ======================================================================
# Install headers and modules
# ======================================================================
noahmp_includes:
	@mkdir -p $(NOAHMP_INC_DIR)
	@echo "Installing module (.mod) and header (.H) files: $(NOAHMP_INC_DIR)"
	@find $(NOAHMP_OBJ_DIR) -maxdepth 1 -type f -name "*.mod" -exec cp -u {} $(NOAHMP_INC_DIR)/ \;
	@find $(NOAHMP_HOME)/drivers/erf -maxdepth 1 -type f -name "*.H" -exec cp -u {} $(NOAHMP_INC_DIR)/ \;

# ======================================================================
# Cleaning
# ======================================================================
clean_noahmp:
	@echo "Cleaning Noah-MP build artifacts..."
	rm -rf $(NOAHMP_OBJ_DIR) $(NOAHMP_LIB_DIR) $(NOAHMP_INC_DIR)
