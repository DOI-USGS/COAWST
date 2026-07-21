macro(set_fire_defaults)
  if("${DM_PARALLEL}" STREQUAL "")
    set(DM_PARALLEL ON)
  endif()
  if("${NUOPC}" STREQUAL "")
    set(NUOPC OFF)
  endif()
  if("${ESMX}" STREQUAL "")
    set(ESMX OFF)
  endif()
endmacro()

function(set_fire_cache)
  set (DM_PARALLEL ${DM_PARALLEL} CACHE BOOL "Turn on MPI build" FORCE)
  set (NUOPC ${NUOPC} CACHE BOOL "Build NUOPC cap" FORCE)
  set (ESMX ${ESMX} CACHE BOOL "Build ESMX application" FORCE)
  set (OPENMP ${OPENMP} CACHE BOOL "Enable OpenMP parallelization" FORCE)
endfunction()

function(print_fire_settings)
  message (STATUS "FIRE BEHAVIOR BUILD SETTINGS\n"
                  "\tOPENMP: ${OPENMP}\n"
                  "\tDM_PARALLEL: ${DM_PARALLEL}\n"
                  "\tNUOPC: ${NUOPC}\n"
		  "\tESMX:  ${ESMX}")
endfunction()
