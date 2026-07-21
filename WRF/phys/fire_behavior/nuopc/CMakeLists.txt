# add ESMF targets
if (NOT TARGET ESMF)
  find_package(ESMF MODULE REQUIRED)
  if (NOT ESMF_FOUND)
    message(FATAL_ERROR "ESMF library not found. Please set environment variable ESMFMKFILE.")
  endif (NOT ESMF_FOUND)
endif (NOT TARGET ESMF)

# configure fire behavior nuopc cap library
add_library(fire_behavior_nuopc fire_behavior_nuopc.F90)
add_dependencies(fire_behavior_nuopc firelib)
target_link_libraries(fire_behavior_nuopc PUBLIC firelib)
target_link_libraries(fire_behavior_nuopc PUBLIC ESMF)
target_include_directories(fire_behavior_nuopc INTERFACE $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
                                                         $<INSTALL_INTERFACE:mod>)
# configure wrfdata nuopc cap library
add_library(fire_behavior_wrfnuopc wrf_nuopc.F90)
add_dependencies(fire_behavior_wrfnuopc firelib)
target_link_libraries(fire_behavior_wrfnuopc PUBLIC firelib)
target_link_libraries(fire_behavior_wrfnuopc PUBLIC ESMF)

# configure fire behavior esmf executable
add_executable(fire_behavior_esmf esmApp.F90 esm.F90)
add_dependencies(fire_behavior_esmf firelib fire_behavior_nuopc fire_behavior_wrfnuopc)
target_link_libraries(fire_behavior_esmf PUBLIC firelib fire_behavior_nuopc fire_behavior_wrfnuopc)
target_link_libraries(fire_behavior_esmf PUBLIC ESMF)

# configure fire behavior esmx executable
if (ESMX)
  if(ESMF_VERSION VERSION_LESS "8.5.0")
    message(FATAL_ERROR "ESMX requires ESMF version 8.5.0 or greater")
  endif (ESMF_VERSION VERSION_LESS "8.5.0")
  include(ExternalProject)
  ExternalProject_Add(esmx
    BINARY_DIR ${CMAKE_BINARY_DIR}/esmx
    SOURCE_DIR ${ESMF_ESMXDIR}
    CMAKE_ARGS -DESMX_BUILD_FILE=${CMAKE_CURRENT_SOURCE_DIR}/esmxBuild.yaml
               -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
               -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    BUILD_ALWAYS TRUE
    DEPENDS fire_behavior_nuopc fire_behavior_wrfnuopc
    EXCLUDE_FROM_ALL TRUE
  )
endif (ESMX)

# install
install(TARGETS firelib fire_behavior_nuopc fire_behavior_wrfnuopc
        EXPORT fire_behavior
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	INCLUDES DESTINATION mod)
install(EXPORT fire_behavior DESTINATION cmake)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION ${CMAKE_INSTALL_PREFIX})
install(TARGETS fire_behavior_esmf DESTINATION bin)
