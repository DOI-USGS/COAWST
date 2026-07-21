# delete build files and directories

# find all files that we want to remove
file( GLOB_RECURSE F90FILES "${CMAKE_SOURCE_DIR}/src/*.f" "${CMAKE_SOURCE_DIR}/src/*.for" "${CMAKE_SOURCE_DIR}/src/*.f90" "${CMAKE_SOURCE_DIR}/src/hcat/*.f" "${CMAKE_SOURCE_DIR}/src/hcat/*.for" )

# place these files into a list
set( DEL ${F90FILES} )

# include build folder to be remove
file( GLOB BUILD "${CMAKE_SOURCE_DIR}/build" )
set( DEL ${DEL} ${BUILD} )

# loop over list items and delete each one
foreach( D ${DEL} )
  if( EXISTS ${D} )
    file( REMOVE_RECURSE ${D} )
  endif()
endforeach()
