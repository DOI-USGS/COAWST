# check MPI status
if( NOT MPI )
  return()
endif()

# define switches
if( WIN32 )
  set( HSWITCH -dos -cvis )
elseif( UNIX )
  set( HSWITCH -unix )
endif()

# set switches
if( PERL_FOUND )
  execute_process( COMMAND ${PERL_EXECUTABLE} ${CMAKE_SOURCE_DIR}/switch.pl ${HSWITCH} ${SRC}/hcat/swanhcat.ftn )
else()
  message( FATAL_ERROR "-- Cannot build hcat without Perl" )
endif()

# build and install the executable
set( swanhcat ${SRC}/hcat/swanhcat.${EXTF} )
add_executable( ${HEXE} ${swanhcat} )
install( TARGETS ${HEXE} RUNTIME DESTINATION bin )
