# define switches
if( WIN32 )
  set( SWITCHES -dos -cvis )
elseif( UNIX )
  set( SWITCHES -unix )
endif()

# include user-dependent switches
if( MPI OR JAC )
  list( APPEND SWITCHES -impi -mpi )
endif()
#
if( JAC )
  set( MPI ON )
  list( APPEND SWITCHES -jac )
endif()
#
if( NETCDF )
  list( APPEND SWITCHES -netcdf )
endif()

# set switches
if( PERL_FOUND )
  execute_process( COMMAND ${PERL_EXECUTABLE} ${CMAKE_SOURCE_DIR}/switch.pl ${SWITCHES} ${SRC}/*.ftn ${SRC}/*.ftn90 )
else()
  message( FATAL_ERROR "-- Cannot build ${CMAKE_PROJECT_NAME} without Perl" )
endif()

# set compiler flags
if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU" )
  # GNU Fortran
  set( CMAKE_Fortran_FLAGS "-O -w -fno-second-underscore -ffree-line-length-none" CACHE STRING "" FORCE )
  if( ${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 10 )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch" )
  endif()
  if( WIN32 AND ${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 7 )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdec" )
  endif()
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
  # Intel Fortran
  if( UNIX )
    set( CMAKE_Fortran_FLAGS "-O2 -W0 -assume byterecl -traceback -diag-disable 8290 -diag-disable 8291 -diag-disable 8293" CACHE STRING "" FORCE )
  elseif( WIN32 )
    set( CMAKE_Fortran_FLAGS "/O2 /assume:byterecl /traceback /nowarn /nologo /Qdiag-disable:8290 /Qdiag-disable:8291 /Qdiag-disable:8293" CACHE STRING "" FORCE )
  endif()
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "PGI" )
  # Portland Group
  set( CMAKE_Fortran_FLAGS "-fast" CACHE STRING "" FORCE )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Fujitsu" )
  # Lahey Fortran
  set( CMAKE_Fortran_FLAGS "-O1 -staticlink -nwo" CACHE STRING "" FORCE )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "XL" )
  # IBM XL Fortran
  set( CMAKE_Fortran_FLAGS "-O3 -qstrict -qarch=auto -qtune=auto -qcache=auto -qalign=4k -w" CACHE STRING "" FORCE )
else()
  message( FATAL_ERROR "-- Current Fortran compiler ${CMAKE_Fortran_COMPILER} not supported" )
endif()

# COAWST include mct params info
set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I$ENV{MCT_PARAMS_DIR} -I$ENV{MCT_INCDIR}" )

# choose parallelisation paradigm: OpenMP or MPI
if( OPENMP AND MPI )
  message( FATAL_ERROR "-- Cannot use both OpenMP and MPI" )
elseif( OPENMP )
  find_package( OpenMP REQUIRED COMPONENTS Fortran )
  if( OpenMP_Fortran_FOUND )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}" )
  else()
    message( FATAL_ERROR "-- Fortran compiler does not support OpenMP" )
  endif()
  # turn off MPI
  unset( MPI_Fortran_FOUND CACHE )
  unset( MPI_Fortran_COMPILER CACHE )
  unset( MPI_Fortran_LIBRARIES CACHE )
elseif( MPI )
  find_package( MPI REQUIRED COMPONENTS Fortran )
  # turn off OpenMP
  unset( OpenMP_Fortran_FLAGS CACHE )
else()
  # turn off both OpenMP and MPI
  unset( OpenMP_Fortran_FLAGS CACHE )
  unset( MPI_Fortran_FOUND CACHE )
  unset( MPI_Fortran_COMPILER CACHE )
  unset( MPI_Fortran_LIBRARIES CACHE )
endif()

# build library
if( NOT TARGET swan${VERSION} )
  if( NOT NETCDF )
    include( ${SRC}/srclist.cmake )
    add_library( swan${VERSION} STATIC ${swan_src} )
    set_target_properties( swan${VERSION} PROPERTIES LINKER_LANGUAGE Fortran )
  else()
    include( ${SRC}/srclistnc.cmake )
    add_library( swan${VERSION} STATIC ${swan_nc_src} )
    set_target_properties( swan${VERSION} PROPERTIES LINKER_LANGUAGE Fortran )
    target_include_directories( swan${VERSION} PUBLIC ${NetCDF_INCLUDE_DIRS} )
#   target_link_libraries( swan${VERSION} PUBLIC NetCDF::NetCDF_Fortran )
    target_link_libraries( swan${VERSION} PUBLIC NetCDF::NetCDF_Fortran  $ENV{MCT_PARAMS_DIR}/mct_coupler_params.o -L$ENV{MCT_LIBDIR} -lmct -lmpeu)
  endif()
endif()

# build executable
if( NOT TARGET ${EXE} )
# set( swanmain ${SRC}/swanmain.${EXTF} )
# add_executable( ${EXE} ${swanmain} )
  set( swanmain2 ${SRC}/swanmain2.${EXTF} )
  add_executable( ${EXE} ${swanmain2} )
# target_link_libraries( ${EXE} swan${VERSION} ${MPI_Fortran_LIBRARIES} )
  target_link_libraries( ${EXE} swan${VERSION} -L$ENV{MCT_LIBDIR} -lmct -lmpeu ${MPI_Fortran_LIBRARIES} )
endif()

# install the package
install( TARGETS ${EXE} swan${VERSION} RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib )
install( DIRECTORY ${CMAKE_BINARY_DIR}/mod/ DESTINATION mod )
install( DIRECTORY ${CMAKE_SOURCE_DIR}/doc/ DESTINATION doc )
install( DIRECTORY ${CMAKE_SOURCE_DIR}/misc/ DESTINATION misc )
install( DIRECTORY ${CMAKE_SOURCE_DIR}/tools/ DESTINATION tools )
if( WIN32 )
  install( PROGRAMS ${CMAKE_SOURCE_DIR}/bin/swanrun.bat DESTINATION bin )
elseif( UNIX )
  install( PROGRAMS ${CMAKE_SOURCE_DIR}/bin/swanrun DESTINATION bin )
endif()
