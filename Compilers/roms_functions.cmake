# git $Id$
# svn $Id: roms_functions.cmake 1151 2023-02-09 03:08:53Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Functions used in CMake to overcome its rudimentary capabilities.

find_program(
  GREP
  NAMES grep egrep
  DOC "Grep command"
)

find_program(
  PERL
  NAMES perl
  DOC "Perl command"
)

find_program(
  CPP_CLEAN
  NAMES cpp_clean
  HINTS "${CMAKE_SOURCE_DIR}/ROMS/Bin/"
  DOC "ROMS CPP Clean command"
)

find_program(
  CPP_EXECUTABLE
  NAMES cpp
  DOC "C-preprocessor command"
)

###########################################################################
# The "preprocess_fortran" function teaches CMake how to use CPP to
# create .f90 files from ROMS .F source files. It is imperative that
# we generate processed .f90 files BEFORE CMake attempts to determine
# Fortran dependencies because its Fortran dependency tracker is
# inconsistent and easily confused by C-preprocessor nested
# if-directives.
###########################################################################

# First, make the directory for the .f90 files.

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/f90")

function(preprocess_fortran)

   # Determine include directories.

   get_directory_property(incdirs INCLUDE_DIRECTORIES)
   set(incflags)
   foreach(i ${incdirs})
     list(APPEND incflags "-I${i}")
   endforeach()

   # Set-up "defines" (not configuration specific).

   get_directory_property(defines COMPILE_DEFINITIONS)
   set(defflags)
   foreach(d ${defines})
     list(APPEND defflags "-D${d}")
   endforeach()

   # Loop over all arguments.

   set(fsrcs)
   foreach(f ${ARGV})

     # Is it a Fortran file?

     if(f MATCHES "\\.[Ff](9[05])?")
       message(STATUS "Got fortran file: ${f}")

       # Get special definitions for file, if any.

       set(fdefflags)
       if(${f} MATCHES ".*mod_strings.*")
         get_property(fdefines
                      SOURCE ${f}
                      PROPERTY COMPILE_DEFINITIONS
         )
         foreach(fd ${fdefines})
           list(APPEND fdefflags "-D${fd}")
         endforeach()
       endif()

       # Construct output filename.

       if(NOT IS_ABSOLUTE "${f}")
         get_filename_component(f "${f}" ABSOLUTE)
       endif()
       file(RELATIVE_PATH r "${CMAKE_CURRENT_SOURCE_DIR}" "${f}")
       get_filename_component(e "${r}" EXT)
       get_filename_component(n "${r}" NAME_WE)
       get_filename_component(p "${r}" PATH)
       set(of "${CMAKE_CURRENT_BINARY_DIR}/f90/${n}.f90")
       message(STATUS "Output name: ${of}")

       # Preprocess the thing.

       add_custom_command(OUTPUT "${of}"
         COMMAND ${CPP_EXECUTABLE} ${CPPFLAGS} ${defflags} ${fdefflags} ${incflags} "${f}" > "${of}"
         COMMAND ${CPP_CLEAN} "${of}"
         IMPLICIT_DEPENDS Fortran "${f}"
         COMMENT "Preprocessing ${f}"
         COMMAND_EXPAND_LISTS
         VERBATIM
       )
       list(APPEND fsrcs "${of}")
     else()
       list(APPEND fsrcs "${f}")
     endif()
   endforeach()

   # Return the (preprocessed) sources.

   set(f90srcs "${fsrcs}" PARENT_SCOPE)
endfunction()

###########################################################################
# The "get_options" function is used by roms_config.cmake to determine if
# the adjoint, tangent linear, and/or representer model(s) are needed
# and whether to add ARPACK/PARPACK and/or parallel I/O using SCORPIO
# to link into the executable.
###########################################################################

function(get_options roms_header)

  # Set ROMS_HEADER C-Preprocessor flag.

  set( APP_HEADER -DROMS_HEADER="${roms_header}" )

  # If MY_CPP_FLAGS was set and passed, cycle through the list and add the -D.

  set(defflags)
  foreach(d ${ARGN})
    list(APPEND defflags "-D${d}")
  endforeach()

  execute_process(
    COMMAND ${CPP_EXECUTABLE} ${CPPFLAGS} ${APP_HEADER} ${defflags} ${CMAKE_CURRENT_SOURCE_DIR}/Compilers/defs_cmake.h
    COMMAND ${GREP} -v ^[[:space:]]*\$
    COMMAND_ECHO STDOUT
    RESULT_VARIABLE status
    ERROR_VARIABLE err
    OUTPUT_VARIABLE models_libs
  )

  # Check whether the cpp command above produced an error and stop CMake if it did.

  if( status AND NOT status EQUAL 0 )
    message( FATAL_ERROR "get_options error: ${err}" )
  endif()

  # Return needed Models/Libraries.

  set(defs "${models_libs}" PARENT_SCOPE)
endfunction()

###########################################################################
# Determine how to link to NetCDF
###########################################################################

function(link_netcdf)
  if( DEFINED ENV{NF_CONFIG} )
    Message( STATUS "Using $ENV{NF_CONFIG}" )

    # Get link line to break down below

    execute_process( COMMAND $ENV{NF_CONFIG} --flibs
                     OUTPUT_VARIABLE tmp
    )

    if( DEFINED ENV{NETCDF_INCDIR} )
      set( idir "$ENV{NETCDF_INCDIR}" )                          # Set include directory
    else()

      # Retrieve include direcotry

      execute_process( COMMAND $ENV{NF_CONFIG} --includedir
                       OUTPUT_VARIABLE idir
      )
    endif()

    string( STRIP "${tmp}" linkline )                            # Strip trailing whitespace
    string( REGEX MATCHALL "-L[^ \t]*" ldirs ${linkline} )       # Create list of dirs
    string( REGEX MATCHALL "-l[^ \t]*" libs ${linkline} )        # Create list of libs

  elseif( DEFINED ENV{NETCDF_LIBS} )
    string( REGEX MATCHALL "-L[^ \t]*" ldirs $ENV{NETCDF_LIBS} ) # Create list of dirs
    string( REGEX MATCHALL "-l[^ \t]*" libs $ENV{NETCDF_LIBS} )  # Create list of libs
    if( DEFINED ENV{NETCDF_INCDIR} )
      set( idir  "$ENV{NETCDF_INCDIR}" )                         # Set include directory
    elseif( DEFINED ENV{NETCDF} )
      set( idir  "$ENV{NETCDF}/include" )                        # Construct include path
    else()
      Message( FATAL_ERROR "NetCDF includes not found!" )
    endif()
  else()
    if( DEFINED ENV{NETCDF_INCDIR} AND DEFINED ENV{NETCDF_LIBDIR} )
      set( ldirs "$ENV{NETCDF_LIBDIR}" )                         # Set lib directory
      set( libs  "netcdf" )                                      # Set NetCDF3 lib name
      set( idir  "$ENV{NETCDF_INCDIR}" )                         # Set include directory
    else()
      Message( FATAL_ERROR "No NetCDF found!" )
    endif()
  endif()

  list( TRANSFORM libs  REPLACE "-l([^ \t]*)" "\\1" )           # remove "-l" from libs
  list( TRANSFORM ldirs REPLACE "-L([^ \t]*)" "\\1" )           # remove "-L" from dirs

  if( "${libs}" STREQUAL "" )
    message( FATAL_ERROR "NetCDF library name not found" )
  endif()

  set( netcdf_ldirs "${ldirs}" PARENT_SCOPE )
  set( netcdf_libs "${libs}" PARENT_SCOPE )
  set( netcdf_idir "${idir}" PARENT_SCOPE )

  Message( STATUS " netcdf_libs = ${libs}" )
  Message( STATUS "netcdf_ldirs = ${ldirs}" )
  Message( STATUS " netcdf_idir = ${idir}" )
endfunction()
