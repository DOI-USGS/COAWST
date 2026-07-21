cmake_minimum_required(VERSION 3.20)
enable_testing()

add_test(NAME test7
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test7.s
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(NAME test8
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test8.s
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

if(NUOPC)
  add_test(NAME test7esmf
           COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test7esmf.s
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
  add_test(NAME test8esmf
           COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test8esmf.s
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
endif()

if(ESMX)
  add_test(NAME test7esmx
           COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test7esmx.s
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
  add_test(NAME test8esmx
           COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test8esmx.s
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
  add_test(NAME testx
           COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/testx.s
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
endif()
