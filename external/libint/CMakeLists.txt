include(ExternalProject)

find_package(libint QUIET COMPONENTS ${LIBINT_OPT_AM})

if(NOT LIBINT_FOUND)
   message(STATUS 
      "Suitable Libint could not be located, building one instead")
   list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
   message(STATUS "${CMAKE_MODULE_PATH}")
   ExternalProject_Add(libint_external
      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libint-src
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/external/libint
                 -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                 -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
                 -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
                 -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                 -DCMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}
                 -DLIBINT_OPT_AM=${LIBINT_OPT_AM}
     CMAKE_CACHE_ARGS -DCMAKE_MODULE_PATH:STRING=${CMAKE_MODULE_PATH}
     INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install DESTDIR=${CMAKE_BINARY_DIR}/stage
      )
else()
    add_library(libint_external INTERFACE)
endif()


