include(ExternalProject)

find_package(libderiv QUIET COMPONENTS ${LIBINT_OPT_AM})
list(APPEND CMAKE_PREFIX_PATH
     "${CMAKE_BINARY_DIR}/stage${CMAKE_INSTALL_PREFIX}/external")
if(NOT LIBDERIV_FOUND)
   message(STATUS 
      "Suitable Libderiv could not be located, building one instead")
   ExternalProject_Add(libderiv_external
      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libderiv-src
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/external/libderiv
                 -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                 -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
                 -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
                 -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                 -DCMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}
                 -DLIBINT_OPT_AM=${LIBINT_OPT_AM}
                 -DBUILD_FPIC=${BUILD_FPIC}
      CMAKE_CACHE_ARGS -DCMAKE_PREFIX_PATH:STRING=${CMAKE_PREFIX_PATH}
      INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install DESTDIR=${CMAKE_BINARY_DIR}/stage
      )
    add_dependencies(libderiv_external libint_external)
else()
    add_library(libderiv_external INTERFACE)
endif()
