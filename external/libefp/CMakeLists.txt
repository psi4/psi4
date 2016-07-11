find_package(libefp QUIET)

if(NOT LIBEFP_FOUND)
   message(STATUS "Suitable libefp could not be located, building one instead")
   ExternalProject_Add(libefp_external
      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libefp-src
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/external/libefp
                 -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                 -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
                 -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                 -DCMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}
                 -DBUILD_FPIC=${BUILD_FPIC}
      INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install DESTDIR=${CMAKE_BINARY_DIR}/stage
   )
else()
   add_library(libefp_external INTERFACE)
endif()
