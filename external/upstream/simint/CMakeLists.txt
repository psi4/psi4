if(NOT(${ENABLE_SIMINT} STREQUAL "None"))
   include(ExternalProject)

   set(SIMINT_VER "0.5a")
   string(TOLOWER ${ENABLE_SIMINT} SIMINT_TYPE)

   message(STATUS "Using SIMINT ${SIMINT_TYPE} v${SIMINT_VER} for integrals")

   list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

   ExternalProject_Add(simint_external
      URL http://www.bennyp.org/research/simint/download/simint-v${SIMINT_VER}.tar.bz2
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/external
                 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                 -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                 -DSIMINT_C_FLAGS=-fPIC
      INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
      DESTDIR=${CMAKE_BINARY_DIR}/stage
   )
else()
    add_library(simint_external INTERFACE)
endif()


