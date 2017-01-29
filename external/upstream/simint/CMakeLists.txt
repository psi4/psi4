if(${ENABLE_simint})
   include(ExternalProject)

   set(SIMINT_VER "0.6")
   string(TOLOWER ${SIMINT_VECTOR} SIMINT_VECTOR)
   list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

   message(STATUS "Using Simint v${SIMINT_VER} (vectorization = ${SIMINT_VECTOR}) for integrals")

   ExternalProject_Add(simint_external
      URL http://www.bennyp.org/research/simint/download/simint-v${SIMINT_VER}.tar.bz2
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                 -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                 -DSIMINT_C_FLAGS=-fPIC
                 -DSIMINT_VECTOR=${SIMINT_VECTOR}
                 -DSIMINT_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install DESTDIR=${CMAKE_BINARY_DIR}/stage
   )

   set(simint_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/simint CACHE PATH "path to internally built simintConfig.cmake" FORCE)

else()
    add_library(simint_external INTERFACE)
endif()


