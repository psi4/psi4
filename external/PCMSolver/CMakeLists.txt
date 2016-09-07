include(ExternalProject)

find_package(PCMSolver QUIET)
if(NOT PCMSolver_FOUND AND ${ENABLE_PCMSOLVER})
  message(STATUS "PCMSolver not found. The pre-packaged version will be built.")
  set(ENABLE_STATIC (NOT ${BUILD_SHARED_LIBS}))

  ExternalProject_Add(pcmsolver_external
       GIT_REPOSITORY https://github.com/PCMSolver/pcmsolver
       GIT_TAG v1.1.4
       UPDATE_COMMAND ""
       CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                  -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/external/pcmsolver
                  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                  -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
                  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                  -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
                  -DENABLE_CXX11_SUPPORT=ON
                  -DBOOST_ROOT=${BOOST_ROOT}
                  -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
                  -DENABLE_TESTS=OFF
                  -DENABLE_LOGGER=OFF
                  -DENABLE_TIMER=OFF
                  -DBUILD_STANDALONE=OFF
                  -DENABLE_FORTRAN_API=OFF
                  -DSTATIC_LIBRARY_ONLY=${ENABLE_STATIC}
                  -DENABLE_GENERIC=${ENABLE_STATIC}
                  -DPYTHON_INTERPRETER=${PYTHON_EXECUTABLE}
       INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
                      DESTDIR=${CMAKE_BINARY_DIR}/stage
    )
    configure_file(${CMAKE_SOURCE_DIR}/psi4/share/python/pcm_placeholder.py.in
       ${CMAKE_SOURCE_DIR}/share/python/pcm_placeholder.py
    )
    add_dependencies(pcmsolver_external boost_external)
else()
   add_library(pcmsolver_external INTERFACE)
endif()
