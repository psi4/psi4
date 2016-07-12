include(ExternalProject)

find_package(PCMSolver QUIET)
if(NOT PCMSolver_FOUND AND ENABLE_PCMSOLVER)
  message(STATUS "PCMSolver not found. The pre-packaged version will be built.")
  ExternalProject_Add(pcmsolver_external
       GIT_REPOSITORY https://github.com/PCMSolver/pcmsolver
       GIT_TAG v1.1.1
       CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                  -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external/pcmsolver
                  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                  -DEXTRA_C_FLAGS=${PCM_EXTRA_C_FLAGS}
                  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                  -DEXTRA_CXX_FLAGS=${PCM_EXTRA_CXX_FLAGS}
                  -DENABLE_CXX11_SUPPORT=ON
                  -DBOOST_INCLUDEDIR=${Boost_INCLUDE_DIRS}
                  -DBOOST_LIBRARYDIR=${Boost_LIBRARY_DIRS}
                  -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
                  -DENABLE_TESTS=OFF
                  -DENABLE_LOGGER=OFF
                  -DENABLE_TIMER=OFF
                  -DBUILD_STANDALONE=OFF
                  -DENABLE_FORTRAN_API=OFF
                  -DSTATIC_LIBRARY_ONLY=ON
                  -DENABLE_GENERIC=${ENABLE_STATIC_LINKING}
                  -DPYTHON_INTERPRETER=${PYTHON_EXECUTABLE}
       INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install 
                      DESTDIR=${CMAKE_BINARY_DIR}/stage
    )
    add_dependencies(pcmsolver_external boost_external)
else()
   add_library(pcmsolver_external INTERFACE)
endif()