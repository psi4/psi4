find_package(pybind11 CONFIG QUIET)

if(${pybind11_FOUND})
    message(STATUS "Found pybind11: ${pybind11_INCLUDE_DIR} (found version ${pybind11_VERSION})")
    add_library(pybind11_external INTERFACE)  # dummy
else()
    include(ExternalProject)
    message(STATUS "Suitable pybind11 could not be located, building one instead.")
    ExternalProject_Add(pybind11_external
        GIT_REPOSITORY https://github.com/psi4/pybind11
        GIT_TAG oldpb11full
        UPDATE_COMMAND ""
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/external/pybind11
                   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                   -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                   -DPYBIND11_CPP_STANDARD=-std=c++11
                   -DPYBIND11_TEST=OFF
                   -DPYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}
        CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                         -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
        INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
        DESTDIR=${CMAKE_BINARY_DIR}/stage)
endif()

