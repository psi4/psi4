find_package(pybind11 CONFIG QUIET)

if(${pybind11_FOUND})
    get_property(_loc TARGET pybind11::pybind11 PROPERTY LOCATION)
    message(STATUS "Suitable pybind11 located at ${_loc}")
    add_library(pybind11_external INTERFACE)  # dummy
else()
    include(ExternalProject)
    message(STATUS "Suitable pybind11 could not be located, building one instead.")
    ExternalProject_Add(pybind11_external
        GIT_REPOSITORY https://github.com/psi4/pybind11
        #GIT_REPOSITORY https://github.com/pybind/pybind11
        #GIT_TAG v1.8.1
        GIT_TAG psipy11
        UPDATE_COMMAND ""
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                   -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                   -DPYBIND11_CPP_STANDARD=-std=c++11
                   -DPYBIND11_TEST=OFF
                   -DPYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}
                   -DPYTHON_INCLUDE_DIR=${PYTHON_INCLUDE_DIR}
                   -DPYTHON_LIBRARY=${PYTHON_LIBRARY}
        CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                         -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
        INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
        DESTDIR=${CMAKE_BINARY_DIR}/stage)
endif()

