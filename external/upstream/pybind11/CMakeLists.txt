
if(${Python_NumPy_VERSION} GREATER_EQUAL "2")
        find_package(pybind11 2.12.0 CONFIG QUIET)  # edit in codedeps
else()
    if(${Python_VERSION_MINOR} GREATER_EQUAL "11")
        find_package( pybind11 2.10.1 CONFIG QUIET)
    else()
        find_package( pybind11 2.6.2 CONFIG QUIET)
    endif()
endif()

if(${pybind11_FOUND})
    message(STATUS "${Cyan}Found pybind11${ColourReset}: ${pybind11_INCLUDE_DIR} (found version ${pybind11_VERSION})")
    add_library(pybind11_external INTERFACE)  # dummy

    if(${CMAKE_INSIST_FIND_PACKAGE_pybind11})
        message(VERBOSE "Suitable pybind11 located externally as user insists")
    endif()

else()
    if(${CMAKE_INSIST_FIND_PACKAGE_pybind11})
        message(FATAL_ERROR "Suitable pybind11 could not be externally located as user insists")
    endif()

    include(ExternalProject)
    message(STATUS "Suitable pybind11 could not be located, ${Magenta}Building pybind11${ColourReset} instead.")
    ExternalProject_Add(pybind11_external
        URL https://github.com/pybind/pybind11/archive/v2.13.1.tar.gz  # edit in codedeps
        UPDATE_COMMAND ""
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                   -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                   -DCMAKE_CXX_STANDARD=${psi4_CXX_STANDARD}
                   -DPYBIND11_TEST=OFF
                   -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                   -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
        CMAKE_CACHE_ARGS -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                         -DCMAKE_CXX_COMPILER_ARG1:STRING=${CMAKE_CXX_COMPILER_ARG1})

    set(pybind11_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/pybind11 CACHE PATH "path to internally built pybind11Config.cmake" FORCE)
endif()

