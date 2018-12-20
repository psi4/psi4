find_package(gau2grid 1.2 CONFIG QUIET COMPONENTS gaussian)

if(${gau2grid_FOUND})
    get_property(_loc TARGET gau2grid::gg PROPERTY LOCATION)
    message(STATUS "${Cyan}Found gau2grid${ColourReset}: ${_loc} (found version ${gau2grid_VERSION})")
    add_library(gau2grid_external INTERFACE)  # dummy
else()
    if(${CMAKE_INSIST_FIND_PACKAGE_gau2grid})
        message(FATAL_ERROR "Suitable gau2grid could not be externally located as user insists")
    endif()

    include(ExternalProject)
    message(STATUS "Suitable gau2grid could not be located, ${Magenta}Building gau2grid${ColourReset} instead.")
    ExternalProject_Add(gau2grid_external
        DEPENDS pybind11_external
        GIT_REPOSITORY https://github.com/dgasmith/gau2grid
        GIT_TAG v1.3.1
        UPDATE_COMMAND ""
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                   -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                   -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                   -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                    # OpenMP irrelevant
                   -DENABLE_XHOST=${ENABLE_XHOST}
                   -DBUILD_FPIC=${BUILD_FPIC}
                   -DMAX_AM=${MAX_AM_ERI}
                   -DCARTESIAN_ORDER=row
                   -DSPHERICAL_ORDER=gaussian
                   -DENABLE_GENERIC=${ENABLE_GENERIC}
                   -DPYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}
        CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS})

    set(gau2grid_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/gau2grid CACHE PATH "path to internally built gau2gridConfig.cmake" FORCE)
endif()
