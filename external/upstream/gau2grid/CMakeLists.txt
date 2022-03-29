find_package(gau2grid 2.0 CONFIG QUIET)

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
        URL https://github.com/dgasmith/gau2grid/archive/v2.0.7.tar.gz
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
                   -DENABLE_GENERIC=${ENABLE_GENERIC}
                   -DPYTHON_EXECUTABLE=${Python_EXECUTABLE}
                   -DPYTHON_INCLUDE_DIR=${Python_INCLUDE_DIRS}
                   -DPYTHON_LIBRARY=${Python_LIBRARIES}
                   -DPython_EXECUTABLE=${Python_EXECUTABLE}
                   -DPython_INCLUDE_DIR=${Python_INCLUDE_DIRS}
                   -DPython_LIBRARY=${Python_LIBRARIES}
                   -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                   -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
        CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS})

    set(gau2grid_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/gau2grid CACHE PATH "path to internally built gau2gridConfig.cmake" FORCE)
endif()
