set(SIMINT_VER "0.7")

if(${ENABLE_simint})
    find_package(simint CONFIG QUIET COMPONENTS am${MAX_AM_ERI} der0)

    if(${simint_FOUND})
        get_property(_loc TARGET simint::simint PROPERTY LOCATION)
        message(STATUS "${Cyan}Found simint ${simint_MAXAM}${ColourReset}: ${_loc} (found version ${simint_VERSION}; vectorization ${simint_VECTOR})")
        add_library(simint_external INTERFACE)  # dummy
    else()
        include(ExternalProject)
        message(STATUS "Suitable simint could not be located, ${Magenta}Building simint${ColourReset} instead.")
        ExternalProject_Add(simint_external
            URL http://www.bennyp.org/research/simint/download/simint-v${SIMINT_VER}.tar.bz2
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                       # OpenMP relevant but always ON and handled by simint
                       # xHost handled through SIMINT_VECTOR. default sse equiv to XHOST=OFF
                       -DCMAKE_POSITION_INDEPENDENT_CODE=${BUILD_FPIC}
                       -DSIMINT_STANDALONE=${ENABLE_GENERIC}
                       -DENABLE_TESTS=OFF
                       -DSIMINT_VECTOR=${SIMINT_VECTOR}
                       -DSIMINT_MAXAM=${MAX_AM_ERI}
            CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                             -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
            INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
            DESTDIR=${CMAKE_BINARY_DIR}/stage)

        set(simint_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/simint CACHE PATH "path to internally built simintConfig.cmake" FORCE)
    endif()
else()
    add_library(simint_external INTERFACE)  # dummy
endif()
