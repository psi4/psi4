set(SIMINT_VER "0.7")

if(${ENABLE_simint})
    find_package(simint 0.7 COMPONENTS am${MAX_AM_ERI} der0 CONFIG QUIET)  # edit in codedeps

    if(${simint_FOUND})
        get_property(_loc TARGET simint::simint PROPERTY LOCATION)
        message(STATUS "${Cyan}Found simint ${simint_MAXAM}${ColourReset}: ${_loc} (found version ${simint_VERSION}; vectorization ${simint_VECTOR})")
        add_library(simint_external INTERFACE)  # dummy

        if(${CMAKE_INSIST_FIND_PACKAGE_simint})
            message(VERBOSE "Suitable simint located externally as user insists")
        endif()

    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_simint})
            message(FATAL_ERROR "Suitable simint could not be externally located as user insists")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable simint could not be located, ${Magenta}Building simint${ColourReset} instead.")

        if (CMAKE_CXX_COMPILER_ID MATCHES Intel)
            if (NOT (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "18.0.0"))
                message(FATAL_ERROR "ICPC version 2018 builds an incorrect simint library. Use the conda package or build simint separately with GNU or Intel 2017 compilers.")
            endif()
        endif()

        ExternalProject_Add(simint_external
            URL https://www.bennyp.org/research/simint/download/simint-v${SIMINT_VER}.tar.bz2  # edit in codedeps
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
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
                       -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                       -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
                       -DCMAKE_POLICY_VERSION_MINIMUM=3.5
            CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                             -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                             -DCMAKE_C_COMPILER_ARG1:STRING=${CMAKE_C_COMPILER_ARG1}
                             -DCMAKE_CXX_COMPILER_ARG1:STRING=${CMAKE_CXX_COMPILER_ARG1})

        set(simint_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/simint CACHE PATH "path to internally built simintConfig.cmake" FORCE)
    endif()
else()
    add_library(simint_external INTERFACE)  # dummy
endif()
