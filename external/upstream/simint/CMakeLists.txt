#if(${ENABLE_simint})
#   include(ExternalProject)

   set(SIMINT_VER "0.7")
   string(TOLOWER ${SIMINT_VECTOR} SIMINT_VECTOR)
#   list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
#
   message(STATUS "Using Simint v${SIMINT_VER} (vectorization = ${SIMINT_VECTOR}) for integrals")
#
#   ExternalProject_Add(simint_external
#      URL http://www.bennyp.org/research/simint/download/simint-v${SIMINT_VER}.tar.bz2
#      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
#                 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
#                 -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
#                 -DSIMINT_C_FLAGS=-fPIC
#                 -DSIMINT_VECTOR=${SIMINT_VECTOR}
#                 -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
#      INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install DESTDIR=${CMAKE_BINARY_DIR}/stage
#   )
#
#   set(simint_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/simint CACHE PATH "path to internally built simintConfig.cmake" FORCE)
#
#else()
#    add_library(simint_external INTERFACE)
#endif()


if(${ENABLE_simint})
    find_package(simint CONFIG QUIET)

    if(${simint_FOUND})
        get_property(_loc TARGET simint::simint PROPERTY LOCATION)
        message(STATUS "${Cyan}Found simint${ColourReset}: ${_loc} (found version ${simint_VERSION})")
        add_library(simint_external INTERFACE)  # dummy
    else()
        include(ExternalProject)
        message(STATUS "Suitable simint could not be located, ${Magenta}Building simint${ColourReset} instead.")
        ExternalProject_Add(simint_external
            # BPP, temp repo until changes calm down, then will return to tarball
            #URL http://www.bennyp.org/research/simint/download/simint-v${SIMINT_VER}.tar.bz2
            GIT_REPOSITORY https://github.com/loriab/simint
            GIT_TAG config6
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       #-DSIMINT_BUILD_SHARED=${BUILD_SHARED_LIBS}  # BPP, simint docs say SININT_BUILD_SHARED, but source doesn't mention it
                       -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
#                       -DENABLE_OPENMP=${ENABLE_OPENMP}  # relevant
                       #-DENABLE_XHOST=${ENABLE_XHOST}  # irrelevant, as simint def. to sse which is ~base for xHost, too
                       -DBUILD_FPIC=${BUILD_FPIC}
                       -DENABLE_GENERIC=${ENABLE_GENERIC}
                       -DENABLE_TESTS=OFF
                       -DSIMINT_VECTOR=${SIMINT_VECTOR}
            CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}  # BPP, you really want SIMINT_C_FLAGS, not CMAKE_C_FLAGS?
                             -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
            INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
            DESTDIR=${CMAKE_BINARY_DIR}/stage)

        set(simint_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/simint CACHE PATH "path to internally built simintConfig.cmake" FORCE)
    endif()
else()
    add_library(simint_external INTERFACE)  # dummy
endif()
