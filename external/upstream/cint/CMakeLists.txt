if(${ENABLE_cint})
    find_package(cint 6.2.0 CONFIG QUIET)  # edit in codedeps

    if(${cint_FOUND})
        get_property(_loc TARGET cint::cint PROPERTY LOCATION)
        message(STATUS "${Cyan}Found cint${ColourReset}: ${_loc} (found version ${cint_VERSION})")
        add_library(cint_external INTERFACE)  # dummy

        if(${CMAKE_INSIST_FIND_PACKAGE_cint})
            message(VERBOSE "Suitable cint located externally as user insists")
        endif()

    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_cint})
            message(FATAL_ERROR "Suitable cint could not be externally located as user insists")
        endif()

        if(MSVC OR CMAKE_C_SIMULATE_ID STREQUAL "MSVC")
            # libcint's spinor integrals use C99 complex arithmetic (double complex, _Complex_I), which the
            # MSVC runtime's complex.h doesn't provide, so compilers on it (cl, clang-cl, clang targeting
            # MSVC) can't build it. MinGW gcc can.
            message(FATAL_ERROR "Suitable cint could not be located, and it can't be built from source with "
                                "${CMAKE_C_COMPILER_ID} (MSVC-like), whose C runtime lacks C99 complex arithmetic. "
                                "Install a gcc-built libcint >= 6.2.0 (e.g., conda-forge libcint) and point "
                                "cint_DIR to its lib/cmake/cint, or set ENABLE_cint=OFF.")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable cint could not be located, ${Magenta}Building cint${ColourReset} instead.")
        ExternalProject_Add(cint_external
            URL https://github.com/sunqm/libcint/archive/v6.2.0.tar.gz  # edit in codedeps
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                       -DWITH_CINT2_INTERFACE=0
                       -DWITH_FORTRAN=0
                       -DENABLE_TEST=0
                       -DENABLE_EXAMPLE=0
                       -DCMAKE_POSITION_INDEPENDENT_CODE=${BUILD_FPIC}
                       -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                       -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
            CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                             -DCMAKE_C_COMPILER_ARG1:STRING=${CMAKE_C_COMPILER_ARG1})

        set(cint_DIR ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/cmake/cint CACHE PATH "path to internally built cintConfig.cmake" FORCE)
    endif()
else()
    message(STATUS "Optional addon ${Blue}unsought: cint${ColourReset}")
    add_library(cint_external INTERFACE)  # dummy
endif()
