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
