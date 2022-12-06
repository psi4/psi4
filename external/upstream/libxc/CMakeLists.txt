find_package(Libxc 6.0.0 CONFIG QUIET COMPONENTS C)
if(NOT ${Libxc_FOUND})
    # LibxcConfig.cmake encodes searching only within major version, so need two calls to allow v5 or v6
    find_package(Libxc 5.1.2 CONFIG QUIET COMPONENTS C)
endif()

if(${Libxc_FOUND})
    get_property(_loc TARGET Libxc::xc PROPERTY LOCATION)
    message(STATUS "${Cyan}Found Libxc${ColourReset}: ${_loc} (found version ${Libxc_VERSION})")
    add_library(libxc_external INTERFACE)  # dummy
else()
    if(${CMAKE_INSIST_FIND_PACKAGE_Libxc})
        message(FATAL_ERROR "Suitable Libxc could not be externally located as user insists")
    endif()

    include(ExternalProject)
    message(STATUS "Suitable Libxc could not be located, ${Magenta}Building Libxc${ColourReset} instead.")
    ExternalProject_Add(libxc_external
        # Default: use a stable release tarball of libxc. To use the
        # development version of libxc, instead, comment the URL line,
        # and uncomment the GIT lines.
        URL https://gitlab.com/libxc/libxc/-/archive/6.0.0/libxc-6.0.0.tar.gz
        #GIT_REPOSITORY https://gitlab.com/libxc/libxc.git
        #GIT_TAG 5.1.5
        #UPDATE_COMMAND ""
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                   -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                   -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
                   -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                   -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
                   -DNAMESPACE_INSTALL_INCLUDEDIR=/
                   -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                   # OpenMP irrelevant
                   -DENABLE_XHOST=${ENABLE_XHOST}
                   -DBUILD_FPIC=${BUILD_FPIC}
                   -DENABLE_GENERIC=${ENABLE_GENERIC}
                   -DLIBC_INTERJECT=${LIBC_INTERJECT}
                   -DBUILD_TESTING=OFF
                   -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                   -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
        CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS})

    set(Libxc_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/Libxc CACHE PATH "path to internally built LibxcConfig.cmake" FORCE)
endif()
