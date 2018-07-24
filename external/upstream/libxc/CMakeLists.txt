find_package(Libxc 4.0.2 CONFIG QUIET)

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
        GIT_REPOSITORY https://github.com/psi4/libxc
        GIT_TAG libxc4retweaked  # v4.0 (0387b1d1) + !62 + MVS corr + 4.0.2 bump + 2 build commits = c754e28
        UPDATE_COMMAND ""
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                   -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                   -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
                   -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                   -DNAMESPACE_INSTALL_INCLUDEDIR=/libxc
                   -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                   # OpenMP irrelevant
                   -DENABLE_XHOST=${ENABLE_XHOST}
                   -DBUILD_FPIC=${BUILD_FPIC}
                   -DENABLE_GENERIC=${ENABLE_GENERIC}
                   -DLIBC_INTERJECT=${LIBC_INTERJECT}
                   -DBUILD_TESTING=OFF
        CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
        INSTALL_COMMAND DESTDIR=${CMAKE_BINARY_DIR}/stage ${CMAKE_MAKE_PROGRAM} install)

    set(Libxc_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/Libxc CACHE PATH "path to internally built LibxcConfig.cmake" FORCE)
endif()
