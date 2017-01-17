find_package(Libxc CONFIG QUIET)

if(${Libxc_FOUND})
    get_property(_loc TARGET Libxc::xc PROPERTY LOCATION)
    message(STATUS "${Cyan}Found Libxc${ColourReset}: ${_loc} (found version ${Libxc_VERSION})")
    add_library(libxc_external INTERFACE)  # dummy
else()
    include(ExternalProject)
    message(STATUS "Suitable Libxc could not be located, ${Magenta}Building Libxc${ColourReset} instead.")
    ExternalProject_Add(libxc_external
        GIT_REPOSITORY https://github.com/loriab/libxc
        GIT_TAG master
        #GIT_TAG v3.0.0  # TODO tag when stable
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
        INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
        DESTDIR=${CMAKE_BINARY_DIR}/stage)

    set(Libxc_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/Libxc CACHE PATH "path to internally built LibxcConfig.cmake" FORCE)
endif()
