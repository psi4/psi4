find_package(Libint2 2.6 CONFIG QUIET COMPONENTS ${MAX_AM_ERI})

if(${Libint2_FOUND})
    get_property(_loc TARGET Libint2::int2 PROPERTY LOCATION)
    message(STATUS "${Cyan}Found Libint2 ${Libint_MAX_AM_ERI}${ColourReset}: ${_loc} (found version ${Libint_VERSION})")
    add_library(libint2_external INTERFACE)  # dummy
else()
    if(${CMAKE_INSIST_FIND_PACKAGE_Libint2})
        message(FATAL_ERROR "Suitable Libint2 could not be externally located as user insists")
    endif()

    include(ExternalProject)
    message(STATUS "Suitable Libint2 could not be located, ${Magenta}Building Libint2${ColourReset} instead.")
    ExternalProject_Add(libint2_external
        # "git checkout" fails on Windows, because of "*" in filenames (e.g. basis set files)
        #URL https://github.com/evaleev/libint/archive/b13e71d3cf9960460c4019e5ecf2546a5f361c71.zip
        URL https://github.com/evaleev/libint/releases/download/v2.6.0/libint-2.6.0.tgz
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                   -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                   -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                   -DMAX_AM_ERI=${MAX_AM_ERI}
                   -DLIBINT2_SHGAUSS_ORDERING=gaussian
                   -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                   # OpenMP irrelevant
                   -DENABLE_XHOST=${ENABLE_XHOST}
                   -DBUILD_FPIC=${BUILD_FPIC}
                   -DENABLE_GENERIC=${ENABLE_GENERIC}
                   -DLIBC_INTERJECT=${LIBC_INTERJECT}
                   -DMERGE_LIBDERIV_INCLUDEDIR=OFF
        CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                         -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS})

    #set(Libint2_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/Libint2 CACHE PATH "path to internally built Libint2Config.cmake" FORCE)
    set(Libint2_DIR ${STAGED_INSTALL_PREFIX}/lib/cmake/libint2 CACHE PATH "path to internally built libint2-config.cmake" FORCE)

endif()

