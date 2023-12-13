if(${ENABLE_Libint1t})
find_package(Libint 1.2.0 CONFIG QUIET COMPONENTS ${MAX_AM_ERI})

if(${Libint_FOUND})
    get_property(_loc TARGET Libint::int PROPERTY LOCATION)
    message(STATUS "${Cyan}Found Libint ${Libint_MAX_AM_ERI}${ColourReset}: ${_loc} (found version ${Libint_VERSION})")
    add_library(libint_external INTERFACE)  # dummy
else()
    if(${CMAKE_INSIST_FIND_PACKAGE_Libint})
        message(FATAL_ERROR "Suitable Libint could not be externally located as user insists")
    endif()

    include(ExternalProject)
    message(STATUS "Suitable Libint could not be located, ${Magenta}Building Libint${ColourReset} instead.")
    ExternalProject_Add(libint_external
        # "git checkout" fails on Windows, because of "*" in filenames (e.g. basis set files)
        URL https://github.com/loriab/libint/archive/libint_t.tar.gz
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                   -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                   -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                   -DMAX_AM_ERI=${MAX_AM_ERI}
                   -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                   # OpenMP irrelevant
                   -DENABLE_XHOST=${ENABLE_XHOST}
                   -DBUILD_FPIC=${BUILD_FPIC}
                   -DENABLE_GENERIC=${ENABLE_GENERIC}
                   -DLIBC_INTERJECT=${LIBC_INTERJECT}
                   -DMERGE_LIBDERIV_INCLUDEDIR=OFF
        CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                         -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                         -DCMAKE_C_COMPILER_ARG1:STRING=${CMAKE_C_COMPILER_ARG1}
                         -DCMAKE_CXX_COMPILER_ARG1:STRING=${CMAKE_CXX_COMPILER_ARG1})

    set(Libint_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/Libint CACHE PATH "path to internally built LibintConfig.cmake" FORCE)
endif()
else()
    add_library(libint_external INTERFACE)  # dummy
endif()

