if(${ENABLE_gdma})
    find_package(gdma 2.3.3 CONFIG QUIET)

    if(${gdma_FOUND})
        get_property(_loc TARGET gdma::gdma PROPERTY LOCATION)
        message(STATUS "${Cyan}Found gdma${ColourReset}: ${_loc} (found version ${gdma_VERSION})")
        add_library(gdma_external INTERFACE)  # dummy
    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_gdma})
            message(FATAL_ERROR "Suitable gdma could not be externally located as user insists")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable gdma could not be located, ${Magenta}Building gdma${ColourReset} instead.")
        ExternalProject_Add(gdma_external
            URL https://github.com/psi4/gdma/archive/9d607d7.tar.gz  # v2.2.6-2-g9d607d7
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                       -DENABLE_OPENMP=${ENABLE_OPENMP}  # relevant for thread safety
                       -DENABLE_XHOST=${ENABLE_XHOST}
                       -DBUILD_FPIC=${BUILD_FPIC}
                       -DENABLE_GENERIC=${ENABLE_GENERIC}
                       -DLIBC_INTERJECT=${LIBC_INTERJECT}
                       -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                       -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
            CMAKE_CACHE_ARGS -DCMAKE_Fortran_FLAGS:STRING=${CMAKE_Fortran_FLAGS}
                             -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS})

        set(gdma_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/gdma CACHE PATH "path to internally built gdmaConfig.cmake" FORCE)
    endif()
else()
    add_library(gdma_external INTERFACE)  # dummy
endif()

