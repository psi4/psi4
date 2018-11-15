if(${ENABLE_cppe})
    find_package(cppe 0.0.0 CONFIG QUIET COMPONENTS shallow)

    if(${cppe_FOUND})
        get_property(_loc TARGET cppe::cppe PROPERTY LOCATION)
        message(STATUS "${Cyan}Found cppe${ColourReset}: ${_loc} (found version ${cppe_VERSION})")
        add_library(cppe_external INTERFACE)  # dummy
    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_cppe})
            message(FATAL_ERROR "Suitable libefp could not be externally located as user insists")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable cppe could not be located, ${Magenta}Building cppe${ColourReset} instead.")
        ExternalProject_Add(cppe_external
            #DEPENDS armadillo
            GIT_REPOSITORY https://github.com/maxscheurer/cppe
            #GIT_TAG 15cd7ce  # v1.5.0 + 10 (docs and a cmake lapack patch)
            GIT_TAG dev-psi4
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
                       -DCMAKE_INSTALL_DATADIR=${CMAKE_INSTALL_DATADIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                       -DENABLE_OPENMP=ON # Psi4 sometimes reacts poorly to threaded efp
                       #-DENABLE_XHOST=${ENABLE_XHOST}
                       -DBUILD_FPIC=${BUILD_FPIC}
                       -DENABLE_GENERIC=${ENABLE_GENERIC}
                       -DLIBC_INTERJECT=${LIBC_INTERJECT}
                       -DINSTALL_DEVEL_HEADERS=ON
                       -DFRAGLIB_UNDERSCORE_L=OFF
                       #-DFRAGLIB_DEEP=OFF
                       -DOpenMP_LIBRARY_DIRS=${OpenMP_LIBRARY_DIRS}
                       #-DTargetLAPACK_DIR=${STAGED_INSTALL_PREFIX}/share/cmake/TargetLAPACK
            CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                             -DTargetOpenMP_FIND_COMPONENTS:STRING=C)

        set(cppe_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/cppe CACHE PATH "path to internally built cppeConfig.cmake" FORCE)
    endif()
else()
    add_library(cppe_external INTERFACE)  # dummy
endif()
