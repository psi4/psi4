if(${ENABLE_cppe})
    find_package(cppe 0.3.1 CONFIG QUIET)

    if(${cppe_FOUND})
        get_property(_loc TARGET cppe::cppe PROPERTY LOCATION)
        message(STATUS "${Cyan}Found cppe${ColourReset}: ${_loc} (found version ${cppe_VERSION})")
        add_library(cppe_external INTERFACE)  # dummy
    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_cppe})
            message(FATAL_ERROR "Suitable cppe could not be externally located as user insists")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable cppe could not be located, ${Magenta}Building cppe${ColourReset} instead.")
        ExternalProject_Add(cppe_external
            # GIT_REPOSITORY https://github.com/maxscheurer/cppe
            # GIT_TAG master
            URL https://github.com/maxscheurer/cppe/archive/v0.3.1.tar.gz
            # UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
                       -DCMAKE_INSTALL_DATADIR=${CMAKE_INSTALL_DATADIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                       -DENABLE_OPENMP=OFF
                       -DENABLE_PYTHON_INTERFACE=ON
                       -DBUILD_FPIC=${BUILD_FPIC}
                       -DENABLE_GENERIC=${ENABLE_GENERIC}
                       -DLIBC_INTERJECT=${LIBC_INTERJECT}
                       -DINSTALL_DEVEL_HEADERS=OFF
                       -DOpenMP_LIBRARY_DIRS=${OpenMP_LIBRARY_DIRS}
            CMAKE_CACHE_ARGS -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                             -DTargetOpenMP_FIND_COMPONENTS:STRING=CXX)

        set(cppe_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/cppe CACHE PATH "path to internally built cppeConfig.cmake" FORCE)
    endif()
else()
    add_library(cppe_external INTERFACE)  # dummy
endif()
