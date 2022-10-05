if(${ENABLE_mdi})
    find_package(mdi 1.2.3 CONFIG QUIET)

    if(${mdi_FOUND})
        get_property(_loc TARGET mdi::mdi PROPERTY LOCATION)
        message(STATUS "${Cyan}Found mdi${ColourReset}: ${_loc} (found version ${mdi_VERSION})")
        add_library(mdi_external INTERFACE)  # dummy
    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_mdi})
            message(FATAL_ERROR "Suitable mdi could not be externally located as user insists")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable mdi could not be located, ${Magenta}Building mdi${ColourReset} instead.")
        ExternalProject_Add(mdi_external
            URL https://github.com/MolSSI-MDI/MDI_Library/archive/v1.2.3.tar.gz
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                       -DENABLE_OPENMP=${ENABLE_OPENMP}  # relevant for thread safety
                       -DENABLE_XHOST=${ENABLE_XHOST}
                       -DBUILD_FPIC=${BUILD_FPIC}
                       -DENABLE_GENERIC=${ENABLE_GENERIC}
                       -DLIBC_INTERJECT=${LIBC_INTERJECT}
                       -DOpenMP_LIBRARY_DIRS=${OpenMP_LIBRARY_DIRS}
                       -Dlanguage=Python
                       -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                       -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
            CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                             -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                             -DTargetOpenMP_FIND_COMPONENTS:STRING=C;CXX)

        set(mdi_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/mdi CACHE PATH "path to internally built mdiConfig.cmake" FORCE)
    endif()
else()
    add_library(mdi_external INTERFACE)  # dummy
endif()
