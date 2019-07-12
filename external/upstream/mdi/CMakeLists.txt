##if(${ENABLE_mdi})
    find_package(mdi CONFIG QUIET)

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
            URL https://github.com/MolSSI/MDI_Library/archive/v0.6.2.tar.gz
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                       -DENABLE_OPENMP=${ENABLE_OPENMP}  # relevant for thread safety
                       -DENABLE_XHOST=${ENABLE_XHOST}
                       -DBUILD_FPIC=${BUILD_FPIC}
                       -DENABLE_GENERIC=${ENABLE_GENERIC}
                       -DLIBC_INTERJECT=${LIBC_INTERJECT}
		       -Dlanguage=Python
            CMAKE_CACHE_ARGS -DCMAKE_Fortran_FLAGS:STRING=${CMAKE_Fortran_FLAGS}
                             -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                             -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS})

        set(mdi_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/mdi CACHE PATH "path to internally built mdiConfig.cmake" FORCE)
    endif()
##else()
##    add_library(mdi_external INTERFACE)  # dummy
##endif()

