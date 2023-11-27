if(${ENABLE_GTFock})
    find_package(GTFock CONFIG QUIET)

    if(${GTFock_FOUND})
        get_property(_loc TARGET GTFock::gtfock PROPERTY LOCATION)
        message(STATUS "${Cyan}Found GTFock${ColourReset}: ${_loc} (found version ${GTFock_VERSION})")
        message(FATAL_ERROR "GTFock has not been hooked back in after the inversion.")
        add_library(gtfock_external INTERFACE)  # dummy

        if(${CMAKE_INSIST_FIND_PACKAGE_GTFock})
            message(VERBOSE "Suitable GTFock located externally as user insists")
        endif()

    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_GTFock})
            message(FATAL_ERROR "Suitable GTFock could not be externally located as user insists")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable GTFock could not be located, ${Magenta}Building GTFock${ColourReset} instead.")
        message(FATAL_ERROR "GTFock has not been hooked back in after the inversion.")

        ExternalProject_Add(gtfock_external
            GIT_REPOSITORY https://github.com/psi4/GTFock
            #GIT_TAG v1.7.1  # TODO tag when stable
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                       # OpenMP irrelevant
                       -DENABLE_XHOST=${ENABLE_XHOST}
                       -DBUILD_FPIC=${BUILD_FPIC}
                       -DENABLE_GENERIC=${ENABLE_GENERIC}
                       -DLIBC_INTERJECT=${LIBC_INTERJECT}
                       -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                       -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
            CMAKE_CACHE_ARGS -DCMAKE_Fortran_FLAGS:STRING=${CMAKE_Fortran_FLAGS}
                             -DCMAKE_Fortran_COMPILER_ARG1:STRING=${CMAKE_Fortran_COMPILER_ARG1})

        set(GTFock_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/GTFock CACHE PATH "path to internally built GTFockConfig.cmake" FORCE)
    endif()
else()
    add_library(gtfock_external INTERFACE)  # dummy
endif()

