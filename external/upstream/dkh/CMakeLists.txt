if(${ENABLE_dkh})
    find_package(dkh 1.2 CONFIG QUIET)

    if(${dkh_FOUND})
        get_property(_loc TARGET dkh::dkh PROPERTY LOCATION)
        message(STATUS "${Cyan}Found dkh${ColourReset}: ${_loc} (found version ${dkh_VERSION})")
        add_library(dkh_external INTERFACE)  # dummy
    else()
        include(ExternalProject)
        message(STATUS "Suitable dkh could not be located, ${Magenta}Building dkh${ColourReset} instead.")
        ExternalProject_Add(dkh_external
            GIT_REPOSITORY https://github.com/psi4/dkh
            GIT_TAG v1.2
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
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
            CMAKE_CACHE_ARGS -DCMAKE_Fortran_FLAGS:STRING=${CMAKE_Fortran_FLAGS}
            INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
            DESTDIR=${CMAKE_BINARY_DIR}/stage)

        set(dkh_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/dkh CACHE PATH "path to internally built dkhConfig.cmake" FORCE)
    endif()
else()
    add_library(dkh_external INTERFACE)  # dummy
endif()

