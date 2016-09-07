option(ENABLE_BOUNDS_CHECK "Enables bounds check in Fortran" OFF)

if(${ENABLE_liberd})
    find_package(liberd CONFIG QUIET)

    if(${liberd_FOUND})
        get_property(_loc TARGET liberd::erd PROPERTY LOCATION)
        message(STATUS "Found liberd: ${_loc} (found version ${liberd_VERSION})")
        add_library(liberd_external INTERFACE)  # dummy
    else()
        include(ExternalProject)
        message(STATUS "Suitable liberd could not be located, building one instead.")
        ExternalProject_Add(liberd_external
            GIT_REPOSITORY https://github.com/psi4/liberd
            #GIT_TAG v1.7.1  # TODO tag when stable
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/external/liberd
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       #-DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                       -DBUILD_SHARED_LIBS=ON  # TODO switch to above line when addons all capable
                       # OpenMP irrelevant
                       -DENABLE_XHOST=${ENABLE_XHOST}
                       -DBUILD_FPIC=${BUILD_FPIC}
                       -DENABLE_GENERIC=${ENABLE_GENERIC}
                       -DLIBC_INTERJECT=${LIBC_INTERJECT}
            CMAKE_CACHE_ARGS -DCMAKE_Fortran_FLAGS:STRING=${CMAKE_Fortran_FLAGS}
            INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
            DESTDIR=${CMAKE_BINARY_DIR}/stage)
    endif()
else()
    add_library(liberd_external INTERFACE)  # dummy
endif()

