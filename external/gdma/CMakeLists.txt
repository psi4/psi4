option(ENABLE_BOUNDS_CHECK "Enables bounds check in Fortran" OFF)

if(${ENABLE_gdma})
    find_package(gdma CONFIG QUIET)

    if(${gdma_FOUND})
        get_property(_loc TARGET gdma::gdma PROPERTY LOCATION)
        message(STATUS "Suitable gdma located at ${_loc}")
        add_library(gdma_external INTERFACE)  # dummy
    else()
        include(ExternalProject)
        message(STATUS "Suitable gdma could not be located, building one instead.")
        ExternalProject_Add(gdma_external
            GIT_REPOSITORY https://github.com/loriab/gdma  # TODO switch origin to psi4 when ready
            #GIT_TAG v1.7.1  # TODO tag when stable
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/external/gdma
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=lib
                       #-DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                       -DBUILD_SHARED_LIBS=ON  # TODO switch to above line when addons all capable
                       -DENABLE_XHOST=${ENABLE_XHOST}
                       -DBUILD_FPIC=${BUILD_FPIC}
                       -DENABLE_GENERIC=${ENABLE_GENERIC}
                       -DLIBC_INTERJECT=${LIBC_INTERJECT}
            CMAKE_CACHE_ARGS -DCMAKE_Fortran_FLAGS:STRING=${CMAKE_Fortran_FLAGS} 
            INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
            DESTDIR=${CMAKE_BINARY_DIR}/stage)
    endif()
else()
    add_library(gdma_external INTERFACE)  # dummy
endif()

