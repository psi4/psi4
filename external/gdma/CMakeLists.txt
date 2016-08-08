include(ExternalProject)
option(ENABLE_BOUNDS_CHECK "Enables bounds check in Fortran" OFF)
find_package(gdma CONFIG QUIET)

if(NOT gdma_FOUND AND ${ENABLE_gdma})
    message(STATUS "Suitable gdma could not be located, building one instead.")
    ExternalProject_Add(gdma_external
        GIT_REPOSITORY https://github.com/loriab/gdma
        #GIT_TAG v1.7.1  # TODO tag when stable
        UPDATE_COMMAND ""
        CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                   -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                   -DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
                   -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/external/gdma
                   #-DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                   -DBUILD_SHARED_LIBS=ON
                   -DENABLE_GENERIC=${ENABLE_GENERIC}
                   -DENABLE_XHOST=${ENABLE_XHOST}
                   -DBUILD_FPIC=${BUILD_FPIC}
                   -DLIBC_INTERJECT=${LIBC_INTERJECT}
                   -DCMAKE_INSTALL_LIBDIR=lib
        INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
                        DESTDIR=${CMAKE_BINARY_DIR}/stage
    )
#    message("in ext pree: ${CMAKE_PREFIX_PATH}")
##    list(APPEND CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR}/stage${CMAKE_INSTALL_PREFIX}/external/gdma)
#    #set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "${CMAKE_BINARY_DIR}/stage${CMAKE_INSTALL_PREFIX}/external/gdma" PARENT_SCOPE)
#    set(CMAKE_PREFIX_PATH "${CMAKE_BINARY_DIR}/stage${CMAKE_INSTALL_PREFIX}/external/gdma" PARENT_SCOPE)
#    message("in ext post: ${CMAKE_PREFIX_PATH}")
else()
    add_library(gdma_external INTERFACE)  # dummy
endif()


