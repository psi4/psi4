find_package(libint CONFIG QUIET COMPONENTS ${MAX_AM_ERI})

if(${libint_FOUND})
    get_property(_loc TARGET libint::int PROPERTY LOCATION)
    message(STATUS "Found libint: ${_loc} (found version ${libint_VERSION})")
    add_library(libint_external INTERFACE)  # dummy
else()
    include(ExternalProject)
    message(STATUS "Suitable libint could not be located, building one instead.")
    ExternalProject_Add(libint_external
        GIT_REPOSITORY https://github.com/psi4/libint
        #GIT_TAG v1.7.1  # TODO tag when stable
        UPDATE_COMMAND ""
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/external/libint
                   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                   -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                   -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                   -DMAX_AM_ERI=${MAX_AM_ERI}
                   -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                   #-DBUILD_SHARED_LIBS=ON  # TODO switch to above line when addons all capable
                   # OpenMP irrelevant
                   -DENABLE_XHOST=${ENABLE_XHOST}
                   -DBUILD_FPIC=${BUILD_FPIC}
                   -DENABLE_GENERIC=${ENABLE_GENERIC}
                   -DLIBC_INTERJECT=${LIBC_INTERJECT}
        CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                         -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
        INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
        DESTDIR=${CMAKE_BINARY_DIR}/stage)
endif()

# prob needed for libm detect
#   list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
#     CMAKE_CACHE_ARGS -DCMAKE_MODULE_PATH:STRING=${CMAKE_MODULE_PATH}
#                      -DCMAKE_PREFIX_PATH:STRING=${CMAKE_PREFIX_PATH}

