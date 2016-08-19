find_package(libefp CONFIG QUIET)

if(${libefp_FOUND})
    get_property(_loc TARGET libefp::efp PROPERTY LOCATION)
    message(STATUS "Suitable libefp located at ${_loc}")
    add_library(libefp_external INTERFACE)  # dummy
else()
    include(ExternalProject)
    message(STATUS "Suitable libefp could not be located, building one instead.")
    ExternalProject_Add(libefp_external
        GIT_REPOSITORY https://github.com/psi4/libefp
        #GIT_TAG v1.7.1  # TODO tag when stable
        UPDATE_COMMAND ""
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/external/libefp
                   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                   -DCMAKE_INSTALL_LIBDIR=lib
                   -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                   #-DBUILD_SHARED_LIBS=ON  # TODO switch to above line when addons all capable
                   -DENABLE_OPENMP=${ENABLE_OPENMP}
                   -DENABLE_XHOST=${ENABLE_XHOST}
                   -DBUILD_FPIC=${BUILD_FPIC}
                   -DENABLE_GENERIC=${ENABLE_GENERIC}
                   -DLIBC_INTERJECT=${LIBC_INTERJECT}
                   -DFRAGLIB_UNDERSCORE_L=OFF
        CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS} 
        INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
        DESTDIR=${CMAKE_BINARY_DIR}/stage)
endif()

