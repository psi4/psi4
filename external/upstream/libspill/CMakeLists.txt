find_package(libspill 0.1 CONFIG QUIET)  # edit in codedeps

if(${libspill_FOUND})
    get_property(_loc TARGET libspill::spill PROPERTY LOCATION)
    message(STATUS "${Cyan}Found libspill${ColourReset}: ${_loc} (found version ${libspill_VERSION})")
    add_library(libspill_external INTERFACE)  # dummy

    if(${CMAKE_INSIST_FIND_PACKAGE_libspill})
        message(VERBOSE "Suitable libspill located externally as user insists")
    endif()

else()
    if(${CMAKE_INSIST_FIND_PACKAGE_libspill})
        message(FATAL_ERROR "Suitable libspill could not be externally located as user insists")
    endif()

    include(ExternalProject)
    message(STATUS "Suitable libspill could not be located, ${Magenta}Building libspill${ColourReset} instead.")
    ExternalProject_Add(libspill_external
        URL https://github.com/ReLibQC/libspill/archive/29eb3ba.tar.gz  # v0.1.0+  # edit in codedeps
        UPDATE_COMMAND ""
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                   -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                   -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                   -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
                   -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                   -DCMAKE_POSITION_INDEPENDENT_CODE=${BUILD_FPIC}
                   # I/O core only: no HDF5 backend, no Fortran module, no tests
                   -DLIBSPILL_WITH_HDF5=OFF
                   -DLIBSPILL_BUILD_FORTRAN=OFF
                   -DLIBSPILL_BUILD_TESTS=OFF
                   -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                   -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
        CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                         -DCMAKE_C_COMPILER_ARG1:STRING=${CMAKE_C_COMPILER_ARG1})

    set(libspill_DIR ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/cmake/libspill CACHE PATH "path to internally built libspillConfig.cmake" FORCE)
endif()
