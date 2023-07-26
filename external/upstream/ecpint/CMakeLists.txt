if(${ENABLE_ecpint})
    find_package(ecpint 1.0.7 CONFIG QUIET)

    if(${ecpint_FOUND})
        get_property(_loc TARGET ECPINT::ecpint PROPERTY LOCATION)
        message(STATUS "${Cyan}Found ecpint${ColourReset}: ${_loc} (found version ${ecpint_VERSION})")
        add_library(ecpint_external INTERFACE)  # dummy
    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_ecpint})
            message(FATAL_ERROR "Suitable ecpint could not be externally located as user insists")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable ecpint could not be located, ${Magenta}Building ecpint${ColourReset} instead.")
        ExternalProject_Add(ecpint_external
            URL https://github.com/robashaw/libecpint/archive/v1.0.7.tar.gz
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                       -DLIBECPINT_USE_PUGIXML=OFF
                       -DLIBECPINT_BUILD_TESTS=OFF
                       -DLIBECPINT_BUILD_DOCS=OFF
                       -DLIBECPINT_MAX_UNROL=1
                       -DLIBECPINT_MAX_L=5
                       -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                       -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
            CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                             -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                             -DCMAKE_C_COMPILER_ARG1:STRING=${CMAKE_C_COMPILER_ARG1}
                             -DCMAKE_CXX_COMPILER_ARG1:STRING=${CMAKE_CXX_COMPILER_ARG1})

        set(ecpint_DIR ${STAGED_INSTALL_PREFIX}/lib/cmake/ecpint CACHE PATH "path to internally built ecpint-config.cmake" FORCE)
    endif()
else()
    add_library(ecpint_external INTERFACE)  # dummy
endif()

