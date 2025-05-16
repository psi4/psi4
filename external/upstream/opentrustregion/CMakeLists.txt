if(${ENABLE_OpenTrustRegion})
    find_package(OpenTrustRegion REQUIRED CONFIG)

    if(TARGET OpenTrustRegion::OpenTrustRegion)
        get_property(_loc TARGET OpenTrustRegion::OpenTrustRegion PROPERTY LOCATION)
        message(STATUS "${Cyan}Found OpenTrustRegion${ColourReset}: ${_loc} (found version ${OpenTrustRegion_VERSION})")
        add_library(opentrustregion_external INTERFACE)  # dummy

        if(${CMAKE_INSIST_FIND_PACKAGE_OpenTrustRegion})
            message(VERBOSE "Suitable OpenTrustRegion located externally as user insists")
        endif()

    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_OpenTrustRegion})
            message(FATAL_ERROR "Suitable OpenTrustRegion could not be externally located as user insists")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable OpenTrustRegion could not be located, ${Magenta}Building OpenTrustRegion${ColourReset} instead.")

        set(_otr_dir "share/cmake/OpenTrustRegion")

        ExternalProject_Add(opentrustregion_external
            DEPENDS lapack_external
            URL https://github.com/eriksen-lab/OpenTrustRegion/archive/main.tar.gz
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DOpenTrustRegion_INSTALL_CMAKEDIR=${_otr_dir}
                       -DOpenTrustRegion_BUILD_TESTING=OFF
                       -DTargetLAPACK_DIR=${TargetLAPACK_DIR}
                       -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                       -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
            CMAKE_CACHE_ARGS -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                             -DCMAKE_CXX_COMPILER_ARG1:STRING=${CMAKE_CXX_COMPILER_ARG1})

        set(OpenTrustRegion_DIR ${STAGED_INSTALL_PREFIX}/${_otr_dir} CACHE PATH "path to internally built OpenTrustRegionConfig.cmake" FORCE)
    endif()
else()
    add_library(opentrustregion_external INTERFACE)  # dummy
endif()
