if(${ENABLE_OpenTrustRegion})
    find_package(OpenTrustRegion 3.0.0 COMPONENTS int32 CONFIG)  # edit in codedeps

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

        set(_otr_dir "lib/cmake/OpenTrustRegion")

        # Link the BLAS and LAPACK Psi4 links rather than whatever OpenTrustRegion's own
        # detection finds. Handing libraries over requires stating their integer size, LP64
        # for Psi4. If the LAPACK list carries the BLAS too, leaving the BLAS list empty, give
        # it as both.
        get_property(_otr_lapack TARGET lapk PROPERTY INTERFACE_LINK_LIBRARIES)
        get_property(_otr_blas TARGET blas PROPERTY INTERFACE_LINK_LIBRARIES)
        if(NOT _otr_blas)
            set(_otr_blas ${_otr_lapack})
        endif()

        ExternalProject_Add(opentrustregion_external
            DEPENDS lapack_external
            URL https://github.com/eriksen-lab/opentrustregion/archive/v3.0.0.tar.gz  # edit in codedeps
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                       -DOpenTrustRegion_ENABLE_XHOST=${ENABLE_XHOST}
                       -DOpenTrustRegion_INSTALL_CMAKEDIR=${_otr_dir}
                       -DOpenTrustRegion_BUILD_TESTING=OFF
                       -DINTEGER_SIZE=4
                       -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                       -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
            CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                             -DCMAKE_Fortran_FLAGS:STRING=${CMAKE_Fortran_FLAGS}
                             -DCMAKE_C_COMPILER_ARG1:STRING=${CMAKE_C_COMPILER_ARG1}
                             -DCMAKE_Fortran_COMPILER_ARG1:STRING=${CMAKE_Fortran_COMPILER_ARG1}
                             -DBLAS_LIBRARIES:STRING=${_otr_blas}
                             -DLAPACK_LIBRARIES:STRING=${_otr_lapack})

        set(OpenTrustRegion_DIR ${STAGED_INSTALL_PREFIX}/${_otr_dir} CACHE PATH "path to internally built OpenTrustRegionConfig.cmake" FORCE)
    endif()
else()
    message(STATUS "Optional addon ${Blue}unsought: OpenTrustRegion${ColourReset}")
    add_library(opentrustregion_external INTERFACE)  # dummy
endif()
