if(${ENABLE_gauxc})
    #== GauXC requires Psi4 to have standard integral ordering rather than gaussian ==#
    if(${psi4_SHGAUSS_ORDERING} STREQUAL gaussian)
        message(FATAL_ERROR "Psi4 with the GauXC interface must be compiled using psi4_SHGAUSS_ORDERING=standard.")
    endif()
  
    find_package(gauxc 0.0.1 CONFIG) 

    if(${gauxc_FOUND})
        get_property(_loc TARGET gauxc::gauxc PROPERTY IMPORTED_LOCATION_RELWITHDEBINFO)
        message(STATUS "${Cyan}Found gauxc ${ColourReset}: ${_loc} (found version ${gauxc_VERSION})")
        add_library(gauxc_external INTERFACE)  # dummy

        if(${CMAKE_INSIST_FIND_PACKAGE_gauxc})
            message(VERBOSE "Suitable gauxc located externally as user insists")
        endif()

    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_gauxc})
            message(FATAL_ERROR "Suitable gauxc could not be externally located as user insists")
        endif()

        message(STATUS "Suitable gauxc could not be located, ${Magenta}Building gauxc${ColourReset} instead.")

        #include(FetchContent)
        #set( GAUXC_ENABLE_TESTS OFF CACHE BOOL "" )
        #set( GAUXC_ENABLE_MPI OFF CACHE BOOL "" )
        #set( GAUXC_ENABLE_HDF5 OFF CACHE BOOL "" )
        #FetchContent_Declare( gauxc
        #    GIT_REPOSITORY https://github.com/wavefunction91/GauXC.git
        #    GIT_TAG master
        # )
        #FetchContent_MakeAvailable( gauxc )

        include(ExternalProject)
        ExternalProject_Add(gauxc_external
            GIT_REPOSITORY https://github.com/wavefunction91/gauxc.git
            #GIT_TAG master
            GIT_TAG hotfix/small_primpair_max_override 
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                   -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                   -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                   -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                   -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
                   -DGAUXC_ENABLE_CUDA=${gauxc_ENABLE_GPU}
                   -DCMAKE_CUDA_ARCHITECTURES=70
                   -DCMAKE_POSITION_INDEPENDENT_CODE=ON
                   -DBUILD_SHARED_LIBS=OFF
                   -DGAUXC_ENABLE_MPI=OFF
                   -DGAUXC_ENABLE_TESTS=OFF
                   -DGAUXC_ENABLE_HDF5=OFF
                   -DGAUXC_ENABLE_MAGMA=OFF
                   -DGAUXC_SHELLPAIR_PRIMPAIR_MAX=256
            CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                         -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                         -DCMAKE_C_COMPILER_ARG1:STRING=${CMAKE_C_COMPILER_ARG1}
                         -DCMAKE_CXX_COMPILER_ARG1:STRING=${CMAKE_CXX_COMPILER_ARG1})
        set(gauxc_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/gauxc CACHE PATH "path to internally built gauxcConfig.cmake" FORCE)
    endif()
else()
    add_library(gauxc_external INTERFACE)  # dummy
endif()
