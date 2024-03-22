if(${ENABLE_gauxc})
    #== perhaps ExchCXX and IntegratorXX were built alongside GauXC, check for this ==#
    if(NOT "${gauxc_DIR}" STREQUAL "")
        if("${ExchCXX_DIR}" STREQUAL "")
           set(ExchCXX_DIR ${gauxc_DIR}/../ExchCXX)
        endif()
    endif()

    #== now search for package ==#
    find_package(gauxc 0.0.1 CONFIG) 

    if(${gauxc_FOUND})
        #== check if GPU-enabled GauXC is found, if GPU is specified ==#
        # commented stuff is from manually-input exported properties modification, getting GPU info as a property of GauXC     
        
        #get_property(_gpu TARGET gauxc::gauxc PROPERTY GAUXC_GPU_ENABLED)
        
        #if(${gauxc_ENABLE_GPU} AND (NOT ${_gpu}))
        if(${gauxc_ENABLE_GPU} AND (NOT ${GAUXC_ENABLE_CUDA}))
           message(FATAL_ERROR "gauxc_ENABLE_GPU turned on, but selected GauXC install is not built for GPU support!")
        #elseif((NOT ${gauxc_ENABLE_GPU}) AND ${_gpu})
        elseif((NOT ${gauxc_ENABLE_GPU}) AND ${GAUXC_ENABLE_CUDA})
           message(FATAL_ERROR "gauxc_ENABLE_GPU turned off, but selected GauXC install is built for GPU support!")
        endif()

        #== normal print-out and library addition ==#
        get_property(_loc TARGET gauxc::gauxc PROPERTY IMPORTED_LOCATION_RELWITHDEBINFO)
        #get_property(_loc TARGET gauxc::gauxc PROPERTY IMPORTED_LOCATION)
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

        include(ExternalProject)
        ExternalProject_Add(gauxc_external
            DEPENDS integratorxx_external 
            URL https://github.com/wavefunction91/GauXC/archive/issue_56_58/reduce_shell_pair_mem.tar.gz 
            #GIT_REPOSITORY https://github.com/wavefunction91/gauxc.git
            #GIT_TAG issue_56_58/reduce_shell_pair_mem  
            #GIT_TAG issue_93/external_gau2grid 
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
            CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                         -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                         -DCMAKE_C_COMPILER_ARG1:STRING=${CMAKE_C_COMPILER_ARG1}
                         -DCMAKE_CXX_COMPILER_ARG1:STRING=${CMAKE_CXX_COMPILER_ARG1})
        set(gauxc_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/gauxc CACHE PATH "path to internally built gauxcConfig.cmake" FORCE)
    endif()
else()
    add_library(gauxc_external INTERFACE)  # dummy
endif()
