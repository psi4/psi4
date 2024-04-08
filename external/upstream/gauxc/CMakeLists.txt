if(${ENABLE_gauxc})
    if((NOT "${gauxc_DIR}" STREQUAL "") AND "${IntegratorXX_DIR}" STREQUAL "")
       set(IntegratorXX_DIR ${gauxc_DIR}/../IntegratorXX)
    endif()

    find_package(ExchCXX CONFIG REQUIRED PATHS ${gauxc_DIR}/../ExchCXX) 
    find_package(gauxc CONFIG) 

    if(${gauxc_FOUND})
        #== check if GPU-enabled GauXC is found, if GPU is specified ==#
        get_property(_gpu TARGET gauxc::gauxc PROPERTY GAUXC_HAS_DEVICE)
        
        if(${gauxc_ENABLE_GPU} AND (NOT ${_gpu}))
           message(FATAL_ERROR "gauxc_ENABLE_GPU turned on, but selected GauXC install is not built for GPU support!")
        elseif((NOT ${gauxc_ENABLE_GPU}) AND ${_gpu})
           message(FATAL_ERROR "gauxc_ENABLE_GPU turned off, but selected GauXC install is built for GPU support!")
        endif()

        #== check if MAGMA-enabled GauXC is found, if MAGMA is specified ==#
        get_property(_magma TARGET gauxc::gauxc PROPERTY GAUXC_HAS_MAGMA)
        
        if(${gauxc_ENABLE_MAGMA} AND (NOT ${GAUXC_ENABLE_GPU}))
           message(FATAL_ERROR "gauxc_ENABLE_MAGMA turned on, but gauxc_ENABLE_GPU is turned off! Magma only works with device-enabled GauXC instances.")
        endif()
        
        if(${gauxc_ENABLE_MAGMA} AND (NOT ${_magma}))
           message(FATAL_ERROR "gauxc_ENABLE_MAGMA turned on, but selected GauXC install is not built for MAGMA support!")
        elseif((NOT ${gauxc_ENABLE_MAGMA}) AND ${_magma})
           message(FATAL_ERROR "gauxc_ENABLE_MAGMA turned off, but selected GauXC install is built for MAGMA support!")
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
            #URL https://github.com/wavefunction91/GauXC/archive/issue_56_58/reduce_shell_pair_mem.tar.gz 
            URL https://github.com/wavefunction91/GauXC/archive/7746c2c.tar.gz 
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
                   -DCMAKE_CUDA_ARCHITECTURES=${CMAKE_CUDA_ARCHITECTURES}
                   -DCMAKE_POSITION_INDEPENDENT_CODE=${BUILD_FPIC}
                   -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                   -DGAUXC_ENABLE_MPI=OFF
                   -DGAUXC_ENABLE_TESTS=OFF
                   -DGAUXC_ENABLE_HDF5=OFF
                   -DGAUXC_ENABLE_MAGMA=${gauxc_ENABLE_MAGMA}
            CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                         -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                         -DCMAKE_C_COMPILER_ARG1:STRING=${CMAKE_C_COMPILER_ARG1}
                         -DCMAKE_CXX_COMPILER_ARG1:STRING=${CMAKE_CXX_COMPILER_ARG1})
        set(gauxc_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/gauxc CACHE PATH "path to internally built gauxcConfig.cmake" FORCE)
    endif()
else()
    add_library(gauxc_external INTERFACE)  # dummy
endif()
