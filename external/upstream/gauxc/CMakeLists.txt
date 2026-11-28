if(${ENABLE_gauxc})
    if(NOT "${gauxc_DIR}" STREQUAL "")
        # a GauXC built standalone installs its ExchCXX/IntegratorXX alongside; default to those
        #   (find_package ignores a <Pkg>_DIR holding no config, so a wrong guess just searches)
        if("${IntegratorXX_DIR}" STREQUAL "")
            set(IntegratorXX_DIR ${gauxc_DIR}/../IntegratorXX)
        endif()
        if("${ExchCXX_DIR}" STREQUAL "")
            set(ExchCXX_DIR ${gauxc_DIR}/../ExchCXX)
        endif()
    endif()

    get_property(BLAS_LIBRARIES TARGET lapk PROPERTY INTERFACE_LINK_LIBRARIES)

    find_package(gauxc CONFIG)  # edit in codedeps

    if(${gauxc_FOUND})
        #== GauXC variant (cpu|gpu) must be the one gauxc_ENABLE_GPU asks for ==#
        get_property(_gpu TARGET gauxc::gauxc PROPERTY GAUXC_HAS_DEVICE)

        # GauXC v1.1 exports GAUXC_HAS_DEVICE as TRUE or FALSE; its config version still reads 1.0.0, so check this instead
        if("${_gpu}" STREQUAL "")
            message(FATAL_ERROR "The GauXC at ${gauxc_DIR} predates v1.1 (its gauxc::gauxc target lacks GAUXC_HAS_DEVICE). "
                                "Psi4 requires GauXC v1.1 or newer.")
        elseif(gauxc_ENABLE_GPU AND NOT _gpu)
            message(FATAL_ERROR "gauxc_ENABLE_GPU=ON, but the GauXC at ${gauxc_DIR} is built without GPU support. "
                                "Point gauxc_DIR at a GPU-enabled GauXC, or set gauxc_ENABLE_GPU=OFF.")
        elseif(NOT gauxc_ENABLE_GPU AND _gpu)
            message(FATAL_ERROR "gauxc_ENABLE_GPU=OFF, but the GauXC at ${gauxc_DIR} is built with GPU support. "
                                "Set gauxc_ENABLE_GPU=ON (plus CUDA hints like CUDAToolkit_ROOT as needed), or point gauxc_DIR at a CPU-only GauXC.")
        endif()

        #== check if MAGMA-enabled GauXC is found, if MAGMA is specified ==#
        get_property(_magma TARGET gauxc::gauxc PROPERTY GAUXC_HAS_MAGMA)
        
        if(${gauxc_ENABLE_MAGMA} AND (NOT ${gauxc_ENABLE_GPU}))
           message(FATAL_ERROR "gauxc_ENABLE_MAGMA turned on, but gauxc_ENABLE_GPU is turned off! Magma only works with device-enabled GauXC instances.")
        endif()
        
        if(${gauxc_ENABLE_MAGMA} AND (NOT ${_magma}))
           message(FATAL_ERROR "gauxc_ENABLE_MAGMA turned on, but selected GauXC install is not built for Magma support!")
        elseif((NOT ${gauxc_ENABLE_MAGMA}) AND ${_magma})
           message(WARNING "gauxc_ENABLE_MAGMA turned off, but selected GauXC install is built for Magma support! If you want Magma support for the Psi4/GauXC interface, make sure to explicitly turn on gauxc_ENABLE_MAGMA!")
        endif()

        #== normal print-out and library addition ==#
        get_property(_loc TARGET gauxc::gauxc PROPERTY LOCATION)
        message(STATUS "${Cyan}Found gauxc ${ColourReset}: ${_loc} (found version ${gauxc_VERSION})")
        add_library(gauxc_external INTERFACE)  # dummy

        if(${CMAKE_INSIST_FIND_PACKAGE_gauxc})
            message(VERBOSE "Suitable gauxc located externally as user insists")
        endif()

    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_gauxc})
            message(FATAL_ERROR "Suitable gauxc could not be externally located as user insists. ${gauxc_NOT_FOUND_MESSAGE}\n"
                                "A GPU-enabled GauXC also needs its CUDA dependencies found: set gauxc_ENABLE_GPU=ON and, as needed, CUDAToolkit_ROOT.")
        endif()

        set(CMAKE_CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}" CACHE STRING "CUDA architectures")

        if(${gauxc_ENABLE_GPU} AND (NOT DEFINED CMAKE_CUDA_ARCHITECTURES OR "${CMAKE_CUDA_ARCHITECTURES}" STREQUAL ""))
           message(FATAL_ERROR "Internal GPU-enabled GauXC construction specified, but CMAKE_CUDA_ARCHITECTURES is not defined!")
        endif()

        message(STATUS "Suitable gauxc could not be located, ${Magenta}Building gauxc${ColourReset} instead.")

        include(ExternalProject)
        ExternalProject_Add(gauxc_external
            DEPENDS lapack_external 
                    libxc_external
                    gau2grid_external
            URL https://github.com/wavefunction91/GauXC/archive/d9cdad3.tar.gz  # edit in codedeps
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
                   -DLibxc_DIR=${Libxc_DIR}
                   -Dgau2grid_DIR=${gau2grid_DIR}
                   -DCMAKE_POSITION_INDEPENDENT_CODE=${BUILD_FPIC}
                   -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                   -DOpenMP_C_FLAGS:STRING=${OpenMP_C_FLAGS}
                   -DOpenMP_CXX_FLAGS:STRING=${OpenMP_CXX_FLAGS}
                   -DOpenMP_omp_LIBRARY=${OpenMP_omp_LIBRARY}
                   -DOpenMP_gomp_LIBRARY=${OpenMP_gomp_LIBRARY}
                   -DOpenMP_iomp5_LIBRARY=${OpenMP_iomp5_LIBRARY}
                   -DOpenMP_pthread_LIBRARY=${OpenMP_pthread_LIBRARY}
                   -DGAUXC_ENABLE_MPI=OFF
                   -DGAUXC_ENABLE_TESTS=OFF
                   -DGAUXC_ENABLE_HDF5=OFF
                   -DGAUXC_ENABLE_MAGMA=${gauxc_ENABLE_MAGMA}
                   -DCMAKE_POLICY_VERSION_MINIMUM=3.5
            CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                         -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                         -DCMAKE_C_COMPILER_ARG1:STRING=${CMAKE_C_COMPILER_ARG1}
                         -DCMAKE_CXX_COMPILER_ARG1:STRING=${CMAKE_CXX_COMPILER_ARG1}
                         -DOpenMP_C_LIB_NAMES:STRING=${OpenMP_C_LIB_NAMES}
                         -DOpenMP_CXX_LIB_NAMES:STRING=${OpenMP_CXX_LIB_NAMES}
                         -DBLAS_LIBRARIES:STRING=${BLAS_LIBRARIES}
                         -DCMAKE_CUDA_ARCHITECTURES:STRING=${CMAKE_CUDA_ARCHITECTURES}
                         # CMAKE_PREFIX_PATH propagates to GauXC's FetchContent builds (ExchCXX)
                         # allowing ExchCXX to find Psi4's LibXC without rebuilding from source
                         -DCMAKE_PREFIX_PATH:STRING=${STAGED_INSTALL_PREFIX};${CMAKE_PREFIX_PATH})
        set(gauxc_DIR ${STAGED_INSTALL_PREFIX}/lib/cmake/gauxc CACHE PATH "path to internally built gauxcConfig.cmake" FORCE)
        # GauXC's FetchContent builds of its dependencies stage alongside it
        set(ExchCXX_DIR ${STAGED_INSTALL_PREFIX}/lib/cmake/ExchCXX CACHE PATH "path to internally built ExchCXXConfig.cmake" FORCE)
        if("${IntegratorXX_DIR}" STREQUAL "")
            set(IntegratorXX_DIR ${STAGED_INSTALL_PREFIX}/lib/cmake/IntegratorXX CACHE PATH "path to internally built IntegratorXXConfig.cmake" FORCE)
        endif()
    endif()
else()
    message(STATUS "Optional addon ${Blue}unsought: gauxc${ColourReset}")
    add_library(gauxc_external INTERFACE)  # dummy
endif()
