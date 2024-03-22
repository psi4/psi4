if(${ENABLE_ExchCXX})
    #find_package(ExchCXX 0.0.1 CONFIG PATHS "${GAUXC_DIR}/../ExchCXX") 
    find_package(ExchCXX 0.0.1 CONFIG) 

    if(${ExchCXX_FOUND})
        get_property(_loc TARGET ExchCXX::ExchCXX PROPERTY IMPORTED_LOCATION_RELWITHDEBINFO)
        message(STATUS "${Cyan}Found ExchCXX ${ColourReset}: ${_loc} (found version ${ExchCXX_VERSION})")
        add_library(exchcxx_external INTERFACE)  # dummy

        if(${CMAKE_INSIST_FIND_PACKAGE_ExchCXX})
            message(VERBOSE "Suitable ExchCXX located externally as user insists")
        endif()

    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_ExchCXX})
            message(FATAL_ERROR "Suitable ExchCXX could not be externally located as user insists")
        endif()

        message(STATUS "Suitable ExchCXX could not be located, ${Magenta}Building ExchCXX${ColourReset} instead.")

        include(ExternalProject)
        ExternalProject_Add(exchcxx_external
            GIT_REPOSITORY https://github.com/wavefunction91/ExchCXX.git
            GIT_TAG master
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                       -DCMAKE_POSITION_INDEPENDENT_CODE=${BUILD_FPIC}
                       -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                       -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
                       -DGAUXC_ENABLE_CUDA=${gauxc_ENABLE_GPU}
                       -DGAUXC_ENABLE_MPI=OFF
                       -DGAUXC_ENABLE_TESTS=OFF
                       -DGAUXC_ENABLE_HDF5=OFF
            CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                             -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                             -DCMAKE_C_COMPILER_ARG1:STRING=${CMAKE_C_COMPILER_ARG1}
                             -DCMAKE_CXX_COMPILER_ARG1:STRING=${CMAKE_CXX_COMPILER_ARG1})

        set(ExchCXX_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/ExchCXX CACHE PATH "path to internally built ExchCXXConfig.cmake" FORCE)
    endif()
else()
    add_library(exchcxx_external INTERFACE)  # dummy
endif()
