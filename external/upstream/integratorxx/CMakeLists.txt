
if(${ENABLE_IntegratorXX})
    find_package(IntegratorXX CONFIG)

    if(${IntegratorXX_FOUND})
        get_property(_loc TARGET IntegratorXX::IntegratorXX PROPERTY LOCATION)
        message(STATUS "${Cyan}Found IntegratorXX${ColourReset}: ${_loc} (found version ${IntegratorXX_VERSION})")
        add_library(integratorxx_external INTERFACE)  # dummy
    else()

        if(${CMAKE_INSIST_FIND_PACKAGE_IntegratorXX})
            message(FATAL_ERROR "Suitable IntegratorXX could not be externally located as user insists")
        endif()
    
        include(ExternalProject)
        message(STATUS "Suitable IntegratorXX could not be located, ${Magenta}Building IntegratorXX${ColourReset} instead.")
        ExternalProject_Add(integratorxx_external
            URL https://github.com/wavefunction91/IntegratorXX/archive/master.tar.gz
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DCMAKE_CXX_STANDARD=${psi4_CXX_STANDARD}
                       -DINTEGRATORXX_ENABLE_TESTS=OFF
                       -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                       -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
            CMAKE_CACHE_ARGS -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                             -DCMAKE_CXX_COMPILER_ARG1:STRING=${CMAKE_CXX_COMPILER_ARG1})
    
        set(IntegratorXX_DIR ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/cmake/IntegratorXX CACHE PATH "path to internally built IntegratorXXConfig.cmake" FORCE)
    endif()
else()
    add_library(integratorxx_external INTERFACE)  # dummy
endif()

