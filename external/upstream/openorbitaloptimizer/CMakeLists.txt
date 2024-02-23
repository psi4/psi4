if(${ENABLE_OpenOrbitalOptimizer})
    find_package(OpenOrbitalOptimizer 0.1 CONFIG)

    if(TARGET OpenOrbitalOptimizer::OpenOrbitalOptimizer)
        get_property(_loc TARGET OpenOrbitalOptimizer::OpenOrbitalOptimizer PROPERTY LOCATION)
        message(STATUS "${Cyan}Found OpenOrbitalOptimizer${ColourReset}: ${_loc} (found version ${OpenOrbitalOptimizer_VERSION})")
        add_library(openorbitaloptimizer_external INTERFACE)  # dummy

        if(${CMAKE_INSIST_FIND_PACKAGE_OpenOrbitalOptimizer})
            message(VERBOSE "Suitable OpenOrbitalOptimizer located externally as user insists")
        endif()

    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_OpenOrbitalOptimizer})
            message(FATAL_ERROR "Suitable OpenOrbitalOptimizer could not be externally located as user insists")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable OpenOrbitalOptimizer could not be located, ${Magenta}Building OpenOrbitalOptimizer${ColourReset} instead.")

        set(_ooo_dir "share/cmake/OpenOrbitalOptimizer")

        ExternalProject_Add(openorbitaloptimizer_external
            DEPENDS lapack_external
            #URL https://github.com/SusiLehtola/OpenOrbitalOptimizer/archive/master.tar.gz
            URL https://github.com/loriab/OpenOrbitalOptimizer/archive/patchup.tar.gz
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DOpenOrbitalOptimizer_INSTALL_CMAKEDIR=${_ooo_dir}
                       -DOpenOrbitalOptimizer_BUILD_TESTING=OFF
                       -DTargetLAPACK_DIR=${TargetLAPACK_DIR}
                       -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                       -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
            CMAKE_CACHE_ARGS -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                             -DCMAKE_CXX_COMPILER_ARG1:STRING=${CMAKE_CXX_COMPILER_ARG1})

        set(OpenOrbitalOptimizer_DIR ${STAGED_INSTALL_PREFIX}/${_ooo_dir} CACHE PATH "path to internally built OpenOrbitalOptimizerConfig.cmake" FORCE)
    endif()
else()
    add_library(openorbitaloptimizer_external INTERFACE)  # dummy
endif()
