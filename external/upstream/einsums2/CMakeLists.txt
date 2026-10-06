if(${ENABLE_Einsums2})
    if(einsums2_DIR AND NOT Einsums2_DIR)
        set(Einsums2_DIR ${einsums2_DIR})
    endif()

    if(Einsums2_DIR)
        set(Einsums_DIR ${Einsums2_DIR})
    endif()

    find_package(Einsums 2.0.0 CONFIG)  # edit in codedeps

    if(TARGET Einsums::Einsums)
        get_property(_loc TARGET Einsums::Einsums PROPERTY LOCATION)
        message(STATUS "${Cyan}Found Einsums2${ColourReset}: ${_loc} (found version ${Einsums_VERSION})")
        add_library(einsums2_external INTERFACE)  # dummy

        # reset below necessary as find_package clears it. better solution sought
        set(TargetHDF5_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/TargetHDF5 CACHE PATH "path to externally detected TargetHDF5Config.cmake" FORCE)
        set(Einsums2_DIR ${Einsums_DIR} CACHE PATH "path to externally detected Einsums v2 Config.cmake" FORCE)

        if(${CMAKE_INSIST_FIND_PACKAGE_Einsums2})
            message(VERBOSE "Suitable Einsums2 located externally as user insists")
        endif()

    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_Einsums2})
            message(FATAL_ERROR "Suitable Einsums2 could not be externally located as user insists")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable Einsums2 could not be located, ${Magenta}Building Einsums2${ColourReset} instead.")

        if(${BUILD_SHARED_LIBS})
            set(_a_only  OFF)
        else()
            set(_a_only  ON)
        endif()

        set(_einsums_dir "share/cmake/Einsums")

        ExternalProject_Add(einsums2_external
            DEPENDS lapack_external
                    hdf5_external
            URL https://github.com/Einsums/Einsums/archive/79305cc.tar.gz
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DEINSUMS_STATIC_BUILD=${_a_only}
                       -DEXTRA_C_FLAGS=${CMAKE_C_FLAGS}
                       -DEXTRA_CXX_FLAGS=${CMAKE_CXX_FLAGS}
                       -DEINSUMS_INSTALL_CMAKEDIR=${_einsums_dir}
#                      -DENABLE_XHOST=${ENABLE_XHOST}
                       -DEINSUMS_ENABLE_TESTING=OFF
                       -DEINSUMS_USE_HPTT=ON
#                      -DOpenMP_LIBRARY_DIRS=${OpenMP_LIBRARY_DIRS}
                       -DTargetLAPACK_DIR=${TargetLAPACK_DIR}
                       -DTargetHDF5_DIR=${TargetHDF5_DIR}
                       -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                       -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
                       -DEINSUMS_WITH_TESTS=OFF
                       -DEINSUMS_WITH_DOCUMENTATION=OFF
                       -DEINSUMS_BUILD_PYTHON=ON
                       -DEINSUMS_INSTALL_PYMODDIR=${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}
                       -DEINSUMS_WITH_EXAMPLES=OFF
            CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                             -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                             -DCMAKE_C_COMPILER_ARG1:STRING=${CMAKE_C_COMPILER_ARG1}
                             -DCMAKE_CXX_COMPILER_ARG1:STRING=${CMAKE_CXX_COMPILER_ARG1}
                             -DTargetOpenMP_FIND_COMPONENTS:STRING=C;CXX)

        set(Einsums2_DIR ${STAGED_INSTALL_PREFIX}/${_einsums_dir} CACHE PATH "path to internally built Einsums v2 Config.cmake" FORCE)
    endif()
else()
    message(STATUS "Optional addon ${Blue}unsought: einsums2${ColourReset}")
    add_library(einsums2_external INTERFACE)  # dummy
endif()
