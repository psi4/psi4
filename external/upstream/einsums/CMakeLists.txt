if(${ENABLE_Einsums})
    find_package(Einsums 0.1 CONFIG QUIET)

    if(TARGET Einsums::einsums)
        get_property(_loc TARGET Einsums::einsums PROPERTY LOCATION)
        message(STATUS "${Cyan}Found Einsums${ColourReset}: ${_loc} (found version ${Einsums_VERSION})")
        add_library(einsums_external INTERFACE)  # dummy

        # reset below necessary as find_package clears it. better solution sought
        set(TargetHDF5_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/TargetHDF5 CACHE PATH "path to externally detected TargetHDF5Config.cmake" FORCE)
    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_Einsums})
            message(FATAL_ERROR "Suitable Einsums could not be externally located as user insists")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable Einsums could not be located, ${Magenta}Building Einsums${ColourReset} instead.")

        if(${BUILD_SHARED_LIBS})
            set(_a_only  OFF)
        else()
            set(_a_only  ON)
        endif()

        set(_einsums_dir "share/cmake/Einsums")

        ExternalProject_Add(einsums_external
            DEPENDS lapack_external
                    hdf5_external
            URL https://github.com/Einsums/Einsums/archive/main.tar.gz
            URL /psi/gits/einsums
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
#                       -DENABLE_XHOST=${ENABLE_XHOST}
                       -DEINSUMS_ENABLE_TESTING=OFF
                       -DEINSUMS_USE_HPTT=ON
#                       -DOpenMP_LIBRARY_DIRS=${OpenMP_LIBRARY_DIRS}
                       -DTargetLAPACK_DIR=${TargetLAPACK_DIR}
                       -DTargetHDF5_DIR=${TargetHDF5_DIR}
                       -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                       -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
            CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                             -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                             -DTargetOpenMP_FIND_COMPONENTS:STRING=C;CXX)

        set(Einsums_DIR ${STAGED_INSTALL_PREFIX}/${_einsums_dir} CACHE PATH "path to internally built EinsumsConfig.cmake" FORCE)
    endif()
else()
    add_library(einsums_external INTERFACE)  # dummy
endif()
