if(${ENABLE_ambit})
    find_package(ambit 0.6 CONFIG QUIET)

    if(${ambit_FOUND})
        get_property(_loc TARGET ambit::ambit PROPERTY LOCATION)
        message(STATUS "${Cyan}Found ambit${ColourReset}: ${_loc} (found version ${ambit_VERSION})")
        add_library(ambit_external INTERFACE)  # dummy

        # reset below necessary as find_package clears it. better solution sought
        set(TargetHDF5_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/TargetHDF5 CACHE PATH "path to externally detected ambitConfig.cmake" FORCE)
    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_ambit})
            message(FATAL_ERROR "Suitable ambit could not be externally located as user insists")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable ambit could not be located, ${Magenta}Building ambit${ColourReset} instead.")

        if(${BUILD_SHARED_LIBS})
            set(_a_only  OFF)
            set(_so_only ON)
        else()
            set(_a_only  ON)
            set(_so_only OFF)
        endif()

        set(_ambit_dir "share/cmake/ambit")

        ExternalProject_Add(ambit_external
            DEPENDS lapack_external
                    hdf5_external
                    pybind11_external
            URL https://github.com/jturney/ambit/archive/733c529.tar.gz
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DPYMOD_INSTALL_LIBDIR=${PYMOD_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DSTATIC_ONLY=${_a_only}
                       -DSHARED_ONLY=${_so_only}
                       -DENABLE_OPENMP=${ENABLE_OPENMP}  # relevant
                       #-DENABLE_VECTORIZATION=${ENABLE_VECTORIZATION}
                       #-DEXTRA_Fortran_FLAGS=${CMAKE_EXTRA_Fortran_FLAGS}
                       -DEXTRA_C_FLAGS=${CMAKE_C_FLAGS}
                       -DEXTRA_CXX_FLAGS=${CMAKE_CXX_FLAGS}
                       -Dambit_INSTALL_CMAKEDIR=${_ambit_dir}
                       -DPython_EXECUTABLE=${Python_EXECUTABLE}
                       -DPython_INCLUDE_DIR=${Python_INCLUDE_DIRS}
                       -DPython_LIBRARY=${Python_LIBRARIES}                       
                       -DENABLE_XHOST=${ENABLE_XHOST}
                       -DBUILD_FPIC=${BUILD_FPIC}
                       -DENABLE_GENERIC=${ENABLE_GENERIC}
                       -DLIBC_INTERJECT=${LIBC_INTERJECT}
                       -DENABLE_TESTS=OFF
                       -DOpenMP_LIBRARY_DIRS=${OpenMP_LIBRARY_DIRS}
                       -DTargetLAPACK_DIR=${TargetLAPACK_DIR}
                       -DTargetHDF5_DIR=${TargetHDF5_DIR}
                       -Dpybind11_DIR=${pybind11_DIR}
                       -DCMAKE_INSIST_FIND_PACKAGE_pybind11=ON
                       -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                       -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
            CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                             -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                             -DCMAKE_C_COMPILER_ARG1:STRING=${CMAKE_C_COMPILER_ARG1}
                             -DCMAKE_CXX_COMPILER_ARG1:STRING=${CMAKE_CXX_COMPILER_ARG1}
                             -DTargetOpenMP_FIND_COMPONENTS:STRING=C;CXX)

        set(ambit_DIR ${STAGED_INSTALL_PREFIX}/${_ambit_dir} CACHE PATH "path to internally built ambitConfig.cmake" FORCE)
    endif()
else()
    add_library(ambit_external INTERFACE)  # dummy
endif()
