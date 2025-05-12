if(${ENABLE_v2rdm_casscf})
    find_package(v2rdm_casscf 0.9 CONFIG QUIET)

    if(${v2rdm_casscf_FOUND})
        get_property(_loc TARGET v2rdm_casscf::v2rdm_casscf PROPERTY LOCATION)
        message(STATUS "${Cyan}Found v2rdm_casscf${ColourReset}: ${_loc} (found version ${v2rdm_casscf_VERSION})")
        add_library(v2rdm_casscf_external INTERFACE)  # dummy

    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_v2rdm_casscf})
            message(FATAL_ERROR "Suitable v2rdm_casscf could not be externally located as user insists")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable v2rdm_casscf could not be located, ${Magenta}Building v2rdm_casscf${ColourReset} instead.")

        ExternalProject_Add(v2rdm_casscf_external
            DEPENDS psi4-core
            URL https://github.com/loriab/v2rdm_casscf/archive/c2f435f.tar.gz
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DPYMOD_INSTALL_LIBDIR=${PYMOD_INSTALL_LIBDIR}
                       -DENABLE_XHOST=${ENABLE_XHOST}
                       -DENABLE_GENERIC=${ENABLE_GENERIC}
                       -DOpenMP_LIBRARY_DIRS=${OpenMP_LIBRARY_DIRS}
                       -Dpsi4_DIR=${STAGED_INSTALL_PREFIX}/share/cmake/psi4
                       -Dpybind11_DIR=${pybind11_DIR}
                       -DTargetLAPACK_DIR=${TargetLAPACK_DIR}
                       -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                       -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
            CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                             -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                             -DCMAKE_Fortran_FLAGS:STRING=${CMAKE_Fortran_FLAGS}
                             -DCMAKE_C_COMPILER_ARG1:STRING=${CMAKE_C_COMPILER_ARG1}
                             -DCMAKE_CXX_COMPILER_ARG1:STRING=${CMAKE_CXX_COMPILER_ARG1}
                             -DCMAKE_Fortran_COMPILER_ARG1:STRING=${CMAKE_Fortran_COMPILER_ARG1}
                             -DCMAKE_PREFIX_PATH:PATH="${STAGED_INSTALL_PREFIX};${CMAKE_PREFIX_PATH}")

        set(v2rdm_casscf_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/v2rdm_casscf CACHE PATH "path to internally built v2rdm_casscfConfig.cmake" FORCE)
    endif()
else()
    add_library(v2rdm_casscf_external INTERFACE)  # dummy
endif()

