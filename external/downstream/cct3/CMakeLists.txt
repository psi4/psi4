if(${ENABLE_cct3})
    find_package(cct3 0.1 CONFIG QUIET)

    if(${cct3_FOUND})
        get_property(_loc TARGET cct3::cct3 PROPERTY LOCATION)
        message(STATUS "${Cyan}Found cct3${ColourReset}: ${_loc} (found version ${cct3_VERSION})")
        add_library(cct3_external INTERFACE)  # dummy

    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_cct3})
            message(FATAL_ERROR "Suitable cct3 could not be externally located as user insists")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable cct3 could not be located, ${Magenta}Building cct3${ColourReset} instead.")

        ExternalProject_Add(cct3_external
            DEPENDS psi4-core
            URL https://github.com/piecuch-group/cct3/archive/bff480f.tar.gz
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
            CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                             -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                             -DCMAKE_Fortran_FLAGS:STRING=${CMAKE_Fortran_FLAGS}
                             -DCMAKE_PREFIX_PATH:PATH=${STAGED_INSTALL_PREFIX})

        set(cct3_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/cct3 CACHE PATH "path to internally built cct3Config.cmake" FORCE)
    endif()
else()
    add_library(cct3_external INTERFACE)  # dummy
endif()
