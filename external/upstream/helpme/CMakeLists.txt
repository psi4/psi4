if(${ENABLE_helpme})
    find_package(helpme 0.1 CONFIG QUIET)

    if(${helpme_FOUND})
        get_property(_loc TARGET helpme::helpmecore PROPERTY LOCATION)
        message(STATUS "${Cyan}Found helpme${ColourReset}: ${_loc} (found version ${helpme_VERSION})")
        add_library(helpme_external INTERFACE)  # dummy

        # reset below necessary as find_package clears it. better solution sought
        #set(TargetHDF5_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/TargetHDF5 CACHE PATH "path to externally detected pylibefpConfig.cmake" FORCE)
    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_helpme})
            message(FATAL_ERROR "Suitable helpme could not be externally located as user insists")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable helpme could not be located, ${Magenta}Building helpme${ColourReset} instead.")

        ExternalProject_Add(helpme_external
            DEPENDS lapack_external
                    pybind11_external
            GIT_REPOSITORY https://github.com/andysim/helpme
            GIT_TAG master #v0.4
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}  # unused, but needs working compiler
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_Fortran_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_CXX_STANDARD=${psi4_CXX_STANDARD}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
                       -DCMAKE_INSTALL_DATADIR=${CMAKE_INSTALL_DATADIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DPYMOD_INSTALL_LIBDIR=${PYMOD_INSTALL_LIBDIR}
                       -DENABLE_OpenMP=${ENABLE_OpenMP}
                       -DENABLE_XHOST=${ENABLE_XHOST}
                       #-DENABLE_GENERIC=${ENABLE_GENERIC}
                       #-DLIBC_INTERJECT=${LIBC_INTERJECT}
                       -DPYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}
                       -DPYTHON_INCLUDE_DIR=${PYTHON_INCLUDE_DIR}
                       -DPYTHON_LIBRARY=${PYTHON_LIBRARY}
                       -DTargetLAPACK_DIR=${STAGED_INSTALL_PREFIX}/share/cmake/TargetLAPACK
                       -DFFTW_ROOT=${FFTW_ROOT}
                       -Dpybind11_DIR=${pybind11_DIR}
            CMAKE_CACHE_ARGS -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                             -DCMAKE_Fortran_FLAGS:STRING=${CMAKE_Fortran_FLAGS})

        set(helpme_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/helpme CACHE PATH "path to internally built helpmeConfig.cmake" FORCE)
    endif()
else()
    add_library(helpme_external INTERFACE)  # dummy
endif()

