if(${ENABLE_libefp})
    find_package(pylibefp 0.3 CONFIG QUIET)

    if(${pylibefp_FOUND})
        get_property(_loc TARGET pylibefp::core PROPERTY LOCATION)
        message(STATUS "${Cyan}Found pylibefp${ColourReset}: ${_loc} (found version ${pylibefp_VERSION})")
        add_library(pylibefp_external INTERFACE)  # dummy

        # reset below necessary as find_package clears it. better solution sought
        #set(TargetHDF5_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/TargetHDF5 CACHE PATH "path to externally detected pylibefpConfig.cmake" FORCE)
    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_pylibefp})
            message(FATAL_ERROR "Suitable pylibefp could not be externally located as user insists")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable pylibefp could not be located, ${Magenta}Building pylibefp${ColourReset} instead.")

        ExternalProject_Add(pylibefp_external
            DEPENDS lapack_external
                    libefp_external
                    pybind11_external
            GIT_REPOSITORY https://github.com/loriab/pylibefp
            GIT_TAG 07ec7e9  # v0.3+7  # add tag when stable
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}  # unused, but needs working compiler
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_CXX_STANDARD=${psi4_CXX_STANDARD}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
                       -DCMAKE_INSTALL_DATADIR=${CMAKE_INSTALL_DATADIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DPYMOD_INSTALL_LIBDIR=${PYMOD_INSTALL_LIBDIR}
                       # OpenMP irrelevant
                       -DENABLE_XHOST=${ENABLE_XHOST}
                       -DENABLE_GENERIC=${ENABLE_GENERIC}
                       #-DLIBC_INTERJECT=${LIBC_INTERJECT}
                       -DPYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}
                       -DPYTHON_INCLUDE_DIR=${PYTHON_INCLUDE_DIR}
                       -DPYTHON_LIBRARY=${PYTHON_LIBRARY}
                       -DTargetLAPACK_DIR=${STAGED_INSTALL_PREFIX}/share/cmake/TargetLAPACK
                       -Dlibefp_DIR=${libefp_DIR}
                       -Dpybind11_DIR=${pybind11_DIR}
            CMAKE_CACHE_ARGS -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
            INSTALL_COMMAND DESTDIR=${CMAKE_BINARY_DIR}/stage ${CMAKE_MAKE_PROGRAM} install)

        set(pylibefp_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/pylibefp CACHE PATH "path to internally built pylibefpConfig.cmake" FORCE)
    endif()
else()
    add_library(pylibefp_external INTERFACE)  # dummy
endif()

