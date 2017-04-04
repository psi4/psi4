if(${ENABLE_PCMSolver})
    find_package(PCMSolver 1.1.10 CONFIG QUIET)

    if(${PCMSolver_FOUND})
        get_property(_loc TARGET PCMSolver::pcm PROPERTY LOCATION)
        message(STATUS "${Cyan}Found PCMSolver${ColourReset}: ${_loc} (found version ${PCMSolver_VERSION})")
        add_library(pcmsolver_external INTERFACE)  # dummy
    else()
        include(ExternalProject)
        message(STATUS "Suitable PCMSolver could not be located, ${Magenta}Building PCMSolver${ColourReset} instead.")

        if(${BUILD_SHARED_LIBS})
            set(_a_only  OFF)
            set(_so_only ON)
        else()
            set(_a_only  ON)
            set(_so_only OFF)
        endif()

        ExternalProject_Add(pcmsolver_external
            GIT_REPOSITORY https://github.com/PCMSolver/pcmsolver
            GIT_TAG v1.1.10
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
                       -DPYMOD_INSTALL_LIBDIR=${PYMOD_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_DATADIR=${CMAKE_INSTALL_DATADIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DPYTHON_INTERPRETER=${PYTHON_EXECUTABLE}
                       -DSTATIC_LIBRARY_ONLY=${_a_only}
                       -DSHARED_LIBRARY_ONLY=${_so_only}
                       -DENABLE_OPENMP=${ENABLE_OPENMP}  # relevant for thread safety
                       #-DENABLE_XHOST=${ENABLE_XHOST}  # option missing
                       # always fpic'd
                       -DENABLE_GENERIC=${ENABLE_GENERIC}
                       -DENABLE_DOCS=OFF
                       -DENABLE_TESTS=OFF
                       -DENABLE_TIMER=OFF
                       -DENABLE_LOGGER=OFF
                       -DBUILD_STANDALONE=OFF
                       -DENABLE_FORTRAN_API=OFF
                       -DENABLE_CXX11_SUPPORT=ON
                       -DENABLE_64BIT_INTEGERS=OFF
            CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                             -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                             -DCMAKE_Fortran_FLAGS:STRING=${CMAKE_Fortran_FLAGS}
            INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
            DESTDIR=${CMAKE_BINARY_DIR}/stage)

        set(PCMSolver_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/PCMSolver CACHE PATH "path to internally built PCMSolverConfig.cmake" FORCE)
    endif()
else()
    add_library(pcmsolver_external INTERFACE)  # dummy
endif()
