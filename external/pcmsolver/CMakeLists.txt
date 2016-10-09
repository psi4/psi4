if(${ENABLE_PCMSolver})
    find_package(PCMSolver QUIET)  # todo CONFIG

    if(${PCMSolver_FOUND})
        #get_property(_loc TARGET PCMSolver::pcm PROPERTY LOCATION)  # todo
        get_property(_loc TARGET PCMSolver::PCMSolver PROPERTY LOCATION)
        message(STATUS "Found PCMSolver: ${_loc} (found version ${PCMSolver_VERSION})")
        add_library(pcmsolver_external INTERFACE)  # dummy
    else()
        include(ExternalProject)
        message(STATUS "Suitable PCMSolver could not be located, building one instead.")

        if(${BUILD_SHARED_LIBS})
#            set(_a_only  OFF)
#            set(_so_only ON)
        else()
            set(_a_only  ON)
#            set(_so_only OFF)
        endif()

        ExternalProject_Add(pcmsolver_external
            GIT_REPOSITORY https://github.com/PCMSolver/pcmsolver
            GIT_TAG v1.1.4  # TODO tag when stable
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/external/pcmsolver
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
                       -DCMAKE_INSTALL_DATADIR=${CMAKE_INSTALL_DATADIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DSTATIC_LIBRARY_ONLY=${_a_only}
#                       -DENABLE_OPENMP=???
                       #-DENABLE_XHOST=???  # option missing
                       #-DENABLE_XHOST=${ENABLE_XHOST}  # option missing
                       #-DBUILD_FPIC=${BUILD_FPIC}  # always fpic'd
                       -DENABLE_GENERIC=${ENABLE_GENERIC}


#                  -DENABLE_CXX11_SUPPORT=ON
#                  -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
#                  -DENABLE_TESTS=OFF
#                  -DENABLE_LOGGER=OFF
#                  -DENABLE_TIMER=OFF
#                  -DBUILD_STANDALONE=OFF
#                  -DENABLE_FORTRAN_API=OFF
#                  -DSTATIC_LIBRARY_ONLY=${ENABLE_STATIC}
#                  -DPYTHON_INTERPRETER=${PYTHON_EXECUTABLE}

            CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                             -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                             -DCMAKE_Fortran_FLAGS:STRING=${CMAKE_Fortran_FLAGS}
            INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
            DESTDIR=${CMAKE_BINARY_DIR}/stage)
#    configure_file(${CMAKE_SOURCE_DIR}/psi4/share/python/pcm_placeholder.py.in
#       ${CMAKE_SOURCE_DIR}/share/python/pcm_placeholder.py
#    )
#    add_dependencies(pcmsolver_external boost_external)
    endif()
else()
    add_library(pcmsolver_external INTERFACE)  # dummy
endif()

####



#                  -DBOOST_ROOT=${BOOST_ROOT}

