if(${ENABLE_libecpint})
    find_package(libecpint 0.9 CONFIG QUIET)

    if(${libecpint_FOUND})
        get_property(_loc TARGET libecpint::ecpint PROPERTY LOCATION)
        message(STATUS "${Cyan}Found libecpint${ColourReset}: ${_loc} (found version ${Lib_VERSION})")
        add_library(libecpint_external INTERFACE)  # dummy
    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_libecpint})
            message(FATAL_ERROR "Suitable Libecpint could not be externally located as user insists")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable Libecpint could not be located, ${Magenta}Building Libecpint${ColourReset} instead.")

        if(${BUILD_SHARED_LIBS})
            set(_a_only  OFF)
            set(_so_only ON)
        else()
            set(_a_only  ON)
            set(_so_only OFF)
        endif()

        ExternalProject_Add(libecpint_external
            GIT_REPOSITORY https://github.com/andysim/libecpint
            GIT_TAG 609b52d665b6595bb2bcad61b58809e7c90b37a7
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
                       -DPYMOD_INSTALL_LIBDIR=${PYMOD_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_DATADIR=${CMAKE_INSTALL_DATADIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DPYTHON_INTERPRETER=${PYTHON_EXECUTABLE}
                       -DSTATIC_LIBRARY_ONLY=${_a_only}
                       -DSHARED_LIBRARY_ONLY=${_so_only}
                       -DENABLE_OPENMP=OFF
                       #-DENABLE_XHOST=${ENABLE_XHOST}  # option missing
                       # always fpic'd
                       -DENABLE_GENERIC=${ENABLE_GENERIC}
                       -DENABLE_TESTS=OFF
                       -DENABLE_TIMER=OFF
                       -DENABLE_LOGGER=OFF
                       -DBUILD_STANDALONE=OFF
                       -DENABLE_CXX11_SUPPORT=ON
                       -DENABLE_64BIT_INTEGERS=OFF
            CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                             -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                             -DCMAKE_Fortran_FLAGS:STRING=${CMAKE_Fortran_FLAGS})

        set(libecpint_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/libecpint CACHE PATH "path to internally built libecpint.cmake" FORCE)
    endif()
else()
    add_library(libecpint_external INTERFACE)  # dummy
endif()
