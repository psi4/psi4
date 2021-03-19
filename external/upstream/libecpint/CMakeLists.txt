find_package(ECPINT 1.0.5 CONFIG QUIET)

if(${LibECPInt_FOUND})
    get_property(_loc TARGET ECPINT::ecpint PROPERTY LOCATION)
    message(STATUS "${Cyan}Found LibECPInt${ColourReset}: ${_loc} (found version ${ECPINT_VERSION})")
    add_library(libecpint_external INTERFACE)  # dummy
else()
    include(ExternalProject)
    message(STATUS "Suitable LibECPInt could not be located, ${Magenta}Building LipECPInt${ColourReset} instead.")

    if(${BUILD_SHARED_LIBS})
        set(_a_only  OFF)
        set(_so_only ON)
    else()
        set(_a_only  ON)
        set(_so_only OFF)
    endif()

    ExternalProject_Add(libecpint_external
        URL https://github.com/robashaw/libecpint/archive/v1.0.5.tar.gz
        UPDATE_COMMAND ""
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                   -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                   -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                   -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
                   -DPYMOD_INSTALL_LIBDIR=${PYMOD_INSTALL_LIBDIR}
                   -DCMAKE_INSTALL_DATADIR=${CMAKE_INSTALL_DATADIR}
                   -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                   -DPYTHON_INTERPRETER=${Python_EXECUTABLE}
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
                   -DTEST_Fortran_API=OFF
                   -DENABLE_CXX11_SUPPORT=ON
                   -DENABLE_64BIT_INTEGERS=OFF
        CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                         -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                         -DCMAKE_Fortran_FLAGS:STRING=${CMAKE_Fortran_FLAGS})

    set(ECPINT_DIR ${STAGED_INSTALL_PREFIX}/lib/cmake/ecpint CACHE PATH "path to internally built ecpint-config.cmake" FORCE)
endif()
