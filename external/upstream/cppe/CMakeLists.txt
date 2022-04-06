if(${ENABLE_cppe})
    if(NOT (${CMAKE_DISABLE_FIND_PACKAGE_cppe}))
        include(FindPythonModule)
        find_python_module(cppe ATLEAST 0.3.1 QUIET)
    endif()

    if(${cppe_FOUND})
        message(STATUS "${Cyan}Found cppe${ColourReset}: ${PY_cppe} (found version ${cppe_VERSION})")
        add_library(cppe_external INTERFACE)  # dummy
    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_cppe})
            message(FATAL_ERROR "Suitable cppe could not be externally located as user insists")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable cppe could not be located, ${Magenta}Building cppe${ColourReset} instead.")

        ExternalProject_Add(cppe_external
            DEPENDS pybind11_external
            URL https://github.com/maxscheurer/cppe/archive/v0.3.1.tar.gz
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                        # runtime segfault if Debug passed from Psi4
                       -DCMAKE_BUILD_TYPE=Release
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
                       -DCMAKE_INSTALL_DATADIR=${CMAKE_INSTALL_DATADIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DPYMOD_INSTALL_LIBDIR=${PYMOD_INSTALL_LIBDIR}
                       -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                       -DENABLE_OPENMP=OFF
                       -DENABLE_PYTHON_INTERFACE=ON
                       -DPython_EXECUTABLE=${Python_EXECUTABLE}
                       -DPython_INCLUDE_DIR=${Python_INCLUDE_DIRS}
                       -DPython_LIBRARY=${Python_LIBRARIES}
                       -DBUILD_FPIC=${BUILD_FPIC}
                       -DENABLE_GENERIC=${ENABLE_GENERIC}
                       -DLIBC_INTERJECT=${LIBC_INTERJECT}
                       -DINSTALL_DEVEL_HEADERS=OFF
                       -DOpenMP_LIBRARY_DIRS=${OpenMP_LIBRARY_DIRS}
                       -Dpybind11_DIR=${pybind11_DIR}
                       -DCMAKE_INSIST_FIND_PACKAGE_pybind11=ON
                       -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                       -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
            CMAKE_CACHE_ARGS -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                             -DTargetOpenMP_FIND_COMPONENTS:STRING=CXX)

    endif()
else()
    add_library(cppe_external INTERFACE)  # dummy
endif()
