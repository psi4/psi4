if(${ENABLE_gdma})
    find_package(gdma 2.3.3 CONFIG QUIET COMPONENTS Python)

    if(${gdma_FOUND})
        get_property(_loc TARGET gdma::pygdma PROPERTY LOCATION)
        message(STATUS "${Cyan}Found gdma${ColourReset}: ${_loc} (found version ${gdma_VERSION})")
        add_library(gdma_external INTERFACE)  # dummy

        if(${CMAKE_INSIST_FIND_PACKAGE_gdma})
            message(VERBOSE "Suitable gdma located externally as user insists")
        endif()

    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_gdma})
            message(FATAL_ERROR "Suitable gdma could not be externally located as user insists")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable gdma could not be located, ${Magenta}Building gdma${ColourReset} instead.")

        set(_gdma_dir "share/cmake/gdma")

        ExternalProject_Add(gdma_external
            DEPENDS pybind11_external
            URL https://github.com/psi4/gdma/archive/v2.3.3.tar.gz  # Stone's upstream c2e0b548 plus commits from Andy Simmonett for lib, Holger Kruse for I/O, Lori Burns for CMake & Python
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                       -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DPYMOD_INSTALL_LIBDIR=${PYMOD_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                       -Dgdma_INSTALL_CMAKEDIR=${_gdma_dir}
                       -Dgdma_ENABLE_PYTHON=ON
                       -DPython_EXECUTABLE=${Python_EXECUTABLE}
                       -DENABLE_OPENMP=${ENABLE_OPENMP}  # relevant for thread safety
                       -DENABLE_XHOST=${ENABLE_XHOST}
                       -DBUILD_FPIC=${BUILD_FPIC}
                       -DENABLE_GENERIC=${ENABLE_GENERIC}
                       -DLIBC_INTERJECT=${LIBC_INTERJECT}
                       -Dpybind11_DIR=${pybind11_DIR}
                       -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                       -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
            CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                             -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                             -DCMAKE_Fortran_FLAGS:STRING=${CMAKE_Fortran_FLAGS}
                             -DCMAKE_C_COMPILER_ARG1:STRING=${CMAKE_C_COMPILER_ARG1}
                             -DCMAKE_CXX_COMPILER_ARG1:STRING=${CMAKE_CXX_COMPILER_ARG1}
                             -DCMAKE_Fortran_COMPILER_ARG1:STRING=${CMAKE_Fortran_COMPILER_ARG1})

        set(gdma_DIR ${STAGED_INSTALL_PREFIX}/${_gdma_dir} CACHE PATH "path to internally built gdmaConfig.cmake" FORCE)
    endif()
else()
    add_library(gdma_external INTERFACE)  # dummy
endif()

