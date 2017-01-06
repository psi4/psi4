include(ExternalProject)
find_package(Ambit QUIET)

if(NOT Ambit_FOUND AND ${ENABLE_AMBIT})
    message(STATUS "Ambit not found. The pre-packaged version will be build.")
    set(ENABLE_STATIC (NOT ${BUILD_SHARED_LIBS}))
    ExternalProject_Add(ambit_external
       GIT_REPOSITORY https://github.com/jturney/ambit
       GIT_TAG v0.2
       CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                  -DENABLE_OMP=${ENABLE_OPENMP}
                  -DENABLE_VECTORIZATION=${ENABLE_VECTORIZATION}
                  -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                  -DEXTRA_Fortran_FLAGS=${CMAKE_EXTRA_Fortran_FLAGS}
                  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                  -DEXTRA_C_FLAGS=${CMAKE_EXTRA_C_FLAGS}
                  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                  -DEXTRA_CXX_FLAGS=${CMAKE_EXTRA_CXX_FLAGS}
                  -DBOOST_ROOT=${BOOST_ROOT}
                  -DPYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}
                  -DPYTHON_INCLUDE_DIR=${PYTHON_INCLUDE_DIR}
                  -DPYTHON_LIBRARY=${PYTHON_LIBRARY}
                  -DPYTHON_INTERPRETER=${PYTHON_EXECUTABLE}
                  -DENABLE_STATIC=${ENABLE_STATIC}
                  -DENABLE_PSI4=ON
                  -DPSI4_SOURCE_DIR=${PROJECT_SOURCE_DIR}/psi4
                  -DPSI4_BINARY_DIR=${CMAKE_BINARY_DIR}
                  -DPSI4_INCLUDE_DIRS=${PYTHON_INCLUDE_DIR}
                  -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                  -DCMAKE_INSTALL_LIBDIR=lib
        CMAKE_CACHE_ARGS -DCMAKE_PREFIX_PATH:STRING=${CMAKE_PREFIX_PATH}
        INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
                    DESTDIR ${CMAKE_BINARY_DIR}/stage
    )
    add_dependencies(ambit_external boost_external)
else()
   add_library(ambit_external INTERFACE)
endif ()


