include(ExternalProject)
find_package(Ambit)

if(NOT Ambit_FOUND AND ENABLE_AMBIT)
    message(STATUS "Ambit not found. The pre-packaged version will be build.")

    ExternalProject_Add(ambit_external
       GIT_REPOSITORY https://github.com/jturney/ambit
       GIT_TAG v0.1.1-alpha
       CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                  -DENABLE_OMP=${OpenMP_FOUND}
                  -DENABLE_VECTORIZATION=${ENABLE_VECTORIZATION}
                  -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                  -DEXTRA_Fortran_FLAGS=${CMAKE_EXTRA_Fortran_FLAGS}
                  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                  -DEXTRA_C_FLAGS=${CMAKE_EXTRA_C_FLAGS}
                  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                  -DEXTRA_CXX_FLAGS=${CMAKE_EXTRA_CXX_FLAGS}
                  #-DBOOST_INCLUDEDIR=${Boost_INCLUDE_DIRS}
                  #-DBOOST_LIBRARYDIR=${Boost_LIBRARY_DIR}
                  -DPYTHON_INTERPRETER=${PYTHON_EXECUTABLE}
                  -DENABLE_STATIC=ON
                  -DENABLE_PSI4=ON
                  -DPSI4_SOURCE_DIR=${PROJECT_SOURCE_DIR}/psi4
                  -DPSI4_BINARY_DIR=${CMAKE_BINARY_DIR}
                  -DPSI4_INCLUDE_DIRS=${PYTHON_INCLUDE_DIR}
                  -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/external/ambit
                  -DCMAKE_INSTALL_LIBDIR=lib
        INSTALL_DIR ${CMAKE_MAKE_PROGRAM} install 
                    DESTDIR ${CMAKE_BINARY_DIR}/stage
    )
    add_dependencies(external_ambit external_int external_deriv external_boost)
else()
   add_library(ambit_external INTERFACE)
endif ()


