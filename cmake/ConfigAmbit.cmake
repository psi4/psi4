
# User signal to try pre-built Ambit (and HDF5)
if (AMBIT_ROOT)
    find_package(Ambit)
    add_definitions(-DHAVE_AMBIT)
endif ()

# Build Ambit as external package if pre-built failed or not signaled
if (NOT Ambit_FOUND)
    message(STATUS "Ambit not found. The pre-packaged version will be build.")

    set(CUSTOM_Ambit_LOCATION ${PROJECT_BINARY_DIR}/interfaces/ambit)
    include(ExternalProject)

    find_package(HDF5 REQUIRED)
    if (NOT HDF5_FOUND)
        message(FATAL_ERROR "No HDF5, no Ambit. Build against existing with -DAMBIT_ROOT=/path/to/ambit or skip with -DENABLE_AMBIT=OFF")
    endif()

    set(Ambit_OPENMP OFF)
    if (ENABLE_OPENMP)
        set(Ambit_OPENMP ON)
    endif ()

    set(Ambit_VECTORIZATION OFF)
    if (ENABLE_VECTORIZATION)
        set(Ambit_VECTORIZATION ON)
    endif ()

    list(APPEND AmbitCMakeArgs
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DENABLE_MPI=OFF
            -DENABLE_OMP=${Ambit_OPENMP}
            -DENABLE_VECTORIZATION=${Ambit_VECTORIZATION}
            -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
            -DEXTRA_Fortran_FLAGS=${CMAKE_EXTRA_Fortran_FLAGS}
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            -DEXTRA_C_FLAGS=${CMAKE_EXTRA_C_FLAGS}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DEXTRA_CXX_FLAGS=${CMAKE_EXTRA_CXX_FLAGS}
            -DBOOST_INCLUDEDIR=${Boost_INCLUDE_DIRS}
            -DBOOST_LIBRARYDIR=${Boost_LIBRARY_DIR}
            -DHDF5_LIBRARIES=${HDF5_LIBRARIES}
            -DHDF5_INCLUDE_DIRS=${HDF5_INCLUDE_DIRS}
            -DPYTHON_INTERPRETER=${PYTHON_EXECUTABLE}
            -DENABLE_STATIC=ON
            -DENABLE_PSI4=ON
            -DPSI4_SOURCE_DIR=${PROJECT_SOURCE_DIR}
            -DPSI4_BINARY_DIR=${PROJECT_BINARY_DIR}
            -DPSI4_INCLUDE_DIRS=${PYTHON_INCLUDE_DIR}
            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
            -DCMAKE_INSTALL_LIBDIR=lib
            )

    ExternalProject_Add(interface_ambit
            PREFIX ${CUSTOM_Ambit_LOCATION}
            GIT_REPOSITORY https://github.com/jturney/ambit
            GIT_TAG v0.1.1-alpha
            CMAKE_ARGS "${AmbitCMakeArgs}"
            INSTALL_DIR "${CUSTOM_Ambit_LOCATION}/install"
            )
    add_dependencies(interface_ambit int deriv)
    if (BUILD_CUSTOM_BOOST)
        add_dependencies(interface_ambit custom_boost)
    endif ()

    ExternalProject_Get_Property(interface_ambit INSTALL_DIR)
    set(Ambit_LIBRARY "${INSTALL_DIR}/lib/libambit.a")
    file(MAKE_DIRECTORY ${INSTALL_DIR}/include/ambit)
    set(Ambit_INCLUDE_DIRS "${INSTALL_DIR}/include" ${HDF5_INCLUDE_DIRS})
    set(Ambit_LIBRARIES "${Ambit_LIBRARY}" "${Boost_LIBRARIES}"
            "${LAPACK_LIBRARIES}" "${BLAS_LIBRARIES}"
            ${HDF5_LIBRARIES})

    add_library(Ambit::Ambit STATIC IMPORTED GLOBAL)
    add_dependencies(Ambit::Ambit interface_ambit)
    set_target_properties(Ambit::Ambit PROPERTIES
            IMPORTED_LOCATION "${Ambit_LIBRARY}"
            INTERFACE_LINK_LIBRARIES "${Ambit_LIBRARIES}"
            INTERFACE_INCLUDE_DIRECTORIES "${Ambit_INCLUDE_DIRS}"
            )

    include_directories(SYSTEM "${Ambit_INCLUDE_DIRS}")
    add_definitions(-DHAVE_AMBIT)
endif ()

