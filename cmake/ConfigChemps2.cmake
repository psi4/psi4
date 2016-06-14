
# User signal to try pre-built CheMPS2 (and HDF5)
if (CHEMPS2_ROOT)
    find_package(CHEMPS2)
endif()

# Build CheMPS2 as external package if pre-built failed or not signaled
if(NOT CHEMPS2_FOUND)
    message(STATUS "CheMPS2 not found. The pre-packaged version will be built.")

    find_package (HDF5 QUIET)
    if (NOT HDF5_FOUND)
        message(FATAL_ERROR "No HDF5, no CheMPS2. Build against an existing CheMPS2 with -DCHEMPS2_ROOT=/path/to/chemps2 or skip with -DENABLE_CHEMPS2=OFF")
    endif()

    set(CUSTOM_CHEMPS2_LOCATION ${PROJECT_BINARY_DIR}/interfaces/chemps2)
    include(ExternalProject)
    ExternalProject_Add(interface_chemps2
        PREFIX ${CUSTOM_CHEMPS2_LOCATION}
        GIT_REPOSITORY https://github.com/SebWouters/CheMPS2
        GIT_TAG v1.7
        CMAKE_ARGS 
                   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                   -DEXTRA_C_FLAGS=${CMAKE_EXTRA_C_FLAGS}
                   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                   -DEXTRA_CXX_FLAGS=${CMAKE_EXTRA_CXX_FLAGS}
                   -DSTATIC_ONLY=ON
                   -DENABLE_TESTS=OFF
                   -DENABLE_GENERIC=${ENABLE_STATIC_LINKING}
                   -DENABLE_XHOST=${ENABLE_XHOST}
                   -DHDF5_LIBRARIES=${HDF5_LIBRARIES}
                   -DHDF5_INCLUDE_DIRS=${HDF5_INCLUDE_DIRS}
                   -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
                   -DCMAKE_INSTALL_LIBDIR=lib
        INSTALL_DIR "${CUSTOM_CHEMPS2_LOCATION}/install")

    # Set also variables usually set by find_package
    ExternalProject_Get_Property(interface_chemps2 INSTALL_DIR)
    set(CHEMPS2_LIBRARY "${INSTALL_DIR}/lib/libchemps2.a")
    file(MAKE_DIRECTORY ${INSTALL_DIR}/include/chemps2)  # note [1] below
    set(CHEMPS2_INCLUDE_DIRS "${INSTALL_DIR}/include" ${HDF5_INCLUDE_DIRS})
    set(CHEMPS2_LIBRARIES  "${CHEMPS2_LIBRARY}" "${LAPACK_LIBRARIES}"
                           "${BLAS_LIBRARIES}" ${HDF5_LIBRARIES})

    # Set target for dmrg and psi4 to depend upon as set by find_package
    add_library(CHEMPS2::CHEMPS2 STATIC IMPORTED GLOBAL)
    add_dependencies(CHEMPS2::CHEMPS2 interface_chemps2)
    set_target_properties(CHEMPS2::CHEMPS2 PROPERTIES
        IMPORTED_LOCATION "${CHEMPS2_LIBRARY}"
        INTERFACE_LINK_LIBRARIES "${CHEMPS2_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${CHEMPS2_INCLUDE_DIRS}")

endif()

add_definitions(-DENABLE_CHEMPS2=1)

# [1] It's nice to have a full CHEMPS2::CHEMPS2 target that has embedded
# the library, the linking library paths with dependencies, and the include
# paths with dependencies, just like FindCHEMPS2 supplies, especially
# since src/bin/dmrg needs the headers for compilation and src/bin/psi4
# needs all the linking libraries. Problem is that conventional target
# derived from ExternalProject_Add is a highly sought but not quite
# certified cmake pattern. Hence INTERFACE_INCLUDE_DIRECTORIES complains
# that the directories don't exist at configure time. Hence the hack to 
# create an empty directory.
