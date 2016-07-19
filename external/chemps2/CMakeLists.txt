
# User signal to try pre-built CheMPS2 (and HDF5)
if (CHEMPS2_ROOT)
    find_package(CHEMPS2)
endif()

# Build CheMPS2 as external package if pre-built failed or not signaled
if((NOT ${CHEMPS2_FOUND}) AND ${ENABLE_CHEMPS2})
    message(STATUS "CheMPS2 not found. The pre-packaged version will be built.")
    include(ExternalProject)
    ExternalProject_Add(chemps2_external
        GIT_REPOSITORY https://github.com/SebWouters/CheMPS2
        GIT_TAG v1.7.1
        CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                   -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
                   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                   -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
                   -DSTATIC_ONLY=ON
                   -DENABLE_TESTS=OFF
                   -DENABLE_GENERIC=${ENABLE_STATIC_LINKING}
                   -DENABLE_XHOST=${ENABLE_XHOST}
                   -DCMAKE_RANLIB=${CMAKE_RANLIB}
                   -DCMAKE_AR=${CMAKE_AR}
                   -DCMAKE_NM=${CMAKE_NM}
                   -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                   -DCMAKE_INSTALL_LIBDIR=lib
        INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install 
                        DESTDIR=${CMAKE_BINARY_DIR}/stage
)

else()
   add_library(chemps2_external INTERFACE)
endif()

