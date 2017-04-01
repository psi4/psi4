if(${ENABLE_CheMPS2})
    find_package(CheMPS2 1.8.3 CONFIG QUIET)

    if(${CheMPS2_FOUND})
        get_property(_loc TARGET CheMPS2::chemps2 PROPERTY LOCATION)
        message(STATUS "${Cyan}Found CheMPS2${ColourReset}: ${_loc} (found version ${CheMPS2_VERSION})")
        add_library(chemps2_external INTERFACE)  # dummy

        # reset below necessary as find_package clears it. better solution sought
        set(TargetHDF5_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/TargetHDF5 CACHE PATH "path to externally detected CheMPS2Config.cmake" FORCE)
    else()
        include(ExternalProject)
        message(STATUS "Suitable CheMPS2 could not be located, ${Magenta}Building CheMPS2${ColourReset} instead.")

        if(${BUILD_SHARED_LIBS})
            set(_a_only  OFF)
            set(_so_only ON)
        else()
            set(_a_only  ON)
            set(_so_only OFF)
        endif()

        ExternalProject_Add(chemps2_external
            DEPENDS lapack_external
                    hdf5_external
            GIT_REPOSITORY https://github.com/SebWouters/CheMPS2
            GIT_TAG 55d6b38  # v1.8.3-12
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                       -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
                       -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                       -DSTATIC_ONLY=${_a_only}
                       -DSHARED_ONLY=${_so_only}
                       -DENABLE_OPENMP=${ENABLE_OPENMP}  # relevant
                       -DENABLE_XHOST=${ENABLE_XHOST}
                       -DBUILD_FPIC=${BUILD_FPIC}
                       -DENABLE_GENERIC=${ENABLE_GENERIC}
                       -DCMAKE_RANLIB=${CMAKE_RANLIB}
                       -DCMAKE_AR=${CMAKE_AR}
                       -DCMAKE_NM=${CMAKE_NM}
                       -DENABLE_TESTS=OFF
                       -DTargetLAPACK_DIR=${STAGED_INSTALL_PREFIX}/share/cmake/TargetLAPACK
                       -DTargetHDF5_DIR=${STAGED_INSTALL_PREFIX}/share/cmake/TargetHDF5
            CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                             -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
            INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
            DESTDIR=${CMAKE_BINARY_DIR}/stage)

        set(CheMPS2_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/CheMPS2 CACHE PATH "path to internally built CheMPS2Config.cmake" FORCE)
    endif()
else()
    add_library(chemps2_external INTERFACE)  # dummy
endif()

