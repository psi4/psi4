if(${ENABLE_CheMPS2})
    find_package(CheMPS2 CONFIG QUIET)

    if(${CheMPS2_FOUND})
        get_property(_loc TARGET CheMPS2::chemps2 PROPERTY LOCATION)
        message(STATUS "Found CheMPS2: ${_loc} (found version ${CheMPS2_VERSION})")
        add_library(chemps2_external INTERFACE)  # dummy
    else()
        include(ExternalProject)
        message(STATUS "Suitable CheMPS2 could not be located, building one instead.")

        if(${BUILD_SHARED_LIBS})
            set(_a_only  OFF)
            set(_so_only ON)
        else()
            set(_a_only  ON)
            set(_so_only OFF)
        endif()

        ExternalProject_Add(chemps2_external
            GIT_REPOSITORY https://github.com/loriab/CheMPS2
            GIT_TAG cmtarget
            #GIT_REPOSITORY https://github.com/SebWouters/CheMPS2  #TODO
            #GIT_TAG v1.7.1  # TODO tag when stable
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/external/chemps2
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
            CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                             -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
            INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
            DESTDIR=${CMAKE_BINARY_DIR}/stage)
    endif()
else()
    add_library(chemps2_external INTERFACE)  # dummy
endif()

