if(${ENABLE_libcint})
    # libcint (Sun, J. Comput. Chem. 36, 1664, 2015) provides an optional
    # alternative two-electron integral engine, selected at runtime with
    # INTEGRAL_PACKAGE LIBCINT.
    find_package(libcint CONFIG QUIET)

    if(TARGET libcint::libcint)
        get_property(_loc TARGET libcint::libcint PROPERTY LOCATION)
        message(STATUS "${Cyan}Found libcint${ColourReset}: ${_loc} (found version ${libcint_VERSION})")
        add_library(libcint_external INTERFACE)  # dummy

        if(${CMAKE_INSIST_FIND_PACKAGE_libcint})
            message(VERBOSE "Suitable libcint located externally as user insists")
        endif()
    else()
        # No libcintConfig.cmake -- common for distro packages (e.g. Fedora).
        # Locate the shared library and header and synthesize a config package
        # so the core build's find_package(libcint CONFIG) succeeds.
        find_library(LIBCINT_LIBRARY NAMES cint)
        find_path(LIBCINT_INCLUDE_DIR NAMES cint.h)

        if(LIBCINT_LIBRARY AND LIBCINT_INCLUDE_DIR)
            message(STATUS "${Cyan}Found system libcint${ColourReset}: ${LIBCINT_LIBRARY}")
            set(_gen ${CMAKE_CURRENT_BINARY_DIR}/libcint-config)
            file(MAKE_DIRECTORY ${_gen})
            file(WRITE ${_gen}/libcintConfig.cmake
"if(NOT TARGET libcint::libcint)
  add_library(libcint::libcint UNKNOWN IMPORTED)
  set_target_properties(libcint::libcint PROPERTIES
    IMPORTED_LOCATION \"${LIBCINT_LIBRARY}\"
    INTERFACE_INCLUDE_DIRECTORIES \"${LIBCINT_INCLUDE_DIR}\")
endif()
set(libcint_VERSION \"system\")
")
            file(WRITE ${_gen}/libcintConfigVersion.cmake
"set(PACKAGE_VERSION \"system\")
set(PACKAGE_VERSION_COMPATIBLE TRUE)
set(PACKAGE_VERSION_EXACT TRUE)
")
            set(libcint_DIR ${_gen} CACHE PATH "path to synthesized libcintConfig.cmake" FORCE)
            add_library(libcint_external INTERFACE)  # dummy
        else()
            if(${CMAKE_INSIST_FIND_PACKAGE_libcint})
                message(FATAL_ERROR "Suitable libcint could not be externally located as user insists")
            endif()

            include(ExternalProject)
            message(STATUS "Suitable libcint could not be located, ${Magenta}Building libcint${ColourReset} instead.")
            ExternalProject_Add(libcint_external
                GIT_REPOSITORY https://github.com/sunqm/libcint
                GIT_TAG v6.1.2
                UPDATE_COMMAND ""
                CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                           -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                           -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                           -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                           -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                           -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                           -DWITH_RANGE_COULOMB=1
                           -DWITH_FORTRAN=0
                           -DENABLE_TEST=0
                           -DENABLE_EXAMPLE=0
                           -DCMAKE_POSITION_INDEPENDENT_CODE=${BUILD_FPIC}
                           -DCMAKE_POLICY_VERSION_MINIMUM=3.5
                CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS})
            set(libcint_DIR ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/cmake/libcint CACHE PATH "path to internally built libcintConfig.cmake" FORCE)
        endif()
    endif()
else()
    message(STATUS "Optional addon ${Blue}unsought: libcint${ColourReset}")
    add_library(libcint_external INTERFACE)  # dummy
endif()
