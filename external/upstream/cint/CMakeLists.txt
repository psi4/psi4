if(${ENABLE_cint})
    # libcint >= 6.2.0 installs a CMake config package, cint
    find_package(cint 6.2.0 CONFIG QUIET)

    if(${cint_FOUND})
        get_property(_loc TARGET cint::cint PROPERTY LOCATION)
        message(STATUS "${Cyan}Found cint${ColourReset}: ${_loc} (found version ${cint_VERSION})")
        add_library(cint_external INTERFACE)  # dummy

        if(${CMAKE_INSIST_FIND_PACKAGE_cint})
            message(VERBOSE "Suitable cint located externally as user insists")
        endif()

    else()
        # libcint before v6.2.0 installs no CMake config (conda-forge 6.1.x, distro packages). Locate the
        # library and header, or take them from LIBCINT_LIBRARY and LIBCINT_INCLUDE_DIR, and synthesize a
        # config package so the core build's find_package(cint CONFIG) succeeds.
        find_library(LIBCINT_LIBRARY NAMES cint)
        find_path(LIBCINT_INCLUDE_DIR NAMES cint_funcs.h)

        if(LIBCINT_LIBRARY AND LIBCINT_INCLUDE_DIR)
            file(STRINGS ${LIBCINT_INCLUDE_DIR}/cint.h _ver REGEX "^#define CINT_VERSION ")
            string(REGEX REPLACE ".*\"([0-9.]+)\".*" "\\1" _ver "${_ver}")
            if(NOT _ver MATCHES "^[0-9]+\\.[0-9]+")
                message(FATAL_ERROR "Can't read CINT_VERSION from ${LIBCINT_INCLUDE_DIR}/cint.h")
            endif()
            if(_ver VERSION_LESS 6.1.2)  # edit in codedeps
                message(FATAL_ERROR "Found libcint ${_ver} (${LIBCINT_LIBRARY}), but psi4 needs >= 6.1.2. Point "
                                    "LIBCINT_LIBRARY and LIBCINT_INCLUDE_DIR to a newer one, or unset them to build it.")
            endif()
            message(STATUS "${Cyan}Found system libcint${ColourReset}: ${LIBCINT_LIBRARY} (found version ${_ver})")
            set(_gen ${CMAKE_CURRENT_BINARY_DIR}/cint-config)
            file(MAKE_DIRECTORY ${_gen})
            file(WRITE ${_gen}/cintConfig.cmake
"if(NOT TARGET cint::cint)
  add_library(cint::cint UNKNOWN IMPORTED)
  set_target_properties(cint::cint PROPERTIES
    IMPORTED_LOCATION \"${LIBCINT_LIBRARY}\"
    INTERFACE_INCLUDE_DIRECTORIES \"${LIBCINT_INCLUDE_DIR}\")
endif()
set(cint_VERSION \"${_ver}\")
")
            string(REGEX MATCH "^[0-9]+" _major "${_ver}")
            file(WRITE ${_gen}/cintConfigVersion.cmake
"set(PACKAGE_VERSION \"${_ver}\")
if(PACKAGE_FIND_VERSION VERSION_GREATER PACKAGE_VERSION OR NOT PACKAGE_FIND_VERSION_MAJOR EQUAL ${_major})
  set(PACKAGE_VERSION_COMPATIBLE FALSE)
else()
  set(PACKAGE_VERSION_COMPATIBLE TRUE)
  if(PACKAGE_FIND_VERSION VERSION_EQUAL PACKAGE_VERSION)
    set(PACKAGE_VERSION_EXACT TRUE)
  endif()
endif()
")
            set(cint_DIR ${_gen} CACHE PATH "path to synthesized cintConfig.cmake" FORCE)
            add_library(cint_external INTERFACE)  # dummy

        else()
            if(${CMAKE_INSIST_FIND_PACKAGE_cint})
                message(FATAL_ERROR "Suitable cint could not be externally located as user insists")
            endif()

            if(MSVC OR CMAKE_C_SIMULATE_ID STREQUAL "MSVC")
                # libcint's spinor integrals use C99 complex arithmetic (double complex, _Complex_I), which the
                # MSVC runtime's complex.h doesn't provide, so compilers on it (cl, clang-cl, clang targeting
                # MSVC) can't build it. MinGW gcc can.
                message(FATAL_ERROR "Suitable cint could not be located, and it can't be built from source with "
                                    "${CMAKE_C_COMPILER_ID} (MSVC-like), whose C runtime lacks C99 complex arithmetic. "
                                    "Install a gcc-built libcint >= 6.2.0 (e.g., conda-forge libcint) and point "
                                    "cint_DIR to its lib/cmake/cint, or set ENABLE_cint=OFF.")
            endif()

            include(ExternalProject)
            message(STATUS "Suitable cint could not be located, ${Magenta}Building cint${ColourReset} instead.")
            ExternalProject_Add(cint_external
                URL https://github.com/sunqm/libcint/archive/v6.2.0.tar.gz  # edit in codedeps
                UPDATE_COMMAND ""
                CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
                           -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                           -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                           -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                           -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
                           -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                           -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                           -DWITH_CINT2_INTERFACE=0
                           -DWITH_FORTRAN=0
                           -DENABLE_TEST=0
                           -DENABLE_EXAMPLE=0
                           -DCMAKE_POSITION_INDEPENDENT_CODE=${BUILD_FPIC}
                           -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
                           -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
                CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                                 -DCMAKE_C_COMPILER_ARG1:STRING=${CMAKE_C_COMPILER_ARG1})

            set(cint_DIR ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/cmake/cint CACHE PATH "path to internally built cintConfig.cmake" FORCE)
        endif()
    endif()
else()
    message(STATUS "Optional addon ${Blue}unsought: cint${ColourReset}")
    add_library(cint_external INTERFACE)  # dummy
endif()
