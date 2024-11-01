cmake_policy(PUSH)
cmake_policy(SET CMP0057 NEW)  # support IN_LISTS

# We require C++20 support from the compiler and standard library.
list(APPEND _allowed_cxx_standards 20)
if(NOT psi4_CXX_STANDARD IN_LIST _allowed_cxx_standards)
  message(FATAL_ERROR "Psi4 requires C++20 at least")
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES GNU)
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10.0)
        message(FATAL_ERROR "GCC version must be at least 10.0!")
    endif()

elseif (CMAKE_CXX_COMPILER_ID MATCHES Intel)
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "19.0.0")
        message(FATAL_ERROR "ICPC version must be at least 2019.0.0 to work with pybind11 2.1!")  # v1.2
    endif()

    set(_testfl ${CMAKE_BINARY_DIR}/test_gcc_version.cc)
    file(WRITE  ${_testfl} "
    #include <stdio.h>

    int main() {
        #ifdef __clang__
        printf(\"%d.%d.%d\", __clang_major__, __clang_minor__, __clang_patchlevel__);
        #else
        printf(\"%d.%d.%d\", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
        #endif
        return 0;
    }
    ")
    try_run(GCCV_COMPILES
            GCCV_RUNS
            ${CMAKE_BINARY_DIR} ${_testfl}
            RUN_OUTPUT_VARIABLE CPLR_VERSION)
    message(STATUS "Found base compiler version ${CPLR_VERSION}")
    file(REMOVE ${_testfl})

    if (APPLE)
        if ("${CPLR_VERSION}" VERSION_LESS 3.6)
            message(FATAL_ERROR "${BoldYellow}Intel ICPC makes use of CLANG (detected: ${CPLR_VERSION}; required for C++11: Clang 3.6 or AppleClang 6.1) so this build won't work without CLANG intervention.\n${ColourReset}")
        endif()
    else ()
        if ("${CPLR_VERSION}" VERSION_LESS 4.9)
            message(FATAL_ERROR "${BoldYellow}Intel ICPC makes use of GCC (detected: ${CPLR_VERSION}; required for C++11: 4.9) so this build won't work without GCC intervention: http://psicode.org/psi4manual/master/build_planning.html#faq-modgcc\n${ColourReset}")
        endif()
    endif()

elseif (CMAKE_CXX_COMPILER_ID MATCHES AppleClang)
    # underlying Clang 3.5 known to not work, so AppleClang value from
    #   https://gist.github.com/yamaya/2924292
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.1)
        message(FATAL_ERROR "APPLECLANG version must be at least 6.1!")
    endif()

elseif (CMAKE_CXX_COMPILER_ID MATCHES Clang)
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11.0)
        message(FATAL_ERROR "CLANG version must be at least 11.0!")
    endif()

elseif (CMAKE_CXX_COMPILER_ID MATCHES MSVC)
    # As for MSVC 14.0, it is not possible to set anything below C++14
    if(MSVC_TOOLSET_VERSION LESS 140)
        message(FATAL_ERROR "MSVC toolset version must be at least 14.0!")
    endif()

else()
    message(WARNING "Please add a check in custom_cxxstandard.cmake for ${CMAKE_CXX_COMPILER_ID}.")
endif()

cmake_policy(POP)
