# We require C++11 support from the compiler and standard library.

if (CMAKE_CXX_COMPILER_ID MATCHES GNU)
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
        message(FATAL_ERROR "GCC version must be at least 4.9!")
    endif()

elseif (CMAKE_CXX_COMPILER_ID MATCHES Intel)
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "16.0.2")
        message(FATAL_ERROR "ICPC version must be at least 2016.0.2!")
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
            RUN_OUTPUT_VARIABLE GCC_VERSION)
    message(STATUS "Found base compiler version ${GCC_VERSION}")
    file(REMOVE ${_testfl})

    if (APPLE)
        if (${GCC_VERSION} VERSION_LESS 6.1)
            message(FATAL_ERROR "${BoldYellow}Intel ICPC makes use of CLANG (detected: ${GCC_VERSION}; required for C++11: 6.1) so this build won't work without CLANG intervention: http://psicode.org/psi4manual/master/build_planning.html\n${ColourReset}")
        endif()
    else ()
        if (${GCC_VERSION} VERSION_LESS 4.9)
            message(FATAL_ERROR "${BoldYellow}Intel ICPC makes use of GCC (detected: ${GCC_VERSION}; required for C++11: 4.9) so this build won't work without GCC intervention: http://psicode.org/psi4manual/master/build_planning.html#faq-modgcc\n${ColourReset}")
        endif()
    endif()

elseif (CMAKE_CXX_COMPILER_ID MATCHES AppleClang)
    # underlying Clang 3.5 known to not work, so AppleClang value from
    #   https://gist.github.com/yamaya/2924292
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.1)
        message(FATAL_ERROR "APPLECLANG version must be at least 6.1!")
    endif()

elseif (CMAKE_CXX_COMPILER_ID MATCHES Clang)
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.6)
        message(FATAL_ERROR "CLANG version must be at least 3.6!")
    endif()

else()
    message(WARNING "Please add a check in CheckCompilerVersion.cmake for ${CMAKE_CXX_COMPILER_ID}.")
endif()

