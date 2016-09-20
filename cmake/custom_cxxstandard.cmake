# We require C++11 support from the compiler and standard library.

if (CMAKE_CXX_COMPILER_ID MATCHES GNU)
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
        message(FATAL_ERROR "GCC version must be at least 4.9!")
    endif()

elseif (CMAKE_CXX_COMPILER_ID MATCHES Intel)
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 16)
        message(FATAL_ERROR "ICPC version must be at least 2016!")
    endif()

    execute_process(COMMAND gcc -dumpversion
                    OUTPUT_VARIABLE GCC_VERSION
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    if (${GCC_VERSION} VERSION_LESS 4.9)
        message(WARNING "Intel ICPC makes use of GCC (detected: ${GCC_VERSION}; required for C++11: 4.9) so this build won't work without GCC intervention: https://github.com/psi4/psi4/wiki/8_FAQ_Contents#modgcc")
    endif()

    if("${CMAKE_VERSION}" VERSION_LESS "3.6")
        # add flag by hand if CMake does not support ICPC
        if("${CMAKE_CXX_STANDARD}" MATCHES "11|14")
            #add_compile_options(-std=c++${CMAKE_CXX_STANDARD})
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++${CMAKE_CXX_STANDARD}")
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

