# We require C++11 support from the compiler and standard library.

if (CMAKE_CXX_COMPILER_ID MATCHES GNU)
    # Require at least g++ 4.9
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
        message(FATAL_ERROR "GCC version must be at least 4.9!")
    endif()
else()
    message(WARNING "Please add a check in CheckCompilerVersion.cmake for ${CMAKE_CXX_COMPILER_ID}.")
endif()

