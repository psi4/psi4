# Downloaded from
#   https://github.com/coderefinery/autocmake/blob/master/modules/safeguards.cmake
# * changed text of in-source message

#.rst:
#
# Provides safeguards against in-source builds and bad build types.
#
# Variables used::
#
#   PROJECT_SOURCE_DIR
#   PROJECT_BINARY_DIR
#   CMAKE_BUILD_TYPE

if(${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please run CMake from top directory and specify a build directory (e.g., cmake -H. -Bbuild).")
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
string(TOUPPER "${CMAKE_BUILD_TYPE}" cmake_build_type_toupper)

if(NOT cmake_build_type_tolower STREQUAL "debug" AND
   NOT cmake_build_type_tolower STREQUAL "release" AND
   NOT cmake_build_type_tolower STREQUAL "relwithdebinfo")
    message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\". Allowed values are Debug, Release, RelWithDebInfo (case-insensitive).")
endif()
