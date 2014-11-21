
# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.")
endif()

# guard against bad build-type strings
if(NOT CMAKE_BUILD_TYPE)
    message("INFO: build type was not defined, using type \"Debug\".")
    set(CMAKE_BUILD_TYPE "Debug")
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
string(TOUPPER "${CMAKE_BUILD_TYPE}" cmake_build_type_toupper)
if(    NOT cmake_build_type_tolower STREQUAL "debug"
   AND NOT cmake_build_type_tolower STREQUAL "release"
   AND NOT cmake_build_type_tolower STREQUAL "profile"
   AND NOT cmake_build_type_tolower STREQUAL "relwithdebinfo")
      message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\". Allowed values are Debug, Release, Profile, RelWithDebInfo (case-insensitive).")
endif()

# guard against math-less build
if(MKL_FLAG_SET)
   # MKL_FLAG_SET is set to ON by using the --mkl setup flag
   # also set BLAS_FOUND and LAPACK_FOUND to TRUE
   set(BLAS_FOUND TRUE)
   set(LAPACK_FOUND TRUE)
endif()	
if(NOT BLAS_FOUND OR NOT LAPACK_FOUND) 
   if(NOT EXPLICIT_LIBS)
      message(FATAL_ERROR "No BLAS/LAPACK implementation found and no explicit libraries specified")
   else()
      message(STATUS "No BLAS/LAPACK implementation found, but explicit libraries specified")
   endif()
endif()

# guard against sanitizer builds
if(ENABLE_ASAN)
   message(STATUS "Building with address sanitizer only works with GCC 4.9 or Clang!!")
elseif(ENABLE_MSAN)
   message(STATUS "Building with memory sanitizer only works with GCC 4.9 or Clang!!")
elseif(ENABLE_TSAN)
   message(STATUS "Building with thread sanitizer only works with GCC 4.9 or Clang!!")
elseif(ENABLE_UBSAN)
   message(STATUS "Building with undefined behaviour sanitizer only works with GCC 4.9 or Clang!!")
endif()
