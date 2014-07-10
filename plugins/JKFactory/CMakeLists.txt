#Root CMakeLists file for the JK Factory Library

#JKFactory variables that should be set before calling this script:
#
#JK_GA_INC directory where global array include files are
#JK_GA_LIB directory where global array libraries are
#
#If one or both of these two global array variables is not set then
#we will build global arrays for you
#
#JK_ERD_LIB The location of ERD, if not set will be built for you
#
if(NOT DEFINED JK_GA_INC OR NOT DEFINED JK_GA_LIB)
  message(STATUS "No Global Arrays library was specified. One will be built")
  set(JK_GA_INC "NONE")
  set(JK_GA_LIB "NONE")
endif()
if(NOT DEFINED JK_ERD)
   message(STATUS "No electron repulsion integral direct (ERD) was specified. One"
            "will be built")
   set(JK_ERD "NONE")
endif()
if(NOT DEFINED JK_BLAS_INC)
   message(STATUS "No BLAS was given.  Will attempt to find one.")
   set(JK_BLAS_INC "NONE")
endif()

cmake_minimum_required(VERSION 2.8.8)
project(JK_FACTORY C CXX)
set(JK_FACTORY_VERSION_MAJOR 1)
set(JK_FACTORY_VERSION_MINOR 0)

enable_language(Fortran)

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a build 
   directory and run CMake from there. You may need to remove CMakeCache.txt.")
endif()

#Find all of our requirements MPI,OpenMP, Boost, and BLAS
include(cmake/FindReqs.txt)

#If we need to build GA do it, also adds the necesary flags to our compilers
include(cmake/BuildGA.txt)

#Now add GTFock
add_subdirectory(GTFock)

#Finally we build the JKFactory library
set(JKROOT ${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(src)
include_directories(src)
add_library(JKFactory STATIC $<TARGET_OBJECTS:APIS>
                             $<TARGET_OBJECTS:SRCOBJS>
)
add_dependencies(JKFactory cint pfock)
install(TARGETS JKFactory DESTINATION ${CMAKE_BINARY_DIR}/lib)