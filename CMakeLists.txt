#       This is the top-level CMakeLists.txt for Psi4.  
# 
# This is version 2.0.0 of the build system.  Version 1 was a swell first
# pass, but quickly became too complicated.  For projects of this size, with
# lots of little (soon to be somewhat) independent pieces the general CMake
# strategy is known as a superbuild.  In such a build all pieces of a package
# are treated as external projects, including the main project itself.  What
# this does is allow us to use the same mechanism (namely find_package) to
# find dependencies and configure each part.
#
#
# Useful CMake variables you should consider setting:
#
#
#

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(psi4 CXX C)

if(${CMAKE_BINARY_DIR}==${CMAKE_SOURCE_DIR})
   message(WARNING "In source builds are prohibited.  Setting build directory "
                   "to build")
   set(CMAKE_BINARY_DIR build)
endif()

option(ENABLE_MPI "Enables modules that depend on MPI (currently on GTFock" OFF)

include(ExternalProject)

find_package(PythonInterp)
add_subdirectory(external)
