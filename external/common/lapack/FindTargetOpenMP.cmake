# FindTargetOpenMP.cmake
# ----------------------
#
# OpenMP cmake module to wrap FindOpenMP.cmake in a target.
# It is designed to imitate or supplement the official Kitware module c. 3.10.
#
# This module sets the following variables in your project: ::
#
#   TargetOpenMP_FOUND - true if OpenMP found on the system
#   TargetOpenMP_MESSAGE - status message with OpenMP library path list
#   OpenMP_CXX_LIB_NAMES - ";"-separated list of CXX libraries for OpenMP (no paths)
#
# * this is a hack until we update min cmake (to 3.10?)
# * it is specific to Psi4 in that:
#   * we only care about C++ OpenMP
#   * minimum Intel in 2017, so no need to add logic for lower
#   * we never use OPENMP_FOUND
# * this is run by the TargetLAPACK module and tgt::MathOpenMP is supplemented by
#   math libraries detected there
#
# This module defines the following imported targets ::
#
#   OpenMP::OpenMP_CXX - if enabled and available for compiler for C++.
#                        may be found by newer FindOpenMP or constructed here.
#

set(PN TargetOpenMP)

# 1st precedence - libraries passed in through -DOpenMP_LIBRARIES
if (OpenMP_LIBRARIES AND OpenMP_FLAGS AND OpenMP_CXX_LIB_NAMES)
    if (NOT ${PN}_FIND_QUIETLY)
        message (STATUS "OpenMP detection suppressed.")
    endif()

    add_library (OpenMP::OpenMP_CXX INTERFACE IMPORTED)
    set_property (TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_FLAGS})
    set_property (TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_LIBRARIES})
else()
    # 2nd precedence - target from modern FindOpenMP.cmake
    find_package (OpenMP QUIET MODULE COMPONENTS CXX)

    if (NOT TARGET OpenMP::OpenMP_CXX)
        # 3rd precedence - construct a target

        if (CMAKE_CXX_COMPILER_ID MATCHES GNU)
            set(OpenMP_CXX_FLAGS "-fopenmp")
            set(OpenMP_CXX_LIB_NAMES "gomp;pthread")
        elseif (CMAKE_CXX_COMPILER_ID MATCHES Intel)
            set(OpenMP_CXX_FLAGS "-qopenmp")
            set(OpenMP_CXX_LIB_NAMES "iomp5;pthread")
        elseif (CMAKE_CXX_COMPILER_ID STREQUAL Clang)
            if (OpenMP_CXX_FLAG MATCHES "-fopenmp=libiomp5")
                set(OpenMP_CXX_FLAGS "-fopenmp=libiomp5")
                set(OpenMP_CXX_LIB_NAMES "iomp5")
            elseif (OpenMP_CXX_FLAG MATCHES "-fopenmp=libgomp")
                set(OpenMP_CXX_FLAGS "-fopenmp=libgomp")
                set(OpenMP_CXX_LIB_NAMES "gomp")
            else()
                set(OpenMP_CXX_FLAGS "-fopenmp=libomp")
                set(OpenMP_CXX_LIB_NAMES "omp")
            endif()
        endif()

        find_omp_libs("TargetOpenMP" ${OpenMP_CXX_LIB_NAMES})
        if (TargetOpenMP_LIBRARIES)
            add_library(OpenMP::OpenMP_CXX INTERFACE IMPORTED)
            set_property (TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_COMPILE_OPTIONS "$<$<COMPILE_LANGUAGE:CXX>:${OpenMP_CXX_FLAGS}>")
            set_property (TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_LINK_LIBRARIES "${TargetOpenMP_LIBRARIES}")
            if (NOT ${PN}_FIND_QUIETLY)
                message (STATUS "OpenMP target constructed.")
            endif()
        endif()
    endif()
endif()

if (TARGET OpenMP::OpenMP_CXX)
    get_property (_ill TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_LINK_LIBRARIES)
    set (${PN}_MESSAGE "Found TargetOpenMP: ${_ill}")
endif()

#include(CMakePrintHelpers)
#cmake_print_properties(TARGETS OpenMP::OpenMP_CXX
#                       PROPERTIES INTERFACE_COMPILE_DEFINITIONS INTERFACE_COMPILE_OPTIONS INTERFACE_INCLUDE_DIRS INTERFACE_LINK_LIBRARIES)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args (${PN} DEFAULT_MSG ${PN}_MESSAGE)
