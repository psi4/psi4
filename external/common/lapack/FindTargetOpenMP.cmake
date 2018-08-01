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
#   OpenMP_<lang>_LIB_NAMES - ";"-separated list of <lang> libraries for OpenMP (no paths)
#                              not present if bypassed by OpenMP_FLAGS & OpenMP_LIBRARIES
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
#   OpenMP::OpenMP_<lang> - if enabled and available for compiler for <lang>.
#                           may be found by newer FindOpenMP or constructed here.
#
#   OpenMP::OpenMP - concatenation of requested language targets.
#
# Inputs ::
#
#   TargetOpenMP_FIND_COMPONENTS - set to languages (C, CXX, Fortran) to find OpenMP,
#                                  default is all enabled.
#   OpenMP_LIBRARY_DIRS - list of directories where OpenMP libraries may be found,
#                         in preference to DEFAULT_PATHS
#
cmake_policy(PUSH)
cmake_policy(SET CMP0057 NEW)  # support IN_LISTS


set(PN TargetOpenMP)

separate_arguments(${PN}_FIND_COMPONENTS)
if(NOT ${PN}_FIND_COMPONENTS)
    set(${PN}_FINDLIST C CXX Fortran)
else()
    set(${PN}_FINDLIST ${${PN}_FIND_COMPONENTS})
endif()
separate_arguments(${PN}_FINDLIST)

if (NOT ${PN}_FIND_QUIETLY)
    message(STATUS "Detecting MathOpenMP -- ?OpenMP=${ENABLE_OPENMP}, ?MKL=${isMKL}, LANG=${${PN}_FINDLIST}, C/CXX/Fortran=${CMAKE_C_COMPILER_ID}/${CMAKE_CXX_COMPILER_ID}/${CMAKE_Fortran_COMPILER_ID}")
endif()

# 1st precedence - libraries passed in through -DOpenMP_LIBRARIES
if (OpenMP_LIBRARIES AND OpenMP_FLAGS)
    if (NOT ${PN}_FIND_QUIETLY)
        message (STATUS "OpenMP detection suppressed.")
    endif()

    add_library (OpenMP::OpenMP INTERFACE IMPORTED)
    set_property (TARGET OpenMP::OpenMP PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_FLAGS})
    set_property (TARGET OpenMP::OpenMP PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_LIBRARIES})
else()
    # 2nd precedence - target from modern FindOpenMP.cmake
    find_package (OpenMP QUIET MODULE COMPONENTS ${${PN}_FINDLIST})

    foreach(_lang ${${PN}_FINDLIST})
        if (NOT TARGET OpenMP::OpenMP_${_lang})
            # 3rd precedence - construct a target

            if (CMAKE_${_lang}_COMPILER_ID MATCHES GNU)
                set(OpenMP_${_lang}_FLAGS "-fopenmp")
                set(OpenMP_${_lang}_LIB_NAMES "gomp;pthread")

            elseif (CMAKE_${_lang}_COMPILER_ID MATCHES Intel)
                set(OpenMP_${_lang}_FLAGS "-qopenmp")
                set(OpenMP_${_lang}_LIB_NAMES "iomp5;pthread")

            elseif (CMAKE_${_lang}_COMPILER_ID STREQUAL Clang)
                if (OpenMP_${_lang}_FLAG MATCHES "-fopenmp=libiomp5")
                    set(OpenMP_${_lang}_FLAGS "-fopenmp=libiomp5")
                    set(OpenMP_${_lang}_LIB_NAMES "iomp5")
                elseif (OpenMP_${_lang}_FLAG MATCHES "-fopenmp=libgomp")
                    set(OpenMP_${_lang}_FLAGS "-fopenmp=libgomp")
                    set(OpenMP_${_lang}_LIB_NAMES "gomp")
                else()
                    set(OpenMP_${_lang}_FLAGS "-fopenmp=libomp")
                    set(OpenMP_${_lang}_LIB_NAMES "omp")
                endif()

            endif()

            find_omp_libs("TargetOpenMP_${_lang}" ${OpenMP_${_lang}_LIB_NAMES})
            if (TargetOpenMP_${_lang}_LIBRARIES)
                add_library(OpenMP::OpenMP_${_lang} INTERFACE IMPORTED)
                set_property (TARGET OpenMP::OpenMP_${_lang} PROPERTY INTERFACE_COMPILE_OPTIONS "$<$<COMPILE_LANGUAGE:${_lang}>:${OpenMP_${_lang}_FLAGS}>")
                set_property (TARGET OpenMP::OpenMP_${_lang} PROPERTY INTERFACE_LINK_LIBRARIES "${TargetOpenMP_${_lang}_LIBRARIES}")
                if (NOT ${PN}_FIND_QUIETLY)
                    message (STATUS "OpenMP_${_lang} target constructed.")
                endif()
            endif()
        endif()
    endforeach()
endif()

unset(_omp_target_for_sought_langs)
add_library (OpenMP::OpenMP INTERFACE IMPORTED)
foreach(_lang IN ITEMS C CXX Fortran)
    if((NOT ${PN}_FIND_COMPONENTS AND CMAKE_${_lang}_COMPILER_LOADED) OR _lang IN_LIST ${PN}_FIND_COMPONENTS)
        if (TARGET OpenMP::OpenMP_${_lang})
            set(${PN}_${_lang}_FOUND 1)
            set_property(TARGET OpenMP::OpenMP APPEND PROPERTY INTERFACE_LINK_LIBRARIES "OpenMP::OpenMP_${_lang}")
        endif()
        list(APPEND _omp_target_for_sought_langs "${PN}_${_lang}_FOUND")
        set (${PN}_MESSAGE "Found TargetOpenMP ${${PN}_FINDLIST}: ${_ill}")  # only last lang
    endif()
endforeach()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(${PN}
                                  REQUIRED_VARS ${_omp_target_for_sought_langs}
                                  HANDLE_COMPONENTS)

#include(CMakePrintHelpers)
#message("Targets after find_package(TargetOpenMP)")
#cmake_print_properties(TARGETS OpenMP::OpenMP_C OpenMP::OpenMP_CXX OpenMP::OpenMP_Fortran OpenMP::OpenMP
#                       PROPERTIES INTERFACE_COMPILE_DEFINITIONS INTERFACE_COMPILE_OPTIONS INTERFACE_INCLUDE_DIRS INTERFACE_LINK_LIBRARIES)
cmake_policy(POP)
