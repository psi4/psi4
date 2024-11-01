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

set(_TargetOpenMP_PN ${PN})
set(PN TargetOpenMP)

separate_arguments(${PN}_FIND_COMPONENTS)
if(${PN}_FIND_COMPONENTS)
    set(_${PN}_FIND_LIST ${${PN}_FIND_COMPONENTS})
else()
    set(_${PN}_FIND_LIST C CXX Fortran)
endif()
separate_arguments(_${PN}_FIND_LIST)

if (NOT ${PN}_FIND_QUIETLY)
    message(STATUS "Detecting MathOpenMP -- ?OpenMP=${ENABLE_OPENMP}, ?MKL=${isMKL}, LANG=${_${PN}_FIND_LIST}, C/CXX/Fortran=${CMAKE_C_COMPILER_ID}/${CMAKE_CXX_COMPILER_ID}/${CMAKE_Fortran_COMPILER_ID}")
endif()

# 1st precedence - libraries passed in through -DOpenMP_LIBRARIES
if (OpenMP_LIBRARIES AND OpenMP_FLAGS)
    if (NOT ${PN}_FIND_QUIETLY)
        message (STATUS "OpenMP detection suppressed.")
    endif()

    add_library (OpenMP::OpenMP INTERFACE IMPORTED)
    target_compile_options(OpenMP::OpenMP INTERFACE ${OpenMP_FLAGS})
    target_link_libraries(OpenMP::OpenMP INTERFACE ${OpenMP_LIBRARIES})
else()
    # 2nd precedence - target from modern FindOpenMP.cmake
    find_package (OpenMP MODULE COMPONENTS ${_${PN}_FIND_LIST})

    if(NOT OpenMP_FOUND)
        message(STATUS "CMake FindOpenMP failed! Trying a custom OpenMP configuration...")
    endif()

    foreach(_lang IN LISTS _${PN}_FIND_LIST)
        if (NOT TARGET OpenMP::OpenMP_${_lang})
            # 3rd precedence - construct a target

            if(WIN32)
                # OpenMP config for clang-cl on Windonws
                # Note: FindOpenMP doesn't yet support clang-cl, so the config has to be done manually.

                # Check if clang-cl is used
                if(CMAKE_${_lang}_COMPILER_ID STREQUAL Clang)
                    unset(USE_CLANG_CL_${_lang} CACHE)
                    if(${_lang} STREQUAL C)
                        check_c_compiler_flag("-Xclang -fopenmp" USE_CLANG_CL_${_lang})
                    elseif(${_lang} STREQUAL CXX)
                        check_cxx_compiler_flag("-Xclang -fopenmp" USE_CLANG_CL_${_lang})
                    else()
                        message(FATAL_ERROR "clang-cl (MSVC compatability wrapper) does not support ${_lang}")
                    endif()
                    if(NOT USE_CLANG_CL_${_lang})
                        message(FATAL_ERROR "clang-cl (MSVC compatability wrapper) is required, if Clang is used")
                    endif()
                else()
                    message(FATAL_ERROR "Clang compiler is required, if ENABLE_OPENMP=ON")
                endif()

                # Set OpenMP compiler options and run-time library
                set(OpenMP_${_lang}_FLAGS "-Xclang;-fopenmp") # This has to be preceed by "-Xclang"
                set(OpenMP_${_lang}_LIB_NAMES "libiomp5md")

            else()
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

            endif()

            find_omp_libs("${PN}_${_lang}" ${OpenMP_${_lang}_LIB_NAMES})
            if (${PN}_${_lang}_LIBRARIES)
                add_library(OpenMP::OpenMP_${_lang} INTERFACE IMPORTED)
                target_compile_options(OpenMP::OpenMP_${_lang} INTERFACE $<$<COMPILE_LANGUAGE:${_lang}>:${OpenMP_${_lang}_FLAGS}>)
                target_link_libraries(OpenMP::OpenMP_${_lang} INTERFACE ${${PN}_${_lang}_LIBRARIES})
                if (NOT ${PN}_FIND_QUIETLY)
                    message (STATUS "OpenMP::OpenMP_${_lang} target constructed (${${PN}_${_lang}_LIBRARIES})")
                endif()
            endif()
        endif()
    endforeach()
endif()

unset(_omp_target_for_sought_langs)
add_library(OpenMP::OpenMP INTERFACE IMPORTED)
foreach(_lang C CXX Fortran)
    if((NOT ${PN}_FIND_COMPONENTS AND CMAKE_${_lang}_COMPILER_LOADED) OR _lang IN_LIST ${PN}_FIND_COMPONENTS)
        if (TARGET OpenMP::OpenMP_${_lang})
            set(${PN}_${_lang}_FOUND 1)
            set_property(TARGET OpenMP::OpenMP APPEND PROPERTY INTERFACE_LINK_LIBRARIES OpenMP::OpenMP_${_lang})
        endif()
        list(APPEND _omp_target_for_sought_langs "${PN}_${_lang}_FOUND")
    endif()
endforeach()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(${PN}
                                  REQUIRED_VARS ${_omp_target_for_sought_langs}
                                  HANDLE_COMPONENTS)

unset(_${PN}_FIND_LIST)

set(PN ${_TargetOpenMP_PN})
unset(_TargetOpenMP_PN)

#include(CMakePrintHelpers)
#message("Targets after find_package(TargetOpenMP)")
#cmake_print_properties(TARGETS OpenMP::OpenMP_C OpenMP::OpenMP_CXX OpenMP::OpenMP_Fortran OpenMP::OpenMP
#                       PROPERTIES INTERFACE_COMPILE_DEFINITIONS INTERFACE_COMPILE_OPTIONS INTERFACE_INCLUDE_DIRS INTERFACE_LINK_LIBRARIES)
cmake_policy(POP)
