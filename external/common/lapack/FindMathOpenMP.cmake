# FindMathOpenMP.cmake
# --------------------
#
# MathOpenMP cmake module to bridge the gap between requirements
# of independently found BLAS/LAPACK and OpenMP and form a practical
# MathOpenMP target.
#
# This module sets the following variables in your project: ::
#
#   MathOpenMP_FOUND - always true
#   MathOpenMP_MESSAGE - status message with MathOpenMP library path list
#
# * this is a hack until we update min cmake (to 3.10?)
# * it is specific to Psi4 in that:
#   * we only care about C++ OpenMP
#   * minimum Intel in 2017, so no need to add logic for lower
#   * we never use OPENMP_FOUND
# * this is run by the TargetLAPACK module and tgt::lapack is supplemented by
#   math libraries detected there
#
# This module defines the following imported targets ::
#
#   tgt::MathOpenMP - always defined, though it may be a dummy target
#                     if OpenMP disabled or compiler incapable. If OpenMP
#                     enabled and found, linked to OpenMP::OpenMP_CXX as
#                     well as carries any additional flags or libraries
#                     needed by BLAS/LAPACK.
#

macro(find_omp_libs _service)
    # Expands any name-only libraries
    # * adapted from: https://github.com/coderefinery/autocmake/blob/22ef2b17da6a26813fc3c8b1240bcabfd9317295/modules/math_libs.cmake#L304-L337
    set(_lib)
    set(_libs)
    set(_ls "${ARGN}")
    foreach(_l IN LISTS _ls)
        get_filename_component(_fullspec "${OpenMP_${_l}_LIBRARY}" DIRECTORY)

        #message("${_service} ${_l} 0 ${_lib} ${_fullspec}")
        set(_stat "")
        find_library(_lib
                     NAMES ${_l}
                     HINTS ${_fullspec}
                           ${OpenMP_CXX_LIBRARY_DIRS}
                           ${OpenMP_ROOT}
                     DOC "Path to the ${_l} library for OpenMP"
                     NO_DEFAULT_PATH)
        #message("${_service} ${_l} A ${_lib}")
        find_library(_lib
                     NAMES ${_l})
        #message("${_service} ${_l} B ${_lib}")
        if(_lib)
            set(_libs ${_libs} ${_lib})
        elseif(${_l} MATCHES "-Wl,")
            set(_libs ${_libs} ${_l})
        else()
            ## remainder could be implicit libs handlable by the linker
            #set(_libs ${_libs} ${_l})
            set(_libs ${_service}_LIBRARIES-NOTFOUND)
            break()
        endif()
        unset(_lib CACHE)
    endforeach()
    set(${_service}_LIBRARIES ${_libs})
    unset(_lib CACHE)
    unset(_libs CACHE)
    #message("    ${_service}_LIBRARIES ${${_service}_LIBRARIES}")
    unset(_service)
endmacro()


add_library(tgt::MathOpenMP INTERFACE IMPORTED)

if (${isMKL} MATCHES "MKL")
    if (CMAKE_CXX_COMPILER_ID STREQUAL GNU)
        set(_MathOpenMP_LIB_NAMES "iomp5;-Wl,--as-needed")
        find_omp_libs("MathOpenMP" ${_MathOpenMP_LIB_NAMES})
        set_property(TARGET tgt::MathOpenMP PROPERTY INTERFACE_LINK_LIBRARIES "${MathOpenMP_LIBRARIES}")
    endif()

    if (CMAKE_CXX_COMPILER_ID STREQUAL Clang)
        if (NOT DEFINED OpenMP_CXX_FLAG)
            set (OpenMP_CXX_FLAG "-fopenmp=libiomp5")
        endif()
    endif()
endif()

if (ENABLE_OPENMP)
    find_package(TargetOpenMP)
endif()

set(PN MathOpenMP)
get_property (_ill TARGET tgt::MathOpenMP PROPERTY INTERFACE_LINK_LIBRARIES)
set (${PN}_MESSAGE "Found MathOpenMP: ${_ill}")
    
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args (${PN} DEFAULT_MSG ${PN}_MESSAGE)

if (TARGET OpenMP::OpenMP_CXX)
    set_property(TARGET tgt::MathOpenMP APPEND PROPERTY INTERFACE_LINK_LIBRARIES "OpenMP::OpenMP_CXX")
endif()

#include(CMakePrintHelpers)
#cmake_print_properties(TARGETS tgt::MathOpenMP
#                       PROPERTIES INTERFACE_COMPILE_DEFINITIONS INTERFACE_COMPILE_OPTIONS INTERFACE_INCLUDE_DIRS INTERFACE_LINK_LIBRARIES)
