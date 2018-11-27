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
#   * we only care about C++ OpenMP. updated for C/CXX/Fortran
#   * minimum Intel in 2017, so no need to add logic for lower
#   * we never use OPENMP_FOUND
# * this is run by the TargetLAPACK module and tgt::lapack is supplemented by
#   math libraries detected there
#
# This module defines the following imported targets ::
#
#   tgt::MathOpenMP - always defined, though it may be a dummy target
#                     if OpenMP disabled or compiler incapable. If OpenMP
#                     enabled and found, linked to OpenMP::OpenMP as
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

        #message("${_service} ${_l} 0 ${_lib} ${_fullspec} ${${_service}_LIBRARIES}")
        set(_stat "")
        find_library(_lib
                     NAMES ${_l}
                     HINTS ${_fullspec}
                           ${OpenMP_LIBRARY_DIRS}
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
            if (NOT ${_service}_FIND_QUIETLY)
                message(STATUS "Could NOT find ${_l} -- consider adding full directory path to OpenMP_LIBRARY_DIRS")
            endif()
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

set(PN MathOpenMP)
add_library(tgt::${PN} INTERFACE IMPORTED)

if (${isMKL} MATCHES "MKL")
    if ((CMAKE_C_COMPILER_ID STREQUAL GNU) OR
        (CMAKE_CXX_COMPILER_ID STREQUAL GNU) OR
        (CMAKE_Fortran_COMPILER_ID STREQUAL GNU))
        if (APPLE)
            set(${PN}_LIB_NAMES "iomp5")
        else()
            # https://stackoverflow.com/questions/25986091/telling-gcc-to-not-link-libgomp-so-it-links-libiomp5-instead
            set(${PN}_LIB_NAMES "iomp5;-Wl,--as-needed")
        endif()
        find_omp_libs("${PN}" ${${PN}_LIB_NAMES})
        set_property(TARGET tgt::${PN} PROPERTY INTERFACE_LINK_LIBRARIES ${${PN}_LIBRARIES})
    endif()

    if (CMAKE_C_COMPILER_ID STREQUAL Clang)
        if (NOT DEFINED OpenMP_C_FLAG)
            set (OpenMP_C_FLAG "-fopenmp=libiomp5")
        endif()
    endif()
    if (CMAKE_CXX_COMPILER_ID STREQUAL Clang)
        if (NOT DEFINED OpenMP_CXX_FLAG)
            set (OpenMP_CXX_FLAG "-fopenmp=libiomp5")
        endif()
    endif()
endif()

if (ENABLE_OPENMP)
    # *not* REQUIRED b/c some compilers don't support OpenMP and -DENABLE_OPENMP isn't a build-or-die-trying
    find_package(TargetOpenMP COMPONENTS ${TargetOpenMP_FIND_COMPONENTS})
    if (TargetOpenMP_FOUND)
        set_property(TARGET tgt::${PN} APPEND PROPERTY INTERFACE_LINK_LIBRARIES OpenMP::OpenMP)
    else()
        message(WARNING "${PN} configuration failed! The code will be built without OpenMP.")
    endif()
endif()

set(_${PN}_REQUIRED 1)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args (${PN} DEFAULT_MSG _${PN}_REQUIRED)
