

# -----------------------------------------------------------------------------
# Copyright 2011-2013 Jonas Juselius <firstname.lastname at uit.no>
#                     Radovan Bast   <lastname at kth.se>
#
# Distributed under the GNU Lesser General Public License.
# See accompanying file LICENSE for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# -----------------------------------------------------------------------------

if(ENABLE_64BIT_INTEGERS)
    set(MATH_LIB_SEARCH_ORDER MKL ACML)
elseif(APPLE)
    set(MATH_LIB_SEARCH_ORDER MKL SYSTEM_NATIVE OPENBLAS ESSL ATLAS ACML)  # prefer Accelerate unless MKL
else()
    set(MATH_LIB_SEARCH_ORDER MKL OPENBLAS ESSL ATLAS ACML SYSTEM_NATIVE)
endif()

if(NOT DEFINED MKL_FLAG)
    if(NOT ENABLE_BUILTIN_BLAS AND NOT ENABLE_BUILTIN_LAPACK)
        if(NOT HAVE_BLAS AND NOT HAVE_LAPACK)
            if(ENABLE_64BIT_INTEGERS)
                message(STATUS "Since you specified 64bit integers the math lib search order is (only) ${MATH_LIB_SEARCH_ORDER}")
                message(STATUS "This is because apart from MKL and ACML default math library installations are built for 32bit integers")
                message(STATUS "If you know that the library you want to use provides 64bit integers, you can select the library")
                message(STATUS "with -D BLAS_TYPE=X or -D LAPACK_TYPE X (X: MKL ESSL ATLAS ACML SYSTEM_NATIVE)")
                message(STATUS "or by redefining MATH_LIB_SEARCH_ORDER")
            else()
                message(STATUS "Math lib search order is ${MATH_LIB_SEARCH_ORDER}")
                message(STATUS "You can select a specific type by defining for instance -D BLAS_TYPE=ATLAS or -D LAPACK_TYPE=ACML")
                message(STATUS "or by redefining MATH_LIB_SEARCH_ORDER")
            endif()
        endif()
    endif()
endif()

include(MathLibs)
include(MathLibsFunctions)

foreach(_service BLAS LAPACK)
    if(ENABLE_BUILTIN_${_service})
        if(ENABLE_AUTO_${_service})
            message(FATAL_ERROR "ENABLE_BUILTIN_${_service} and ENABLE_AUTO_${_service} together makes no sense")
        endif()
        set(USE_BUILTIN_${_service} TRUE)
        message("-- ${_service}: using builtin implementation (slow)")
    else()
        set(USE_BUILTIN_${_service} FALSE)
    endif()

    #if(DEFINED EXPLICIT_${_service}_LIB)
    #    if(ENABLE_AUTO_${_service})
    #        message(FATAL_ERROR "EXPLICIT_${_service}_LIB and ENABLE_AUTO_${_service} together makes no sense")
    #    endif()
    #    if(ENABLE_BUILTIN_${_service})
    #        message(FATAL_ERROR "EXPLICIT_${_service}_LIB and ENABLE_BUILTIN_${_service} together makes no sense")
    #    endif()
    #    set(EXTERNAL_LIBS
    #        ${EXTERNAL_LIBS}
    #        ${EXPLICIT_${_service}_LIB}
    #        )
    #    message("-- ${_service}: using explicit library (${EXPLICIT_${_service}_LIB})")
    #endif()
endforeach()

if(ENABLE_CRAY_WRAPPERS)
    foreach(_service BLAS LAPACK)
        if(ENABLE_BUILTIN_${_service})
            message(FATAL_ERROR "ENABLE_CRAY_WRAPPERS and ENABLE_BUILTIN_${_service} together makes no sense")
        endif()
        if(ENABLE_AUTO_${_service})
            message(FATAL_ERROR "ENABLE_CRAY_WRAPPERS and ENABLE_AUTO_${_service} together makes no sense")
        endif()
    endforeach()
    message("-- Use CRAY wrappers; this disables math detection and builtin math libraries")
endif()

# EXTERNAL_LIBS no longer supported, and -mkl won't get written to LAPACKTarget correctly
#
#if(DEFINED MKL_FLAG)
#    foreach(_service BLAS LAPACK)
#        if(ENABLE_BUILTIN_${_service})
#            message(FATAL_ERROR "MKL_FLAG and ENABLE_BUILTIN_${_service} together makes no sense")
#        endif()
#        if(ENABLE_AUTO_${_service})
#            message(FATAL_ERROR "MKL_FLAG and ENABLE_AUTO_${_service} together makes no sense")
#        endif()
#    endforeach()
#    set(EXTERNAL_LIBS
#        ${EXTERNAL_LIBS}
#        ${MKL_FLAG}
#        )
#    message("-- User set explicit MKL flag which is passed to the compiler and linker: ${MKL_FLAG}")
#    message("-- This disables math detection and builtin math libraries")
#    message("-- Setting -DHAVE_MKL_BLAS and -DHAVE_MKL_LAPACK")
#    add_definitions(-DHAVE_MKL_BLAS)
#    add_definitions(-DHAVE_MKL_LAPACK)
#endif()

foreach(_service BLAS LAPACK)
    if(ENABLE_AUTO_${_service})
        config_math_service(${_service})
    endif()
    if(${_service}_FOUND)
        include_directories(${${_service}_INCLUDE_DIRS})
    endif()
endforeach()

