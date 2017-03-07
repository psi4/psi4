

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


include(FindPackageHandleStandardArgs)
include(FindPackageMessage)

foreach(_service BLAS LAPACK)
    if(NOT ${_service}_LANG)
        set(${_service}_LANG C)
    elseif(${_service}_LANG STREQUAL "C" OR ${_service}_LANG STREQUAL "CXX")
        set(${_service}_LANG C)
    elseif(NOT ${_service}_LANG STREQUAL "Fortran")
        message(FATAL_ERROR "Invalid ${_service} library linker language: ${${_service}_LANG}")
    endif()
endforeach()

macro(find_math_header _service _header)
    string(TOUPPER ${_service} _SERVICE)
    find_path(${_SERVICE}_INCLUDE_DIRS
        NAMES ${_header}
        PATHS ${${_SERVICE}_ROOT}
        HINTS ${${_SERVICE}_ROOT}/include
        PATH_SUFFIXES ${MATH_INCLUDE_PATH_SUFFIXES}
        NO_DEFAULT_PATH
        )
    # the following is needed for Atlas' clapack.h
    # this whole code needs major cleanup soon (2014-10-31)
    find_path(${_SERVICE}_INCLUDE_DIRS
        NAMES ${_header}
        PATH_SUFFIXES ${MATH_INCLUDE_PATH_SUFFIXES}
        )
    find_path(${_SERVICE}_INCLUDE_DIRS
        NAMES ${_header}
        PATH_SUFFIXES include
        )
    set(${_SERVICE}_H ${_header})
    unset(_SERVICE)
endmacro()

macro(find_math_libs _service)
    string(TOUPPER ${_service} _SERVICE)
    if(${_SERVICE}_FOUND)
        return()
    endif()
    set(_lib)
    set(_libs)
    set(NARGN ${ARGN})
    foreach(l ${NARGN})
        set(_stat "")
        if(ENABLE_GENERIC)
            IF((NOT "${l}" STREQUAL "iomp5") AND (NOT "${l}" STREQUAL "pthread") AND
               (NOT "${l}" STREQUAL "dl") AND (NOT "${l}" STREQUAL "m"))
                set(_stat "lib${l}.a")
            endif()
        endif()
        find_library(_lib
            NAMES ${_stat} ${l}
            PATHS ${${_SERVICE}_ROOT}
            HINTS ${${_SERVICE}_ROOT}/lib64 ${${_SERVICE}_ROOT}/lib
            PATH_SUFFIXES ${MATH_LIBRARY_PATH_SUFFIXES}
            NO_DEFAULT_PATH
            )
        find_library(_lib
            NAMES ${_stat} ${l}
            PATH_SUFFIXES ${MATH_LIBRARY_PATH_SUFFIXES}
            )
        if(_lib)
            set(_libs ${_libs} ${_lib})
        elseif((${l} STREQUAL "-Wl,--start-group") OR (${l} STREQUAL "-Wl,--end-group"))
            set(_libs ${_libs} ${l})
        else()
            set(${_SERVICE}_LIBRARIES ${_SERVICE}_LIBRARIES-NOTFOUND)
            set(_libs ${_SERVICE}_LIBRARIES-NOTFOUND)
            break()
        endif()
        unset(_lib CACHE)
    endforeach()
    set(${_SERVICE}_LIBRARIES ${_libs})
    unset(_lib CACHE)
    unset(_libs CACHE)
    unset(_SERVICE)
    unset(l)
endmacro()

macro(cache_math_result _service MATH_TYPE)
    string(TOUPPER ${_service} _SERVICE)
    set(${_SERVICE}_FIND_QUIETLY TRUE)
    if(DEFINED ${_SERVICE}_INCLUDE_DIRS)
        find_package_handle_standard_args(
            ${_SERVICE}
            "Could NOT find ${MATH_TYPE} ${_SERVICE}"
            ${_SERVICE}_LIBRARIES
            ${_SERVICE}_INCLUDE_DIRS
            )
    else()
        find_package_handle_standard_args(
            ${_SERVICE}
            "Could NOT find ${MATH_TYPE} ${_SERVICE}"
            ${_SERVICE}_LIBRARIES
            )
    endif()

    if(${_SERVICE}_FOUND)
        set(${_SERVICE}_TYPE ${MATH_TYPE} CACHE STRING
            "${_SERVICE} type")
        mark_as_advanced(${_SERVICE}_TYPE)

        add_definitions(-DHAVE_${MATH_TYPE}_${_SERVICE})
        set(HAVE_${_SERVICE} ON CACHE INTERNAL
            "Defined if ${_SERVICE} is available"
            )
        set(HAVE_${MATH_TYPE}_${_SERVICE} ON CACHE INTERNAL
            "Defined if ${MATH_TYPE}_${_SERVICE} is available"
            )
        set(${_SERVICE}_LIBRARIES ${${_SERVICE}_LIBRARIES} CACHE STRING
            "${_SERVICE} libraries"
            )
        mark_as_advanced(${_SERVICE}_LIBRARIES)
        if(DEFINED ${_SERVICE}_INCLUDE_DIRS)
            set(${_SERVICE}_H ${${_SERVICE}_H} CACHE STRING
                "${_SERVICE} header file")
            mark_as_advanced(${_SERVICE}_H)
            set(${_SERVICE}_INCLUDE_DIRS ${${_SERVICE}_INCLUDE_DIRS}
                CACHE STRING "${_SERVICE} include directory"
                )
            mark_as_advanced(${_SERVICE}_INCLUDE_DIRS)
        endif()
    else()
        set(${_SERVICE}_LIBRARIES ${_SERVICE}_LIBRARIES-NOTFOUND)
        if(DEFINED ${_SERVICE}_H)
            set(${_SERVICE}_INCLUDE_DIRS ${_SERVICE}_INCLUDE_DIRS-NOTFOUND)
            unset(${_SERVICE}_H)
        endif()
    endif()
    set(${_SERVICE}_FOUND ${${_SERVICE}_FOUND} PARENT_SCOPE)
    unset(MATH_TYPE)
    unset(_SERVICE)
endmacro()

macro(config_math_service _SERVICE)
    set(ENABLE_AUTO_${_SERVICE}
        ENABLE_AUTO_${_SERVICE}
        CACHE BOOL
        "Enable ${_SERVICE}"
        )
    set(${_SERVICE}_FOUND FALSE)
    if(ENABLE_AUTO_${_SERVICE})
        if(EXISTS $ENV{MATH_ROOT})
            if("${${_SERVICE}_ROOT}" STREQUAL "")
                set(${_SERVICE}_ROOT $ENV{MATH_ROOT})
                message("-- ${_SERVICE} will be searched for based on MATH_ROOT=${${_SERVICE}_ROOT} ")
            endif()
        endif()

        if(EXISTS $ENV{${_SERVICE}_ROOT})
            if("${${_SERVICE}_ROOT}" STREQUAL "")
                set(${_SERVICE}_ROOT $ENV{${_SERVICE}_ROOT})
                message("-- ${_SERVICE} will be searched for based on ${_SERVICE}_ROOT=${${_SERVICE}_ROOT}")
            endif()
        endif()

        if(EXISTS $ENV{MKL_ROOT})
            if("${${_SERVICE}_ROOT}" STREQUAL "")
                set(${_SERVICE}_ROOT $ENV{MKL_ROOT})
                message("-- ${_SERVICE} will be searched for based on MKL_ROOT=${${_SERVICE}_ROOT}")
            endif()
        endif()

        if(EXISTS $ENV{MKLROOT})
            if("${${_SERVICE}_ROOT}" STREQUAL "")
                set(${_SERVICE}_ROOT $ENV{MKLROOT})
                message("-- ${_SERVICE} will be searched for based on MKLROOT=${${_SERVICE}_ROOT}")
            endif()
        endif()

        if(${_SERVICE}_INCLUDE_DIRS AND ${_SERVICE}_LIBRARIES)
            set(${_SERVICE}_FIND_QUIETLY TRUE)
        endif()

        if(NOT ${_SERVICE}_FIND_COMPONENTS)
            if(DEFINED ${_SERVICE}_TYPE)
                set(${_SERVICE}_FIND_COMPONENTS ${${_SERVICE}_TYPE})
            else()
                set(${_SERVICE}_FIND_COMPONENTS ${MATH_LIB_SEARCH_ORDER})
            endif()
        endif()

        find_service(${_SERVICE})
    endif()

    if(${_SERVICE}_FOUND)

        # this codeblock is dead
        # take care of omp flags
        if(ENABLE_THREADED_MKL)
            set(_omp_flag)
            if(HAVE_MKL_BLAS OR HAVE_MKL_LAPACK)
                if(MKL_COMPILER_BINDINGS MATCHES Intel)
                    set(_omp_flag -qopenmp)
                endif()
                if(MKL_COMPILER_BINDINGS MATCHES GNU)
                    set(_omp_flag -fopenmp)
                endif()
                if(MKL_COMPILER_BINDINGS MATCHES PGI)
                    set(_omp_flag -mp)
                endif()
            endif()
            if(HAVE_MKL_${_SERVICE})
                if(APPLE)
                    set(${_SERVICE}_LIBRARIES ${${_SERVICE}_LIBRARIES} ${_omp_flag})
                else()
                    set(${_SERVICE}_LIBRARIES -Wl,--start-group ${${_SERVICE}_LIBRARIES} ${_omp_flag} -Wl,--end-group)
                endif()
            endif()
            unset(_omp_flag)
        endif()

        find_package_message(${_SERVICE}
            "Found ${_SERVICE}: ${${_SERVICE}_TYPE} (${${_SERVICE}_LIBRARIES})"
            "[${${_SERVICE}_LIBRARIES}]"
            )

        set(EXTERNAL_LIBS
            ${EXTERNAL_LIBS}
            ${${_SERVICE}_LIBRARIES}
            )
    else()
        if(ENABLE_AUTO_${_SERVICE})
            message(FATAL_ERROR "-- No external ${_SERVICE} library found (have you set the MATH_ROOT environment variable?)")
        endif()
        add_definitions(-DUSE_BUILTIN_${_SERVICE})
        set(USE_BUILTIN_${_SERVICE} TRUE)
    endif()
endmacro()

macro(find_math_library _myservice _mytype)
    set(MATH_INCLUDE_PATH_SUFFIXES ${${_mytype}_${_myservice}_INCLUDE_PATH_SUFFIXES})
    #if(${_myservice}_LANG STREQUAL "C")
    #    find_math_header(${_myservice} ${${_mytype}_${_myservice}_HEADERS})
    #endif()
    # yes really, both unsets are needed
    unset(${_myservice}_INCLUDE_DIRS)
    unset(${_myservice}_INCLUDE_DIRS CACHE)
    if(${_mytype} STREQUAL "MKL")
        # fyi: special-casing is an ugly way of doing this.
        #   better way would be extra variable specifying whether calling code (per dist) needs headers found.
        find_math_header(${_myservice} ${${_mytype}_${_myservice}_HEADERS})
    endif()
    set(MATH_LIBRARY_PATH_SUFFIXES ${${_mytype}_${_myservice}_LIBRARY_PATH_SUFFIXES})

    find_math_libs(${_myservice} ${${_mytype}_${_myservice}_LIBS})
    # try some alternative patterns (if defined) until we find it
    foreach(_i 2 3 4 5 6 7 8 9)
        if(NOT ${_myservice}_LIBRARIES)
            if(DEFINED ${_mytype}_${_myservice}_LIBS${_i})
                find_math_libs(${_myservice} ${${_mytype}_${_myservice}_LIBS${_i}})
            endif()
        endif()
    endforeach()
endmacro()

function(find_service _myservice)
    foreach(_component ${${_myservice}_FIND_COMPONENTS})
        find_math_library(${_myservice} ${_component})
        cache_math_result(${_myservice} ${_component})
        if(${_myservice}_FOUND)
            break()
        endif()
    endforeach()
endfunction()
