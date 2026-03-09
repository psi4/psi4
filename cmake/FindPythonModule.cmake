# Downloaded from
#   http://www.cmake.org/pipermail/cmake/2011-January/041666.html
# * Added CMake stage dir to PYTHONPATH
# * Added FORCE to location var
# * Function to macro so module_FOUND shows up
# * Remove ``if(NOT PY_${module})`` so runs each time, also remove module caps
# * Module name, not path in fphsa
# * Added version handling via parse_version call
# * Switched PYTHON -> Python

#.rst:
#
# Find if a Python module is installed.
# Usage: find_python_module(<package> [[ATLEAST | EXACT] version] [QUIET] [REQUIRED])

macro(find_python_module module)
    cmake_parse_arguments(ARG "QUIET;REQUIRED" "ATLEAST;EXACT" "" ${ARGN})

    if(ARG_QUIET)
        set(${module}_FIND_QUIETLY TRUE)
    endif()

    if(ARG_REQUIRED)
        set(${module}_FIND_REQUIRED TRUE)
    endif()

    if(ARG_ATLEAST AND ARG_EXACT)
        message(FATAL_ERROR "Can't be both ATLEAST and EXACT")
    endif()
    if(ARG_ATLEAST)
        set(_op ">=")
        set(${module}_tgtver ${ARG_ATLEAST})
    elseif(ARG_EXACT)
        set(_op "==")
        set(${module}_tgtver ${ARG_EXACT})
    else()
        # deceive handle_standard_arguments into not caring about version
        set(_${module}_requested_version_found "${Python_EXECUTABLE}")
    endif()

    unset(PY_${module} CACHE)
    unset(${module}_VERSION CACHE)

    # * A module's location is usually a directory, but for binary modules
    #   it's a .so file.
    # * Unsure of the balance btwn submission to user's PYTHONPATH and avoiding
    #   strays in same. So clobbering user for now with `sys.path.insert(0`
    execute_process(COMMAND "${Python_EXECUTABLE}" "-c"
                            "import re, sys; \
                             sys.path.insert(0, '${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}'); \
                             import ${module}; \
                             print(re.compile('/__init__.py.*').sub('', ${module}.__file__))"
        RESULT_VARIABLE _${module}_status
        OUTPUT_VARIABLE _${module}_location
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT ${_${module}_status})
        set(PY_${module} ${_${module}_location} CACHE STRING
            "Location of Python module ${module}" FORCE)

        execute_process(COMMAND "${Python_EXECUTABLE}" "-c"
                                "import sys; \
                                 sys.path.insert(0, '${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}'); \
                                 import ${module}; \
                                 print(${module}.__version__)"
        RESULT_VARIABLE _${module}_ver_status
        OUTPUT_VARIABLE _${module}_version
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(NOT ${_${module}_ver_status})
            set(${module}_VERSION ${_${module}_version} CACHE STRING
                "Version of Python module ${module}" FORCE)

            if(${module}_tgtver)
                execute_process(COMMAND "${Python_EXECUTABLE}" "-c"
                                        "from packaging.version import parse; \
                                         print(parse('${${module}_VERSION}') ${_op} parse('${${module}_tgtver}'))"
                RESULT_VARIABLE _${module}_verenuf_status
                OUTPUT_VARIABLE _${module}_verenuf
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE)
                if(NOT ${_${module}_verenuf_status})
                    if(${_${module}_verenuf} STREQUAL "True")
                        set(_${module}_requested_version_found "${Python_EXECUTABLE}")
                    endif()
                endif()
            endif()
        endif()
    endif()
    find_package_handle_standard_args(${module} DEFAULT_MSG PY_${module} _${module}_requested_version_found)
endmacro()
