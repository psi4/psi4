# Downloaded from
#   http://www.cmake.org/pipermail/cmake/2011-January/041666.html
# * Added CMake stage dir to PYTHONPATH
# * Added FORCE to location var
# * Function to macro so module_FOUND shows up
# * Remove ``if(NOT PY_${module})`` so runs each time, also remove module caps
# * Module name, not path in fphsa

#.rst:
#
# Find if a Python module is installed.
# Usage: find_python_module(PyQt4 REQUIRED)

macro(find_python_module module)
    unset(PY_${module} CACHE)
    if(ARGC GREATER 1 AND ARGV1 STREQUAL "REQUIRED")
        set(${module}_FIND_REQUIRED TRUE)
    endif()
    # * A module's location is usually a directory, but for binary modules
    #   it's a .so file.
    # * Unsure of the balance btwn submission to user's PYTHONPATH and avoiding
    #   strays in same. So clobbering user for now with `sys.path.insert(0`
    execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c" 
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
    endif()
    find_package_handle_standard_args(${module} DEFAULT_MSG PY_${module})
endmacro()
