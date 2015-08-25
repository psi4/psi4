# If the PYTHON_INTERPRETER variable is not empty, then the user passed a custom version
# of the interpreter to be used. If that's the case, set PYTHON_EXECUTABLE accordingly and
# proceed with detection of Python interpreter accordingly.

if("${PYTHON_INTERPRETER}" STREQUAL "")
   find_package(PythonInterp REQUIRED)
else()
    if(NOT EXISTS "${PYTHON_INTERPRETER}")
        find_program(PYTHON_EXECUTABLE NAMES ${PYTHON_INTERPRETER})
        if (NOT EXISTS "${PYTHON_EXECUTABLE}")
            set(PYTHONINTERP_FOUND FALSE)
        endif()
    else()
        set(PYTHONINTERP_FOUND TRUE)
        set(PYTHON_EXECUTABLE "${PYTHON_INTERPRETER}")
    endif()
endif()
find_package(PythonInterp REQUIRED)



# Now find Python libraries and headers of the EXACT SAME VERSION
if(PYTHONINTERP_FOUND)
   # Get Python include path from Python interpreter
   execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c
                         "import distutils.sysconfig, sys; sys.stdout.write(distutils.sysconfig.get_python_inc())"
                  OUTPUT_VARIABLE _PYTHON_INCLUDE_PATH
                  RESULT_VARIABLE _PYTHON_INCLUDE_RESULT)
   # Get Python library path from interpreter
   execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c
                        "from distutils.sysconfig import get_config_var; import sys; sys.stdout.write(get_config_var('LIBDIR'))"
                  OUTPUT_VARIABLE _PYTHON_LIB_PATH
                  RESULT_VARIABLE _PYTHON_LIB_RESULT)


   set(PYTHON_INCLUDE_DIR ${_PYTHON_INCLUDE_PATH} CACHE PATH "Path to a directory")
   set(_PYTHON_VERSION "${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}")
   set(_PYTHON_VERSION_NO_DOTS "${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}")

   find_library(PYTHON_LIBRARY
     NAMES
     python${_PYTHON_VERSION_NO_DOTS}
     python${_PYTHON_VERSION}mu
     python${_PYTHON_VERSION}m
     python${_PYTHON_VERSION}u
     python${_PYTHON_VERSION}
     NO_DEFAULT_PATH
     HINTS
     "${_PYTHON_LIB_PATH}"
     DOC "Path to Python library file."
   )
   if (NOT EXISTS "${PYTHON_LIBRARY}")
       # redo with default paths
       find_library(PYTHON_LIBRARY
         NAMES
         python${_PYTHON_VERSION_NO_DOTS}
         python${_PYTHON_VERSION}mu
         python${_PYTHON_VERSION}m
         python${_PYTHON_VERSION}u
         python${_PYTHON_VERSION}
         HINTS
         "${_PYTHON_LIB_PATH}"
         DOC "Path to Python library file."
       )
   endif()

   mark_as_advanced(CLEAR PYTHON_EXECUTABLE)
   mark_as_advanced(FORCE PYTHON_LIBRARY)
   mark_as_advanced(FORCE PYTHON_INCLUDE_DIR)
endif()

find_package_handle_standard_args(Python
    FOUND_VAR Python_FOUND
    REQUIRED_VARS
        PYTHON_LIBRARY
        PYTHON_INCLUDE_DIR
        PYTHON_EXECUTABLE)

if(NOT Python_FOUND)
    message(FATAL_ERROR "Could NOT find Python")
endif()
