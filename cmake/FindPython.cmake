function(get_Python_version_string pyVersion) 
    execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c
                            "import sys; sys.stdout.write(';'.join([str(x) for x in sys.version_info[:3]]))"
                    OUTPUT_VARIABLE _VERSION
                    RESULT_VARIABLE _PYTHON_VERSION_RESULT
                    ERROR_QUIET)
    if(NOT _PYTHON_VERSION_RESULT)
        string(REPLACE ";" "." pyVersion "${_VERSION}")
        list(GET _VERSION 0 PYTHON_VERSION_MAJOR)
        list(GET _VERSION 1 PYTHON_VERSION_MINOR)
        list(GET _VERSION 2 PYTHON_VERSION_PATCH)
        if(PYTHON_VERSION_PATCH EQUAL 0)
            # it's called "Python 2.7", not "2.7.0"
            string(REGEX REPLACE "\\.0$" "" pyVersion "${pyVersion}")
        endif()
    else()
        # sys.version predates sys.version_info, so use that
        execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c "import sys; sys.stdout.write(sys.version)"
                        OUTPUT_VARIABLE _VERSION
                        RESULT_VARIABLE _PYTHON_VERSION_RESULT
                        ERROR_QUIET)
        if(NOT _PYTHON_VERSION_RESULT)
            string(REGEX REPLACE " .*" "" pyVersion "${_VERSION}")
            string(REGEX REPLACE "^([0-9]+)\\.[0-9]+.*" "\\1" PYTHON_VERSION_MAJOR "${pyVersion}")
            string(REGEX REPLACE "^[0-9]+\\.([0-9])+.*" "\\1" PYTHON_VERSION_MINOR "${pyVersion}")
            if(pyVersion MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+.*")
                string(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" PYTHON_VERSION_PATCH "${pyVersion}")
            else()
                set(PYTHON_VERSION_PATCH "0")
            endif()
        else()
            # sys.version was first documented for Python 1.5, so assume
            # this is older.
            set(pyVersion "1.4")
            set(PYTHON_VERSION_MAJOR "1")
            set(PYTHON_VERSION_MAJOR "4")
            set(PYTHON_VERSION_MAJOR "0")
        endif()
    endif()
    set(pyVersion "${pyVersion}" PARENT_SCOPE)
    unset(_PYTHON_VERSION_RESULT)
    unset(_VERSION)
endfunction()

function(check_Python_compiles pyCompiles)
   set(_bindir  "${PROJECT_BINARY_DIR}/pyCompiles")	
   set(_srcfile "${PROJECT_SOURCE_DIR}/cmake/checkPython.cpp")
   
   try_compile(pyCompiles "${_bindir}" "${_srcfile}" 
             CMAKE_FLAGS 
	     "-DINCLUDE_DIRECTORIES=${PYTHON_INCLUDE_DIRS}"
             "-DLINK_LIBRARIES=${PYTHON_LIBRARIES}")

   set(pyCompiles "${pyCompiles}" PARENT_SCOPE)	
endfunction()


# 1. a Python interpreter of the suitable version exists
if("${PYTHON_INTERPRETER}" STREQUAL "")
   find_package(PythonInterp 2.6 REQUIRED)
else()
   set(PYTHONINTERP_FOUND TRUE)
   set(PYTHON_EXECUTABLE "${PYTHON_INTERPRETER}")
   get_Python_version_string(pyVersion)
   set(PYTHON_VERSION_STRING "${pyVersion}")
   message(STATUS "Passed PythonInterp: ${PYTHON_INTERPRETER} (found suitable version \"${PYTHON_VERSION_STRING}\", minimum required is \"2.6\")")
endif()	

# 2. Python libraries and headers of the same version exist
# Set variables to help find Python library that is compatible with interpreter
# Copied from https://bitbucket.org/fenics-project/dolfin
if(PYTHONINTERP_FOUND)
   # Get Python include path from Python interpretter
   execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c
                         "import distutils.sysconfig, sys; sys.stdout.write(distutils.sysconfig.get_python_inc())"
                  OUTPUT_VARIABLE _PYTHON_INCLUDE_PATH
                  RESULT_VARIABLE _PYTHON_INCLUDE_RESULT)
   # Get Python library path from interpreter
   execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c
                        "import os, sys, inspect; sys.stdout.write(os.path.split(os.path.split(inspect.getfile(inspect))[0])[0])"
                  OUTPUT_VARIABLE _PYTHON_LIB_PATH
                  RESULT_VARIABLE _PYTHON_LIB_RESULT)
   # Set include path, if returned by interpreter
   if("${_PYTHON_INCLUDE_RESULT}" STREQUAL "0")
      set(PYTHON_INCLUDE_DIR ${_PYTHON_INCLUDE_PATH} CACHE PATH "Path to a file")
   endif()
   # Add a search path for Python library based on output from
   # interpreter
   set(CMAKE_LIBRARY_PATH_SAVE ${CMAKE_LIBRARY_PATH})
   if("${_PYTHON_LIB_RESULT}" STREQUAL "0")
      set(CMAKE_LIBRARY_PATH ${_PYTHON_LIB_PATH})
   endif()
   # Find Pythons libs
   find_package(PythonLibs ${PYTHON_VERSION_STRING} EXACT REQUIRED)
   # Check headers
   set(CMAKE_REQUIRED_INCLUDES "${PYTHON_INCLUDE_DIRS}")
   check_include_files(Python.h HAS_PYTHON_H)
   check_include_files(pyconfig.h HAS_PYCONFIG_H)
   if(NOT HAS_PYTHON_H)
      message(STATUS "No Python.h header found!!")	   
   endif()	   
   if(NOT HAS_PYCONFIG_H)
      message(STATUS "No pyconfig.h header found!!")	   
   endif()	   
   unset(CMAKE_REQUIRED_INCLUDES)
   # 3. that we can link against those libraries
   check_Python_compiles(pyCompiles)
   if(NOT pyCompiles)
      message(STATUS "Cannot link against Python!!")
   endif()	 
endif()

# Iff the answer is "YES" to all of the above, we set
# EMBEDDED_PYTHON to TRUE enabling the embedded Python code
if(PYTHONINTERP_FOUND AND PYTHONLIBS_FOUND 
   AND HAS_PYTHON_H AND HAS_PYCONFIG_H AND pyCompiles)
   message(STATUS "Python ${PYTHON_VERSION_STRING} FOUND")
else()
   # If detection of interpreter or libs goes wrong CMake stops.
   # We need an additional check here to be sure that all conditions are met
   message(FATAL_ERROR "Something wrong with Python... Check CMake logs")
endif()
