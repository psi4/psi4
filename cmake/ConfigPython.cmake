# If the PYTHON_INTERPRETER variable is not empty, then the user passed a custom version
# of the interpreter to be used. If that's the case, set PYTHON_EXECUTABLE accordingly and
# proceed with detection of Python interpreter accordingly.
if("${PYTHON_INTERPRETER}" STREQUAL "")
   find_package(PythonInterp 2.6 REQUIRED)
else()
   set(PYTHONINTERP_FOUND TRUE)
   set(PYTHON_EXECUTABLE "${PYTHON_INTERPRETER}")
   find_package(PythonInterp 2.6 REQUIRED)
endif()	

# Now find Python libraries and headers of the EXACT SAME VERSION
# Set variables to help find Python library with the EXACT SAME VERSION as the interpreter
# Copied from https://bitbucket.org/fenics-project/dolfin
if(PYTHONINTERP_FOUND)
   # Get Python include path from Python interpreter
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
   # Finally, find Python libs using the standard module
   find_package(PythonLibs ${PYTHON_VERSION_STRING} EXACT REQUIRED)
endif()
