# - Find the NumPy libraries
# This module finds if NumPy is installed, and creates an interface library.
# For the moment the interface library only carries include directories.

#============================================================================
# Copyright 2012 Continuum Analytics, Inc.
#
# MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to permit
# persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
#
#============================================================================

include(FindPythonModule)

# Finding NumPy involves calling the Python interpreter
if(NumPy_FIND_REQUIRED)
    find_python_module(numpy REQUIRED QUIET)
else()
    find_python_module(numpy QUIET)
endif()

# Get include directories and prepare a target
if(numpy_FOUND)
  execute_process(
    COMMAND
      "${PYTHON_EXECUTABLE}" "-c"
      "import numpy as n; print(n.__version__); print(n.get_include());"
    RESULT_VARIABLE
      _NUMPY_SEARCH_SUCCESS
    OUTPUT_VARIABLE
      _NUMPY_VALUES
    ERROR_VARIABLE
      _NUMPY_ERROR_VALUE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

  # Convert the process output into a list
  string(REGEX REPLACE ";" "\\\\;" _NUMPY_VALUES ${_NUMPY_VALUES})
  string(REGEX REPLACE "\n" ";" _NUMPY_VALUES ${_NUMPY_VALUES})
  list(GET _NUMPY_VALUES 0 NUMPY_VERSION)
  list(GET _NUMPY_VALUES 1 NUMPY_INCLUDE_DIRS)

  # Make sure all directory separators are '/'
  string(REGEX REPLACE "\\\\" "/" NUMPY_INCLUDE_DIRS ${NUMPY_INCLUDE_DIRS})

  # Get the major and minor version numbers
  string(REGEX REPLACE "\\." ";" _NUMPY_VERSION_LIST ${NUMPY_VERSION})
  list(GET _NUMPY_VERSION_LIST 0 NUMPY_VERSION_MAJOR)
  list(GET _NUMPY_VERSION_LIST 1 NUMPY_VERSION_MINOR)
  list(GET _NUMPY_VERSION_LIST 2 NUMPY_VERSION_PATCH)
  string(REGEX MATCH "[0-9]*" NUMPY_VERSION_PATCH ${NUMPY_VERSION_PATCH})
  math(EXPR NUMPY_VERSION_DECIMAL
      "(${NUMPY_VERSION_MAJOR} * 10000) + (${NUMPY_VERSION_MINOR} * 100) + ${NUMPY_VERSION_PATCH}")

  add_library(NumPy INTERFACE)
  target_include_directories(NumPy
    INTERFACE
      ${NUMPY_INCLUDE_DIRS}
    )

  set(NUMPY_FOUND TRUE)
endif()
