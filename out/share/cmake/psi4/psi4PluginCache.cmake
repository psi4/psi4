# psi4PluginCache.cmake
# ---------------------
#
# This module sets some likely variable values to initialize the CMake cache in your plugin.
# See ``psi4 --plugin-compile`` for use.
#

set(CMAKE_C_COMPILER          "/usr/bin/cc" CACHE STRING "")
set(CMAKE_C_FLAGS             " -march=native" CACHE STRING "")
set(CMAKE_CXX_COMPILER        "/usr/bin/c++" CACHE STRING "")
set(CMAKE_CXX_FLAGS           " -march=native" CACHE STRING "")
set(CMAKE_Fortran_COMPILER    "" CACHE STRING "")
set(CMAKE_Fortran_FLAGS       "" CACHE STRING "")

#set(CMAKE_INSTALL_PREFIX      "/home/connor/git/psi4/build/stage" CACHE PATH "")
set(CMAKE_INSTALL_LIBDIR      "lib" CACHE STRING "")
set(CMAKE_INSTALL_BINDIR      "bin" CACHE STRING "")
set(CMAKE_INSTALL_DATADIR     "share" CACHE STRING "")
set(CMAKE_INSTALL_INCLUDEDIR  "include" CACHE STRING "")
set(PYMOD_INSTALL_LIBDIR      "/" CACHE STRING "")

set(CMAKE_INSTALL_MESSAGE     "LAZY" CACHE STRING "")
set(pybind11_DIR              "/home/connor/anaconda3/envs/p4dev/share/cmake/pybind11" CACHE PATH "")

set(PYTHON_VERSION_MAJORMINOR "3.7" CACHE STRING "")
set(Python_VERSION_MAJORMINOR "3.7" CACHE STRING "")
set(PYTHON_EXECUTABLE         "/home/connor/anaconda3/envs/p4dev/bin/python3.7" CACHE STRING "")
set(Python_EXECUTABLE         "/home/connor/anaconda3/envs/p4dev/bin/python3.7" CACHE STRING "")

