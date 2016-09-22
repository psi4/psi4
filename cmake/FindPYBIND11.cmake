#.rst:
# FindPYBIND11
# --------
#
# Find the native PYBIND11 includes and dependencies.
# Modeled on FindZLIB.cmake
#
# IMPORTED Targets
# ^^^^^^^^^^^^^^^^
#
# This module defines :prop_tgt:`IMPORTED` target ``PYBIND11::PYBIND11``, if
# PYBIND11 has been found.
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module defines the following variables:
#
# ::
#
#   PYBIND11_INCLUDE_DIRS   - where to find pybind11/pybind11.h, etc.
#   PYBIND11_FOUND          - True if Pybind11 found.
#
# Hints
# ^^^^^
#
# A user may set ``PYBIND11_ROOT`` to a PyBind11 installation root to tell
# this module where to look.

# defines CMAKE_INSTALL_LIBDIR for Debian package
include (GNUInstallDirs)

set(_PYBIND11_SEARCHES)

find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)

# Search PYBIND11_ROOT first if it is set.
if(PYBIND11_ROOT)
  set(_PYBIND11_SEARCH_ROOT PATHS ${PYBIND11_ROOT} NO_DEFAULT_PATH)
  list(APPEND _PYBIND11_SEARCHES _PYBIND11_SEARCH_ROOT)
endif()

# Normal search.
set(_PYBIND11_SEARCH_NORMAL
  PATHS "[HKEY_LOCAL_MACHINE\\SOFTWARE\\GnuWin32\\PyBind11;InstallPath]"
        "$ENV{PROGRAMFILES}/pybind11"
        "/usr"
  )
list(APPEND _PYBIND11_SEARCHES _PYBIND11_SEARCH_NORMAL)

set(PYBIND11_NAMES pybind11)

# Try each search configuration.
unset(PYBIND11_INCLUDE_DIR CACHE)
foreach(search ${_PYBIND11_SEARCHES})
    find_path(PYBIND11_INCLUDE_DIR 
        NAMES pybind11/pybind11.h
        ${${search}}
        PATH_SUFFIXES include ${CMAKE_INSTALL_INCLUDEDIR})
endforeach()

mark_as_advanced(PYBIND11_INCLUDE_DIR)

# handle the QUIETLY and REQUIRED arguments and set PYBIND11_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PYBIND11 
    REQUIRED_VARS PYBIND11_INCLUDE_DIR)

if(PYBIND11_FOUND)
    set(PYBIND11_INCLUDE_DIRS ${PYBIND11_INCLUDE_DIR} ${PYTHON_INCLUDE_DIR})

    if(NOT TARGET PYBIND11::PYBIND11)
        add_library(PYBIND11::PYBIND11 INTERFACE IMPORTED)
        set_target_properties(PYBIND11::PYBIND11 PROPERTIES
            INTERFACE_LINK_LIBRARIES "${PYTHON_LIBRARIES}"
            INTERFACE_INCLUDE_DIRECTORIES "${PYBIND11_INCLUDE_DIRS}")
    endif()
endif()
