#.rst:
# FindPCMSolver
# --------
#
# Find the native PCMSolver includes and library and dependencies.
# Modeled on FindCHEMPS2.cmake
#
# IMPORTED Targets
# ^^^^^^^^^^^^^^^^
#
# This module defines :prop_tgt:`IMPORTED` target ``PCMSolver::PCMSolver``, if
# PCMSolver has been found.
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module defines the following variables:
#
# ::
#
#   PCMSolver_INCLUDE_DIRS   - where to find pcmsolver/pcmsolver.h, etc.
#   PCMSolver_LIBRARIES      - List of libraries when using PCMSolver.
#   PCMSolver_FOUND          - True if PCMSolver found.
#
# Hints
# ^^^^^
#
# A user may set ``PCMSOLVER_ROOT`` to a PCMSolver installation root to tell
# this module where to look.

# defines CMAKE_INSTALL_LIBDIR for Debian package
include (GNUInstallDirs)

find_package(ZLIB QUIET)

set(_PCMSolver_SEARCHES)

# Search PCMSOLVER_ROOT first if it is set.
if(PCMSOLVER_ROOT)
  set(_PCMSolver_SEARCH_ROOT PATHS ${PCMSOLVER_ROOT} NO_DEFAULT_PATH)
  list(APPEND _PCMSolver_SEARCHES _PCMSolver_SEARCH_ROOT)
endif()

# Normal search.
set(_PCMSolver_SEARCH_NORMAL
  PATHS "[HKEY_LOCAL_MACHINE\\SOFTWARE\\GnuWin32\\PCMSolver;InstallPath]"
        "$ENV{PROGRAMFILES}/pcmsolver"
        "/usr"
  )
list(APPEND _PCMSolver_SEARCHES _PCMSolver_SEARCH_NORMAL)

set(PCMSolver_NAMES pcm)
#set(_hold_suffix ${CMAKE_FIND_LIBRARY_SUFFIXES})
#set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")

# Try each search configuration.
foreach(search ${_PCMSolver_SEARCHES})
    find_path(PCMSolver_INCLUDE_DIR
        NAMES PCMSolver/pcmsolver.h
        ${${search}}
        PATH_SUFFIXES include)
      find_path(PCMSolver_PARSE_DIR
        NAMES pcmsolver.py
        ${${search}}
        PATH_SUFFIXES bin)
    find_library(PCMSolver_LIBRARY
        NAMES ${PCMSolver_NAMES}
        ${${search}}
        PATH_SUFFIXES lib lib64 ${CMAKE_INSTALL_LIBDIR})
endforeach()

mark_as_advanced(PCMSolver_LIBRARY PCMSolver_INCLUDE_DIR PCMSolver_PARSE_DIR)
#set(CMAKE_FIND_LIBRARY_SUFFIXES ${_hold_suffix})
if(DETECT_EXTERNAL_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif()

# handle the QUIETLY and REQUIRED arguments and set PCMSolver_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PCMSolver
    FOUND_VAR PCMSolver_FOUND
    REQUIRED_VARS PCMSolver_LIBRARY PCMSolver_INCLUDE_DIR PCMSolver_PARSE_DIR
                  ZLIB_LIBRARIES ZLIB_INCLUDE_DIRS)

if(PCMSolver_FOUND)
  set(PCMSolver_INCLUDE_DIRS ${PCMSolver_INCLUDE_DIR} ${ZLIB_INCLUDE_DIRS})
  set(PCMSolver_LIBRARIES ${PCMSolver_LIBRARY} ${ZLIB_LIBRARIES})

  if(NOT TARGET PCMSolver::PCMSolver)
    add_library(PCMSolver::PCMSolver UNKNOWN IMPORTED)
    set_target_properties(PCMSolver::PCMSolver PROPERTIES
      IMPORTED_LOCATION ${PCMSolver_LIBRARY}
      #INTERFACE_LINK_LIBRARIES "${PCMSolver_LIBRARIES}"
      INTERFACE_LINK_LIBRARIES ${ZLIB_LIBRARIES}
      INTERFACE_COMPILE_DEFINITIONS USING_PCMSolver
      INTERFACE_INCLUDE_DIRECTORIES "${PCMSolver_INCLUDE_DIRS}")
  endif()
  include_directories(SYSTEM "${PCMSolver_INCLUDE_DIRS}")
endif()
