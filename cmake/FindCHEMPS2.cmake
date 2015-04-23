#.rst:
# FindCHEMPS2
# --------
#
# Find the native CHEMPS2 includes and library and dependencies.
# Modeled on FindZLIB.cmake
#
# IMPORTED Targets
# ^^^^^^^^^^^^^^^^
#
# This module defines :prop_tgt:`IMPORTED` target ``CHEMPS2::CHEMPS2``, if
# CHEMPS2 has been found.
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module defines the following variables:
#
# ::
#
#   CHEMPS2_INCLUDE_DIRS   - where to find chemps2/DMRG.h, etc.
#   CHEMPS2_LIBRARIES      - List of libraries when using chemps2.
#   CHEMPS2_FOUND          - True if CheMPS2 found.
#
# Hints
# ^^^^^
#
# A user may set ``CHEMPS2_ROOT`` to a CheMPS2 installation root to tell
# this module where to look.

find_package (GSL QUIET)
find_package (HDF5 QUIET)

set(_CHEMPS2_SEARCHES)

# Search CHEMPS2_ROOT first if it is set.
if(CHEMPS2_ROOT)
  set(_CHEMPS2_SEARCH_ROOT PATHS ${CHEMPS2_ROOT} NO_DEFAULT_PATH)
  list(APPEND _CHEMPS2_SEARCHES _CHEMPS2_SEARCH_ROOT)
endif()

# Normal search.
set(_CHEMPS2_SEARCH_NORMAL
  PATHS "[HKEY_LOCAL_MACHINE\\SOFTWARE\\GnuWin32\\CheMPS2;InstallPath]"
        "$ENV{PROGRAMFILES}/chemps2"
  )
list(APPEND _CHEMPS2_SEARCHES _CHEMPS2_SEARCH_NORMAL)

set(CHEMPS2_NAMES chemps2)
#set(_hold_suffix ${CMAKE_FIND_LIBRARY_SUFFIXES})
#set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")

# Try each search configuration.
foreach(search ${_CHEMPS2_SEARCHES})
    find_path(CHEMPS2_INCLUDE_DIR 
        NAMES chemps2/DMRG.h
        ${${search}}
        PATH_SUFFIXES include include/chemps2)
    find_library(CHEMPS2_LIBRARY  
        NAMES ${CHEMPS2_NAMES} 
        ${${search}} 
        PATH_SUFFIXES lib lib64)
endforeach()

mark_as_advanced(CHEMPS2_LIBRARY CHEMPS2_INCLUDE_DIR)
#set(CMAKE_FIND_LIBRARY_SUFFIXES ${_hold_suffix})

# handle the QUIETLY and REQUIRED arguments and set CHEMPS2_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CHEMPS2 
    REQUIRED_VARS CHEMPS2_LIBRARY CHEMPS2_INCLUDE_DIR
                  GSL_LIBRARIES GSL_INCLUDE_DIRS
                  HDF5_LIBRARIES HDF5_INCLUDE_DIRS)

if(CHEMPS2_FOUND)
    set(CHEMPS2_INCLUDE_DIRS ${CHEMPS2_INCLUDE_DIR} ${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
    set(CHEMPS2_LIBRARIES ${CHEMPS2_LIBRARY} ${GSL_LIBRARIES} ${HDF5_LIBRARIES})

    if(NOT TARGET CHEMPS2::CHEMPS2)
        add_library(CHEMPS2::CHEMPS2 UNKNOWN IMPORTED)
        set_target_properties(CHEMPS2::CHEMPS2 PROPERTIES
            IMPORTED_LOCATION ${CHEMPS2_LIBRARY}
            INTERFACE_LINK_LIBRARIES "${CHEMPS2_LIBRARIES}"
            INTERFACE_INCLUDE_DIRECTORIES "${CHEMPS2_INCLUDE_DIRS}")
    endif()
endif()
