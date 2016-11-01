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

# defines CMAKE_INSTALL_LIBDIR for Debian package
include (GNUInstallDirs)

#find_package (HDF5 QUIET)
set(LIBXC_ROOT "/Users/daniel/Gits/psixc/objdir/stage/usr/local/external")

set(_LIBXC_SEARCHES)

# Search CHEMPS2_ROOT first if it is set.
if(LIBXC_ROOT)
  set(_LIBXC_SEARCH_ROOT PATHS ${LIBXC_ROOT} NO_DEFAULT_PATH)
  list(APPEND _LIBXC_SEARCHES _LIBXC_SEARCH_ROOT)
endif()

# Normal search.
set(_LIBXC_SEARCH_NORMAL
  PATHS "[HKEY_LOCAL_MACHINE\\SOFTWARE\\GnuWin32\\LIBXC;InstallPath]"
        "$ENV{PROGRAMFILES}/libxc"
        "/usr"
  )
list(APPEND _LIBXC_SEARCHES _LIBXC_SEARCH_NORMAL)
#list(APPEND _LIBXC_SEARCHES "/Users/daniel/Gits/psixc/objdir/stage/usr/local")

set(LIBXC_NAMES xc)
#set(_hold_suffix ${CMAKE_FIND_LIBRARY_SUFFIXES})
if(DETECT_EXTERNAL_STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif()

# Try each search configuration.
foreach(search ${_LIBXC_SEARCHES})
    find_path(LIBXC_INCLUDE_DIR 
        NAMES libxc/xc.h
        ${${search}}
        PATH_SUFFIXES include external/include)
    find_library(LIBXC_LIBRARY  
        NAMES ${LIBXC_NAMES} 
        ${${search}} 
        PATH_SUFFIXES lib lib64 ${CMAKE_INSTALL_LIBDIR})
endforeach()

mark_as_advanced(LIBXC_LIBRARY LIBXC_INCLUDE_DIR)
#set(CMAKE_FIND_LIBRARY_SUFFIXES ${_hold_suffix})

# handle the QUIETLY and REQUIRED arguments and set CHEMPS2_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LIBXC
    REQUIRED_VARS LIBXC_LIBRARY LIBXC_INCLUDE_DIR)

if(LIBXC_FOUND)
    set(LIBXC_INCLUDE_DIRS ${LIBXC_INCLUDE_DIR})
    set(LIBXC_LIBRARIES ${LIBXC_LIBRARY})

    if(NOT TARGET libxc::xc)
        add_library(libxc::xc STATIC IMPORTED)
        set_target_properties(libxc::xc PROPERTIES
            IMPORTED_LOCATION ${LIBXC_LIBRARY}
            INTERFACE_LINK_LIBRARIES "${LIBXC_LIBRARIES}"
            INTERFACE_INCLUDE_DIRECTORIES "${LIBXC_INCLUDE_DIRS}")
    endif()
endif()
