# FindcuEST.cmake
# ---------------
# Finds cuEST headers and libraries.
#
# Result variables:
#   cuEST_FOUND
#   cuEST_INCLUDE_DIR
#   cuEST_LIBRARY            (shared, if found)
#   cuEST_STATIC_LIBRARY     (static, if found)
#
# Imported targets:
#   cuEST::cuEST         - shared library target if available, else static
#
# Hints:
#   -DcuEST_ROOT=/path/to/prefix
#   or set environment variable cuEST_ROOT

cmake_minimum_required(VERSION 3.15)

include(FindPackageHandleStandardArgs)

# ---- user hints ----
set(_cuEST_ROOT_HINTS "")
if(DEFINED cuEST_ROOT)
  list(APPEND _cuEST_ROOT_HINTS "${cuEST_ROOT}")
endif()
if(DEFINED ENV{cuEST_ROOT})
  list(APPEND _cuEST_ROOT_HINTS "$ENV{cuEST_ROOT}")
endif()

# ---- include dir ----
find_path(cuEST_INCLUDE_DIR
  NAMES cuest.h
  HINTS ${_cuEST_ROOT_HINTS}
  PATH_SUFFIXES include
)

# ---- library dir suffixes ----
set(_cuEST_LIB_SUFFIXES
  lib
  lib/12
  lib/13
)

# Shared
find_library(cuEST_LIBRARY
  NAMES cuest
  HINTS ${_cuEST_ROOT_HINTS}
  PATH_SUFFIXES ${_cuEST_LIB_SUFFIXES}
)

# Static
find_library(cuEST_STATIC_LIBRARY
  NAMES cuest_static
  HINTS ${_cuEST_ROOT_HINTS}
  PATH_SUFFIXES ${_cuEST_LIB_SUFFIXES}
)

# ---- found? ----
# Treat the package as found if we have headers and at least one library.
find_package_handle_standard_args(cuEST
  REQUIRED_VARS cuEST_INCLUDE_DIR
  HANDLE_COMPONENTS
)

set(cuEST_FOUND FALSE)
if(cuEST_INCLUDE_DIR AND (cuEST_LIBRARY OR cuEST_STATIC_LIBRARY))
  set(cuEST_FOUND TRUE)
endif()

mark_as_advanced(cuEST_INCLUDE_DIR cuEST_LIBRARY cuEST_STATIC_LIBRARY)

# ---- imported targets ----
if(cuEST_FOUND)
  # Convenience target that picks shared if available.
  if(NOT TARGET cuEST::cuEST)
    if(cuEST_LIBRARY)
    add_library(cuEST::cuEST SHARED IMPORTED GLOBAL)
      set_target_properties(cuEST::cuEST PROPERTIES
        IMPORTED_LOCATION "${cuEST_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${cuEST_INCLUDE_DIR}"
      )
    elseif(cuEST_STATIC_LIBRARY)
      add_library(cuEST::cuEST STATIC IMPORTED GLOBAL)
      set_target_properties(cuEST::cuEST PROPERTIES
        IMPORTED_LOCATION "${cuEST_STATIC_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${cuEST_INCLUDE_DIR}"
      )
    endif()
  endif()
endif()

