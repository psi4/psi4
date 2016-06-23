# PSI4: an ab initio quantum chemistry software package
#
# Copyright (c) 2007-2013 The PSI4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# - Find librt
# Find the native librt includes and library
#
#  LIBRT_INCLUDE_DIR - where to find time.h
#  LIBRT_LIBRARIES   - List of libraries when using librt.
#  LIBRT_FOUND       - True if librt found.

if(LIBRT_INCLUDE_DIR)
  # Already in cache, be silent
  set(LIBRT_FIND_QUIETLY TRUE)
endif(LIBRT_INCLUDE_DIR)

find_path(LIBRT_INCLUDE_DIR time.h
	  PATHS
	  /usr/include
	  /usr/local/include
	  )

set(LIBRT_NAMES rt librt)
find_library(LIBRT_LIBRARY 
	     NAMES ${LIBRT_NAMES}
	     PATHS 
	     /usr/lib
	     /usr/local/lib
	     /lib
	     )

# handle the QUIETLY and REQUIRED arguments and set LIBRT_FOUND to TRUE if 
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(librt DEFAULT_MSG LIBRT_LIBRARY LIBRT_INCLUDE_DIR)

if(LIBRT_FOUND)
   set(LIBRT_LIBRARIES ${LIBRT_LIBRARY})
else(LIBRT_FOUND)
   set(LIBRT_LIBRARIES)
endif(LIBRT_FOUND)

mark_as_advanced(LIBRT_LIBRARY LIBRT_INCLUDE_DIR)
