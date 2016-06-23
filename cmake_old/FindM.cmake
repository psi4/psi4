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
# - Find libm
# Find the native libm includes and library
#
#  LIBM_INCLUDE_DIR - where to find math.h
#  LIBM_LIBRARIES   - List of libraries when using libm.
#  LIBM_FOUND       - True if libm found.

if(LIBM_INCLUDE_DIR)
  # Already in cache, be silent
  set(LIBM_FIND_QUIETLY TRUE)
endif(LIBM_INCLUDE_DIR)

find_path(LIBM_INCLUDE_DIR math.h
	  PATHS
	  /usr/include
	  /usr/local/include
	  )

set(LIBM_NAMES m libm)
find_library(LIBM_LIBRARY 
	     NAMES ${LIBM_NAMES}
	     PATHS 
	     /usr/lib
	     /usr/local/lib
	     /lib
	     )

# handle the QUIETLY and REQUIRED arguments and set LIBM_FOUND to TRUE if 
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(libm DEFAULT_MSG LIBM_LIBRARY LIBM_INCLUDE_DIR)

if(LIBM_FOUND)
   set(LIBM_LIBRARIES ${LIBM_LIBRARY})
else(LIBM_FOUND)
   set(LIBM_LIBRARIES)
endif(LIBM_FOUND)

mark_as_advanced(LIBM_LIBRARY LIBM_INCLUDE_DIR)
