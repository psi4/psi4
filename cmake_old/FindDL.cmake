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
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# - Find libdl
# Find the native libdl includes and library
#
#  LIBDL_INCLUDE_DIR - where to find dlfcn.h
#  LIBDL_LIBRARIES   - List of libraries when using libdl.
#  LIBDL_FOUND       - True if libdl found.

if(LIBDL_INCLUDE_DIR)
  # Already in cache, be silent
  set(LIBDL_FIND_QUIETLY TRUE)
endif(LIBDL_INCLUDE_DIR)

find_path(LIBDL_INCLUDE_DIR dlfcn.h
	  PATHS
	  /usr/include
	  /usr/local/include
	  )

set(LIBDL_NAMES dl libdl ltdl libltdl)
find_library(LIBDL_LIBRARY 
	     NAMES ${LIBDL_NAMES}
	     PATHS
	     /usr/lib
	     /usr/local/lib
	     /lib
	     )

# handle the QUIETLY and REQUIRED arguments and set LIBDL_FOUND to TRUE if 
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(libdl DEFAULT_MSG LIBDL_LIBRARY LIBDL_INCLUDE_DIR)

if(LIBDL_FOUND)
  set(LIBDL_LIBRARIES ${LIBDL_LIBRARY})
else(LIBDL_FOUND)
  set(LIBDL_LIBRARIES)
endif(LIBDL_FOUND)

mark_as_advanced(LIBDL_LIBRARY LIBDL_INCLUDE_DIR)
