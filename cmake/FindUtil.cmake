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
# - Find libutil
# Find the native libutil includes and library
#
#  LIBUTIL_INCLUDE_DIR - where to find pty.h and utmp.h (UNIX) or util.h (Mac OS X)
#  LIBUTIL_LIBRARIES   - List of libraries when using libutil.
#  LIBUTIL_FOUND       - True if libutil found.
find_lib_x(util util.h pty.h)