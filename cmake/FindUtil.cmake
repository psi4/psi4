# - Find libutil
# Find the native libutil includes and library
#
#  LIBUTIL_INCLUDE_DIR - where to find pty.h and utmp.h (UNIX) or util.h (Mac OS X)
#  LIBUTIL_LIBRARIES   - List of libraries when using libutil.
#  LIBUTIL_FOUND       - True if libutil found.
find_lib_x(util util.h pty.h)