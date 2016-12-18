# - Find libm
# Find the native libm includes and library
#
#  LIBM_INCLUDE_DIR - where to find math.h
#  LIBM_LIBRARIES   - List of libraries when using libm.
#  LIBM_FOUND       - True if libm found.
include(FindLibX)
find_lib_x(m math.h math.h)
