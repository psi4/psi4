# - Find librt
# Find the native librt includes and library
#
#  LIBRT_INCLUDE_DIR - where to find time.h
#  LIBRT_LIBRARIES   - List of libraries when using librt.
#  LIBRT_FOUND       - True if librt found.
find_lib_x(rt time.h time.h)

