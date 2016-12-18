# - Find libdl
# Find the native libdl includes and library
#
#  LIBDL_INCLUDE_DIR - where to find dlfcn.h
#  LIBDL_LIBRARIES   - List of libraries when using libdl.
#  LIBDL_FOUND       - True if libdl found.
include(FindLibX)
find_lib_x(dl dlfcn.h dlfcn.h ltdl)
