cmake_minimum_required(VERSION 2.6)

# The name of our project is "PSI4". CMakeLists files in this project can
# refer to the root source directory of the project as ${PSI4_SOURCE_DIR}
# and to the root binary directory of the project as ${PSI4_BINARY_DIR}
project(PSI4)

# We need BLAS and LAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Detect Python and its libraries
find_package(PythonInterp 2.4 REQUIRED)
find_package(PythonLibs 2.4 REQUIRED)

# Some Perl is used in the tests
find_package(Perl REQUIRED)

# Determine settings for OpenMP
find_package(OpenMP)

# MPI
find_package(MPI)

# Check for C++ header files
#  This requires a cmake module to do these checks; include it.
include(CheckIncludeFileCXX)
check_include_file_cxx(cmath HAVE_CMATH)
check_include_file_cxx(cstdio HAVE_CSTDIO)
check_include_file_cxx(cstdlib HAVE_CSTDLIB)
check_include_file_cxx(cstring HAVE_CSTRING)
check_include_file_cxx(cstddef HAVE_CSTDDEF)
check_include_file_cxx(dlfcn.h HAVE_DLFCN_H)

# Check for specific symbols in libraries
include(CheckLibraryExists)
check_library_exists(m isinf "" HAVE_FUNC_ISINF)
check_library_exists(m finite "" HAVE_FUNC_FINITE)

# There are headers in include that we need to be able to find.
include_directories(include ${PSI4_BINARY_DIR}/include)

# Additional header search path as a consequence of the checks done above
include_directories(${PYTHON_INCLUDE_DIRS})

# Prepare to handle the psiconfig.h header file
set(PSI_TOP_SRCDIR ${PSI4_SOURCE_DIR})
configure_file(${PSI4_SOURCE_DIR}/include/psiconfig.h.cmake ${PSI4_BINARY_DIR}/include/psiconfig.h)

# Recurse into the src subdirectory. This does not actually cause another
# cmake executable to run. The same process will walk through the project's
# entire directory structure.
add_subdirectory(src)

