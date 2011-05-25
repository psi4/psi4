cmake_minimum_required(VERSION 2.6)

# The name of our project is "PSI4". CMakeLists files in this project can
# refer to the root source directory of the project as ${PSI4_SOURCE_DIR}
# and to the root binary directory of the project as ${PSI4_BINARY_DIR}
project(PSI4)

# Set version
set(PSI_VERSION 4.0.0-alpha)
set(PSI_BUILDID alpha)
set(PSI_BUGREPORT psicode@users.sourceforge.net)

# Any executables/libraries we create go here:
set(EXECUTABLE_OUTPUT_PATH ${PSI4_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PSI4_BINARY_DIR}/lib)

# We need BLAS and LAPACK
#find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)   # LAPACK requires BLAS for us.

# Determine Fortran name mangling for BLAS/LAPACK
#    Interesting fact: Apple's BLAS/LAPACK provides all possible manglings.
include(CheckFunctionExists)
set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARIES})
check_function_exists(dgemm FORTRAN_LC)
check_function_exists(dgemm_ FORTRAN_LCU)
check_function_exists(DGEMM FORTRAN_UC)
check_function_exists(DGEMM_ FORTRAN_UCU)

if (FORTRAN_LC)
    add_definitions(-DFC_SYMBOL=1)
elseif (FORTRAN_LCU)
    add_definitions(-DFC_SYMBOL=2)
elseif (FC_SYMBOL_UC)
    add_definitions(-DFC_SYMBOL=3)
elseif (FC_SYMBOL_UCU)
    add_definitions(-DFC_SYMBOL=4)
else(FORTRAN_LC)
    message(FATAL_ERROR "Unable to determine Fortran name mangling for BLAS/LAPACK.")
endif(FORTRAN_LC)

# Detect Python and its libraries
find_package(PythonInterp 2.4 REQUIRED)
find_package(PythonLibs 2.4 REQUIRED)

# Some Perl is used in the tests
find_package(Perl REQUIRED)

# Determine compiler flags for OpenMP
find_package(OpenMP)

# MPI (currently not linked in, modify src/bin/psi4/CMakeList.txt to add it)
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
include_directories(include ${PSI4_BINARY_DIR}/include src/lib/libint src/lib/libderiv)

# Additional header search path as a consequence of the checks done above
include_directories(${PYTHON_INCLUDE_DIRS})

# Prepare to handle the psiconfig.h header file
set(PSI_TOP_SRCDIR ${PSI4_SOURCE_DIR})
set(INSTALLEDPSIDATADIR /usr/local/psi4/lib)
configure_file(${PSI4_SOURCE_DIR}/include/psiconfig.h.cmake ${PSI4_BINARY_DIR}/include/psiconfig.h)

# Recurse into the src subdirectory. This does not actually cause another
# cmake executable to run. The same process will walk through the project's
# entire directory structure.
add_subdirectory(src)

