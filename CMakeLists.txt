cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# Declare project name and programming languages
project(psi4 CXX C)

# Custom CMake modules location
list(APPEND CMAKE_MODULE_PATH
     ${CMAKE_SOURCE_DIR}/cmake
     ${CMAKE_SOURCE_DIR}/cmake/compilers
     ${CMAKE_SOURCE_DIR}/cmake/math
     ${CMAKE_SOURCE_DIR}/cmake/testing
     )

#  Options
option(ENABLE_BOUNDS_CHECK   "Enable bounds check"                         OFF)
option(ENABLE_MPI            "Enable MPI parallelization"                  OFF)
option(ENABLE_OMP            "Enable OpenMP parallelization"               ON)
option(ENABLE_UNIT_TESTS     "Enable compilation of unit test suite"       OFF)
option(ENABLE_AUTO_BLAS      "Enable CMake to autodetect BLAS"             ON)
option(ENABLE_AUTO_LAPACK    "Enable CMake to autodetect LAPACK"           ON)
option(ENABLE_ACCELERATE     "Enable use of Mac OS X Accelerate Framework" OFF)
option(ENABLE_CSR            "Enable MKL compressed sparse row"            OFF)
option(ENABLE_SCALAPACK      "Enable SCALAPACK"                            OFF)
option(ENABLE_SCALASCA       "Enable scalasca profiler mode"               OFF)
option(ENABLE_STATIC_LINKING "Enable mostly static libraries linking"      OFF)
option(ENABLE_XHOST          "Enable processor-specific optimizations"     ON)
option(ENABLE_GPU_DFCC       "Enable GPU-DFCC plugin"                      OFF)
option(ENABLE_PLUGINS        "Enable compilation of example plugins"       OFF)
option(ENABLE_DUMMY_PLUGIN   "Enable dummy plugin"                         OFF)
option(ENABLE_CHEMPS2        "Enable CheMPS2 for DMRG"                     OFF)
option(ENABLE_CXX11_SUPPORT  "Enable C++11 compiler support"               ON)
# Modules requiring Fortran
if(ENABLE_LIBERD OR ENABLE_PCMSOLVER OR ENABLE_JKFACTORY OR ENABLE_DKH)
   if(NOT CMAKE_Fortran_COMPILER)
      message(FATAL_ERROR "The Fortran compiler has to be explicitly specified!")
   endif()
endif()
include(CMakeDependentOption)
CMAKE_DEPENDENT_OPTION(
       ENABLE_DKH            "Enable DKH in libmints"                      ON
                             "CMAKE_Fortran_COMPILER"                      OFF)
CMAKE_DEPENDENT_OPTION(
       ENABLE_PCMSOLVER      "Enable PCMSolver library"                    OFF
                             "CMAKE_Fortran_COMPILER"                      OFF)
CMAKE_DEPENDENT_OPTION(
       ENABLE_LIBERD         "Enable use of LibERD instead of LibInts"     OFF
                             "CMAKE_Fortran_COMPILER"                      OFF)
# Options for dashboard builds
option(ENABLE_CODE_COVERAGE  "Enable code coverage"                        OFF)
option(ENABLE_MEMCHECK       "Enable Valgrind memory check"                OFF)
option(ENABLE_ASAN           "Enable address sanitizer"                    OFF)
option(ENABLE_MSAN           "Enable memory sanitizer"                     OFF)
option(ENABLE_TSAN           "Enable thread sanitizer"                     OFF)
option(ENABLE_UBSAN          "Enable undefined behaviour sanitizer"        OFF)
option(ENABLE_JKFACTORY      "Enable the distributed JKFactory (requires MPI)"            OFF)

set(EXTERNAL_LIBS)

# Decide if Fortran is needed
if(ENABLE_LIBERD OR ENABLE_PCMSOLVER OR ENABLE_JKFACTORY OR ENABLE_DKH)
   include(FortranEnabler)
   enable_language(Fortran)
   # add_definitions(-DHAVE_FORTRAN)  # unused
   # This is to use the CMake generated macros and not those based on FC_SYMBOL
   add_definitions(-DUSE_FCMANGLE_H)
   set(FORTRAN_ENABLED TRUE)
   message(STATUS "Fortran ENABLED")
   fortran_enabler()
endif()

# Include CMake modules as needed
include(ConfigVersion)
include(ConfigArchitecture)
include(ConfigCompilerFlags)
include(PerTargetCompilerFlags) # This module has to be included **after** ConfigCompilerFlags.cmake
include(ConfigOMP)

# Math: find BLAS and LAPACK
set(BLAS_LANG   "CXX")
set(BLAS_FOUND FALSE)
set(LAPACK_LANG "CXX")
set(LAPACK_FOUND FALSE)
set(MKL_COMPILER_BINDINGS "${CMAKE_CXX_COMPILER_ID}")
include(ConfigMath)

include(ConfigMPI)
include(ConfigExplicitLibs)
include(ConfigSafeGuards)
include(GenericMacros)
#include(BinaryInfo)
#include(mergestaticlibs)
include(CheckIncludeFiles)
include(CheckFunctionExists)
# Here we look for Perl, Sphinx, Doxygen and LaTeX
include(ConfigDocumentation)

if(NOT FORTRAN_ENABLED)
   # For linking to math subroutines use the FC_SYMBOL strategy, if Fortran was not enabled
   set(FC_SYMBOL 2)
   add_definitions(-DFC_SYMBOL=${FC_SYMBOL})
endif()

if(CMAKE_INSTALL_PREFIX STREQUAL "/usr/local/psi4")
   set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT TRUE CACHE INTERNAL "")
endif()

# This section attempts to find some system libraries:
#  a. libutil : needed for linking against static libpython (Linux only)
#  b. libm    : needed for lib{int,deriv}_compiler (Linux only)
#  c. librt   : needed for the timers (Linux only)
#  d. libdl   : needed for dynamic linking (Linux and Mac OS X)
#  e. Threads   (Linux and Mac OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
   find_package(Util REQUIRED)
   include_directories(SYSTEM "${LIBUTIL_INCLUDE_DIR}")
   link_libraries("${LIBUTIL_LIBRARIES}")

   find_package(M REQUIRED)
   include_directories(SYSTEM "${LIBM_INCLUDE_DIR}")
   link_libraries("${LIBM_LIBRARIES}")

   find_package(RT REQUIRED)
   include_directories(SYSTEM "${LIBRT_INCLUDE_DIR}")
   link_libraries("${LIBRT_LIBRARIES}")
endif()

find_package(DL REQUIRED)
if(LIBDL_FOUND)
   set(HAVE_DLFCN_H "1")
   include_directories(SYSTEM "${LIBDL_INCLUDE_DIR}")
   link_libraries("${LIBDL_LIBRARIES}")
endif()

find_package(Threads REQUIRED)
link_libraries("${CMAKE_THREAD_LIBS_INIT}")

# Python Detection
include(ConfigPython)
link_libraries("${PYTHON_LIBRARY}")
include_directories(SYSTEM "${PYTHON_INCLUDE_DIR}")

# Boost Detection
# We need Boost.Python, so this has to come _after_ Python detection
include(ConfigBoost)
link_directories("${Boost_LIBRARY_DIRS}")
include_directories(SYSTEM "${Boost_INCLUDE_DIRS}")

# CheMPS2 Detection
if(ENABLE_CHEMPS2)
   include(ConfigChemps2)
endif()

# Append the suffix given from input to all generated executables
if(EXECUTABLE_SUFFIX)
   set(CMAKE_EXECUTABLE_SUFFIX "${EXECUTABLE_SUFFIX}")
   message(STATUS "Suffix ${CMAKE_EXECUTABLE_SUFFIX} will be appended to executables")
endif()

if(LDFLAGS)
   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LDFLAGS}")
endif()

if(ENABLE_PLUGINS)
   message(STATUS "Plugins support ENABLED")
endif()

# MKL_Free_Buffers NEEDS TO BE TESTED
check_function_exists(MKL_Free_Buffers HAVE_MKL)

# Check for mm_malloc.h for the _mm_malloc and _mm_free aligned memory allocation
check_include_files(mm_malloc.h HAVE_MM_MALLOC_H)
if(HAVE_MM_MALLOC_H)
    add_definitions(-DHAVE_MM_MALLOC_H)
endif()

# Detection of erf, __builtin_expect, __builtin_prefetch, __builtin_constant_p
# This step is skipped for GCC and Clang: they are all defined but for some reasons not detected by CMake
if(CMAKE_CXX_COMPILER_ID MATCHES GNU OR CMAKE_CXX_COMPILER_ID MATCHES Clang)
   set(HAVE_BUILTIN_EXPECT     TRUE)
   set(HAVE_BUILTIN_PREFETCH   TRUE)
   set(HAVE_BUILTIN_CONSTANT_P TRUE)
   set(HAVE_FUNC_ERF           TRUE)
else()
   # Check that __builtin_expect is available
   check_function_exists(__builtin_expect HAVE_BUILTIN_EXPECT)
   # Check that __builtin_prefetch is available
   check_function_exists(__builtin_prefetch HAVE_BUILTIN_PREFETCH)
   # Check that __builtin_constant_p is available
   check_function_exists(__builtin_constant_p HAVE_BUILTIN_CONSTANT_P)
   # Error function
   check_function_exists(erf HAVE_FUNC_ERF)
endif()

if(NOT HAVE_FUNC_ERF)
   message(FATAL_ERROR "ERF was not found")
endif()

if(ENABLE_LIBERD)
   message(STATUS "Building LibERD in addition to LibInts! When able LibERD will be used by default.")
   add_definitions(-DHAVE_ERD)
endif()

if(ENABLE_DKH)
   message(STATUS "DKH in LibMints ENABLED")
   add_definitions(-DHAVE_DKH)
endif()

# Libint configuration
# Set default for LIBINT_OPT_AM if not set by setup script.
if(NOT LIBINT_OPT_AM)
   set(LIBINT_OPT_AM 5)
endif()

# Psi4-specific includes and libraries
set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(
    ${PROJECT_BINARY_DIR}/include
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src/lib
    ${PROJECT_BINARY_DIR}/src/lib
    ${PROJECT_BINARY_DIR}/interfaces/include
)


###set(headers_list "")
###set(headers_list include/compiler.h include/psi4-dec.h include/exception.h include/process.h)
###install_list_FILES("${headers_list}" include)


#
# The location of compiled libraries and executables
#
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set_property(GLOBAL PROPERTY PSILIB)

# PCMSolver
if(ENABLE_PCMSOLVER)
   find_package(ZLIB)
    if (NOT ZLIB_FOUND)
        message(FATAL_ERROR "No Zlib, no PCMSolver. Skip with -DENABLE_PCMSOLVER=OFF")
        # TODO ability to build psi4 against pre-built pcmsolver library like chemps2
        #message(FATAL_ERROR "No Zlib, no PCMSolver. Build against existing with -DPCMSOLVER_ROOT=/path/to/pcmsolver or skip with -DENABLE_PCMSOLVER=OFF")
   else()
      add_definitions(-DHAVE_PCMSOLVER)
      message(STATUS "Polarizable Continuum Model via PCMSolver ENABLED")
      link_directories(${PROJECT_BINARY_DIR}/interfaces/lib)
   endif()
endif()

#If we have MPI we may want to also build JKFactory which makes J and K's
#in distributed parallel
set(BUILD_JK_FACTORY FALSE)
if(MPI_FOUND AND ENABLE_JKFACTORY)
   set(BUILD_JK_FACTORY TRUE)
endif()
if(BUILD_JK_FACTORY)
   set(JK_BLAS_LIB "")
   set(JK_BLAS_INC "")
   if(EXPLICIT_BLAS_LIB)
     set(JK_BLAS_LIB ${EXPLICIT_BLAS_LIB})
     set(JK_BLAS_INC ${CMAKE_CXX_FLAGS})
   endif()
   message(STATUS "LibJK BLAS Paths: ${JK_BLAS_INC} ${JK_BLAS_LIB}")
   set(HAVE_GLOBAL_ARRAYS FALSE)
   add_definitions(-DHAVE_JK_FACTORY)
   set(JKROOT ${CMAKE_CURRENT_SOURCE_DIR}/src/lib/libJKFactory)
endif()

# Add the documentation subdirectory
add_subdirectory(doc)

# Recursively add interface directories (some containing dependencies of src)
add_subdirectory(interfaces)

# Recursively add source directories
add_subdirectory(src)

# Add the library directory, for install purposes
add_subdirectory(lib)

# Handle creation of Makefile for plugins
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES) ###
get_property(defs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY COMPILE_DEFINITIONS)
###set(dirs "${CMAKE_INSTALL_PREFIX}/include" "${PYTHON_INCLUDE_DIR}"  "${Boost_INCLUDE_DIRS}")
foreach(dir ${dirs})
    set(PLUGIN_INCLUDE_DIRECTORIES "-I${dir} ${PLUGIN_INCLUDE_DIRECTORIES}")
endforeach()
foreach(def ${defs})
    set(PLUGIN_DEFINES "-D${def} ${PLUGIN_DEFINES}")
endforeach()
set(CXX_FLAGS_PLUGIN "")
set(PLUGIN_LDFLAGS "${CMAKE_CXX_LINK_FLAGS} ${LIBC_INTERJECT}")
include(ConfigPlugins) # This module has to be included **after** ConfigCompilerFlags.cmake
string(STRIP ${CXX_FLAGS_PLUGIN} CXX_FLAGS_PLUGIN)
string(STRIP ${PLUGIN_INCLUDE_DIRECTORIES} PLUGIN_INCLUDE_DIRECTORIES)
string(STRIP ${PLUGIN_DEFINES} PLUGIN_DEFINES)
string(STRIP ${PLUGIN_LDFLAGS} PLUGIN_LDFLAGS)
configure_file(include/psiconfig.h.cmake.in include/psiconfig.h)
###install_list_FILES("${CMAKE_CURRENT_BINARY_DIR}/include/psiconfig.h" include/)

# Make known that the setup command given and the corresponding CMake
# line are available in the file setup_command in the build directory
message(STATUS "Use the commands in setup_command to reproduce this build")

#If we are making Boost add it as a dependency
if(BUILD_CUSTOM_BOOST)
   add_dependencies(psi4 custom_boost)
endif()

#FILE(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/bin)
INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/bin/psi4 DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

# Install samples directory
INSTALL(
  DIRECTORY "samples"
  DESTINATION ${CMAKE_INSTALL_PREFIX}/share/psi
  USE_SOURCE_PERMISSIONS
  PATTERN "example_psi4rc_file" EXCLUDE
)

# Configure some scripts
configure_files()
# Configure testing
# This must come after ConfigDocumentation, as it needs Perl detection
# It also needs Python so it must go _after_ Python detection!
include(ConfigTesting)
# This must come after ConfigTesting to register plugins tests
if(ENABLE_PLUGINS)
   add_subdirectory(plugins)
endif()
# This has to be the very last CMake module to be included
include(ConfigInfo)
