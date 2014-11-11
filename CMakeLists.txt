cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# Declare project name and programming languages
project(psi4 CXX C Fortran)

# Custom CMake modules location
list(APPEND CMAKE_MODULE_PATH
     ${CMAKE_SOURCE_DIR}/cmake
     ${CMAKE_SOURCE_DIR}/cmake/compilers
     ${CMAKE_SOURCE_DIR}/cmake/math
     ${CMAKE_SOURCE_DIR}/cmake/testing
     )

#  Options  
option(ENABLE_BOUNDS_CHECK   "Enable bounds check"                         OFF)
option(ENABLE_CODE_COVERAGE  "Enable code coverage"                        OFF)
option(ENABLE_MPI            "Enable MPI parallelization"                  OFF)
option(ENABLE_OMP            "Enable OpenMP parallelization"               ON)
option(ENABLE_UNIT_TESTS     "Enable compilation of unit test suite"       OFF)
option(ENABLE_AUTO_BLAS      "Enable CMake to autodetect BLAS"             ON)
option(ENABLE_AUTO_LAPACK    "Enable CMake to autodetect LAPACK"           ON)
option(ENABLE_ACCELERATE     "Enable use of Mac OS X Accelerate Framework" OFF)
option(ENABLE_CSR            "Enable MKL compressed sparse row"            OFF) 
option(ENABLE_SCALAPACK      "Enable SCALAPACK"                            OFF)
option(ENABLE_SCALASCA       "Enable scalasca profiler mode"               OFF)
option(ENABLE_LIBERD         "Enable use of LibERD instead of LibInts"     OFF)
option(ENABLE_STATIC_LINKING "Enable static libraries linking"             OFF)
option(ENABLE_GPU_DFCC       "Enable GPU-DFCC plugin"                      OFF)
option(ENABLE_PLUGINS        "Enable plugins"                              ON)
option(ENABLE_DUMMY_PLUGIN   "Enable dummy plugin"                         OFF)
option(ENABLE_CXX11_SUPPORT  "Enable C++11 compiler support"               ON)

set(EXTERNAL_LIBS)

# Include CMake modules as needed
include(ConfigVersion)
include(ConfigArchitecture)
include(ConfigCompilerFlags)
include(ConfigExternal)
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

# Determine Fortran name mangling, used for external linking
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)
init_FCMangle()
# This is to use the CMake generated macros and not those based on FC_SYMBOL
add_definitions(-DUSE_FCMANGLE_H)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
   set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR}/external)
   set(EXTERNAL_PROJECT_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
   set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT ${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT} CACHE INTERNAL "")
else()
   set(EXTERNAL_PROJECT_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
   set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT ${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT} CACHE INTERNAL "")
endif()

# Python Detection 
find_package(PythonInterp 2.6 REQUIRED)
find_package(PythonLibs 2.6 REQUIRED)
link_directories("${PYTHON_LIBRARIES}")
include_directories(SYSTEM "${PYTHON_INCLUDE_DIRS}")

# Boost Detection 
# We need Boost.Python, so this has to come _after_ Python detection
include(ConfigBoost)
link_directories("${Boost_LIBRARY_DIRS}")
include_directories(SYSTEM "${Boost_INCLUDE_DIRS}")


if(UNIX)
   link_libraries(m)
endif(UNIX)

#Not setting ldflags was causing a problem
if(LDFLAGS)
   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LDFLAGS}")
endif()

if(ENABLE_PLUGINS)
   message(STATUS "Plugins support ENABLED")
endif()   

# Dynamic linking header
check_include_files(dlfcn.h HAVE_DLFCN_H)
# MKL_Free_Buffers NEEDS TO BE TESTED
check_function_exists(MKL_Free_Buffers HAVE_MKL) 

# Detection of erf, __builtin_expect, __builtin_prefetch, __builtin_constant_p
# This step is skipped for GCC and Clang: they are all defined for but for some reasons not detected by CMake
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
   message(WARNING "Building/Using LibERD instead of LibInts!!!!")
   add_definitions(-DHAVE_ERD)
endif()

# Is this needed for the script?
# TODO: Fix gitversion.py to work without setting top_srcdir...maybe use 
# PROJECT_SOURCE_DIR?
set(top_srcdir ${PROJECT_SOURCE_DIR}) # The quotes are added in the source, in this case
configure_file(src/bin/psi4/gitversion.py.in src/bin/psi4/gitversion.py)

#
# Libint configuration
# Perhaps move to libint's cmakelists.txt file

math(EXPR LIBDERIV_OPT_AM1 ${LIBINT_OPT_AM}-1) # A.M. level for 1st derivative ERIs
math(EXPR LIBDERIV_OPT_AM2 ${LIBINT_OPT_AM}-2) # A.M. level for 2nd derivative ERIs
math(EXPR LIBINT_NEW_AM ${LIBINT_OPT_AM}*2)
math(EXPR LIBDERIV_NEW_AM1 ${LIBDERIV_OPT_AM1}*2)
math(EXPR LIBDERIV_NEW_AM2 ${LIBDERIV_OPT_AM2}*2)
math(EXPR LIBDERIV_NEW_AM12 ${LIBDERIV_OPT_AM2}*2)
configure_file(src/lib/libint/libint_config.h.in src/lib/libint/libint_config.h)
configure_file(src/lib/libderiv/libderiv_config.h.in src/lib/libderiv/libderiv_config.h)

# Psi4-specific includes and libraries
set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(
    ${PROJECT_BINARY_DIR}/include 
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src/lib
    ${PROJECT_BINARY_DIR}/src/lib
)

# PThreads
find_package(Threads)

#
# The location of compiled libraries and executables
#
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set(PSILIB mints_wrapper dfmp2 dfocc scf ccenergy ccsort psimrcc transqt2 
           cctriples scf_solver fock dcft lmp2 mcscf sapt dftsapt sapt_solver
           cchbar cclambda ccdensity transqt ccresponse detci detcas occ mrcc fnocc
           cceom adc thermo functional disp thce 3index deriv_wrapper optking
           findif mints trans dpd chkpt iwl psio qt ciomr options moinfo util
           stable deriv scfgrad int util diis plugin parallel frag)


#RMR---This is our plugin workaround
#Here we start a list of plugins
set(PLUGINS "")
add_subdirectory(plugins)
#Add plugins to the libraries Psi depends on
set(PSILIB ${PSILIB} ${PLUGINS})


if(ENABLE_LIBERD)
    set(PSILIB ${PSILIB} erd)
    set(JK_ERD ${CMAKE_BINARY_DIR}/lib/liberd.a)
endif()

#If we have MPI also build the JKFactory library for awesome SCF performance
if(MPI_FOUND)
   set(BUILD_JK_FACTORY "TRUE")
endif()

if(BUILD_JK_FACTORY)
    add_definitions(-DHAVE_JK_FACTORY)
    set(JKROOT ${CMAKE_CURRENT_SOURCE_DIR}/src/lib/libJKFactory)
    set(JK_BLAS_INC ${LAPACK_INCLUDE_DIR})
    set(PSILIB ${PSILIB} JKFactory)
    set(PSILIB ${PSILIB} pfock)
    set(PSILIB ${PSILIB} cint)
    set(PSILIB ${PSILIB} ${CMAKE_BINARY_DIR}/lib/libga.a)
    set(PSILIB ${PSILIB} ${CMAKE_BINARY_DIR}/lib/libarmci.a)
    set(PSILIB ${PSILIB} oed)
    if(NOT DEFINED ENABLE_LIBERD)
       set(PSILIB ${PSILIB} erd)
    endif()
    set(PSILIB ${PSILIB} -lmkl_scalapack_lp64 -lmkl_blacs_lp64)
endif()

#Add the documentation subdirectory
add_subdirectory(doc)

# Recursively add source directories
add_subdirectory(src)

# Add the library directory, for install purposes
add_subdirectory(lib)

# Handle creation of Makefile for plugins
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
get_property(defs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY COMPILE_DEFINITIONS)
foreach(dir ${dirs})
    set(PLUGIN_INCLUDE_DIRECTORIES "-I${dir} ${PLUGIN_INCLUDE_DIRECTORIES}")
endforeach()
foreach(def ${defs})
    set(PLUGIN_DEFINES "-D${def} ${PLUGIN_DEFINES}")
endforeach()
string(STRIP ${CMAKE_CXX_FLAGS} CMAKE_CXX_FLAGS)
string(STRIP ${PLUGIN_INCLUDE_DIRECTORIES} PLUGIN_INCLUDE_DIRECTORIES)
string(STRIP ${PLUGIN_DEFINES} PLUGIN_DEFINES)
configure_file(include/psiconfig.h.cmake.in include/psiconfig.h)

# Make known that the setup command given and the corresponding CMake
# line are available in the file setup_command in the build directory
message(STATUS "Use the commands in setup_command to reproduce this build")

#If we are making Boost add it as a dependency
if(BUILD_CUSTOM_BOOST)
   add_dependencies(psi4 custom_boost)
endif()

#FILE(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/bin)
INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/bin/psi4 DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

add_subdirectory(external)

# Configure some scripts
configure_files()
# Configure testing
# This must come after ConfigDocumentation, as it needs Perl detection 
# It also needs Python so it must go _after_ Python detection!
include(ConfigTesting)
# This has to be the very last CMake module to be included
include(ConfigInfo)
