cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
# Prefer custom CMake modules
cmake_policy(SET CMP0017 OLD)

# Declare project name and programming languages
project(psi4 CXX C Fortran)

# Custom CMake modules location
list(APPEND CMAKE_MODULE_PATH
     ${CMAKE_SOURCE_DIR}/cmake
     ${CMAKE_SOURCE_DIR}/cmake/compilers
     ${CMAKE_SOURCE_DIR}/cmake/math
     )

#  Options  
option(ENABLE_BOUNDS_CHECK   "Enable bounds check"                              OFF)
option(ENABLE_CODE_COVERAGE  "Enable code coverage"                             OFF)
option(ENABLE_MPI            "Enable MPI parallelization"                       OFF)
option(ENABLE_OMP            "Enable OpenMP parallelization"                    ON)
option(ENABLE_UNIT_TESTS     "Enable compilation of unit test suite"            OFF)
option(ENABLE_AUTO_BLAS      "Enable CMake to autodetect BLAS"                  ON)
option(ENABLE_AUTO_LAPACK    "Enable CMake to autodetect LAPACK"                ON)
option(ENABLE_CSR            "Enable MKL compressed sparse row"                 OFF)
option(ENABLE_SCALAPACK      "Enable SCALAPACK"                                 OFF)
option(ENABLE_SCALASCA       "Enable scalasca profiler mode"                    OFF)
option(ENABLE_LIBERD         "Enable use of LibERD instead of LibInts"          OFF) 
option(ENABLE_STATIC_LINKING "Enable static libraries linking"                  OFF)
option(ENABLE_GPU_DFCC       "Enable GPU-DFCC plugin"                           OFF)
option(ENABLE_PLUGINS        "Enable plugins"                                   OFF)
option(ENABLE_DUMMY_PLUGIN   "Enable dummy plugin"                              OFF)
option(ENABLE_CXX11_SUPPORT  "Enable C++11 compiler support"                    ON)

set(EXTERNAL_LIBS)

# Include CMake modules as needed
include(ConfigVersion)
include(ConfigArchitecture)
include(ConfigCompilerFlags)
include(ConfigExternal)
include(ConfigOMP)

# Math: find BLAS and LAPACK
set(BLAS_LANG   "CXX")
set(LAPACK_LANG "CXX")
set(MKL_COMPILER_BINDINGS "${CMAKE_CXX_COMPILER_ID}")
include(ConfigMath)

include(ConfigMPI)
include(ConfigExplicitLibs)
include(ConfigSafeGuards)
#include(ConfigTesting)
include(GenericMacros)
#include(BinaryInfo)
#include(mergestaticlibs)
include(CheckIncludeFiles)
include(CheckFunctionExists)
include(ConfigDocumentation)

# Determine Fortran name mangling, used for external linking
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)
init_FCMangle()
# The interface files for BLAS/LAPACK should be rewritten
# to use the FCMangle.h header produced by CMake
get_fc_symbol(FC_SYMBOL)
add_definitions(-DFC_SYMBOL=${FC_SYMBOL})

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
   set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR}/external)
   set(EXTERNAL_PROJECT_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
   set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT ${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT} CACHE INTERNAL "")
else()
   set(EXTERNAL_PROJECT_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
   set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT ${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT} CACHE INTERNAL "")
endif()

if(UNIX)
   link_libraries(m)
endif(UNIX)

#Not setting ldflags was causing a problem
if(LDFLAGS)
   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LDFLAGS}")
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

# Python Detection 
find_package(Python 2.6 REQUIRED)
link_directories("${PYTHON_LIBRARIES}")
include_directories(SYSTEM "${PYTHON_INCLUDE_DIRS}")

# Boost Detection 
# Just change the Boost version number here
set(BOOSTVER 1.55.0)
set(BUILD_CUSTOM_BOOST FALSE)
# List all components needed (except mpi and unit_test_framework) here.
# mpi and unit_test_framework will be added afterwards, if needed.
list(APPEND needed_components filesystem python regex serialization system thread)
set(Boost_USE_STATIC_LIBS    ON)
set(Boost_USE_MULTITHREADED  ON)
set(Boost_USE_STATIC_RUNTIME OFF)
if(ENABLE_UNIT_TESTS)
   list(APPEND needed_components unit_test_framework)
   find_package(Boost ${BOOSTVER} COMPONENTS "${needed_components}")
else()
   find_package(Boost ${BOOSTVER} COMPONENTS "${needed_components}")
endif()
if(NOT Boost_FOUND)
   # Set also variables usually set by find_package
   message(STATUS "Boost ${BOOSTVER} not found. The pre-packaged version will be built.")
   set(BUILD_CUSTOM_BOOST TRUE)
   set(CUSTOM_BOOST_LOCATION ${PROJECT_BINARY_DIR}/boost)
   string(REGEX REPLACE "\\." "0" Boost_VERSION ${BOOSTVER})
   math(EXPR Boost_MAJOR_VERSION "${Boost_VERSION} / 100000")
   math(EXPR Boost_MINOR_VERSION "${Boost_VERSION} / 100 % 1000")
   math(EXPR Boost_SUBMINOR_VERSION "${Boost_VERSION} % 100")
   set(Boost_LIB_VERSION ${Boost_MAJOR_VERSION}_${Boost_MINOR_VERSION})
   add_subdirectory(boost)
   set(Boost_FOUND TRUE)
   set(Boost_LIBRARIES "")
   # Read documentation in FindBoost.cmake for the difference between the singular and plural forms
   set(Boost_INCLUDE_DIR  ${CUSTOM_BOOST_LOCATION}/include)
   set(Boost_INCLUDE_DIRS ${CUSTOM_BOOST_LOCATION}/include) 
   set(Boost_LIBRARY_DIR  ${CUSTOM_BOOST_LOCATION}/lib)
   set(Boost_LIBRARY_DIRS ${CUSTOM_BOOST_LOCATION}/lib)
   # This is the one that was in use in the PSI4 cmake files and is included
   # to maintain the scripts in some working order, but the former should
   # be preferred
   set(BOOSTLIBDIR ${Boost_LIBRARY_DIR})
   # We will link statically, so just set the Boost_<C>_LIBRARY for the static library 
   foreach(_component ${needed_components})
   	string(TOUPPER ${_component} _COMP)
   	set(Boost_${_COMP}_FOUND TRUE)
   	set(Boost_${_COMP}_LIBRARY libboost_${_component}-${Boost_LIB_VERSION}.a)
   	list(APPEND Boost_LIBRARIES ${Boost_${_COMP}_LIBRARY})
   endforeach()
   if(CMAKE_SYSTEM_NAME MATCHES "Linux")
      list(APPEND Boost_LIBRARIES rt)
   endif()
endif()
link_directories("${Boost_LIBRARY_DIRS}")
include_directories(SYSTEM "${Boost_INCLUDE_DIRS}")

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

enable_testing()
add_subdirectory(tests)

# Configure some scripts
configure_files()
# This has to be the very last CMake module to be included
include(ConfigInfo)
