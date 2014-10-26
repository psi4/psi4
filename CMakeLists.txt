cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

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

# Include CMake modules as needed
#include(ConfigVersion)
include(ConfigArchitecture)
include(ConfigCompilerFlags)
include(CheckCXX11Features)
#include(ConfigExternal)
include(ConfigOMP)
include(ConfigMPI)

set(BLAS_LANG   "CXX")
set(LAPACK_LANG "CXX")
set(MKL_COMPILER_BINDINGS ${CMAKE_CXX_COMPILER_ID})
include(ConfigMath)

#include(ConfigExplicitLibs)
#include(ConfigSafeGuards)
#include(ConfigTesting)
include(GenericMacros)
#include(BinaryInfo)
#include(mergestaticlibs)
include(CheckIncludeFiles)
include(CheckFunctionExists)

# Determine Fortran name mangling, used for external linking
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)
init_FCMangle()
# The interface files for BLAS/LAPACK should be rewritten
# to use the FCMangle.h header produced by CMake
get_fc_symbol(FC_SYMBOL)
add_definitions(-DFC_SYMBOL=${FC_SYMBOL})

#Not setting ldflags was causing a problem
if(LDFLAGS)
   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LDFLAGS}")
endif()

#On Ubuntu -lm needs to be last...
if(UNIX)
   link_libraries(m)
endif(UNIX)

#
# Psi4 configuration
#
set(PSI_VERSION \"4.0\")
set(PACKAGE_VERSION "4.0git")
set(PACKAGE_NAME "Psi")
set(PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")
set(PACKAGE_BUGREPORT "psicode@users.sourceforge.net")

# Dynamic linking header
check_include_files(dlfcn.h HAVE_DLFCN)
# MKL_Free_Buffers NEEDS TO BE TESTED
check_function_exists(MKL_Free_Buffers HAVE_MKL) 
# Error function
check_function_exists(erf HAVE_FUNC_ERF)
#GCC thinks it's funny "building in" types including erf
if(CMAKE_C_COMPILER_ID MATCHES GNU)
   set(HAVE_FUNC_ERF TRUE)
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

#
# BLAS and LAPACK
#
#if(NOT ${useroptLAPACKLIBS} STREQUAL ${EMPTY})
#    #RMR only have a library check at this time
#    # and then set the variables that would be set by find_package
#    set(LAPACK_LIBRARIES ${useroptLAPACKLIBS})
#    set(LAPACK_INCLUDE_DIR ${useroptLAPACKINCS})
#    message(WARNING "Skipping LAPACK detection and using these linker flags: ${useroptLAPACKLIBS}")
#    set(LAPACK_FOUND TRUE)
#else()
#    find_package(LAPACK REQUIRED)
#endif()


# Python Detection 
find_package(Python 2.6 REQUIRED)
#if(NOT useroptPYTHON STREQUAL EMPTY)
#    set(PYTHON_EXECUTABLE ${useroptPYTHON})
#    message(WARNING "Skipping Python detection and using this version: ${useroptPYTHON}")
#else()
#    find_package(PythonInterp 2.6 REQUIRED)
#endif()
#set(PYCONFIG ${PYTHON_EXECUTABLE}-config)
#execute_process(COMMAND ${PYCONFIG} --ldflags  OUTPUT_VARIABLE PYTHON_LIBRARIES OUTPUT_STRIP_TRAILING_WHITESPACE)
#execute_process(COMMAND ${PYCONFIG} --includes  OUTPUT_VARIABLE PYTHON_INCLUDE_DIRS OUTPUT_STRIP_TRAILING_WHITESPACE)
#if (("${PYTHON_LIBRARIES}" STREQUAL "") OR ("${PYTHON_INCLUDE_DIRS}" STREQUAL ""))
#    message(FATAL_ERROR "Python development libraries not detected. Make sure that the development version is installed and that 'python-config' exists in the same directory as the python executable.")
#endif()
## Strip leading whitespace & turn the flags into proper lists
#string(STRIP ${PYTHON_INCLUDE_DIRS} PYTHON_INCLUDE_DIRS)
#string(STRIP ${PYTHON_LIBRARIES} PYTHON_LIBRARIES)
#string(REPLACE "-I" "" PYTHON_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS})
#string(REPLACE " " ";" PYTHON_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS})

include_directories(${PYTHON_INCLUDE_DIRS})

#
# Boost Detection 
#
# N.B. Update boost/CMakeLists.txt if this list of components changes
set(Boost_USE_STATIC_LIBS ON)
if(MPI_FOUND)
  set(Boost_COMPONENTS filesystem python mpi regex serialization system thread )
else()
  set(Boost_COMPONENTS filesystem python regex serialization system thread )
endif()
if(NOT useroptUSE_SYSTEM_BOOST)
    message(WARNING "Skipping Boost detection, and forcing Psi to build the bundled version.")
else()
    find_package(Boost 1.55.0 COMPONENTS ${Boost_COMPONENTS})
endif()
if(NOT Boost_FOUND)
    # Add the Boost version bundled with Psi4 as a dependency
    message(WARNING "Boost not found. The pre-packaged version will be built.")
    # The version bundled with Psi4
    set(BOOSTVER boost_1_55_0)
    # Link and include information

	#This is the official cmake boost library macro
    set(BOOST_LIBRARYDIR ${PROJECT_BINARY_DIR}/boost/${BOOSTVER}/stage/lib)
    #This is the one that was in use in the PSI4 cmake files and is included
    #to maintain the scripts in some working order, but the former should
    #be preferred
    set(BOOSTLIBDIR ${BOOST_LIBRARYDIR})
    
    #Same as for the lib dir above,this is the "real" one
    set(Boost_INCLUDE_DIR ${PROJECT_BINARY_DIR}/boost/${BOOSTVER})
    set(Boost_INCLUDE_DIRS ${Boost_INCLUDE_DIR})
    
    foreach(component ${Boost_COMPONENTS})
        set(Boost_LIBRARIES ${Boost_LIBRARIES} "${BOOSTLIBDIR}/libboost_${component}.a")
    endforeach()
    if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        set(Boost_LIBRARIES ${Boost_LIBRARIES} -lrt)
    endif()
    add_subdirectory(boost)
    set(Boost_FOUND TRUE)
endif()
include_directories(${Boost_INCLUDE_DIRS})

#
# PThreads
#
find_package(Threads)
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

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


if(DEFINED USEERD)
    set(PSILIB ${PSILIB} erd)
    set(JK_ERD ${CMAKE_BINARY_DIR}/lib/liberd.a)
endif()

#If we have MPI also build the JKFactory library for awesome SCF performance
if(MPI_FOUND)
   set(BUILD_JK_FACTORY "FALSE")
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
    if(NOT DEFINED USEERD)
       set(PSILIB ${PSILIB} erd)
    endif()
    set(PSILIB ${PSILIB} -lmkl_scalapack_lp64 -lmkl_blacs_lp64)
endif()

#Goodbye python warning messages and boost warning messages...
#all remaining warnings should be addressed as they are actually caused by
#poor psi4 coding...
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd47 -wd2196 -wd411 -wd279")
endif()


#Add the doxygen documentation
add_subdirectory(doc/doxygen)

# Add the sphinx documentation
add_subdirectory(doc/sphinxman)

#
# Recursively add source directories
#
add_subdirectory(src)

#
# Add the library directory, for install purposes
#
add_subdirectory(lib)

#
# Handle creation of Makefile for plugins
#
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

#
# Echo the directory information back
#
set(RECONFIGURE_PSI "${PROJECT_BINARY_DIR}/reconfigure.sh")
message("The following command can be used to exactly reproduce this build. It can also be found in ${RECONFIGURE_PSI}:\n")
message("${CONFIGURE_CMD}\n")
configure_file(cmake/reconfigure.sh.in reconfigure.sh)

#If we are making boost add it as a dependency
if(NOT Boost_FOUND)
   add_dependencies(psi4 ${Boost_LIBRARIES})
endif()

#FILE(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/bin)
INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/bin/psi4 DESTINATION
${CMAKE_INSTALL_PREFIX}/bin)

add_subdirectory(external)

enable_testing()
add_subdirectory(tests)

