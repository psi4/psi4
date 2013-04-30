project(psi4)
cmake_minimum_required(VERSION 2.8)

# cmake modules are relative to this file
get_filename_component(ROOT_SRC_DIR ${CMAKE_CURRENT_LIST_FILE} PATH)
set(CMAKE_MODULE_PATH "${ROOT_SRC_DIR}/cmake/")

#
# Psi4 configuration
#
set(PSI_VERSION \"4.0\")
set(psi4datadir \"${CMAKE_INSTALL_PREFIX}\")
set(top_srcdir \"${PROJECT_SOURCE_DIR}\")
set(top_objdir \"${PROJECT_BINARY_DIR}\")

#set(PSI_VERSION "4.0")
set(PACKAGE_VERSION "4.0git")
set(PACKAGE_NAME "Psi")
set(PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")
set(PACKAGE_BUGREPORT "psicode@users.sourceforge.net")

if( CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR AND NOT MSVC_IDE )
  message(FATAL_ERROR "In-source builds are not allowed.
Please create a directory and run cmake from there, passing the path
to this source directory as the last argument.
This process created the file `CMakeCache.txt' and the directory `CMakeFiles'.
Please delete them.")
endif()

#
# OpenMP
#
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include(CheckIncludeFiles)
check_include_files(dlfcn.h HAVE_DLFCN)

include(CheckFunctionExists)
check_function_exists(MKL_Free_Buffers HAVE_MKL) # NEEDS TO BE TESTED

include(CheckSymbolExists)
check_symbol_exists(MPI_Finalize mpi.h HAVE_MPI) # NEEDS TO BE TESTED

include(cmake/Restrict.cmake)
test_restrict(restrict)

configure_file(include/psiconfig.h.cmake.in include/psiconfig.h)

# Is this needed for the script?
# TODO: Fix gitversion.py to work without setting top_srcdir...maybe use PROJECT_SOURCE_DIR
set(top_srcdir ${PROJECT_SOURCE_DIR}) # The quotes are added in the source, in this case
configure_file(src/bin/psi4/gitversion.py.in src/bin/psi4/gitversion.py)

#
# Libint configuration
# Perhaps move to libint's cmakelists.txt file
#
set(LIBINT_TMPDIR "") # I haven't quite figured out the scratch stuff yet
set(LIBINT_OPT_AM 3)    # A.M. level for ERIs
set(LIBDERIV_OPT_AM1 2) # A.M. level for 1st derivative ERIs
set(LIBDERIV_OPT_AM2 1) # A.M. level for 2nd derivative ERIs
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
#find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)   # "Finding LAPACK implies BLAS was found too"
add_definitions(-DFC_SYMBOL=2)
#
# Python Detection 
#  We don't really need the interpreter. Could we just use the PythonLibs method?

find_package(PythonInterp REQUIRED)
set(PYCONFIG ${PYTHON_EXECUTABLE}-config)
execute_process(COMMAND ${PYCONFIG} --libs     OUTPUT_VARIABLE PYTHON_LIBRARIES OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${PYCONFIG} --include  OUTPUT_VARIABLE PYTHON_INCLUDE_DIRS OUTPUT_STRIP_TRAILING_WHITESPACE)
# Turn the flags into proper lists
string(REPLACE "-I" "" PYTHON_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS})
string(REPLACE " " ";" PYTHON_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS})

##if(PYLIB STREQUAL "" OR PYINC STREQUAL "")
##   message(FATAL_ERROR "Python development libraries were not detected. Make sure that the development version "
##                       "is installed and that python-config exists in the same directory as the python executable.")
##endif()
#
#find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

#
# Boost Detection 
#
# N.B. Update boost/CMakeLists.txt if this list of components changes
set(Boost_USE_STATIC_LIBS ON)
set(Boost_COMPONENTS filesystem python regex serialization system thread )
find_package(Boost 1.47.0 COMPONENTS ${Boost_COMPONENTS})
if(NOT Boost_FOUND)
    # Add the Boost version bundled with Psi4 as a dependency
    message(WARNING "Boost not found. The pre-packaged version will be built.")
    # The version bundled with Psi4
    set(BOOSTVER boost_1_53_0)
    # Link and include information
    set(BOOSTLIBDIR ${PROJECT_BINARY_DIR}/boost/${BOOSTVER}/stage/lib)
    set(Boost_INCLUDE_DIRS ${PROJECT_BINARY_DIR}/boost/${BOOSTVER})
    foreach(component ${Boost_COMPONENTS})
        set(Boost_LIBRARIES ${Boost_LIBRARIES} "${BOOSTLIBDIR}/libboost_${component}.a")
    endforeach()
    add_subdirectory(boost)
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

set(PSILIB 3index chkpt ciomr deriv diis disp dpd fock functional iwl int mints moinfo options
           parallel plugin psio qt sapt_solver scf_solver thce trans util adc ccdensity
           ccenergy cceom cchbar cclambda ccresponse ccsort cctriples dcft deriv_wrapper detci
           dfmp2 dftsapt findif fnocc mcscf mints_wrapper mp2 mrcc occ optking psimrcc
           sapt scf scfgrad stable thermo transqt transqt2
)

add_subdirectory(src)

#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
