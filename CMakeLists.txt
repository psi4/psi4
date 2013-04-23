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
set(PACKAGE_NAME PSI)
set(PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")
set(PACKAGE_BUGREPORT "psicode@users.sourceforge.net")

if( CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR AND NOT MSVC_IDE )
  message(FATAL_ERROR "In-source builds are not allowed.
Please create a directory and run cmake from there, passing the path
to this source directory as the last argument.
This process created the file `CMakeCache.txt' and the directory `CMakeFiles'.
Please delete them.")
endif()

#set(top_srcdir "${PROJECT_SOURCE_DIR}\")
#set(top_objdir "${PROJECT_BINARY_DIR}\")

find_package(OpenMP)

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
set(PSIINC ${PROJECT_BINARY_DIR}/include 
           ${PROJECT_SOURCE_DIR}/include
           ${PROJECT_SOURCE_DIR}/src/lib
           ${PROJECT_BINARY_DIR}/src/lib)


#
# BLAS and LAPACK
#
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
add_definitions(-DFC_SYMBOL=2)
#
# Python Detection 
#  We don't really need the interpreter. Could we just use the PythonLibs method?
#
find_package(PythonInterp REQUIRED)
set(PYCONFIG ${PYTHON_EXECUTABLE}-config)
execute_process(COMMAND ${PYCONFIG} --libs     OUTPUT_VARIABLE PYLIB OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${PYCONFIG} --include  OUTPUT_VARIABLE PYINC OUTPUT_STRIP_TRAILING_WHITESPACE)
# Turn the flags into proper lists
string(REPLACE "-I" "" PYINC ${PYINC})
string(REPLACE " " ";" PYINC ${PYINC})

if(PYLIB STREQUAL "" OR PYINC STREQUAL "")
   message(FATAL_ERROR "Python development libraries were not detected. Make sure that the development version "
                       "is installed and that python-config exists in the same directory as the python executable.")
endif()

#
# Boost Detection 
#
# N.B. Update boost/CMakeLists.txt if this list of components changes
set(Boost_USE_STATIC_LIBS   ON)
set(BOOST_COMPONENTS filesystem python regex thread system serialization)
find_package(Boost 1.47.0 COMPONENTS ${BOOST_COMPONENTS})
if(Boost_FOUND)
    # Use the Boost installation detected on this system
    set(BOOSTLIB ${Boost_LIBRARIES})
    set(BOOSTINC ${Boost_INCLUDE_DIRS})
else()
    # Add the Boost version bundled with Psi4 as a dependency
    message(WARNING "Boost not found. The pre-packaged version will be built.")
    # The version bundled with Psi4
    set(BOOSTVER boost_1_53_0)
    # Link and include information
    set(BOOSTLIBDIR ${PROJECT_BINARY_DIR}/boost/${BOOSTVER}/stage/lib)
    set(BOOSTINC ${PROJECT_BINARY_DIR}/boost)
    foreach(component ${BOOST_COMPONENTS})
        set(BOOSTLIB ${BOOSTLIB} "${BOOSTLIBDIR}/libboost_${component}.a")
    endforeach()
    add_subdirectory(boost)
endif()


#
# The location of compiled libraries and executables
#
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set(PSILIB 3index chkpt ciomr deriv diis disp dpd fock functional iwl int mints moinfo options
           parallel plugin psio qt sapt_solver scf_solver trans util adc ccdensity
           ccenergy cceom cchbar cclambda ccresponse ccsort cctriples dcft deriv_wrapper detci
           dfmp2 dftsapt findif fnocc mcscf mints_wrapper mp2 mrcc occ optking psimrcc
           sapt scf scfgrad stable thermo transqt transqt2)

set(CMAKE_INCLUDE_CURRENT_DIR true)
set(INCDIRS ${BOOSTINC};${PSIINC};${PYINC})
include_directories(${INCDIRS})

add_subdirectory(src)


#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
