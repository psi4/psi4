# Useful CMake variables you should consider setting:
#    - CMAKE_C_COMPILER The C compiler
#    - CMAKE_CXX_COMPILER The C++ compiler
#    - Fortran_COMPILER The Fortran compiler (required for some external projs)
#    - CMAKE_C_FLAGS Any extra C flags you want to pass
#    - CMAKE_CXX_FLAGS Any extra C++ flags you want to pass
#    - Fortran_FLAGS Any extra Fortran flags you want to pass
#    - CMAKE_BINARY_DIR The build directory (default: build)
#    - CMAKE_BUILD_TYPE Is this Debug, Release, or Profiling
#    - BUILD_SHARED_LIBS If ON Psi4 libraries will be shared
#    - CMAKE_PREFIX_PATH List of places find_package should look for packages
#    - PYTHON_EXECUTABLE Path to the python executable you want to use
#    - PYTHON_LIBRARY Path to the python library that goes with the exe
#    - PYTHON_INCLUDE_DIR Path to the python include files    
#
#    CMake by default installs things to:
#
#      CMAKE_INSTALL_PREFIX/CMAKE_INSTALL_XDIR
#
#    where X is BIN, LIB, or INCLUDE for binaries, libraries, and include files
#    respectively.  The defaults are bin, lib, and include respectively.
#

#Need 3.1 for build_always
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(psi4 CXX C)
#set(CMAKE_CXX_STANDARD 11) Not working with Intel compilers at moment

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake
                              ${CMAKE_SOURCE_DIR}/cmake/math
)
include(CheckCompilerVersion)
include(ExternalProject)

############################# Options: Build What? #############################
option(ENABLE_AMBIT "Enables the Ambit tensor library" OFF)
option(ENABLE_CheMPS2 "Enables CheMPS2 for DMRG" OFF)
option(ENABLE_dkh "Enables DKH integrals (requires Fortran)" OFF)
option(ENABLE_liberd "Enables use of LibERD instead of LibInt (requires Fortran)" OFF)
option(ENABLE_gdma "Enables Stone's gdma multipole code (requires Fortran)" OFF)
option(ENABLE_GTFOCK "Enables GTFock (requires MPI)" OFF)
option(ENABLE_PCMSOLVER "Enables PCMSolver library (requires Fortran)" OFF)

if(ENABLE_gdma OR ENABLE_dkh OR ENABLE_liberd) # OR ENABLE_PCMSOLVER OR ENABLE_JKFACTORY
    enable_language(Fortran)
    set(Fortran_ENABLED ON)  # communicate required languages with psi4-core
    message(STATUS "Enabling Fortran")
endif()

if(ENABLE_liberd)
    message(FATAL_ERROR "LibERD will build, link, and run in Psi4 just fine, but it will give *wrong* answers for density-fitting (default) calculations. To continue with ERD, accept that you have been warned, and comment out this line in the main CMakeLists.txt")
endif()

############################# Options: Build How? ##############################
include(Psi4Macros)
option_with_print(BUILD_SHARED_LIBS "Build internally built Psi4 add-on libraries as shared, not static" OFF)
option(ENABLE_AUTO_BLAS "Enables CMake to auto-detect BLAS" ON)
option(ENABLE_AUTO_LAPACK "Enables CMake to auto-detect LAPACK" ON)
option_with_flags(ENABLE_XHOST "Enables processor-specific optimization" ON
                  "-xHost" "-march=native")
option_with_flags(ENABLE_CODE_COVERAGE "Enables details on code coverage" OFF
                  "-ftest-coverage")
option_with_flags(ENABLE_ASAN "Enables address sanitizer" OFF
                  "-fsanitize=address" "-fno-omit-frame-pointer")
option_with_flags(ENABLE_TSAN "Enables thread sanitizer" OFF
                  "-fsanitize=thread" "-fno-omit-frame-pointer -pie")
option_with_flags(ENABLE_UBSAN "Enables undefined behavior sanitizer" OFF
                  "-fsanitize=undefined" "-fno-omit-frame-pointer")
option_with_default(MAX_AM_ERI "Maximum angular momentum for integrals" 5)
option_with_default(CMAKE_BUILD_TYPE "The type of our build" Release)
option_with_default(FC_SYMBOL "The type of Fortran name mangling" 2)
option_with_default(BUILD_FPIC
          "Libraries will be compiled with position independent code" TRUE)
option_with_print(ENABLE_OPENMP "Enable OpenMP parallelization" ON)
option_with_default(CMAKE_INSTALL_LIBDIR "Directory to which libraries installed" lib)
set(CMAKE_INSTALL_MESSAGE LAZY)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include(GNUInstallDirs)
set(STAGED_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/stage${CMAKE_INSTALL_PREFIX})

######################### Process & Validate Options ###########################
include(autocmake_safeguards)
include(autocmake_omp)
#include(autocmake_static_library)

################################# Dependencies #################################
find_package(PythonInterp REQUIRED)
find_package(PythonLibs ${PYTHON_VERSION_MAJOR} REQUIRED)

set(MKL_COMPILER_BINDINGS "${CMAKE_CXX_COMPILER_ID}")
set(BLAS_LANG "CXX")
set(LAPACK_LANG "CXX")
include(ConfigMath)
#These packages don't believe in configuration files...
foreach(pkg in pybind11 pcmsolver)
   string(TOUPPER "${pkg}" PKG)
   if(NOT ${PKG}_ROOT)
      set(${PKG}_ROOT ${STAGED_INSTALL_PREFIX}/external/${pkg})
   endif()
endforeach()



add_subdirectory(external)


list(APPEND CMAKE_PREFIX_PATH  ${STAGED_INSTALL_PREFIX}/external)
################################# Main Project #################################

ExternalProject_Add(psi4_core
   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/psi4
   CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
              -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
              -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
              -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
              -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
              -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
              -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
              -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
              -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
              -DCMAKE_INSTALL_MESSAGE=${CMAKE_INSTALL_MESSAGE}
              -DMAX_AM_ERI=${MAX_AM_ERI}
              -DBUILD_FPIC=${BUILD_FPIC}
              -DPYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}
              -DPYTHON_INCLUDE_DIR=${PYTHON_INCLUDE_DIR}
              -DPYTHON_LIBRARY=${PYTHON_LIBRARY}
              -DPSI4_ROOT=${CMAKE_CURRENT_SOURCE_DIR}
              -DPCMSOLVER_ROOT=${PCMSOLVER_ROOT}
              -DPYBIND11_ROOT=${PYBIND11_ROOT}
              -DENABLE_LIBERD=${ENABLE_LIBERD}
              -DENABLE_AMBIT=${ENABLE_AMBIT}
              -DENABLE_CheMPS2=${ENABLE_CheMPS2}
              -DENABLE_dkh=${ENABLE_dkh}
              -DENABLE_liberd=${ENABLE_liberd}
              -DENABLE_gdma=${ENABLE_gdma}
              -DENABLE_GTFOCK=${ENABLE_GTFOCK}
              -DENABLE_PCMSOLVER=${ENABLE_PCMSOLVER}
              -DFortran_ENABLED=${Fortran_ENABLED}
              -DRESTRICT_KEYWORD=${RESTRICT_KEYWORD}
              -DFC_SYMBOL=${FC_SYMBOL}
              -DDESTDIR=${CMAKE_BINARY_DIR}/stage
   CMAKE_CACHE_ARGS -DCMAKE_PREFIX_PATH:PATH=${CMAKE_PREFIX_PATH}
                    -DLAPACK_LIBRARIES:LIST=${LAPACK_LIBRARIES}
                    -DLAPACK_INCLUDE_DIRS:LIST=${LAPACK_INCLUDE_DIRS}
                    -DBLAS_LIBRARIES:LIST=${BLAS_LIBRARIES}
                    -DBLAS_INCLUDE_DIRS:LIST=${BLAS_INCLUDE_DIRS}
   BUILD_ALWAYS 1
   INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
                   DESTDIR=${CMAKE_BINARY_DIR}/stage
)
add_dependencies(psi4_core ambit_external
                           chemps2_external
                           dkh_external
                           libefp_external
                           liberd_external
                           gdma_external
                           libint_external
                           pcmsolver_external
                           pybind11_external
)
add_subdirectory(doc)
enable_testing()
add_subdirectory(tests)
install(DIRECTORY ${STAGED_INSTALL_PREFIX}/
        DESTINATION ${CMAKE_INSTALL_PREFIX}
        USE_SOURCE_PERMISSIONS)
#install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/cmake
#        DESTINATION ${CMAKE_INSTALL_PREFIX} USE_SOURCE_PERMISSIONS)
