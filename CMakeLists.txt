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
include(GNUInstallDirs)
#set(CMAKE_CXX_STANDARD 11) Not working with Intel compilers at moment
include(ExternalProject)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake
                              ${CMAKE_SOURCE_DIR}/cmake/math
)
include(Psi4Macros)
no_in_source()

############################## Options #########################################
option(ENABLE_PCMSOLVER "Enables the PCMSolver module" OFF)
option(ENABLE_GTFOCK "Enables GTFock (requires MPI)" OFF)
option(ENABLE_CHEMPS2 "Enables CheMPS2 for DMRG" OFF)
option(ENABLE_DKH "Enables DKH integrals" OFF)
option(ENABLE_LIBERD "Enables ERD instead of libint" OFF)
option(ENABLE_GDMA "Enables Stone's GDMA multipole code" OFF)
option(ENABLE_AMBIT "Enables the Ambit tensor library" OFF)
option(BUILD_SHARED_LIBS "All Psi4 libraries will be shared objects" OFF)
option(ENABLE_PYTHON_INTERFACE "Final library will be shared" ON)
option(ENABLE_AUTO_BLAS "Enables CMake to auto-detect BLAS" ON)
option(ENABLE_AUTO_LAPACK "Enables CMake to auto-detect LAPACK" ON)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
   option(ENABLE_ACCELERATE "Forces use of accelerate framework on Mac" ON)
endif()
option_with_flags(ENABLE_XHOST "Enables processor-specific optimization" ON 
                  "-xHost")
option_with_flags(ENABLE_CODE_COVERAGE "Enables details on code coverage" OFF
                  "-ftest-coverage")
option_with_flags(ENABLE_ASAN "Enables address sanitizer" OFF 
                  "-fsanitize=address" "-fno-omit-frame-pointer")
option_with_flags(ENABLE_TSAN "Enables thread sanitizer" OFF
                  "-fsanitize=thread" "-fno-omit-frame-pointer -pie")
option_with_flags(ENABLE_UBSAN "Enables undefined behavior sanitizer" OFF
                  "-fsanitize=undefined" "-fno-omit-frame-pointer")
option_with_default(LIBINT_OPT_AM "Maximum angular momentum for integrals" 5)
option_with_default(CMAKE_BUILD_TYPE "The type of our build" Release)
option_with_default(FC_SYMBOL "The type of Fortran name mangling" 2)
option_with_default(BUILD_FPIC
          "Libraries will be compiled with position independent code" TRUE)

set(PSI4_STAGE_EX ${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/external)

############################## Compiler Features ###############################
find_package(OpenMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++11")

############################## Dependencies ####################################
find_package(PythonInterp REQUIRED)
find_package(PythonLibs ${PYTHON_VERSION_MAJOR} REQUIRED)

set(MKL_COMPILER_BINDINGS "${CMAKE_CXX_COMPILER_ID}")
set(BLAS_LANG "CXX")
set(LAPACK_LANG "CXX")
include(ConfigMath)
#These packages don't believe in configuration files...
foreach(pkg in boost pcmsolver chemps2)
   string(TOUPPER "${pkg}" PKG)
   if(NOT ${PKG}_ROOT)
      set(${PKG}_ROOT ${PSI4_STAGE_EX}/${pkg})
   endif()
endforeach()



add_subdirectory(external)


list(APPEND CMAKE_PREFIX_PATH  ${PSI4_STAGE_EX})

############################## Main Project ####################################

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
              -DLIBINT_OPT_AM=${LIBINT_OPT_AM}
              -DBUILD_FPIC=${BUILD_FPIC}
              -DPYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}
              -DPYTHON_INCLUDE_DIR=${PYTHON_INCLUDE_DIR}
              -DPYTHON_LIBRARY=${PYTHON_LIBRARY}
              -DLAPACK_LIBRARIES=${LAPACK_LIBRARIES}
              -DLAPACK_INCLUDE_DIRS=${LAPACK_INCLUDE_DIRS}
              -DPSI4_ROOT=${CMAKE_CURRENT_SOURCE_DIR}
              -DBOOST_ROOT=${BOOST_ROOT}
              -DPCMSOLVER_ROOT=${PCMSOLVER_ROOT}
              -DCHEMPS2_ROOT=${CHEMPS2_ROOT}
              -DENABLE_LIBERD=${ENABLE_LIBERD}
              -DENABLE_AMBIT=${ENABLE_AMBIT}
              -DENABLE_PCMSOLVER=${ENABLE_PCMSOLVER}
              -DENABLE_GTFOCK=${ENABLE_GTFOCK}
              -DENABLE_GDMA=${ENABLE_GDMA}
              -DENABLE_DKH=${ENABLE_DKH}
              -DENABLE_CHEMPS2=${ENABLE_CHEMPS2}
              -DENABLE_PYTHON_INTERFACE=${ENABLE_PYTHON_INTERFACE}
              -DRESTRICT_KEYWORD=${RESTRICT_KEYWORD}
              -DFC_SYMBOL=${FC_SYMBOL}
              -DDESTDIR=${CMAKE_BINARY_DIR}/stage
   CMAKE_CACHE_ARGS -DCMAKE_PREFIX_PATH:STRING=${CMAKE_PREFIX_PATH}
   BUILD_ALWAYS 1
   INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install 
                   DESTDIR=${CMAKE_BINARY_DIR}/stage
)
add_dependencies(psi4_core boost_external 
                           libint_external
                           libderiv_external 
                           libefp_external
                           pcmsolver_external
                           liberd_external
                           ambit_external
                           dkh_external
                           gdma_external
                           chemps2_external
)
add_subdirectory(doc)
enable_testing()
add_subdirectory(tests)
install(DIRECTORY ${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/
        DESTINATION ${CMAKE_INSTALL_PREFIX}
        USE_SOURCE_PERMISSIONS)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/cmake
        DESTINATION ${CMAKE_INSTALL_PREFIX} USE_SOURCE_PERMISSIONS)
