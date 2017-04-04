#############################  Superbuild Project  #############################
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)  # 3.3: EP_add(DEPENDS interfacelib, CheckFortranCompilerFlag
                                                 # 3.2: continue()
                                                 # 3.1: CMAKE_CXX_STANDARD, BUILD_ALWAYS
project(psi4
        LANGUAGES C CXX)
set(psi4_URL "http://www.psicode.org/")
set(psi4_EMAIL "psi4aiqc+cmake@gmail.com")
set(psi4_LICENSE "GNU General Public License, version 2")
set(psi4_DESCRIPTION "Open-Source Quantum Chemistry")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

################  Options: Overview and Not Otherwise Mentioned  ###############

#  <<<  CMake build overview  >>>
#
#    >>> ls
#    external/  LICENSE  psi4/  tests/  ...
#    >>> cmake -H. -Bobjdir -DCMAKE_INSTALL_PREFIX=/path/to/install-psi4 ...
#    ...
#    -- Generating done
#    -- Build files have been written to: /current/dir/objdir
#    >>> cd objdir && make -j`getconf _NPROCESSORS_ONLN`
#    >>> make install

#  These three "###  Options  ###" sections contain useful CMake variables for build configuration.

#  <<<  Compilers and flags  >>>
#
#    - CMAKE_C_COMPILER "C compiler"
#    - CMAKE_C_FLAGS "Additional C flags"
#    - CMAKE_CXX_COMPILER "C++ compiler"
#    - CMAKE_CXX_FLAGS "Additional C++ flags"
#    - CMAKE_Fortran_COMPILER "Fortran compiler (required for some add-ons)"
#    - CMAKE_Fortran_FLAGS "Additional Fortran flags"

#  <<<  Detecting dependencies and add-ons  >>>
#
#    - PYTHON_EXECUTABLE "Python interpreter to use (e.g., /path/to/bin/python2.7)"
#    - PYTHON_LIBRARY "Python library that goes with the interpreter (e.g., /path/to/lib/python2.7.so)"
#    - PYTHON_INCLUDE_DIR "Path to the python include files (e.g., /path/to/include/python2.7)"
#    - SPHINX_ROOT "Root directory for Sphinx: 'bin/sphinx-build' (or similar) should be in this dir."
#
#    For any ${AddOn} of: CheMPS2, dkh, libefp, erd, gdma, libint, PCMSolver, pybind11, simint
#    - CMAKE_PREFIX_PATH "Set to list of root directories to look for externally built add-ons and dependencies
#      (e.g., /path/to/install-libint;/path/to/install-gdma where exists /path/to/install-libint/lib/libderiv.a)"
#    - ${AddOn}_DIR "Set to directory containing ${AddOn}Config.cmake file to facilitate detection of external build"
#    - CMAKE_DISABLE_FIND_PACKAGE_${AddON} "Set to OFF to force internal build"

#  <<<  Detecting BLAS/LAPACK  >>>
#
#    - ENV(MATH_ROOT) "Root directory where BLAS/LAPACK libraries should be detected (e.g., ${MATH_ROOT}/lib/libblas.so)"
#    - BLAS_TYPE "Target BLAS distribution for math detection
#                 (default: search order MKL>OPENBLAS>ESSL>ATLAS>ACML>SYSTEM_NATIVE on Linux; MKL>Accelerate>... on Mac)"
#    - LAPACK_TYPE "Target LAPACK distribution for math detection
#                 (default: search order MKL>OPENBLAS>ESSL>ATLAS>ACML>SYSTEM_NATIVE on Linux; MKL>Accelerate>... on Mac)"
#    - LAPACK_LIBRARIES "Location of BLAS/LAPACK libraries as ";"-separated list of full paths, bypassing math detection"
#    - LAPACK_INCLUDE_DIRS "Location of BLAS/LAPACK headers (only needed for MKL), bypassing math detection"

#  <<<  Install  >>>
#
#    - CMAKE_INSTALL_PREFIX "Location to which Psi4 and internally built add-ons are installed (default: /usr/local/psi4)"
#    - CMAKE_INSTALL_BINDIR "Location within CMAKE_INSTALL_PREFIX to which executables are installed (default: bin)"
#    - CMAKE_INSTALL_LIBDIR "Location within CMAKE_INSTALL_PREFIX to which libraries are installed (default: lib)"
#    - CMAKE_INSTALL_DATADIR "Location within CMAKE_INSTALL_PREFIX to which resources are installed (default: share)"
#    - CMAKE_INSTALL_INCLUDEDIR "Location within CMAKE_INSTALL_PREFIX to which headers are installed (default: include)"
#    - PYMOD_INSTALL_LIBDIR "Location within CMAKE_INSTALL_LIBDIR to which python modules are installed (default: /)
#                            Must start with /"

############################  Options: Build What?  ############################
option(ENABLE_AMBIT "Enables the Ambit tensor library" OFF)
option(ENABLE_CheMPS2 "Enables CheMPS2 for DMRG (requires HDF5)" OFF)
option(ENABLE_dkh "Enables DKH integrals (requires Fortran)" OFF)
option(ENABLE_libefp "Enables LIBEFP for fragments" OFF)
option(ENABLE_erd "Enables use of ERD instead of libint (requires Fortran)" OFF)
option(ENABLE_simint "Enables use of SIMINT two-electron integral library" OFF)
option(ENABLE_gdma "Enables Stone's GDMA multipole code (requires Fortran)" OFF)
option(ENABLE_PCMSolver "Enables PCMSolver library (requires Fortran)" OFF)
# These options are relevant to pasture
option(ENABLE_ccsort "Enables ccsort plugin installed from psi4pasture" OFF)
option(ENABLE_transqt2 "Enables transqt2 plugin installed from psi4pasture" OFF)

# Append modules added to pasture as needed
if(ENABLE_ccsort OR ENABLE_transqt2)
  set(ENABLE_pasture ON)
  message(STATUS "Enabling pasture plugins")
endif()

if(ENABLE_gdma OR ENABLE_dkh OR ENABLE_erd OR ENABLE_PCMSolver)
    enable_language(Fortran)
    set(Fortran_ENABLED ON)  # communicate required languages with psi4-core
    message(STATUS "Enabling Fortran")
endif()

if(ENABLE_erd)
    message(WARNING "ERD will build, link, and run in Psi4 just fine. However, it has not been hooked into Psi4 in all roles, notably gradients, LRC DFT energies, and ESP. So upon activating through ``set integral_package erd``, known failures will be caught and halted, but perhaps other types not tested and identified will give *wrong* answers. Consider this your warning.")
endif()

############################  Options: Build How?  #############################
include(psi4OptionsTools)
option_with_print(BUILD_SHARED_LIBS "Build internally built Psi4 add-on libraries as shared, not static" OFF)
option_with_print(ENABLE_OPENMP "Enables OpenMP parallelization" ON)
option_with_print(ENABLE_AUTO_BLAS "Enables CMake to auto-detect BLAS" ON)
option_with_print(ENABLE_AUTO_LAPACK "Enables CMake to auto-detect LAPACK" ON)
option_with_print(ENABLE_PLUGIN_TESTING "Test the plugin templates build and run" OFF)
option_with_flags(ENABLE_XHOST "Enables processor-specific optimization" ON
                  "-xHost" "-march=native")
option_with_flags(ENABLE_CODE_COVERAGE "Enables details on code coverage" OFF
                  "-ftest-coverage")
option_with_flags(ENABLE_BOUNDS_CHECK "Enables bounds check in Fortran" OFF
                  "-ftrapuv -check all -fpstkchk" "-fcheck=all" "-fbounds-check -fcheck-array-temporaries")
option_with_flags(ENABLE_ASAN "Enables address sanitizer" OFF
                  "-fsanitize=address" "-fno-omit-frame-pointer")
option_with_flags(ENABLE_TSAN "Enables thread sanitizer" OFF
                  "-fsanitize=thread" "-fno-omit-frame-pointer -pie")
option_with_flags(ENABLE_UBSAN "Enables undefined behavior sanitizer" OFF
                  "-fsanitize=undefined" "-fno-omit-frame-pointer")
option_with_default(MAX_AM_ERI "Maximum angular momentum for integrals" 5)
option_with_default(CMAKE_BUILD_TYPE "Build type (Release or Debug)" Release)
option_with_default(FC_SYMBOL "The type of Fortran name mangling" 2)
option_with_default(BUILD_FPIC "Compile static libraries with position independent code" ON)
option_with_default(CMAKE_INSTALL_LIBDIR "Directory to which libraries installed" lib)
option_with_default(PYMOD_INSTALL_LIBDIR "Location within CMAKE_INSTALL_LIBDIR to which python modules are installed" /)
option_with_default(ENABLE_GENERIC "Enables mostly static linking of system and math libraries for shared library" OFF)
option_with_default(CMAKE_INSTALL_MESSAGE "Specify verbosity of installation messages" LAZY)
option_with_default(psi4_CXX_STANDARD "Specify C++ standard for core Psi4" 11)
option_with_default(SIMINT_VECTOR "Vectorization type to use for simint (scalar sse avx avxfma micavx512)" avx)

########################  Process & Validate Options  ##########################
include(GNUInstallDirs)
include(autocmake_safeguards)
include(custom_color_messages)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "/usr/local/psi4" CACHE PATH "Install path" FORCE)
endif()
message(STATUS "Psi4 install: ${CMAKE_INSTALL_PREFIX}")

# Python use in psi4 (not incl. external/)
# * interpreter: run the versioner in psi4-core. run tests and build docs in psi4
# * headers: Python.h for pybind11 to build against in psi4-core
# * library: not at all

set(Python_ADDITIONAL_VERSIONS 3.7 3.6 3.5)  # adjust with CMake minimum FindPythonInterp
find_package(PythonLibsNew REQUIRED)
message(STATUS "${Cyan}Found Python ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}${ColourReset}: ${PYTHON_EXECUTABLE} (found version ${PYTHON_VERSION_STRING})")

set(STAGED_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/stage${CMAKE_INSTALL_PREFIX})
add_subdirectory(external/common)
add_subdirectory(external/upstream)

# external projects manage their own OpenMP and c++YY flags, so only add to CXX_FLAGS for psi4-core
include(autocmake_omp)
include(custom_cxxstandard)
include(custom_static_library)

################################  Main Project  ################################
include(ExternalProject)
ExternalProject_Add(psi4-core
   DEPENDS lapack_external
           hdf5_external
           ambit_external
           chemps2_external
           dkh_external
           libefp_external
           erd_external
           gdma_external
           libint_external
           pcmsolver_external
           pybind11_external
           simint_external
   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/psi4
   CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
              -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
              -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
              -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
              -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
              -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
              -DCMAKE_CXX_STANDARD=${psi4_CXX_STANDARD}
              -DCMAKE_CXX_STANDARD_REQUIRED=ON
              -DCMAKE_CXX_EXTENSIONS=OFF
              -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
              -DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
              -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
              -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
              -DCMAKE_INSTALL_DATADIR=${CMAKE_INSTALL_DATADIR}
              -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
              -DCMAKE_INSTALL_MESSAGE=${CMAKE_INSTALL_MESSAGE}
              -DPYMOD_INSTALL_LIBDIR=${PYMOD_INSTALL_LIBDIR}
              -DMAX_AM_ERI=${MAX_AM_ERI}
              -DBUILD_FPIC=${BUILD_FPIC}
              -DPYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}
              -DPYTHON_INCLUDE_DIR=${PYTHON_INCLUDE_DIR}
              -DPYTHON_LIBRARY=${PYTHON_LIBRARY}
              -DPSI4_ROOT=${CMAKE_CURRENT_SOURCE_DIR}
              -DENABLE_AMBIT=${ENABLE_AMBIT}
              -DENABLE_CheMPS2=${ENABLE_CheMPS2}
              -DENABLE_dkh=${ENABLE_dkh}
              -DENABLE_libefp=${ENABLE_libefp}
              -DENABLE_erd=${ENABLE_erd}
              -DENABLE_simint=${ENABLE_simint}
              -DENABLE_gdma=${ENABLE_gdma}
              -DENABLE_PCMSolver=${ENABLE_PCMSolver}
              -DTargetLAPACK_DIR=${TargetLAPACK_DIR}
              -DTargetHDF5_DIR=${TargetHDF5_DIR}
              -DCheMPS2_DIR=${CheMPS2_DIR}
              -Ddkh_DIR=${dkh_DIR}
              -Dlibefp_DIR=${libefp_DIR}
              -Derd_DIR=${erd_DIR}
              -Dgdma_DIR=${gdma_DIR}
              -Dlibint_DIR=${libint_DIR}
              -DPCMSolver_DIR=${PCMSolver_DIR}
              -Dpybind11_DIR=${pybind11_DIR}
              -Dsimint_DIR=${simint_DIR}
              -DFortran_ENABLED=${Fortran_ENABLED}
              -DLIBC_INTERJECT=${LIBC_INTERJECT}
              -DRESTRICT_KEYWORD=${RESTRICT_KEYWORD}
              -DFC_SYMBOL=${FC_SYMBOL}
              -DDESTDIR=${CMAKE_BINARY_DIR}/stage
   CMAKE_CACHE_ARGS -DCMAKE_PREFIX_PATH:PATH=${CMAKE_PREFIX_PATH}
   BUILD_ALWAYS 1
   INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
                   DESTDIR=${CMAKE_BINARY_DIR}/stage)

add_subdirectory(external/downstream)
add_subdirectory(doc)
include(ConfigTesting)

# <<<  Install  >>>

install(DIRECTORY ${STAGED_INSTALL_PREFIX}/
        DESTINATION ${CMAKE_INSTALL_PREFIX}
        USE_SOURCE_PERMISSIONS
        PATTERN "*.pyc" EXCLUDE)

install(DIRECTORY samples
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/psi4
        USE_SOURCE_PERMISSIONS
        PATTERN "example_psi4rc_file" EXCLUDE)
