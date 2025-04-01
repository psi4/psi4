#############################  Superbuild Project  #############################
cmake_minimum_required(VERSION 3.19 FATAL_ERROR)
cmake_policy(SET CMP0094 NEW)
cmake_policy(SET CMP0075 NEW)
if (POLICY CMP0135)
    # Avoid warning about DOWNLOAD_EXTRACT_TIMESTAMP in CMake 3.24:
    cmake_policy(SET CMP0135 NEW)
endif ()
if (POLICY CMP0144)
    cmake_policy(SET CMP0144 NEW)
endif()
if (POLICY CMP0167)
    # L2 not prepared for BoostConfig instead of FindBoost
    cmake_policy(SET CMP0167 OLD)
endif ()
                                              # 3.15: useable FindPython with NumPy component and LOCATION
                                              # 3.8: CXX_STANDARD recognizes C++17
                                              # 3.3: EP_add(DEPENDS interfacelib, CheckFortranCompilerFlag
                                              # 3.2: continue()
                                              # 3.1: CMAKE_CXX_STANDARD, BUILD_ALWAYS

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build.")  # foil Ninja Debug on Windows

project(psi4
        LANGUAGES C CXX)
set(psi4_URL "http://www.psicode.org/")
set(psi4_EMAIL "psi4aiqc+cmake@gmail.com")
set(psi4_LICENSE "GNU Lesser General Public License, version 3 (LGPL-3.0)")
set(psi4_DESCRIPTION "Open-Source Quantum Chemistry")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

message(STATUS "Building using CMake ${CMAKE_VERSION} Generator ${CMAKE_GENERATOR}")

################  Options: Overview and Not Otherwise Mentioned  ###############

#  <<<  CMake build overview  >>>
#
#    >>> ls
#    external/  COPYING  psi4/  tests/  ...
#    >>> cmake -S. -Bobjdir -DCMAKE_INSTALL_PREFIX=/path/to/install-psi4 ...
#    ...
#    -- Generating done
#    -- Build files have been written to: /current/dir/objdir
#    >>> cmake --build objdir --target install -j`getconf _NPROCESSORS_ONLN`

#  <<<  Required build dependencies that Psi4 can't build itself  >>>
#
#    - CMake (e.g., `conda install cmake`)
#    - C++ and C compilers (C++20 compliant)
#    - BLAS/LAPACK (also runtime; e.g., `conda install mkl-devel`)
#    - Python (also runtime; interpreter and headers; e.g., `conda install python`)
#    - NumPy (also runtime; avoidable at buildtime if gau2grid pre-built; e.g., `conda install numpy`)
#    - networkx >=2.4 (runtime only; e.g., `conda install networkx`)
#    - pint >=0.10 (runtime only; e.g., `conda install pint -c conda-forge`)
#    - pydantic >=1.0 (runtime only; e.g., `conda install pydantic -c conda-forge`)
#    - msgpack-python (runtime only; e.g., `conda install msgpack-python`)
#    - Eigen (for Libint2; e.g., `conda install eigen`)
#    - Boost (header-only libs for Libint2; e.g., `conda install boost-cpp`)
#    - SciPy (runtime only; avoidable through keyword setting; e.g., `conda install scipy`)
#    - py-cpuinfo (runtime only; e.g., `conda install py-cpuinfo`)
#    - psutil (runtime only; e.g., `conda install psutil`)

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
#    - Python_EXECUTABLE "Python interpreter to use (e.g., /path/to/bin/python3.7)"
#    - Python_LIBRARY "Python library that goes with the interpreter (e.g., /path/to/lib/python3.7.so)"
#    - Python_INCLUDE_DIR "Path to the python include files (e.g., /path/to/include/python3.7)"
#    - SPHINX_ROOT "Root directory for Sphinx: 'bin/sphinx-build' (or similar) should be in this dir."
#
#    For any ${AddOn} of: ambit, CheMPS2, dkh, ecpint, libefp, Einsums, gau2grid, gauxc, gdma,
#                         IntegratorXX, Libint2, OpenOrbitalOptimizer, optking, PCMSolver, pybind11,
#                         pylibefp, qcelemental, qcengine, simint, Libxc
#    - CMAKE_PREFIX_PATH "Set to list of root directories to look for externally built add-ons and dependencies
#      (e.g., /path/to/install-libint;/path/to/install-gdma where exists /path/to/install-libint/lib/libderiv.a)"
#    - ${AddOn}_DIR "Set to directory containing ${AddOn}Config.cmake file to facilitate detection of external build"
#    - CMAKE_DISABLE_FIND_PACKAGE_${AddON} "Set to ON to force internal build"
#    - CMAKE_INSIST_FIND_PACKAGE_${AddON} "Set to ON to force external detect"

#  <<<  Detecting BLAS/LAPACK  >>>
#
#    - ENV(MATH_ROOT) "Root directory where BLAS/LAPACK libraries should be detected (e.g., ${MATH_ROOT}/lib/libblas.so)"
#    - BLAS_TYPE "Target BLAS distribution for math detection
#                 (default: search order MKL>FLEXIBLAS>OPENBLAS>ESSL>ATLAS>ACML>SYSTEM_NATIVE on Linux; MKL>Accelerate>... on Mac)"
#    - LAPACK_TYPE "Target LAPACK distribution for math detection
#                 (default: search order MKL>FLEXIBLAS>OPENBLAS>ESSL>ATLAS>ACML>SYSTEM_NATIVE on Linux; MKL>Accelerate>... on Mac)"
#    - LAPACK_LIBRARIES "Location of BLAS/LAPACK libraries as ";"-separated list of full paths, bypassing math detection"
#    - LAPACK_INCLUDE_DIRS "Location of BLAS/LAPACK headers (only needed for MKL), bypassing math detection"
#    - OpenMP_LIBRARY_DIRS "Location of OpenMP libraries (iomp5/gomp/omp) as ";"-separated list, hinting OpenMP detection"

#  <<<  Install  >>>
#
#    - CMAKE_INSTALL_PREFIX "Location to which Psi4 and internally built add-ons are installed (default: /usr/local/psi4)"
#    - CMAKE_INSTALL_BINDIR "Location within CMAKE_INSTALL_PREFIX to which executables are installed (default: bin)"
#    - CMAKE_INSTALL_LIBDIR "Location within CMAKE_INSTALL_PREFIX to which libraries are installed (default: lib)"
#    - CMAKE_INSTALL_DATADIR "Location within CMAKE_INSTALL_PREFIX to which resources are installed (default: share)"
#    - CMAKE_INSTALL_INCLUDEDIR "Location within CMAKE_INSTALL_PREFIX to which headers are installed (default: include)"
#    - PYMOD_INSTALL_LIBDIR "Location within CMAKE_INSTALL_LIBDIR to which python modules are installed (default: /)
#                            Must start with: / . Used to imitate python install: /python3.6/site-packages ."

############################  Options: Build What?  ############################
option(ENABLE_ambit "Enables the ambit tensor library" OFF)
option(ENABLE_bse "Enables the basis_set_exchange for requesting basis sets (requires jsonschema, regex; can also be added at runtime)" OFF)
option(ENABLE_CheMPS2 "Enables CheMPS2 for DMRG (requires HDF5)" OFF)
option(ENABLE_cppe "Enables CPPE for Polarizable Embedding" OFF)
option(ENABLE_adcc "Enables adcc for algebraic-diagrammatic construction methods (requires libtensorlight; can also be added at runtime)" OFF)
option(ENABLE_ddx "Enables ddx for domain-decomposition contiuum solvation models (requires Fortran; can also be added at runtime)" OFF)
option(ENABLE_dkh "Enables DKH integrals (requires Fortran)" OFF)
option(ENABLE_ecpint "Enables libecpint for effective core potentials (ECP)" OFF)
option(ENABLE_libefp "Enables LIBEFP and PylibEFP for fragments" OFF)
option(ENABLE_Einsums "Enables the Einsums tensor library" OFF)
option(ENABLE_simint "Enables use of SIMINT two-electron integral library" OFF)
option(ENABLE_gauxc "Enables use of GauXC density functional library" OFF)
option(ENABLE_gdma "Enables Stone's GDMA multipole code (requires Fortran; can also be added at runtime)" OFF)
option(ENABLE_OpenOrbitalOptimizer "Enables the OpenOrbitalOptimizer library (requires Armadillo)" OFF)
option(ENABLE_PCMSolver "Enables PCMSolver library (requires Fortran)" OFF)
option(ENABLE_snsmp2 "Enables SNSMP2 plugin (can also be added at runtime)" OFF)
option(ENABLE_v2rdm_casscf "Enables V2RDM_CASSCF plugin (requires Fortran; can also be added at runtime)" OFF)
option(ENABLE_psi4fockci "Enables Psi4FockCI plugin (can also be added at runtime)" OFF)
option(ENABLE_cct3 "Enables cct3 plugin (requires Fortran; can also be added at runtime)" OFF)
option(ENABLE_gpu_dfcc "Enables GPU_DFCC plugin for gpu-accelerated df-cc (requires CUDA; can also be added at runtime)" OFF)
option(ENABLE_mdi "Enables MolSSI driver interface" OFF)
option(ENABLE_BrianQC "Enables the BrianQC GPU module (requires CUDA; requires separate installation and licensing of the BrianQC module)" OFF)
# These options are relevant to pasture, expert only
option(ENABLE_ccsort "Enables ccsort plugin installed from psi4pasture" OFF)
option(ENABLE_transqt2 "Enables transqt2 plugin installed from psi4pasture" OFF)
option(psi4_SKIP_ENABLE_Fortran "Suppress enabling Fortran language even when Fortran addons enabled. Use when only linking pre-built, not compiling addon code." OFF)

# Append modules added to pasture as needed
if(ENABLE_ccsort OR ENABLE_transqt2)
  set(ENABLE_pasture ON)
  message(STATUS "Enabling pasture plugins")
endif()

if(ENABLE_gdma OR ENABLE_dkh OR ENABLE_PCMSolver)
  if(psi4_SKIP_ENABLE_Fortran)
    message(STATUS "NOT Enabling Fortran")
  else()
    enable_language(Fortran)
    set(Fortran_ENABLED ON)  # communicate required languages with psi4-core
    message(STATUS "Enabling Fortran")
  endif()
endif()

############################  Options: Build How?  #############################
include(psi4OptionsTools)
option_with_print(BUILD_SHARED_LIBS "Build internally built Psi4 add-on libraries as shared, not static" OFF)
option_with_default(CMAKE_BUILD_TYPE "Build type (Release or Debug)" Release)
include(xhost)  # defines: option(ENABLE_XHOST "Enable processor-specific optimization" ON)
# below are uncommon to adjust
option_with_print(ENABLE_OPENMP "Enables OpenMP parallelization" ON)
option_with_print(ENABLE_AUTO_BLAS "Enables CMake to auto-detect BLAS" ON)
option_with_print(ENABLE_AUTO_LAPACK "Enables CMake to auto-detect LAPACK" ON)
option_with_print(ENABLE_PLUGIN_TESTING "Test the plugin templates build and run" OFF)
option_with_print(ENABLE_CYTHONIZE "Compile each python file rather than plaintext (requires cython) !experimental!" OFF)
# sadly, it seems FORCE_PEDANTIC runs into issues with GauXC when option_with_flag is used
# option_with_print + setting -pedantic-errors in psi4/CMakeLists.txt seems to fix this
option_with_print(FORCE_PEDANTIC "Forces code to comply with the language standard." OFF)
option_with_flags(ENABLE_CODE_COVERAGE "Enables details on code coverage" OFF
                  "-ftest-coverage")
option_with_flags(ENABLE_BOUNDS_CHECK "Enables bounds check in Fortran" OFF
                  "-ftrapuv -check all -fpstkchk" "-fcheck=all" "-fbounds-check -fcheck-array-temporaries")
option_with_flags(ENABLE_ASAN "Enables address sanitizer (requires similarly compiled Python and Numpy)" OFF
                  "-fsanitize=address -fno-omit-frame-pointer")
option_with_flags(ENABLE_TSAN "Enables thread sanitizer (requires similarly compiled Python and Numpy)" OFF
                  "-fsanitize=thread -fPIE -pie -fno-omit-frame-pointer")
option_with_flags(ENABLE_UBSAN "Enables undefined behavior sanitizer (requires similarly compiled Python and Numpy)" OFF
                  "-fsanitize=undefined -fno-omit-frame-pointer")
option_with_flags(ENABLE_MSAN "Enables memory sanitizer (requires similarly compiled Python and Numpy)" OFF
                  "-fsanitize=memory -fPIE -pie -fno-omit-frame-pointer")
option_with_default(MAX_AM_ERI "Maximum angular momentum for integrals. The definition changed with Libint2 but still single integer. Default fine for production." 5)
option_with_default(FC_SYMBOL "The type of Fortran name mangling" 2)
option_with_default(BUILD_FPIC "Compile static libraries with position independent code" ON)
option_with_default(CMAKE_INSTALL_LIBDIR "Directory to which libraries installed" lib)
option_with_default(PYMOD_INSTALL_LIBDIR "Location within CMAKE_INSTALL_LIBDIR to which python modules are installed" /)
option_with_default(psi4_INSTALL_CMAKEDIR "Directory to which psi4 CMake config files installed." share/cmake/psi4)
option_with_default(TargetLAPACK_INSTALL_CMAKEDIR "Directory to which psi4 CMake config files installed." share/cmake/TargetLAPACK)
option_with_default(TargetHDF5_INSTALL_CMAKEDIR "Directory to which psi4 CMake config files installed." share/cmake/TargetHDF5)
option_with_default(ENABLE_GENERIC "Enables mostly static linking of language libraries for shared library" OFF)
option_with_default(CMAKE_INSTALL_MESSAGE "Specify verbosity of installation messages" LAZY)
option_with_default(psi4_CXX_STANDARD "Specify C++ standard for core Psi4" 20)
option_with_default(psi4_CXX_STANDARD_Libint2 "Specify C++ standard for Psi4 files including Libint2 headers. Older (e.g., 2021.4) Intel compilers may need 14" 20)
option_with_default(SIMINT_VECTOR "Vectorization type to use for simint (scalar sse avx avxfma micavx512)" avx)
option_with_default(SPHINX_THEME "Theme for Sphinx documentation and extensions" sphinx_psi_theme)
option_with_default(SPHINXMAN_STRICT "Turn warnings into errors for docs target sphinxman" OFF)
option_with_default(BUILD_Libint2_GENERATOR "If prebuilt not found, don't build Libint2 from source tarball but instead build and run generator (several hours for production library). Edit angular momentum and integrals classes directly in `external/upstream/libint2/CMakeLists.txt` (requires Ninja)." OFF)
option_with_default(psi4_SHGAUSS_ORDERING "Ordering for shells of solid harmonic Gaussians: standard -- standard ordering (-l, -l+1 ... l) -or- gaussian -- the Gaussian ordering (0, 1, -1, 2, -2, ... l, -l). Changing this is experimental, expert only, and not a future plan. Internal calculations should be correct and fchk files should be correct, but other exports like Molden and QCSchema.Wavefunction won't be." gaussian)
option_with_default(gauxc_ENABLE_GPU "Enable GPU support for GauXC" OFF)
option_with_default(gauxc_ENABLE_MAGMA "Enable MAGMA support for GauXC. Required GPU-enabled GauXC" OFF)

include(custom_int_orderings)


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

set(Python_ADDITIONAL_VERSIONS 3.12 3.11 3.10 3.9 3.8)  # adjust with CMake minimum FindPythonInterp
find_package(Python 3.8 COMPONENTS Interpreter Development NumPy REQUIRED)
message(STATUS "${Cyan}Found Python ${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}${ColourReset}: ${Python_EXECUTABLE} (found version ${Python_VERSION}; NumPy version: ${Python_NumPy_VERSION})")
set(STAGED_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/stage)
message(STATUS "Psi4 staging: ${STAGED_INSTALL_PREFIX}")
message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
add_subdirectory(external/common)
add_subdirectory(external/upstream)
# external projects manage their own OpenMP and c++YY flags, so only add to CXX_FLAGS for psi4-core
include(custom_cxxstandard)
include(custom_static_library)

################################  Main Project  ################################
include(ExternalProject)
ExternalProject_Add(psi4-core
   DEPENDS lapack_external
           hdf5_external
           ambit_external
           bse_external
           chemps2_external
           cppe_external
           dkh_external
           ecpint_external
           libefp_external
           einsums_external
           gau2grid_external
           gdma_external
           libint2_external
           pcmsolver_external
           ddx_external
           pybind11_external
           pylibefp_external
           qcelemental_external
           qcengine_external
           openorbitaloptimizer_external
           optking_external
           simint_external
           libxc_external
           mdi_external
           brianqc_external
           integratorxx_external
           gauxc_external
 
   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/psi4
   CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
              -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
              -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
              -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
              -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
              -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
              -DCMAKE_CXX_STANDARD=${psi4_CXX_STANDARD}
              -Dpsi4_CXX_STANDARD_Libint2=${psi4_CXX_STANDARD_Libint2}
              -DCMAKE_CXX_STANDARD_REQUIRED=ON
              -DCMAKE_CXX_EXTENSIONS=OFF
              -DFORCE_PEDANTIC=${FORCE_PEDANTIC}
              -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
              -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}
              -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
              -DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
              -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
              -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
              -DCMAKE_INSTALL_DATADIR=${CMAKE_INSTALL_DATADIR}
              -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
              -Dpsi4_INSTALL_CMAKEDIR=${psi4_INSTALL_CMAKEDIR}
              -DCMAKE_INSTALL_MESSAGE=${CMAKE_INSTALL_MESSAGE}
              -DPYMOD_INSTALL_LIBDIR=${PYMOD_INSTALL_LIBDIR}
              -DMAX_AM_ERI=${MAX_AM_ERI}
              -DPython_EXECUTABLE=${Python_EXECUTABLE}
              -DPython_INCLUDE_DIR=${Python_INCLUDE_DIRS}
              -DPython_LIBRARY=${Python_LIBRARIES}
              -DPython_NumPy_INCLUDE_DIR=${Python_NumPy_INCLUDE_DIR}
              -DPSI4_ROOT=${CMAKE_CURRENT_SOURCE_DIR}
              -DENABLE_ambit=${ENABLE_ambit}
              -DENABLE_CheMPS2=${ENABLE_CheMPS2}
              -DENABLE_cppe=${ENABLE_cppe}
              -DENABLE_dkh=${ENABLE_dkh}
              -DENABLE_IntegratorXX=${ENABLE_IntegratorXX}
              -DENABLE_gauxc=${ENABLE_gauxc}
              -DENABLE_ecpint=${ENABLE_ecpint}
              -DENABLE_libefp=${ENABLE_libefp}
              -DENABLE_Einsums=${ENABLE_Einsums}
              -DENABLE_simint=${ENABLE_simint}
              -DENABLE_gdma=${ENABLE_gdma}
              -DENABLE_PCMSolver=${ENABLE_PCMSolver}
              -DENABLE_ddx=${ENABLE_ddx}
              -DENABLE_mdi=${ENABLE_mdi}
              -DENABLE_BrianQC=${ENABLE_BrianQC}
              -DENABLE_OpenOrbitalOptimizer=${ENABLE_OpenOrbitalOptimizer}
              -DENABLE_OPENMP=${ENABLE_OPENMP}
              -DTargetLAPACK_DIR=${TargetLAPACK_DIR}
              -DTargetHDF5_DIR=${TargetHDF5_DIR}
              -DEigen3_DIR=${Eigen3_DIR}
              -Dambit_DIR=${ambit_DIR}
              -DCheMPS2_DIR=${CheMPS2_DIR}
              -Ddkh_DIR=${dkh_DIR}
              -Decpint_DIR=${ecpint_DIR}
              -Dlibefp_DIR=${libefp_DIR}
              -DEinsums_DIR=${Einsums_DIR}
              -Dgau2grid_DIR=${gau2grid_DIR}
              -Dgdma_DIR=${gdma_DIR}
              -DLibint2_DIR=${Libint2_DIR}
              -DOpenOrbitalOptimizer_DIR=${OpenOrbitalOptimizer_DIR}
              -DPCMSolver_DIR=${PCMSolver_DIR}
              -Dmdi_DIR=${mdi_DIR}
              -Dpybind11_DIR=${pybind11_DIR}
              -DIntegratorXX_DIR=${IntegratorXX_DIR}
              -Dgauxc_DIR=${gauxc_DIR}
              -Dpylibefp_DIR=${pylibefp_DIR}
              -Dsimint_DIR=${simint_DIR}
              -DLibxc_DIR=${Libxc_DIR}
              -DBrianQC_DIR=${BrianQC_DIR}
              -DFortran_ENABLED=${Fortran_ENABLED}
              -DLIBC_INTERJECT=${LIBC_INTERJECT}
              -DRESTRICT_KEYWORD=${RESTRICT_KEYWORD}
              -DFC_SYMBOL=${FC_SYMBOL}
              -DOpenMP_LIBRARY_DIRS=${OpenMP_LIBRARY_DIRS}
              -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
              -DENABLE_CYTHONIZE=${ENABLE_CYTHONIZE}
              -Dpsi4_SHGSHELL_ORDERING=${psi4_SHGSHELL_ORDERING}
   CMAKE_CACHE_ARGS -DCMAKE_PREFIX_PATH:PATH=${CMAKE_PREFIX_PATH}
              -DCMAKE_C_COMPILER_ARG1:STRING=${CMAKE_C_COMPILER_ARG1}
              -DCMAKE_CXX_COMPILER_ARG1:STRING=${CMAKE_CXX_COMPILER_ARG1}
              -DCMAKE_Fortran_COMPILER_ARG1:STRING=${CMAKE_Fortran_COMPILER_ARG1}
               # CMAKE_<lang>_COMPILER_ARG1 vars pass args thru when CMAKE_<lang>_COMPILER set as "compiler --vital-args"
   USES_TERMINAL_BUILD 1
   BUILD_ALWAYS 1)

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
