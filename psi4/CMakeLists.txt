cmake_minimum_required(VERSION 3.19 FATAL_ERROR)
cmake_policy(SET CMP0094 NEW)
project(psi4-core
        LANGUAGES C CXX)
# no Fortran in psi4-core proper, but language needs to be declared
#   for CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES to be populated so
#   static Fortran add-ons can be linked
if(${Fortran_ENABLED})
    enable_language(Fortran)
endif()

if (NOT MSVC)
  set(CMAKE_C_FLAGS_DEBUG "-g -O0")
  set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
endif()

list(APPEND CMAKE_MODULE_PATH ${PSI4_ROOT}/cmake)

include(psi4OptionsTools)
include(TestRestrict)
include(TestBigEndian)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(custom_color_messages)
test_restrict(restrict)
test_big_endian(endianness)
if(${endianness})
    set(_endian_ext ".big")
else()
    set(_endian_ext ".little")
endif()
set(CMAKECONFIG_INSTALL_DIR "${psi4_INSTALL_CMAKEDIR}")

# <<<  Marshal Dependencies & Add-ons  >>>

set(_addons)
set(TargetOpenMP_FIND_COMPONENTS "CXX")
find_package(TargetLAPACK CONFIG REQUIRED)
get_property(_ill TARGET tgt::lapk PROPERTY INTERFACE_LINK_LIBRARIES)
list(GET _ill 0 _ill0)
get_property(_cd TARGET tgt::lapack PROPERTY INTERFACE_COMPILE_DEFINITIONS)
if(${_cd} MATCHES "USING_LAPACK_MKL")
    set(_isMKL " MKL")
endif()
get_property(_illb TARGET tgt::blas PROPERTY INTERFACE_LINK_LIBRARIES)
list(APPEND _addons ${_ill} ${_illb})
message(STATUS "${Cyan}Using LAPACK${_isMKL}${ColourReset}: ${_ill0};...")

if(ENABLE_ambit OR ENABLE_CheMPS2 OR ENABLE_Einsums)
    find_package(TargetHDF5 CONFIG REQUIRED)
    get_property(_ill TARGET tgt::hdf5 PROPERTY INTERFACE_LINK_LIBRARIES)
    list(GET _ill 0 _ill0)
    list(APPEND _addons ${_ill})
    message(STATUS "${Cyan}Using HDF5${ColourReset}: ${_ill0};... (found version ${${PN}_VERSION})")
else()
    message(STATUS "Disabled HDF5 ${TargetHDF5_DIR}")
endif()

#  <<  Pybind11 & Python  >>
find_package(Python 3.8 COMPONENTS Interpreter Development NumPy REQUIRED)
find_package(pybind11 2.6.2 CONFIG REQUIRED)
message(STATUS "${Cyan}Using pybind11${ColourReset}: ${pybind11_INCLUDE_DIR} (version ${pybind11_VERSION} for Py${Python_VERSION} and ${CMAKE_CXX_STANDARD})")
message(STATUS "${Cyan}Using Python ${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}${ColourReset}: ${Python_EXECUTABLE}")

find_file(LIBDL_FOUND dlfcn.h)
if(LIBDL_FOUND)
    message(STATUS "Using DL: ${LIBDL_FOUND} ${CMAKE_DL_LIBS}")
    # use ${CMAKE_DL_LIBS} if found to be explicitly needed
else()
    message(STATUS "Disabled DL ${LIBDL_FOUND}")
endif()

if(${ENABLE_ambit})
    find_package(ambit 0.6 CONFIG REQUIRED)
    get_property(_loc TARGET ambit::ambit PROPERTY LOCATION)
    list(APPEND _addons ${_loc})
    message(STATUS "${Cyan}Using ambit${ColourReset}: ${_loc} (version ${ambit_VERSION})")
else()
    message(STATUS "Disabled ambit ${ambit_DIR}")
endif ()

if(${ENABLE_CheMPS2})
    find_package(CheMPS2 1.8.7 CONFIG REQUIRED)
    get_property(_loc TARGET CheMPS2::chemps2 PROPERTY LOCATION)
    list(APPEND _addons ${_loc})
    message(STATUS "${Cyan}Using CheMPS2${ColourReset}: ${_loc} (version ${CheMPS2_VERSION})")
else()
    message(STATUS "Disabled CheMPS2 ${CheMPS2_DIR}")
endif ()

if(${ENABLE_dkh})
    find_package(dkh 1.2 CONFIG REQUIRED)
    get_property(_loc TARGET dkh::dkh PROPERTY LOCATION)
    list(APPEND _addons ${_loc})
    message(STATUS "${Cyan}Using dkh${ColourReset}: ${_loc} (version ${dkh_VERSION})")
else()
    message(STATUS "Disabled dkh ${dkh_DIR}")
endif()

if(${ENABLE_IntegratorXX})
    find_package(IntegratorXX CONFIG REQUIRED)
    get_property(_loc TARGET IntegratorXX::IntegratorXX PROPERTY LOCATION)
    list(APPEND _addons ${_loc})
    message(STATUS "${Cyan}Using IntegratorXX${ColourReset}: ${_loc} (version ${IntegratorXX_VERSION})")
else()
    message(STATUS "Disabled IntegratorXX")
endif()

if(${ENABLE_ecpint})
    find_package(ecpint 1.0.7 CONFIG REQUIRED)
    get_property(_loc TARGET ECPINT::ecpint PROPERTY LOCATION)
    list(APPEND _addons ${_loc})
    message(STATUS "${Cyan}Using ecpint${ColourReset}: ${_loc} (version ${ecpint_VERSION})")
else()
    message(STATUS "Disabled ecpint ${ecpint_DIR}")
endif()

if(${ENABLE_libefp})
    # neither libefp nor pylibefp actually needed c-side; runtime detection would
    #   do just as well. included for completeness and so pylibefp_PYMOD can be
    #   hard-loaded into PYTHONPATH. also incl for finding efp frag files below.
    find_package(libefp 1.5.0 CONFIG REQUIRED COMPONENTS shallow)
    get_property(_loc TARGET libefp::efp PROPERTY LOCATION)
    list(APPEND _addons ${_loc})
    message(STATUS "${Cyan}Using libefp${ColourReset}: ${_loc} (version ${libefp_VERSION})")

    #find_package(pylibefp 0.6.1 CONFIG REQUIRED)
    #get_property(_loc TARGET pylibefp::core PROPERTY LOCATION)
    include(FindPythonModule)
    find_python_module(pylibefp QUIET REQUIRED)
    #list(APPEND _addons ${_loc})
    message(STATUS "${Cyan}Using pylibefp${ColourReset}: ${_loc} (version ${pylibefp_VERSION})")
else()
    message(STATUS "Disabled libefp ${libefp_DIR} ${pylibefp_DIR}")
endif()

if(${ENABLE_Einsums})
    find_package(Einsums 0.6 CONFIG REQUIRED)
    get_property(_loc TARGET Einsums::einsums PROPERTY LOCATION)
    list(APPEND _addons ${_loc})
    message(STATUS "${Cyan}Using Einsums${ColourReset}: ${_loc} (version ${Einsums_VERSION})")
else()
    message(STATUS "Disabled Einsums ${Einsums_DIR}")
endif ()

find_package(gau2grid 2.0 CONFIG REQUIRED)
get_property(_loc TARGET gau2grid::gg PROPERTY LOCATION)
list(APPEND _addons ${_loc})
message(STATUS "${Cyan}Using gau2grid${ColourReset}: ${_loc} (version ${gau2grid_VERSION})")

if(${ENABLE_gauxc})
    if(NOT "${gauxc_DIR}" STREQUAL "")
        if("${IntegratorXX_DIR}" STREQUAL "")
            set(IntegratorXX_DIR ${gauxc_DIR}/../IntegratorXX)
        endif()
        set(ExchCXX_DIR ${gauxc_DIR}/../ExchCXX)
    endif()

    get_property(BLAS_LIBRARIES TARGET tgt::lapk PROPERTY INTERFACE_LINK_LIBRARIES)

    find_package(gauxc CONFIG REQUIRED COMPONENTS IntegratorXX ExchCXX)
    get_property(_loc TARGET gauxc::gauxc PROPERTY LOCATION)
    list(APPEND _addons ${_loc})
    message(STATUS "${Cyan}Using gauxc${ColourReset}: ${_loc} (version ${gauxc_VERSION})")
else()
    message(STATUS "Disabled GauXC ${gauxc_DIR}")
endif()

if(${ENABLE_gdma})
    # neither gdma nor pygdma actually needed c-side; runtime detection would
    #   do just as well. included for completeness and so gdma_PYMOD can be
    #   hard-loaded into PYTHONPATH.
    find_package(gdma 2.3.3 CONFIG REQUIRED COMPONENTS Python)
    get_property(_loc TARGET gdma::pygdma PROPERTY LOCATION)
    message(STATUS "${Cyan}Found gdma${ColourReset}: ${_loc} (version ${gdma_VERSION})")
else()
    message(STATUS "Disabled gdma ${gdma_DIR}")
endif()

find_package(
  Libint2
  2.8.1
  CONFIG
  REQUIRED
#  COMPONENTS
#    sss
#    CXX_ho
#    impure_sh
#    "eri_c4_d0_l${MAX_AM_ERI}" eri_c3_d0_l4 eri_c2_d0_l4 onebody_d0_l4
#     eri_c4_d1_l2              eri_c3_d1_l3 eri_c2_d1_l3 onebody_d1_l3
#                                                         onebody_d2_l3
#  OPTIONAL_COMPONENTS
#     eri_c4_d2_l2              eri_c3_d2_l3 eri_c2_d2_l3
  )
get_property(_loc TARGET Libint2::int2 PROPERTY LOCATION)
# defer until upstream provides full-dress targets get_property(Libint2_VERSION TARGET Libint2::int2 PROPERTY Libint2_VERSION)
# defer until upstream provides full-dress targets get_property(Libint2_MAX_AM_ERI TARGET Libint2::int2 PROPERTY Libint2_MAX_AM_ERI)
list(APPEND _addons ${_loc})
message(STATUS "${Cyan}Using Libint2 ${Libint2_MAX_AM_ERI}${ColourReset}: ${_loc} (version ${Libint2_VERSION})")

if(${ENABLE_mdi})
    find_package(mdi 1.2.3 CONFIG REQUIRED)
    get_property(_loc TARGET mdi::mdi PROPERTY LOCATION)
    list(APPEND _addons ${_loc})
    message(STATUS "${Cyan}Using mdi${ColourReset}: ${_loc} (version ${mdi_VERSION})")
else()
    message(STATUS "Disabled MDI ${mdi_DIR}")
endif ()

if(${ENABLE_OpenOrbitalOptimizer})
    find_package(OpenOrbitalOptimizer 0.1 CONFIG REQUIRED)
    get_property(_loc TARGET OpenOrbitalOptimizer::OpenOrbitalOptimizer PROPERTY LOCATION)
    list(APPEND _addons ${_loc})
    message(STATUS "${Cyan}Using OpenOrbitalOptimizer${ColourReset}: ${_loc} (version ${OpenOrbitalOptimizer_VERSION})")
else()
    message(STATUS "Disabled OpenOrbitalOptimizer ${OpenOrbitalOptimizer_DIR}")
endif ()

if(${ENABLE_PCMSolver})
    find_package(PCMSolver 1.2.1 CONFIG REQUIRED)
    get_property(_loc TARGET PCMSolver::pcm PROPERTY LOCATION)
    list(APPEND _addons ${_loc})
    message(STATUS "${Cyan}Using PCMSolver${ColourReset}: ${_loc} (version ${PCMSolver_VERSION})")
else()
    message(STATUS "Disabled PCMSolver ${PCMSolver_DIR}")
endif ()

if(${ENABLE_simint})
    find_package(simint 0.7 CONFIG REQUIRED COMPONENTS am${MAX_AM_ERI} der0)
    get_property(_loc TARGET simint::simint PROPERTY LOCATION)
    list(APPEND _addons ${_loc})
    message(STATUS "${Cyan}Using simint ${simint_MAXAM}${ColourReset}: ${_loc} (version ${simint_VERSION}; vectorization ${simint_VECTOR})")
else()
    message(STATUS "Disabled simint ${simint_DIR}")
endif()

if(${ENABLE_BrianQC})
    find_package(BrianQC 1.1 CONFIG REQUIRED)
    get_property(_loc TARGET BrianQC::static_wrapper PROPERTY LOCATION)
    list(APPEND _addons ${_loc})
    message(STATUS "${Cyan}Using BrianQC${ColourReset}: ${BrianQC_DIR}")
else()
    message(STATUS "Disabled BrianQC ${BrianQC_DIR}")
endif()

find_package(Libxc 7.0.0 CONFIG QUIET COMPONENTS C)
if(NOT ${Libxc_FOUND})
    find_package(Libxc 6.1.0 CONFIG REQUIRED COMPONENTS C)
endif()
get_property(_loc TARGET Libxc::xc PROPERTY LOCATION)
list(APPEND _addons ${_loc})
message(STATUS "${Cyan}Using Libxc${ColourReset}: ${_loc} (version ${Libxc_VERSION})")

if(APPLE)
    set(PRE_LIBRARY_OPTION -Wl,-all_load)
elseif(UNIX)
    set(PRE_LIBRARY_OPTION -Wl,--whole-archive)
    set(POST_LIBRARY_OPTION -Wl,--no-whole-archive)
endif()

# add -pedantic-errors flag here, if FORCE_PEDANTIC is enabled 
if(${FORCE_PEDANTIC}) 
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pedantic-errors")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic-errors")
endif()

# <<<  Build  >>>

if(MSVC)
    # MSVC does not include <cmath> constants, unless _USE_MATH_DEFINES is defined.
    # _CRT_* to squash some getenv, strdup, strncpy, ctime, fopen warnings
    add_compile_definitions(
      _USE_MATH_DEFINES
      _CRT_NONSTDC_NO_DEPRECATE
      _CRT_NONSTDC_NO_WARNINGS
      _CRT_SECURE_NO_WARNINGS
      _ALLOW_COMPILER_AND_STL_VERSION_MISMATCH
      )
    # Set the exception handling model
    add_compile_options("/EHsc")
endif()

include_directories(include)
include_directories(src)
add_subdirectory(src)
if(ENABLE_CYTHONIZE)
    add_subdirectory(driver)
endif()

# <<<  Version  >>>
# * computes version from metadata.py and git info
# * calls cmake to run write_basic_package_version_file

add_custom_target(update_version ALL
                  COMMAND ${Python_EXECUTABLE} versioner.py --metaout ${CMAKE_CURRENT_BINARY_DIR}/metadata.py
                                                            --cmakeout ${CMAKE_CURRENT_BINARY_DIR}/metadata.cmake
                  COMMAND ${CMAKE_COMMAND} -DWTO="${CMAKE_CURRENT_BINARY_DIR}/${CMAKECONFIG_INSTALL_DIR}"
                                           -DPN="psi4"
                                           -P ${CMAKE_CURRENT_BINARY_DIR}/metadata.cmake
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  COMMENT "Generating version info")
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/metadata.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4)

# <<<  Install  >>>

    # <<<  install bin/  >>>
configure_file(run_psi4.py psi4 @ONLY)
install(PROGRAMS ${CMAKE_BINARY_DIR}/psi4
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

if(WIN32)
    # Generate a batch script, which wraps "bin/psi4" script to "python bin/psi4"
    configure_file(psi4.bat psi4.bat @ONLY)
    install(PROGRAMS ${CMAKE_BINARY_DIR}/psi4.bat
            DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})
endif()

    # <<<  install lib/  >>>
if(NOT ENABLE_CYTHONIZE)
    install(DIRECTORY driver
            DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4
            FILES_MATCHING
              PATTERN "*.py"
              PATTERN pytest EXCLUDE
      )
endif()

install(FILES header.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4)

configure_file(__init__.py __init__.py @ONLY)
configure_file(extras.py extras.py @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/__init__.py
              ${CMAKE_BINARY_DIR}/extras.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4)


install(
  FILES ../pytest.ini
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4
  )
install(
  DIRECTORY
    ../tests/pytests/
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4/tests/
  FILES_MATCHING
    PATTERN "test_*.py"
    PATTERN "__init__.py"
    PATTERN "addons.py"
    PATTERN "conftest.py"
    PATTERN "standard_suite_runner.py"
    PATTERN "standard_suite_ref_local.py"
    PATTERN "adcc_reference_data.json"
    PATTERN "oei_reference_data.json"
    PATTERN "f12_libint1.json"
    PATTERN "tdscf_reference_data.json"
    PATTERN "utils.py"
  )
install(
  DIRECTORY
    ../tests/pytests/test_fchk_writer
    ../tests/pytests/test_gauxc_writer
    ../tests/pytests/test_molden_writer
    ../tests/pytests/test_psi4_qcschema
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4/tests/
  )

install(
  DIRECTORY
    ../tests/
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4/tests/
  FILES_MATCHING
    PATTERN ../tests/pytests EXCLUDE
    PATTERN ../tests/cfour EXCLUDE
    PATTERN "test_input.py"
    PATTERN "input.dat"
    PATTERN "input.py"
    PATTERN "__pycache__" EXCLUDE
    PATTERN "*pyc" EXCLUDE
    # below place extra files the ctests-run-as-pytests need
    PATTERN "adcc/formaldehyde-pe-adc1/fa_6w.pot"
    PATTERN "adcc/formaldehyde-pe-adc2/fa_6w.pot"
    PATTERN "ci-property/grid.dat"
    PATTERN "scf-property/grid.dat"
    PATTERN "mp2-property/grid.dat"
    PATTERN "cubeprop*/*cube.ref"
    PATTERN "fcidump/*INTDUMP.ref"
    PATTERN "psithon2/psiaux1/*"
    PATTERN "psithon2/psiaux1/myplugin1/*"
    PATTERN "psithon2/psiaux2/*"
    PATTERN "dftd3/psithon2/psiaux1/*"
    PATTERN "dftd3/psithon2/psiaux1/myplugin1/*"
    PATTERN "dftd3/psithon2/psiaux2/*"
    PATTERN "mints3/L*dat"
    PATTERN "mints3/P*dat"
  )

    # <<<  install psi4 share/ & include/  >>>

install(DIRECTORY share/psi4/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/psi4
        MESSAGE_NEVER
        PATTERN "*.bin" EXCLUDE
        PATTERN "*.pyc" EXCLUDE)

install(FILES share/psi4/quadratures/1_x/R_avail${_endian_ext}.bin
        RENAME quadratures/1_x/R_avail.bin
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/psi4)

install(FILES share/psi4/quadratures/1_x/error${_endian_ext}.bin
        RENAME quadratures/1_x/error.bin
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/psi4)

install(DIRECTORY include/psi4/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/psi4)

    # <<<  install external's share/ to psi4 share/  >>>

if(TARGET libefp::efp)
    # bring libefp's fraglib to PSIDATADIR's attention
    set(_fraglib_in_psi4_dir ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/psi4/efpfrag)

    add_custom_target(consolidate_psidatadir ALL
                      VERBATIM
                      COMMAND ${CMAKE_COMMAND} -E make_directory ${_fraglib_in_psi4_dir}
                      COMMENT "Symlink external resources into PSIDATADIR")

    file(GLOB _dotefps "${libefp_FRAGLIB_DIRS}/*.efp")
    foreach(_dotefp ${_dotefps})
        get_filename_component(_efpfile ${_dotefp} NAME)
        add_custom_command(TARGET consolidate_psidatadir
                           POST_BUILD
                           COMMAND ${CMAKE_COMMAND} -E create_symlink ${_dotefp} ${_fraglib_in_psi4_dir}/${_efpfile})
    endforeach()
endif()

# <<<  Export Config  >>>

configure_file(psi4PluginCache.cmake.in psi4PluginCache.cmake @ONLY)

install(FILES ${CMAKE_BINARY_DIR}/psi4PluginCache.cmake
              ${PSI4_ROOT}/cmake/psi4OptionsTools.cmake
              ${PSI4_ROOT}/cmake/custom_static_library.cmake
              ${PSI4_ROOT}/cmake/custom_cxxstandard.cmake
              ${PSI4_ROOT}/cmake/xhost.cmake
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKECONFIG_INSTALL_DIR})

    # <<<  Export Config for Plugins  >>>
# Determine relative path from Psi4's include directory to PyBind11's
file(RELATIVE_PATH RELATIVE_PYBIND11_INCLUDE_DIR
    ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
    ${pybind11_INCLUDE_DIR})
configure_package_config_file(
        psi4Config.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/${psi4_INSTALL_CMAKEDIR}/psi4Config.cmake
        INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${psi4_INSTALL_CMAKEDIR}/psi4Config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/${psi4_INSTALL_CMAKEDIR}/psi4ConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKECONFIG_INSTALL_DIR})
