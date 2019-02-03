cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(psi4-core
        LANGUAGES C CXX)
# no Fortran in psi4-core proper, but language needs to be declared
#   for CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES to be populated so
#   static Fortran add-ons can be linked
if(${Fortran_ENABLED})
    enable_language(Fortran)
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
# GNUInstallDirs "DATADIR" wrong here; CMake search path wants "share".
set(CMAKECONFIG_INSTALL_DIR "share/cmake/psi4")

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

if(ENABLE_ambit OR ENABLE_CheMPS2)
    find_package(TargetHDF5 CONFIG REQUIRED)
    get_property(_ill TARGET tgt::hdf5 PROPERTY INTERFACE_LINK_LIBRARIES)
    list(GET _ill 0 _ill0)
    list(APPEND _addons ${_ill})
    message(STATUS "${Cyan}Using HDF5${ColourReset}: ${_ill0};... (found version ${${PN}_VERSION})")
else()
    message(STATUS "Disabled HDF5")
endif()

#  <<  Pybind11 & Python  >>
if(MSVC)
  set(PYBIND11_CPP_STANDARD "/std:c++${CMAKE_CXX_STANDARD}")
else()
  set(PYBIND11_CPP_STANDARD "-std=c++${CMAKE_CXX_STANDARD}")
endif()
find_package(pybind11 2.2.3 CONFIG REQUIRED)
message(STATUS "${Cyan}Using pybind11${ColourReset}: ${pybind11_INCLUDE_DIR} (version ${pybind11_VERSION} for Py${PYTHON_VERSION_STRING} and ${PYBIND11_CPP_STANDARD})")
message(STATUS "${Cyan}Using Python ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}${ColourReset}: ${PYTHON_EXECUTABLE}")

find_package(DL)

if(${ENABLE_ambit})
    find_package(ambit CONFIG REQUIRED)
    get_property(_loc TARGET ambit::ambit PROPERTY LOCATION)
    message(STATUS "${Cyan}Using ambit${ColourReset}: ${_loc} (version ${ambit_VERSION})")
else()
    message(STATUS "Disabled ambit")
endif ()

if(${ENABLE_CheMPS2})
    find_package(CheMPS2 1.8.7 CONFIG REQUIRED)
    get_property(_loc TARGET CheMPS2::chemps2 PROPERTY LOCATION)
    list(APPEND _addons ${_loc})
    message(STATUS "${Cyan}Using CheMPS2${ColourReset}: ${_loc} (version ${CheMPS2_VERSION})")
else()
    message(STATUS "Disabled CheMPS2")
endif ()

if(${ENABLE_dkh})
    find_package(dkh 1.2 CONFIG REQUIRED)
    get_property(_loc TARGET dkh::dkh PROPERTY LOCATION)
    list(APPEND _addons ${_loc})
    message(STATUS "${Cyan}Using dkh${ColourReset}: ${_loc} (version ${dkh_VERSION})")
else()
    message(STATUS "Disabled dkh")
endif()

if(${ENABLE_libefp})
    # neither libefp nor pylibefp actually needed c-side; runtime detection would
    #   do just as well. included for completeness and so pylibefp_PYMOD can be
    #   hard-loaded into PYTHONPATH. also incl for finding efp frag files below.
    find_package(libefp 1.5.0 CONFIG REQUIRED COMPONENTS shallow)
    get_property(_loc TARGET libefp::efp PROPERTY LOCATION)
    list(APPEND _addons ${_loc})
    message(STATUS "${Cyan}Using libefp${ColourReset}: ${_loc} (version ${libefp_VERSION})")

    find_package(pylibefp 0.4 CONFIG REQUIRED)
    get_property(_loc TARGET pylibefp::core PROPERTY LOCATION)
    list(APPEND _addons ${_loc})
    message(STATUS "${Cyan}Using pylibefp${ColourReset}: ${_loc} (version ${pylibefp_VERSION})")
else()
    message(STATUS "Disabled libefp")
endif()

if(${ENABLE_erd})
    find_package(erd 3.0.6 CONFIG REQUIRED)
    get_property(_loc TARGET erd::erd PROPERTY LOCATION)
    list(APPEND _addons ${_loc})
    message(STATUS "${Cyan}Using erd${ColourReset}: ${_loc} (version ${erd_VERSION})")
else()
    message(STATUS "Disabled erd")
endif()

find_package(gau2grid 1.2 CONFIG REQUIRED COMPONENTS gaussian)
get_property(_loc TARGET gau2grid::gg PROPERTY LOCATION)
list(APPEND _addons ${_loc})
message(STATUS "${Cyan}Using gau2grid${ColourReset}: ${_loc} (version ${gau2grid_VERSION})")

if(${ENABLE_gdma})
    find_package(gdma 2.2.6 CONFIG REQUIRED)
    get_property(_loc TARGET gdma::gdma PROPERTY LOCATION)
    list(APPEND _addons ${_loc})
    message(STATUS "${Cyan}Using gdma${ColourReset}: ${_loc} (version ${gdma_VERSION})")
else()
    message(STATUS "Disabled gdma")
endif()

find_package(Libint 1.2.0 CONFIG REQUIRED COMPONENTS ${MAX_AM_ERI})
get_property(_loc TARGET Libint::int PROPERTY LOCATION)
list(APPEND _addons ${_loc})
message(STATUS "${Cyan}Using Libint ${Libint_MAX_AM_ERI}${ColourReset}: ${_loc} (version ${Libint_VERSION})")

if(${ENABLE_PCMSolver})
    find_package(PCMSolver 1.2.1 CONFIG REQUIRED)
    get_property(_loc TARGET PCMSolver::pcm PROPERTY LOCATION)
    list(APPEND _addons ${_loc})
    message(STATUS "${Cyan}Using PCMSolver${ColourReset}: ${_loc} (version ${PCMSolver_VERSION})")
else()
    message(STATUS "Disabled PCMSolver")
endif ()

if(${ENABLE_simint})
    find_package(simint 0.7 CONFIG REQUIRED COMPONENTS am${MAX_AM_ERI} der0)
    get_property(_loc TARGET simint::simint PROPERTY LOCATION)
    list(APPEND _addons ${_loc})
    message(STATUS "${Cyan}Using simint ${simint_MAXAM}${ColourReset}: ${_loc} (version ${simint_VERSION}; vectorization ${simint_VECTOR})")
else()
    message(STATUS "Disabled simint")
endif()

find_package(Libxc 4.3.0 CONFIG REQUIRED)
get_property(_loc TARGET Libxc::xc PROPERTY LOCATION)
list(APPEND _addons ${_loc})
message(STATUS "${Cyan}Using Libxc${ColourReset}: ${_loc} (version ${Libxc_VERSION})")

if(APPLE)
    set(PRE_LIBRARY_OPTION -Wl,-all_load)
elseif(UNIX)
    set(PRE_LIBRARY_OPTION -Wl,--whole-archive)
    set(POST_LIBRARY_OPTION -Wl,--no-whole-archive)
endif()

# <<<  Build  >>>

if(MSVC)
    # MSVC does not include <cmath> constants, unless _USE_MATH_DEFINES is defined.
    add_definitions("/D_USE_MATH_DEFINES")
    # Set the exception handling model
    add_definitions("/EHsc")
endif()

include_directories(include)
include_directories(src)
add_subdirectory(src)

# <<<  Version  >>>
# * computes version from metadata.py and git info
# * calls cmake to run write_basic_package_version_file

add_custom_target(update_version ALL
                  COMMAND ${PYTHON_EXECUTABLE} versioner.py --metaout ${CMAKE_CURRENT_BINARY_DIR}/metadata.py
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
install(DIRECTORY driver
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4
        FILES_MATCHING PATTERN "*.py")

install(FILES header.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4)

configure_file(__init__.py __init__.py @ONLY)
configure_file(extras.py extras.py @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/__init__.py
              ${CMAKE_BINARY_DIR}/extras.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4)

install(FILES ../tests/pytest/conftest.py
              ../tests/pytest/addons.py
              ../tests/pytest/utils.py
              ../tests/pytest/test_psi4.py
              ../tests/pytest/test_addons.py
              ../tests/pytest/test_molecule_to_string.py
              ../tests/pytest/test_matrix.py
              ../tests/pytest/test_vector.py
              ../tests/pytest/test_qcvars.py
              ../tests/pytest/test_dft_benchmarks.py
              ../tests/pytest/test_vibanalysis.py
              ../tests/pytest/test_dftd3.py
              ../tests/pytest/test_fcidump_energy.py
              ../tests/pytest/test_solvers.py
              ../tests/pytest/test_tdscf_products.py
              ../tests/pytest/test_tdscf_excitations.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4/tests/)

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
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})

    # <<<  Export Config for Plugins  >>>
# Determine relative path from Psi4's include directory to PyBind11's
file(RELATIVE_PATH RELATIVE_PYBIND11_INCLUDE_DIR
    ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
    ${pybind11_INCLUDE_DIR})
configure_package_config_file(
        psi4Config.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/share/cmake/psi4/psi4Config.cmake
        INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/share/cmake/psi4/psi4Config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/share/cmake/psi4/psi4ConfigVersion.cmake
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})
