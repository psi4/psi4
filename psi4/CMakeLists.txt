cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(psi4-core LANGUAGES C CXX VERSION 1.0)
# no Fortran in psi4-core proper, but language needs to be declared
#   for CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES to be populated so
#   static Fortran add-ons can be linked
if(${Fortran_ENABLED})
    enable_language(Fortran)
endif()

list(APPEND CMAKE_MODULE_PATH ${PSI4_ROOT}/cmake)

include(Psi4Macros)
include(TestRestrict)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
test_restrict(restrict)

# <<<  Marshal Dependencies & Add-ons  >>>

# using Python in psi4-core
# * interpreter: run the versioner (and elsewhere to run tests, docs)
# * headers: Python.h for pybind11 to build against
# * library: not at all
find_package(PythonLibsNew REQUIRED)
find_package(pybind11 CONFIG REQUIRED)
message(STATUS "Found pybind11: ${pybind11_INCLUDE_DIR} (found version ${pybind11_VERSION} for Py${PYTHON_VERSION_STRING})")
set_property(TARGET pybind11::pybind11 APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${PYTHON_INCLUDE_DIRS})

find_package(Threads)
find_package(DL)

set_property(GLOBAL PROPERTY PSI4_MODULES "")

if(${ENABLE_CheMPS2})
    find_package(CheMPS2 1.8.3 CONFIG REQUIRED)
    get_property(_loc TARGET CheMPS2::chemps2 PROPERTY LOCATION)
    message(STATUS "Found CheMPS2: ${_loc} (found version ${CheMPS2_VERSION})")
endif ()
if(${ENABLE_dkh})
    find_package(dkh CONFIG REQUIRED)
    get_property(_loc TARGET dkh::dkh PROPERTY LOCATION)
    message(STATUS "Found dkh: ${_loc} (found version ${dkh_VERSION})")
endif()
find_package(libefp CONFIG REQUIRED COMPONENTS shallow)
get_property(_loc TARGET libefp::efp PROPERTY LOCATION)
message(STATUS "Found libefp: ${_loc} (found version ${libefp_VERSION})")
if(${ENABLE_erd})
    find_package(erd CONFIG REQUIRED)
    get_property(_loc TARGET erd::erd PROPERTY LOCATION)
    message(STATUS "Found erd: ${_loc} (found version ${erd_VERSION})")
endif()
if(${ENABLE_gdma})
    find_package(gdma CONFIG REQUIRED)
    get_property(_loc TARGET gdma::gdma PROPERTY LOCATION)
    message(STATUS "Found gdma: ${_loc} (found version ${gdma_VERSION})")
endif()
find_package(libint CONFIG REQUIRED COMPONENTS ${MAX_AM_ERI})
get_property(_loc TARGET libint::int PROPERTY LOCATION)
message(STATUS "Found libint: ${_loc} (found version ${libint_VERSION})")
if(${ENABLE_PCMSolver})
    find_package(PCMSolver MODULE REQUIRED)  # todo CONFIG
    #get_property(_loc TARGET PCMSolver::PCMSolver PROPERTY LOCATION)
    #message(STATUS "Found PCMSolver: ${_loc} (found version ${PCMSolver_VERSION})")
endif ()

if(APPLE)
    set(PRE_LIBRARY_OPTION -Wl,-all_load)
elseif(UNIX)
    set(PRE_LIBRARY_OPTION -Wl,--whole-archive)
    set(POST_LIBRARY_OPTION -Wl,--no-whole-archive)
endif()

# <<<  Build  >>>

include_directories(include)
include_directories(src)
add_subdirectory(src)

# <<<  Version  >>>
add_custom_target(update_version ALL
                  COMMAND ${PYTHON_EXECUTABLE} versioner.py --metaout ${CMAKE_CURRENT_BINARY_DIR}/metadata.py
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  COMMENT "Generating version info")
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/metadata.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4)
    # TODO absorb some version info into CM variables, write to package, adj psi4ConfigVersion to sort dev

# <<<  Install  >>>

    # <<<  install bin/  >>>
configure_file(run_psi4.py psi4 @ONLY)
install(PROGRAMS ${CMAKE_BINARY_DIR}/psi4
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

# configure_file("psi4-config.in" "psi4-config.tmp" @ONLY)
# install(PROGRAMS ${CMAKE_BINARY_DIR}/psi4-config
#         DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

    # <<<  install lib/  >>>
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4)

install(DIRECTORY share/psi4/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/psi4
        PATTERN "*.pyc" EXCLUDE)

install(DIRECTORY include/psi4/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/psi4)

install(DIRECTORY driver
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4
        FILES_MATCHING PATTERN "*.py")

install(FILES header.py extras.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4)

configure_file(__init__.py __init__.py @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/__init__.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4)

    # <<<  install psi4 share/ & include/  >>>

install(DIRECTORY share/psi4/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/psi4
        MESSAGE_NEVER
        PATTERN "*.pyc" EXCLUDE)

install(DIRECTORY include/psi4/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/psi4)

configure_file(config.py.in config.py @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/config.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4)

    # <<<  install external's share/ to psi4 share/  >>>

        # bring libefp's fraglib to PSIDATADIR's attention
# bring libefp's fraglib to PSIDATADIR's attention
set(_fraglib_in_psi4_dir ${DESTDIR}${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/psi4/efpfrag)

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

# <<<  Export Config  >>>

    # explicit "share" not "DATADIR" for CMake search path
set(CMAKECONFIG_INSTALL_DIR "share/cmake/psi4")
#configure_package_config_file(psi4Config.cmake.in
#                              "${CMAKE_CURRENT_BINARY_DIR}/psi4Config.cmake"
#                              INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})
#write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/psi4ConfigVersion.cmake
#                                 VERSION ${psi4_VERSION}
#                                 COMPATIBILITY SameMajorVersion)
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/psi4Config.cmake
#              ${CMAKE_CURRENT_BINARY_DIR}/psi4ConfigVersion.cmake
#        DESTINATION ${CMAKECONFIG_INSTALL_DIR})
#install(EXPORT "psi4Targets"
#        NAMESPACE "psi4::"
#        DESTINATION ${CMAKECONFIG_INSTALL_DIR})

    # <<<  Export Config for Plugins  >>>
# Determine relative path from Psi4's include directory to PyBind11's
file(RELATIVE_PATH RELATIVE_PYBIND11_INCLUDE_DIR
    ${DESTDIR}${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
    ${pybind11_INCLUDE_DIR})
set(CMAKECONFIG_INSTALL_DIR "${CMAKE_INSTALL_DATADIR}/cmake/psi4")
configure_package_config_file(
        share/cmake/psi4/psi4Config.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/share/cmake/psi4/psi4Config.cmake
        INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/share/cmake/psi4/psi4ConfigVersion.cmake
        VERSION ${PSI4_VERSION}
        COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/share/cmake/psi4/psi4Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/share/cmake/psi4/psi4ConfigVersion.cmake
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})
