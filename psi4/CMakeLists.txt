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

find_package(PythonLibsNew REQUIRED)
#find_package(PythonInterp REQUIRED)
#find_package(PythonLibs ${PYTHON_VERSION_MAJOR} REQUIRED)
find_package(pybind11 CONFIG REQUIRED)
find_package(Threads)
find_package(DL)

set_property(GLOBAL PROPERTY PSI4_MODULES "")

if(${ENABLE_CheMPS2})
    find_package(CheMPS2 CONFIG REQUIRED)
endif ()
if(${ENABLE_dkh})
    find_package(dkh CONFIG REQUIRED)
endif()
find_package(libefp CONFIG REQUIRED COMPONENTS shallow)
if(${ENABLE_erd})
    find_package(erd CONFIG REQUIRED)
endif()
if(${ENABLE_gdma})
    find_package(gdma CONFIG REQUIRED)
endif()
find_package(libint CONFIG REQUIRED COMPONENTS ${MAX_AM_ERI})
if(${ENABLE_PCMSolver})
    find_package(PCMSolver MODULE REQUIRED)  # todo CONFIG
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
                  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/versioner.py
                                               ${CMAKE_CURRENT_SOURCE_DIR}
                                               ${CMAKE_CURRENT_BINARY_DIR}/metadata.py
                  COMMENT "Generating version info")
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/metadata.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/psi4)
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
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/psi4)

install(DIRECTORY driver
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/psi4
        FILES_MATCHING PATTERN "*.py")

install(FILES header.py extras.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/psi4)

configure_file(__init__.py __init__.py @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/__init__.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/psi4)

    # <<<  install psi4 share/ & include/  >>>

install(DIRECTORY share/psi4/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/psi4
        MESSAGE_NEVER
        PATTERN "*.pyc" EXCLUDE)

install(DIRECTORY include/psi4/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/psi4)

    # <<<  install external's share/ to psi4 share/  >>>

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
#set(CMAKECONFIG_INSTALL_DIR "${CMAKE_INSTALL_DATADIR}/cmake/Psi4Plugins")
configure_package_config_file(
        share/cmake/Psi4Plugins/Psi4PluginsConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/share/cmake/Psi4Plugins/Psi4PluginsConfig.cmake
        INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/share/cmake/Psi4Plugins/Psi4PluginsConfigVersion.cmake
        VERSION ${PSI4_VERSION}
        COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/share/cmake/Psi4Plugins/Psi4PluginsConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/share/cmake/Psi4Plugins/Psi4PluginsConfigVersion.cmake
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})
