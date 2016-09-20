cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
#project(psi4-core LANGUAGES C CXX Fortran VERSION 1.0)
project(psi4-core LANGUAGES C CXX VERSION 1.0)
# no Fortran in psi4-core proper, but language needs to be declared
#   for CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES to be populated so
#   static Fortran add-ons can be linked
if(${Fortran_ENABLED})
    enable_language(Fortran)
endif()

#set(CMAKE_CXX_STANDARD 11) enable if they fix this for Intel
list(APPEND CMAKE_MODULE_PATH ${PSI4_ROOT}/cmake)

include(Psi4Macros)
include(TestRestrict)
include(GNUInstallDirs)
test_restrict(restrict)

### >> Locate Python Interpreter and Library << ###
#find_package(PythonLibsNew REQUIRED)
find_package(PythonInterp REQUIRED)
find_package(PythonLibs ${PYTHON_VERSION_MAJOR} REQUIRED)
find_package(pybind11 CONFIG REQUIRED
             HINTS ${CMAKE_PREFIX_PATH}/pybind11)
find_package(Threads)
find_package(DL)

find_package(libint CONFIG REQUIRED COMPONENTS ${MAX_AM_ERI}
        HINTS ${CMAKE_PREFIX_PATH}/libint)

set_property(GLOBAL PROPERTY PSI4_MODULES "")

#optional_plugin(Ambit)

if(${ENABLE_PCMSolver})
    find_package(PCMSolver MODULE REQUIRED)  # todo CONFIG
                           #HINTS ${CMAKE_PREFIX_PATH}/pcmsolver)
endif ()
if(${ENABLE_CheMPS2})
    find_package(CheMPS2 CONFIG REQUIRED
                     HINTS ${CMAKE_PREFIX_PATH}/chemps2)
endif ()
if(${ENABLE_dkh})
    find_package(dkh CONFIG REQUIRED
                     HINTS ${CMAKE_PREFIX_PATH}/dkh)
endif()
if(${ENABLE_liberd})
    find_package(liberd CONFIG REQUIRED
                        HINTS ${CMAKE_PREFIX_PATH}/liberd)
endif()
if (${ENABLE_gdma})
    find_package(gdma CONFIG REQUIRED
                 HINTS ${CMAKE_PREFIX_PATH}/gdma)
endif()
find_package(libefp CONFIG REQUIRED COMPONENTS shallow
             HINTS ${CMAKE_PREFIX_PATH}/libefp)

include_directories(include)
include_directories(src)

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

add_subdirectory(src)

configure_file("psi4-config.in" "psi4-config.tmp" @ONLY)
install(PROGRAMS ${CMAKE_BINARY_DIR}/psi4-config
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

configure_file("psi4Config.cmake.in" "psi4Config.cmake" @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/psi4Config.cmake
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/cmake/psi4)

install(DIRECTORY share/psi4/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/psi4
        PATTERN fsapt EXCLUDE
        PATTERN "*.pyc" EXCLUDE)

install(DIRECTORY include/psi4/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/psi4)

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

# <<< Export Config for Plugins >>>
include(CMakePackageConfigHelpers)
# Determine relative path from Psi4's include directory to PyBind11's
file(RELATIVE_PATH RELATIVE_PYBIND11_INCLUDE_DIR
    ${DESTDIR}${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR} 
    ${pybind11_INCLUDE_DIR})
set(CMAKECONFIG_INSTALL_DIR "${CMAKE_INSTALL_DATADIR}/cmake/Psi4Plugins")
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
