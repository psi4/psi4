if(ENABLE_Ambit OR ENABLE_CheMPS2)
    project(TargetHDF5)
    set(PN ${PROJECT_NAME})

    # <<<  "Build"  >>>

    find_package(HDF5 QUIET REQUIRED)

    add_library(hdf5 INTERFACE)
    set_property(TARGET hdf5 PROPERTY INTERFACE_LINK_LIBRARIES ${HDF5_LIBRARIES})
    set_property(TARGET hdf5 PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIRS})
    set (${PN}_VERSION ${HDF5_VERSION})

    get_property(_ill TARGET hdf5 PROPERTY INTERFACE_LINK_LIBRARIES)
    get_property(_iid TARGET hdf5 PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
    set(${PN}_MESSAGE "Found HDF5: ${_ill} (found version ${${PN}_VERSION})")
    list(GET _ill 0 _ill0)
    message(STATUS "${Cyan}Found HDF5${ColourReset}: ${_ill0};... (found version ${${PN}_VERSION})")

    install(TARGETS hdf5
            EXPORT "${PN}Targets")

    # <<<  Export Config  >>>

    include(CMakePackageConfigHelpers)

    set(CMAKECONFIG_INSTALL_DIR "share/cmake/${PN}")
    configure_package_config_file(${PN}Config.cmake.in
                                  "${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake"
                                  INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})
    write_basic_package_version_file (${CMAKE_CURRENT_BINARY_DIR}/${PN}ConfigVersion.cmake
                                      VERSION ${${PN}_VERSION}
                                      COMPATIBILITY ExactVersion)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake
                  ${CMAKE_CURRENT_BINARY_DIR}/${PN}ConfigVersion.cmake
            DESTINATION ${CMAKECONFIG_INSTALL_DIR})
    install(EXPORT "${PN}Targets"
            NAMESPACE "tgt::"
            DESTINATION ${CMAKECONFIG_INSTALL_DIR})

    unset(HDF5_FOUND)
    unset(HDF5_VERSION)
    unset(HDF5_LIBRARIES)
    unset(HDF5_INCLUDE_DIRS)

    add_custom_target(hdf5_external
        ${CMAKE_COMMAND} "-DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/stage${CMAKE_INSTALL_PREFIX}" -P "${PROJECT_BINARY_DIR}/cmake_install.cmake"
        COMMENT "Installing Psi4-detected HDF5")

    set(${PN}_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/${PN} CACHE PATH "path to externally detected ${PN}Config.cmake" FORCE)
else()
    add_library(hdf5_external INTERFACE)  # dummy
endif()
