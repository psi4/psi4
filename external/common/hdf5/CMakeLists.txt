if(ENABLE_ambit OR ENABLE_CheMPS2 OR ENABLE_Einsums)
    project(TargetHDF5)
    set(PN ${PROJECT_NAME})

    # <<<  "Build"  >>>

    find_package(HDF5 QUIET REQUIRED)

    # HDF5 1.x's FindHDF5 sets HDF5_LIBRARIES to absolute library file paths,
    # but HDF5 2.x ships its own HDF5Config.cmake which (in config mode) sets
    # HDF5_LIBRARIES to imported-target names like 'hdf5-shared'. Those
    # targets are only visible in the scope of the HDF5 config file, so once
    # we stash them into our INTERFACE library and re-export it for
    # downstream consumers the bare names fall through as link library
    # names ('-lhdf5-shared'), and the link fails. Resolve any imported
    # targets to their actual library file paths before storing.
    set(_HDF5_RESOLVED_LIBS)
    foreach(_lib IN LISTS HDF5_LIBRARIES)
        if(TARGET ${_lib})
            get_target_property(_imp_cfgs ${_lib} IMPORTED_CONFIGURATIONS)
            set(_loc)
            if(_imp_cfgs)
                list(GET _imp_cfgs 0 _first_cfg)
                get_target_property(_loc ${_lib} IMPORTED_LOCATION_${_first_cfg})
            endif()
            if(NOT _loc)
                get_target_property(_loc ${_lib} IMPORTED_LOCATION)
            endif()
            if(_loc)
                list(APPEND _HDF5_RESOLVED_LIBS ${_loc})
            else()
                list(APPEND _HDF5_RESOLVED_LIBS ${_lib})
            endif()
        else()
            list(APPEND _HDF5_RESOLVED_LIBS ${_lib})
        endif()
    endforeach()

    add_library(hdf5 INTERFACE)
    set_property(TARGET hdf5 PROPERTY INTERFACE_LINK_LIBRARIES ${_HDF5_RESOLVED_LIBS})
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

    set(CMAKECONFIG_INSTALL_DIR "${TargetHDF5_INSTALL_CMAKEDIR}")
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
        ${CMAKE_COMMAND} "-DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}" -P "${PROJECT_BINARY_DIR}/cmake_install.cmake"
        COMMENT "Installing Psi4-detected HDF5")

    set(${PN}_DIR ${STAGED_INSTALL_PREFIX}/${TargetHDF5_INSTALL_CMAKEDIR} CACHE PATH "path to externally detected ${PN}Config.cmake" FORCE)
else()
    add_library(hdf5_external INTERFACE)  # dummy
endif()
