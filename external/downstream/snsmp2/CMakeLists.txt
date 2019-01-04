if(${ENABLE_snsmp2})
    if(NOT (${CMAKE_DISABLE_FIND_PACKAGE_snsmp2}))
        include(FindPythonModule)
        find_python_module(snsmp2 ATLEAST 1.0.0 QUIET)
    endif()

    if(${snsmp2_FOUND})
        message(STATUS "${Cyan}Found snsmp2${ColourReset}: ${PY_snsmp2} (found version ${snsmp2_VERSION})")
        add_library(snsmp2_external INTERFACE)  # dummy

    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_snsmp2})
            message(FATAL_ERROR "Suitable snsmp2 could not be externally located as user insists")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable snsmp2 could not be located, ${Magenta}Building snsmp2${ColourReset} instead.")

        ExternalProject_Add(snsmp2_external
            DEPENDS qcelemental_external
            BUILD_ALWAYS 1
            GIT_REPOSITORY https://github.com/DEShawResearch/sns-mp2
            GIT_TAG 99e2a9c  # v1.0 + 16
            CONFIGURE_COMMAND ""
            UPDATE_COMMAND ""
            BUILD_COMMAND ${PYTHON_EXECUTABLE} setup.py build
            BUILD_IN_SOURCE 1
            INSTALL_COMMAND ${PYTHON_EXECUTABLE} setup.py install
                                                          --record=record.txt
                                                          --single-version-externally-managed
                                                          --install-lib=${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR})
    endif()
else()
    add_library(snsmp2_external INTERFACE)  # dummy
endif()
                                                          #--install-scripts=${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}
                                                          #--install-data=${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}
                                                          #--install-headers=${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
