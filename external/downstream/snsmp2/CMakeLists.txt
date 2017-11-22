if(${ENABLE_snsmp2})
    if(NOT (${CMAKE_DISABLE_FIND_PACKAGE_snsmp2}))
        include(FindPythonModule)
        find_python_module(snsmp2 QUIET)
    endif()

    if(${snsmp2_FOUND})
        message(STATUS "${Cyan}Found snsnp2${ColourReset}: ${PY_snsmp2}") # (found version ${snsmp2_VERSION})")
    else()
        include(ExternalProject)
        message(STATUS "Suitable snsmp2 could not be located, ${Magenta}Building snsmp2${ColourReset} instead.")


        ExternalProject_Add(snsmp2_external
            #GIT_REPOSITORY https://github.com/DEShawResearch/sns-mp2
            #GIT_TAG 1.0
            GIT_REPOSITORY https://github.com/loriab/sns-mp2
            GIT_TAG psifixes
            CONFIGURE_COMMAND ""
            UPDATE_COMMAND ""
            BUILD_COMMAND ${PYTHON_EXECUTABLE} setup.py build
            BUILD_IN_SOURCE 1
            INSTALL_COMMAND ${PYTHON_EXECUTABLE} setup.py install
                                                          --record=record.txt
                                                          --single-version-externally-managed
                                                          --install-lib=${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR})
    endif()
endif()

                                                          #--install-scripts=${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}
                                                          #--install-data=${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}
                                                          #--install-headers=${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
