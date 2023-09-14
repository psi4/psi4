if(${ENABLE_adcc})
    if(NOT (${CMAKE_DISABLE_FIND_PACKAGE_adcc}))
        include(FindPythonModule)
        find_python_module(adcc ATLEAST 0.15.16 QUIET)
    endif()

    if(${adcc_FOUND})
        message(STATUS "${Cyan}Found adcc${ColourReset}: ${PY_adcc} (found version ${adcc_VERSION})")
        add_library(adcc_external INTERFACE)  # dummy

    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_adcc})
            message(FATAL_ERROR "Suitable adcc could not be externally located as user insists")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable adcc could not be located, ${Magenta}Building adcc${ColourReset} instead.")

        ExternalProject_Add(adcc_external
            BUILD_ALWAYS 1
            URL https://github.com/adc-connect/adcc/archive/v0.15.16.tar.gz
            CONFIGURE_COMMAND ""
            UPDATE_COMMAND ""
            BUILD_COMMAND ${Python_EXECUTABLE} setup.py build
            BUILD_IN_SOURCE 1
            INSTALL_COMMAND ${Python_EXECUTABLE} setup.py install
                                                          --record=record.txt
                                                          --single-version-externally-managed
                                                          --install-lib=${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR})
    endif()
else()
    add_library(adcc_external INTERFACE)  # dummy
endif()
