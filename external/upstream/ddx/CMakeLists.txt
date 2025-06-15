if(${ENABLE_ddx})
    if(NOT (${CMAKE_DISABLE_FIND_PACKAGE_ddx}))
        include(FindPythonModule)
        find_python_module(pyddx ATLEAST 0.4.2 QUIET)
    endif()

    if(${pyddx_FOUND})
        message(STATUS "${Cyan}Found ddx${ColourReset}: ${PY_pyddx} (found version ${pyddx_VERSION})")
        add_library(ddx_external INTERFACE)  # dummy

        if(${CMAKE_INSIST_FIND_PACKAGE_ddx})
            message(VERBOSE "Suitable ddx located externally as user insists")
        endif()

    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_ddx})
            message(FATAL_ERROR "Suitable ddx could not be externally located as user insists")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable ddx could not be located, ${Magenta}Building ddx${ColourReset} instead.")

        file(TO_NATIVE_PATH "${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}" _install_lib)
        file(TO_NATIVE_PATH "${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}" _install_bin)

        ExternalProject_Add(ddx_external
            BUILD_ALWAYS 1
            URL https://github.com/ddsolvation/ddX/archive/a43f9b4.tar.gz  # edit in codedeps
            DOWNLOAD_NO_EXTRACT 1
            CONFIGURE_COMMAND ""
            UPDATE_COMMAND ""
            BUILD_COMMAND ""
            BUILD_IN_SOURCE 1
            INSTALL_COMMAND ${Python_EXECUTABLE} -m pip install <DOWNLOADED_FILE>
                                                      --target=${_install_lib}
                                                      --no-deps
            )
    endif()
else()
    add_library(ddx_external INTERFACE)  # dummy
endif()
