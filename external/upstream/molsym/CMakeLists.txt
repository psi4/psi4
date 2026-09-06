if(${ENABLE_molsym})
    if(NOT (${CMAKE_DISABLE_FIND_PACKAGE_molsym}))
        include(FindPythonModule)
        find_python_module(molsym ATLEAST 1.1.0 QUIET)
    endif()

    if(${molsym_FOUND})
        message(STATUS "${Cyan}Found molsym${ColourReset}: ${PY_molsym} (found version ${molsym_VERSION})")
        add_library(molsym_external INTERFACE)

        if(${CMAKE_INSIST_FIND_PACKAGE_molsym})
            message(VERBOSE "Suitable molsym located externally as user insists")
        endif()
    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_molsym})
            message(FATAL_ERROR "Suitable molsym could not be externally located as user insists")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable molsym could not be located, ${Magenta}Installing molsym${ColourReset} instead.")

        file(TO_NATIVE_PATH "${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}" _install_lib)

        ExternalProject_Add(molsym_external
            BUILD_ALWAYS 1
            DOWNLOAD_COMMAND ""
            CONFIGURE_COMMAND ""
            UPDATE_COMMAND ""
            BUILD_COMMAND ""
            BUILD_IN_SOURCE 1
            INSTALL_COMMAND ${Python_EXECUTABLE} -m pip install molsym==1.1.0 --target=${_install_lib} --no-deps --upgrade
        )
    endif()
else()
    message(STATUS "Optional addon ${Blue}unsought: molsym${ColourReset} (molsym; deferrable to runtime)")
    add_library(molsym_external INTERFACE)
endif()