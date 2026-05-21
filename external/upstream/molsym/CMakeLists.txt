if(NOT (${CMAKE_DISABLE_FIND_PACKAGE_molsym}))
    include(FindPythonModule)
    find_python_module(molsym QUIET)
endif()

if(${molsym_FOUND})
    message(STATUS "${Cyan}Found molsym${ColourReset}: ${PY_molsym}")
    add_library(molsym_external INTERFACE)

    if(${CMAKE_INSIST_FIND_PACKAGE_molsym})
        message(VERBOSE "Suitable molsym located externally as user insists")
    endif()

else()
    if(${CMAKE_INSIST_FIND_PACKAGE_molsym})
        message(FATAL_ERROR "Suitable molsym could not be externally located as user insists")
    endif()

    include(ExternalProject)
    message(STATUS "Suitable molsym could not be located, ${Magenta}Building molsym${ColourReset} instead.")

    file(TO_NATIVE_PATH "${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}" _install_lib)

    ExternalProject_Add(molsym_external
        BUILD_ALWAYS 1
        GIT_REPOSITORY https://github.com/NASymmetry/MolSym.git
        GIT_TAG f9f93d1
        CONFIGURE_COMMAND ""
        UPDATE_COMMAND ""
        BUILD_COMMAND ""
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND ${Python_EXECUTABLE} -m pip install .
                                                      --target=${_install_lib}
                                                      --no-deps
    )
endif()
