if(NOT (${CMAKE_DISABLE_FIND_PACKAGE_optking}))
    include(FindPythonModule)
    find_python_module(optking ATLEAST 0.4.1 QUIET)  # edit in codedeps
endif()

if(${optking_FOUND})
    message(STATUS "${Cyan}Found optking${ColourReset}: ${PY_optking} (found version ${optking_VERSION})")
    add_library(optking_external INTERFACE)  # dummy

    if(${CMAKE_INSIST_FIND_PACKAGE_optking})
        message(VERBOSE "Suitable optking located externally as user insists")
    endif()

else()
    if(${CMAKE_INSIST_FIND_PACKAGE_optking})
        message(FATAL_ERROR "Suitable optking could not be externally located as user insists")
    endif()

    include(ExternalProject)
    message(STATUS "Suitable optking could not be located, ${Magenta}Building optking${ColourReset} instead.")

    file(TO_NATIVE_PATH "${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}" _install_lib)

    ExternalProject_Add(optking_external
        BUILD_ALWAYS 1
        URL https://github.com/psi-rking/optking/archive/0.4.1.tar.gz  # edit in codedeps
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

