if(NOT (${CMAKE_DISABLE_FIND_PACKAGE_optking}))
    include(FindPythonModule)
    find_python_module(optking ATLEAST 0.2.1 QUIET)
endif()

if(${optking_FOUND})
    message(STATUS "${Cyan}Found optking${ColourReset}: ${PY_optking} (found version ${optking_VERSION})")
    add_library(optking_external INTERFACE)  # dummy

else()
    if(${CMAKE_INSIST_FIND_PACKAGE_optking})
        message(FATAL_ERROR "Suitable optking could not be externally located as user insists")
    endif()

    include(ExternalProject)
    message(STATUS "Suitable optking could not be located, ${Magenta}Building optking${ColourReset} instead.")

    file(TO_NATIVE_PATH "${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}" _install_lib)

    ExternalProject_Add(optking_external
        BUILD_ALWAYS 1
        URL https://github.com/psi-rking/optking/tarball/master
        CONFIGURE_COMMAND ""
        UPDATE_COMMAND ""
        BUILD_COMMAND ${Python_EXECUTABLE} setup.py build
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND ${Python_EXECUTABLE} setup.py install
                                                      --record=record.txt
                                                      --single-version-externally-managed
                                                      --install-lib=${_install_lib})
endif()

