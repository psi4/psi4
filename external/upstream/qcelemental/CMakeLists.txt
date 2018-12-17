if(NOT (${CMAKE_DISABLE_FIND_PACKAGE_qcelemental}))
    include(FindPythonModule)
    find_python_module(qcelemental ATLEAST 0.1.3 QUIET)
endif()

if(${qcelemental_FOUND})
    message(STATUS "${Cyan}Found qcelemental${ColourReset}: ${PY_qcelemental} (found version ${qcelemental_VERSION})")
    add_library(qcelemental_external INTERFACE)  # dummy

else()
    if(${CMAKE_INSIST_FIND_PACKAGE_qcelemental})
        message(FATAL_ERROR "Suitable qcelemental could not be externally located as user insists")
    endif()

    include(ExternalProject)
    message(STATUS "Suitable qcelemental could not be located, ${Magenta}Building qcelemental${ColourReset} instead.")

    file(TO_NATIVE_PATH "${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}" _install_lib)

    ExternalProject_Add(qcelemental_external
        BUILD_ALWAYS 1
        GIT_REPOSITORY https://github.com/MolSSI/QCElemental
        GIT_TAG v0.1.3
        CONFIGURE_COMMAND ""
        UPDATE_COMMAND ""
        BUILD_COMMAND ${PYTHON_EXECUTABLE} setup.py build
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND ${PYTHON_EXECUTABLE} setup.py install
                                                      --record=record.txt
                                                      --single-version-externally-managed
                                                      --install-lib=${_install_lib})
endif()
