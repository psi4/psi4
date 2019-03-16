if(NOT (${CMAKE_DISABLE_FIND_PACKAGE_qcengine}))
    include(FindPythonModule)
    find_python_module(qcengine ATLEAST 0.6.3 QUIET)
endif()

if(${qcengine_FOUND})
    message(STATUS "${Cyan}Found qcengine${ColourReset}: ${PY_qcengine} (found version ${qcengine_VERSION})")
    add_library(qcengine_external INTERFACE)  # dummy

else()
    if(${CMAKE_INSIST_FIND_PACKAGE_qcengine})
        message(FATAL_ERROR "Suitable qcengine could not be externally located as user insists")
    endif()

    include(ExternalProject)
    message(STATUS "Suitable qcengine could not be located, ${Magenta}Building qcengine${ColourReset} instead.")

    file(TO_NATIVE_PATH "${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}" _install_lib)

    ExternalProject_Add(qcengine_external
        DEPENDS qcelemental_external
        BUILD_ALWAYS 1
        URL https://github.com/MolSSI/QCEngine/archive/v0.6.3.tar.gz
        CONFIGURE_COMMAND ""
        UPDATE_COMMAND ""
        BUILD_COMMAND ${PYTHON_EXECUTABLE} setup.py build
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND ${PYTHON_EXECUTABLE} setup.py install
                                                      --record=record.txt
                                                      --single-version-externally-managed
                                                      --install-lib=${_install_lib})
endif()
