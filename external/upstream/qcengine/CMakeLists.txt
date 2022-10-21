if(NOT (${CMAKE_DISABLE_FIND_PACKAGE_qcengine}))
    include(FindPythonModule)
    find_python_module(qcengine ATLEAST 0.24.2 QUIET)
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
    file(TO_NATIVE_PATH "${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}" _install_bin)

    ExternalProject_Add(qcengine_external
        DEPENDS qcelemental_external
        BUILD_ALWAYS 1
        # TODO-dfocc2: reset URL after qcng v0.24.2
        # URL https://github.com/MolSSI/QCEngine/archive/v0.24.2.tar.gz
        URL https://github.com/loriab/QCEngine/archive/stdsuite_refs_for_dfocc_5.tar.gz
        CONFIGURE_COMMAND ""
        UPDATE_COMMAND ""
        BUILD_COMMAND ${Python_EXECUTABLE} setup.py build
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND ${Python_EXECUTABLE} setup.py install
                                                      --record=record.txt
                                                      --single-version-externally-managed
                                                      --install-scripts=${_install_bin}
                                                      --install-lib=${_install_lib})
endif()
