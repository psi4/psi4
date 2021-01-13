# Format borrowed from sns-mp2 CMakeLists

if(${ENABLE_psi4fockci})
    if(NOT (${CMAKE_DISABLE_FIND_PACKAGE_psi4fockci}))
        include(FindPythonModule)
        find_python_module(psi4fockci ATLEAST 0.1.0 QUIET)
    endif()

    if(${psi4fockci_FOUND})
        message(STATUS "${Cyan}Found Psi4FockCI${ColourReset}: ${PY_psi4fockci} (found version ${psi4fockci_VERSION})")
        add_library(psi4fockci_external INTERFACE)  # dummy

    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_psi4fockci})
            message(FATAL_ERROR "Suitable Psi4FockCI could not be externally located as user insists")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable Psi4FockCI could not be located, ${Magenta}Building Psi4FockCI${ColourReset} instead.")

        ExternalProject_Add(psi4fockci_external
            BUILD_ALWAYS 1
            URL https://github.com/shannonhouck/psi4fockci/archive/267e171.tar.gz # v0.1.0
            CONFIGURE_COMMAND ""
            UPDATE_COMMAND ""
            BUILD_COMMAND ${Python_EXECUTABLE} setup.py build
            BUILD_IN_SOURCE 1
            INSTALL_COMMAND ${Python_EXECUTABLE} setup.py install
                                                          --record=record.txt
                                                          --install-lib=${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR})
    endif()
else()
    add_library(psi4fockci_external INTERFACE)  # dummy
endif()
