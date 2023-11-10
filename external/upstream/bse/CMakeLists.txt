if(${ENABLE_bse})
    if(NOT (${CMAKE_DISABLE_FIND_PACKAGE_bse}))
        include(FindPythonModule)
        find_python_module(basis_set_exchange ATLEAST 0.9.0 QUIET)
    endif()
    
    if(${basis_set_exchange_FOUND})
        message(STATUS "${Cyan}Found bse${ColourReset}: ${PY_basis_set_exchange} (found version ${basis_set_exchange_VERSION})")
        add_library(bse_external INTERFACE)  # dummy

        if(${CMAKE_INSIST_FIND_PACKAGE_bse})
            message(VERBOSE "Suitable bse located externally as user insists")
        endif()

    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_bse})
            message(FATAL_ERROR "Suitable bse could not be externally located as user insists")
        endif()
    
        include(ExternalProject)
        message(STATUS "Suitable bse could not be located, ${Magenta}Building bse${ColourReset} instead.")
    
        file(TO_NATIVE_PATH "${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}" _install_lib)
    
        ExternalProject_Add(bse_external
            BUILD_ALWAYS 1
            URL https://github.com/MolSSI-BSE/basis_set_exchange/archive/v0.9.tar.gz
            CONFIGURE_COMMAND ""
            UPDATE_COMMAND ""
            BUILD_COMMAND ${Python_EXECUTABLE} setup.py build
            BUILD_IN_SOURCE 1
            INSTALL_COMMAND ${Python_EXECUTABLE} setup.py install
                                                          --record=record.txt
                                                          --single-version-externally-managed
                                                          --install-lib=${_install_lib})
    endif()
else()
    add_library(bse_external INTERFACE)  # dummy
endif()

