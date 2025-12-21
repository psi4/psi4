if(NOT (${CMAKE_DISABLE_FIND_PACKAGE_qcelemental}))
    include(FindPythonModule)
    if(${Python_NumPy_VERSION} GREATER_EQUAL "2")
        find_python_module(qcelemental ATLEAST 0.28.0 QUIET)  # edit in codedeps
    else()
        find_python_module( qcelemental ATLEAST 0.26.0 QUIET)
    endif()
endif()

if(${qcelemental_FOUND})
    message(STATUS "${Cyan}Found qcelemental${ColourReset}: ${PY_qcelemental} (found version ${qcelemental_VERSION})")
    add_library(qcelemental_external INTERFACE)  # dummy

    if(${CMAKE_INSIST_FIND_PACKAGE_qcelemental})
        message(VERBOSE "Suitable qcelemental located externally as user insists")
    endif()

else()
    if(${CMAKE_INSIST_FIND_PACKAGE_qcelemental})
        message(FATAL_ERROR "Suitable qcelemental could not be externally located as user insists")
    endif()

    include(ExternalProject)
    message(STATUS "Suitable qcelemental could not be located, ${Magenta}Building qcelemental${ColourReset} instead.")

    file(TO_NATIVE_PATH "${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}" _install_lib)

    ExternalProject_Add(qcelemental_external
        BUILD_ALWAYS 1
        DOWNLOAD_COMMAND ""
        CONFIGURE_COMMAND ""
        UPDATE_COMMAND ""
        BUILD_COMMAND ""
        BUILD_IN_SOURCE 1
        # master: INSTALL_COMMAND ${Python_EXECUTABLE} -m pip install git+https://github.com/MolSSI/ QCElemental.git#egg=proj
        #INSTALL_COMMAND ${Python_EXECUTABLE} -m pip install git+https://github.com/MolSSI/QCElemental.git@v0.50.0#egg=proj  # edit in codedeps
        INSTALL_COMMAND ${Python_EXECUTABLE} -m pip install git+https://github.com/loriab/QCElemental.git@py314#egg=proj  # edit in codedeps
                                                      --target=${_install_lib}
                                                      --no-deps
                                                      #-v
        )
endif()
