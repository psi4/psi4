if(NOT (${CMAKE_DISABLE_FIND_PACKAGE_qcmanybody}))
    include(FindPythonModule)
    find_python_module(qcmanybody ATLEAST 0.2.1 QUIET)
endif()

if(${qcmanybody_FOUND})
    message(STATUS "${Cyan}Found qcmanybody${ColourReset}: ${PY_qcmanybody} (found version ${qcmanybody_VERSION})")
    add_library(qcmanybody_external INTERFACE)  # dummy

    if(${CMAKE_INSIST_FIND_PACKAGE_qcmanybody})
        message(VERBOSE "Suitable qcmanybody located externally as user insists")
    endif()

else()
    if(${CMAKE_INSIST_FIND_PACKAGE_qcmanybody})
        message(FATAL_ERROR "Suitable qcmanybody could not be externally located as user insists")
    endif()

    include(ExternalProject)
    message(STATUS "Suitable qcmanybody could not be located, ${Magenta}Building qcmanybody${ColourReset} instead.")

    file(TO_NATIVE_PATH "${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}" _install_lib)

    # Must have specialty egg to compute version correctly, so the URL is only grabbing a dummy file. It's
    #   the pip command that's really selecting the source.
    ExternalProject_Add(qcmanybody_external
        DEPENDS qcelemental_external
                qcengine_external
        BUILD_ALWAYS 1
        #URL https://github.com/MolSSI/QCManyBody/archive/main.tar.gz
        URL https://github.com/loriab/QCManyBody/archive/modeltweaks.tar.gz
        DOWNLOAD_NO_EXTRACT 1
        CONFIGURE_COMMAND ""
        UPDATE_COMMAND ""
        BUILD_COMMAND ""
        BUILD_IN_SOURCE 1
        # master: INSTALL_COMMAND ${Python_EXECUTABLE} -m pip install git+https://github.com/MolSSI/QCManyBody.git#egg=proj
        INSTALL_COMMAND ${Python_EXECUTABLE} -m pip install git+https://github.com/MolSSI/QCManyBody.git@v0.2.1#egg=proj
                                                      --target=${_install_lib}
                                                      --no-deps
                                                      #-v
        )
endif()
