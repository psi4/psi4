if(NOT (${CMAKE_DISABLE_FIND_PACKAGE_qcengine}))
    include(FindPythonModule)
    find_python_module(qcengine ATLEAST 0.28.1 QUIET)
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

    # Note: Without downgrading setuptools or `sed`ing "git_refnames" (_version.py for versioneer-generated),
    #   pip will strenuously resist ("Invalid version: <branch or commit>") arbitrary branch or commit refs from GH.
    #   Make a proper PEP440 tag (e.g., `git tag -a v0.26.0.dev1 -m "v0.26.0.dev1"`), push it, and use in URL below.

    ExternalProject_Add(qcengine_external
        DEPENDS qcelemental_external
        BUILD_ALWAYS 1
        #URL https://github.com/MolSSI/QCEngine/archive/v0.27.0.tar.gz
        URL https://github.com/loriab/QCEngine/archive/v0.28.1.dev2.tar.gz
        DOWNLOAD_NO_EXTRACT 1
        CONFIGURE_COMMAND ""
        UPDATE_COMMAND ""
        BUILD_COMMAND ""
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND ${Python_EXECUTABLE} -m pip install <DOWNLOADED_FILE>
                                                      --target=${_install_lib}
                                                      --no-deps
                                                      #-v
        )
endif()
