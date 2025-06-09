if(NOT (${CMAKE_DISABLE_FIND_PACKAGE_qcengine}))
    include(FindPythonModule)
    find_python_module(qcengine ATLEAST 0.28.1 QUIET)
endif()

if(${qcengine_FOUND})
    message(STATUS "${Cyan}Found qcengine${ColourReset}: ${PY_qcengine} (found version ${qcengine_VERSION})")
    add_library(qcengine_external INTERFACE)  # dummy

    if(${CMAKE_INSIST_FIND_PACKAGE_qcengine})
        message(VERBOSE "Suitable qcengine located externally as user insists")
    endif()

else()
    if(${CMAKE_INSIST_FIND_PACKAGE_qcengine})
        message(FATAL_ERROR "Suitable qcengine could not be externally located as user insists")
    endif()

    include(ExternalProject)
    message(STATUS "Suitable qcengine could not be located, ${Magenta}Building qcengine${ColourReset} instead.")

    file(TO_NATIVE_PATH "${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}" _install_lib)
    file(TO_NATIVE_PATH "${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/" _install_bin)
    file(TO_NATIVE_PATH "${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/bin/" _pip_install_bin)
    if(WIN32)
        set(_qcengineexe "qcengine.exe")
    else()
        set(_qcengineexe "qcengine")
    endif()

    # Note: Without downgrading setuptools or `sed`ing "git_refnames" (_version.py for versioneer-generated),
    #   pip will strenuously resist ("Invalid version: <branch or commit>") arbitrary branch or commit refs from GH.
    #   Make a proper PEP440 tag (e.g., `git tag -a v0.26.0.dev1 -m "v0.26.0.dev1"`), push it, and use in URL below.
    #
    #   With downgrading to (setuptools"<66.0.0" packaging"<23.0" ; former vendors latter), one can use arbitrary
    #   branch or commit tarballs without the tagging procedure described above.

    ExternalProject_Add(qcengine_external
        DEPENDS qcelemental_external
        BUILD_ALWAYS 1
        URL https://github.com/MolSSI/QCEngine/archive/v0.32.0.tar.gz
        DOWNLOAD_NO_EXTRACT 1
        CONFIGURE_COMMAND ""
        UPDATE_COMMAND ""
        BUILD_COMMAND ""
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND ${Python_EXECUTABLE} -m pip install <DOWNLOADED_FILE>
                                                      --target=${_install_lib}
                                                      --no-deps
                                                      --upgrade
                                                      #-v
        )
    ExternalProject_Add_Step(qcengine_external
        fix_qcng_bin
        DEPENDEES install
        COMMAND ${CMAKE_COMMAND} -E copy "${_pip_install_bin}${_qcengineexe}" "${_install_bin}${_qcengineexe}"
        COMMENT "Moves qcengine entry-point into place alongside bin/psi4 for `qcengine info`"
        )
endif()
