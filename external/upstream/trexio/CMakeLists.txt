# TrexIO Python bindings, for psi4.trexio() / trexio_from_file() / trexio_to_wavefunction().
#
# Psi4 needs only the Python module. It is built from the PyPI sdist, whose
# setup.py compiles the TREXIO C sources straight into the extension -- no
# separate libtrexio is produced or needed. (The GitHub tarballs cannot be used:
# the C sources and SWIG wrapper under python/src are generated, not tracked.)
#
# HDF5 is the reason this file is more than the usual pip recipe. psi4.trexio()
# writes the HDF5 back end by default, and trexio's setup.py enables it only if
# H5_CFLAGS/H5_LDFLAGS are set or pkg-config knows an `hdf5` module -- it does
# not fall back to h5cc. On a system without hdf5.pc (Fedora, for one) a plain
# `pip install trexio` therefore builds a TEXT-only module and says nothing, and
# the first HDF5 file then fails with a misleading "Invalid file". PyPI has no
# wheel for newer Pythons (none for 3.14 as of 2.6.1), so that is the ordinary
# outcome there, not an unusual one. Hence:
#   * a found trexio counts only if its HDF5 back end is compiled in, and
#   * when building, the HDF5 flags come from Psi4's own detection
#     (external/common/hdf5) and are passed explicitly.

if(${ENABLE_trexio})
    if(NOT (${CMAKE_DISABLE_FIND_PACKAGE_trexio}))
        include(FindPythonModule)
        find_python_module(trexio ATLEAST 2.5.0 QUIET)  # edit in codedeps
        if(${trexio_FOUND})
            # `import trexio` succeeding says nothing about HDF5, and the public
            # API has no has_backend(); the low-level binding does.
            execute_process(
                COMMAND ${Python_EXECUTABLE} -c "import sys, trexio, pytrexio.pytrexio as p; sys.exit(0 if p.trexio_has_backend(trexio.TREXIO_HDF5) else 1)"
                RESULT_VARIABLE _trexio_has_hdf5
                OUTPUT_QUIET ERROR_QUIET)
            if(NOT _trexio_has_hdf5 EQUAL 0)
                message(STATUS "Found trexio ${trexio_VERSION} at ${PY_trexio}, but it was built without the HDF5 back end; not using it")
                set(trexio_FOUND FALSE)
            endif()
        endif()
    endif()

    if(${trexio_FOUND})
        message(STATUS "${Cyan}Found trexio${ColourReset}: ${PY_trexio} (found version ${trexio_VERSION}, with HDF5)")
        add_library(trexio_external INTERFACE)  # dummy

        if(${CMAKE_INSIST_FIND_PACKAGE_trexio})
            message(VERBOSE "Suitable trexio located externally as user insists")
        endif()

    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_trexio})
            message(FATAL_ERROR "Suitable trexio (>= 2.5.0, with the HDF5 back end) could not be externally located as user insists")
        endif()

        # external/common/hdf5 has run by now (it precedes external/upstream) and
        # left an INTERFACE target whose libraries are already resolved to files,
        # which also covers HDF5 2.x's config-mode imported targets.
        if(NOT TARGET hdf5)
            message(FATAL_ERROR "Building trexio needs HDF5, but Psi4's HDF5 detection did not run")
        endif()
        get_property(_h5_incs TARGET hdf5 PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
        get_property(_h5_libs TARGET hdf5 PROPERTY INTERFACE_LINK_LIBRARIES)
        set(_h5_libdir "")
        foreach(_lib IN LISTS _h5_libs)
            get_filename_component(_name "${_lib}" NAME)
            if(_name MATCHES "^(lib)?hdf5[._]" AND NOT _name MATCHES "_hl")
                get_filename_component(_h5_libdir "${_lib}" DIRECTORY)
                break()
            endif()
        endforeach()
        if(NOT _h5_libdir)
            message(FATAL_ERROR "Could not identify the HDF5 C library among: ${_h5_libs}")
        endif()
        # HDF5 2.x in config mode can leave the include directories empty; the
        # header then sits in the library directory's sibling include/, as in
        # both conda and system layouts.
        set(_h5_incdir "")
        foreach(_inc IN LISTS _h5_incs)
            if(EXISTS "${_inc}/hdf5.h")
                set(_h5_incdir "${_inc}")
                break()
            endif()
        endforeach()
        if(NOT _h5_incdir)
            get_filename_component(_guess "${_h5_libdir}/../include" ABSOLUTE)
            if(EXISTS "${_guess}/hdf5.h")
                set(_h5_incdir "${_guess}")
            else()
                message(FATAL_ERROR "Could not find hdf5.h (looked in: ${_h5_incs}; ${_guess})")
            endif()
        endif()

        include(ExternalProject)
        message(STATUS "Suitable trexio could not be located, ${Magenta}Building trexio${ColourReset} instead (HDF5 from ${_h5_libdir}).")

        file(TO_NATIVE_PATH "${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}" _install_lib)

        # The PyPI URL is content-addressed: its path is the file's BLAKE2b-256
        # digest, so it pins the sdist without a separate URL_HASH.
        #
        # setup.py keeps only the FIRST token of each variable and adds -lhdf5
        # itself, so exactly one -I and one -L. The rpath is ours to add: setup.py
        # sets none, and an HDF5 outside the system paths (conda, modules) would
        # otherwise link but not load.
        ExternalProject_Add(trexio_external
            BUILD_ALWAYS 1
            URL https://files.pythonhosted.org/packages/65/e1/9deb3a2f66d322e0ff3c850463798c5c9197d37a2914e612701b40d7b030/trexio-2.6.1.tar.gz  # edit in codedeps
            DOWNLOAD_NO_EXTRACT 1
            CONFIGURE_COMMAND ""
            UPDATE_COMMAND ""
            BUILD_COMMAND ""
            BUILD_IN_SOURCE 1
            INSTALL_COMMAND ${CMAKE_COMMAND} -E env
                                "H5_CFLAGS=-I${_h5_incdir}"
                                "H5_LDFLAGS=-L${_h5_libdir}"
                                "LDFLAGS=-Wl,-rpath,${_h5_libdir}"
                            ${Python_EXECUTABLE} -m pip install <DOWNLOADED_FILE>
                                                      --target=${_install_lib}
                                                      --no-deps
        )
    endif()
else()
    message(STATUS "Optional addon ${Blue}unsought: trexio${ColourReset} (deferrable to runtime)")
    add_library(trexio_external INTERFACE)  # dummy
endif()
