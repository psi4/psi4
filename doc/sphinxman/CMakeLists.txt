# Settle docs-specific dependencies
find_package(Perl QUIET)
if(NOT PERL_FOUND)
   message(STATUS "No Perl, no docs. Pre-built documentation at http://psicode.org/psi4manual/master/index.html")
endif()

find_package(Sphinx QUIET)
if(NOT SPHINX_FOUND)
   message(STATUS "No Sphinx, no docs. Pre-built documentation at http://psicode.org/psi4manual/master/index.html")
endif()

find_package(LATEX QUIET)
if(NOT (LATEX_COMPILER AND PDFLATEX_COMPILER))
   message(STATUS "No LaTeX (incl. pdflatex), no PDF docs. Pre-built documentation at http://psicode.org/psi4manual/master/index.html")
endif()


# Setup and build Sphinx documentation
if(PERL_FOUND AND SPHINX_FOUND)
    set(CCSD "${CMAKE_CURRENT_SOURCE_DIR}")
    set(CCBD "${CMAKE_CURRENT_BINARY_DIR}")

    # Configure Sphinx build
    configure_file(
        "${CCSD}/source/conf.py.in"
        "${CCBD}/source/conf.py"
        @ONLY)

    # Static-doc reST files
    set(STATICDOC abbr_accents.rst adc.rst appendices.rst basissets.rst basissets_byelement.rst bibliography.rst cbs.rst cbs_eqn.rst contributing.rst nbody.rst customizing.rst db.rst detci.rst dfmp2.rst diatomic.rst documentation.rst energy.rst external.rst fchk.rst freq.rst gdma.rst glossary_psivariables.rst index.rst intercalls.rst interfacing.rst introduction.rst methods.rst notes_c.rst notes_py.rst oeprop.rst opt.rst optking.rst plugins.rst programming.rst prop.rst psithonfunc.rst psithoninput.rst quickaddalias.rst basissets_tables.rst quickadddatabase.rst sapt.rst sowreap.rst testsuite.rst tutorial.rst dft_byfunctional.rst scf.rst mrcc.rst psimrcc.rst cc.rst dft.rst fnocc.rst optionshandling.rst basissets_byfamily.rst psipep.rst pep0001.rst dftd3.rst occ.rst molden.rst cubeprop.rst dcft.rst relativistic.rst dkh.rst proc_py.rst psithonmol.rst cfour.rst mrcc_table_energy.rst cfour_table_energy.rst cfour_table_grad.rst libefp.rst conda.rst pcmsolver.rst thermo.rst chemps2.rst fisapt.rst plugin_v2rdm_casscf.rst psi4api.rst addon.rst numpy.rst build_planning.rst build_faq.rst build_obtaining.rst libint.rst erd.rst simint.rst gcp.rst index_tutorials.rst)
    # * compute relative path btwn top_srcdir and objdir/doc/sphinxman
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
        "import sys; import os; print('/' + os.path.relpath('${psi4_SOURCE_DIR}', '${CCBD}') + '/../')"
        OUTPUT_VARIABLE SFNX_INCLUDE OUTPUT_STRIP_TRAILING_WHITESPACE)
    # * add ipynb doc if build env can handle it
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
        "import nbsphinx; import jupyter_client; import ipykernel"
        OUTPUT_QUIET ERROR_QUIET RESULT_VARIABLE _exec_python_err)
    if(_exec_python_err)
        message(STATUS "Some ipynb docs require extra modules (try: conda install nbsphinx jupyter_client ipykernel -c conda-forge)")
        list(APPEND STATICDOC "psiapi.rst")
        set(_jupy "")
    else()
        list(APPEND STATICDOC "psiapi.ipynb")
        set(_jupy "w/Jupyter ")
    endif()
    # * copy all over to objdir, substituting literalinclude along the way
    foreach(rst ${STATICDOC})
       configure_file(
        ${CCSD}/source/${rst}
        ${CCBD}/source/${rst}
        @ONLY)
    endforeach()
    # * duplicate index.rst as contents.rst for toc link
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CCBD}/source/index.rst ${CCBD}/source/contents.rst)

    # Static resource png and pdf files
    set(STATICRES psi4square.png favicon-psi4.ico AlphaMOs.png gdma-2.2.06.pdf Total.png VA.png prflow.001.jpeg prflow.002.jpeg)
    # * copy image and pdf files all over to objdir (configure_file corrupts file)
    foreach(fl ${STATICRES})
        execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CCSD}/source/resources/${fl} ${CCBD}/source)
    endforeach()

    # Autodoc basis sets
    add_custom_command(
        COMMAND ${PERL_EXECUTABLE} ${CCSD}/document_bases.pl ${CCSD}
        OUTPUT  ${CCBD}/source/autodoc_psi4bases.rst
        DEPENDS ${CCSD}/document_bases.pl
                ${psi4_SOURCE_DIR}/psi4/share/psi4/basis/*.gbs
        COMMENT "Autodoc basis sets")


    # Autodoc binary scratch files
    add_custom_command(
        COMMAND ${PYTHON_EXECUTABLE} ${CCSD}/document_psifiles.py ${CCSD}
        OUTPUT  ${CCBD}/source/autodoc_psifiles.rst
        DEPENDS ${CCSD}/document_psifiles.py
                ${psi4_SOURCE_DIR}/psi4/include/psi4/psifiles.h
        COMMENT "Autodoc binary scratch files")


    # Autodoc physical constants
    add_custom_command(
        COMMAND ${PERL_EXECUTABLE} ${CCSD}/document_physconst.pl ${CCSD}
        OUTPUT  ${CCBD}/source/autodoc_physconst.rst
        DEPENDS ${CCSD}/document_physconst.pl
                ${psi4_SOURCE_DIR}/psi4/include/psi4/physconst.h
        COMMENT "Autodoc physical constants")


    # Autodoc options c-side
    add_custom_command(
        COMMAND ${CMAKE_COMMAND} -E make_directory source/autodir_options_c
        COMMAND ${PERL_EXECUTABLE} ${CCSD}/document_options_c.pl ${CCSD}
        OUTPUT  ${CCBD}/source/autodoc_abbr_options_c.rst
                ${CCBD}/source/autodoc_glossary_options_c.rst
                ${CCBD}/source/autodoc_options_c_bymodule.rst
                ${CCBD}/source/autodir_options_c/module__scf.rst  # representative
                ${CCBD}/source/autodir_options_c/scf__scf_type.rst  # representative
        DEPENDS ${CCSD}/document_options_c.pl
                ${psi4_SOURCE_DIR}/psi4/src/read_options.cc
        COMMENT "Autodoc options c-side")


    # Autodoc PSI variables c-side
    add_custom_command(
        COMMAND ${CMAKE_COMMAND} -E make_directory source/autodir_psivariables
        COMMAND ${PERL_EXECUTABLE} ${CCSD}/document_psivariables.pl ${CCSD}
        OUTPUT  ${CCBD}/source/autodoc_psivariables_bymodule.rst
                ${CCBD}/source/autodir_psivariables/module__scf.rst  # representative
        DEPENDS ${CCSD}/document_psivariables.pl
                ${psi4_SOURCE_DIR}/psi4/src/*/*/*.cc
        COMMENT "Autodoc PSI variables c-side")


    # Autodoc Cfour PSI variables
    add_custom_command(
        COMMAND ${CMAKE_COMMAND} -E make_directory source/autodir_psivariables
        COMMAND ${PYTHON_EXECUTABLE} ${CCSD}/document_cfour.py ${CCSD}
        OUTPUT  ${CCBD}/source/autodir_psivariables/module__cfour.rst
        DEPENDS ${CCSD}/document_cfour.py
                ${psi4_SOURCE_DIR}/psi4/driver/qcdb/cfour.py
                ${CCBD}/source/autodoc_psivariables_bymodule.rst
        COMMENT "Autodoc Cfour PSI variables")


    # Autodoc databases
    add_custom_command(
        COMMAND ${PYTHON_EXECUTABLE} ${CCSD}/document_databases.py ${CCSD}
        OUTPUT  ${CCBD}/source/autodoc_available_databases.rst
        DEPENDS ${CCSD}/document_databases.py
                ${psi4_SOURCE_DIR}/psi4/share/psi4/databases/*.py
        COMMENT "Autodoc databases")


    # Autodoc test cases
    add_custom_command(
        COMMAND ${PERL_EXECUTABLE} ${CCSD}/document_tests.pl ${CCSD}
        OUTPUT  ${CCBD}/source/autodoc_testsuite_corepsi4.rst
                ${CCBD}/source/autodoc_testsuite_dftd3.rst
        DEPENDS ${CCSD}/document_tests.pl
                ${psi4_SOURCE_DIR}/tests/*/input.dat
                ${psi4_SOURCE_DIR}/tests/*/*/input.dat
        COMMENT "Autodoc test cases")


#    # Autodoc driver
#    add_custom_command(
#        COMMAND ${PYTHON_EXECUTABLE} ${CCSD}/document_driver.py ${CCSD}
#        OUTPUT  ${CCBD}/source/autodoc_driver.rst
#        DEPENDS ${psi4_SOURCE_DIR}/psi4/driver/*.py
#                ${CCSD}/document_driver.py
#        COMMENT "Autodoc driver")


#    # Autodoc psi4 module
#    add_custom_command(
#        COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_BINARY_DIR} ctest -R docs-psimod
#        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/tests/docs-psimod/autodoc_psimod.rst ${CCBD}/source
#        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/tests/docs-psimod/autodoc_man_psi4.txt ${CCBD}/source
#        OUTPUT  ${CCBD}/source/autodoc_psimod.rst
#                ${CCBD}/source/autodoc_man_psi4.txt
#        DEPENDS ${psi4_SOURCE_DIR}/psi4/src/*/*/*.cc
#                ${psi4_SOURCE_DIR}/tests/docs-psimod/input.dat
#                ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/psi4
#        COMMENT "Autodoc psi4 module")


    # Autodoc functionals
    add_custom_command(
        COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_BINARY_DIR} ctest -R docs-dft
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/tests/docs-dft/autodoc_dft_energy.rst ${CCBD}/source
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/tests/docs-dft/autodoc_dft_opt.rst ${CCBD}/source
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/tests/docs-dft/autodoc_dft_all.rst ${CCBD}/source
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/tests/docs-dft/autodoc_dft_gga.rst ${CCBD}/source
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/tests/docs-dft/autodoc_dft_meta.rst ${CCBD}/source
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/tests/docs-dft/autodoc_dft_lrc.rst ${CCBD}/source
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/tests/docs-dft/autodoc_dft_hybrid.rst ${CCBD}/source
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/tests/docs-dft/autodoc_dft_dhybrid.rst ${CCBD}/source
        #COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/tests/docs-dft/autodoc_dft_scsdhybrid.rst ${CCBD}/source  # commented until populated again to appease Sphinx
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/tests/docs-dft/autodoc_dft_disp.rst ${CCBD}/source
        OUTPUT  ${CCBD}/source/autodoc_dft_energy.rst
                ${CCBD}/source/autodoc_dft_opt.rst
                ${CCBD}/source/autodoc_dft_all.rst
                ${CCBD}/source/autodoc_dft_gga.rst
                ${CCBD}/source/autodoc_dft_meta.rst
                ${CCBD}/source/autodoc_dft_lrc.rst
                ${CCBD}/source/autodoc_dft_hybrid.rst
                ${CCBD}/source/autodoc_dft_dhybrid.rst
                #${CCBD}/source/autodoc_dft_scsdhybrid.rst  # commented until populated again to appease Sphinx
                ${CCBD}/source/autodoc_dft_disp.rst
        DEPENDS ${psi4_SOURCE_DIR}/psi4/driver/procrouting/dft_functional.py
                ${psi4_SOURCE_DIR}/tests/docs-dft/input.dat
                ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/psi4
        COMMENT "Autodoc functionals")


    # Autodoc auxiliary basis sets
    add_custom_command(
        COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_BINARY_DIR} ctest -R docs-bases
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/tests/docs-bases/autodoc_basissets_byfamily.rst ${CCBD}/source
        OUTPUT ${CCBD}/source/autodoc_basissets_byfamily.rst
        DEPENDS ${psi4_SOURCE_DIR}/psi4/driver/qcdb/basislist.py
                ${psi4_SOURCE_DIR}/psi4/driver/qcdb/basislistdunning.py
                ${psi4_SOURCE_DIR}/psi4/driver/qcdb/basislistother.py
                ${psi4_SOURCE_DIR}/tests/docs-bases/input.dat
                ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/psi4
        COMMENT "Autodoc auxiliary basis sets")


    # Autodoc plugins
    #   Not at all in working order
    add_custom_command(
        COMMAND ${CMAKE_COMMAND} -E make_directory source/autodir_plugins
        COMMAND ${PYTHON_EXECUTABLE} ${CCSD}/document_plugins.py ${CCSD}
        OUTPUT  ${CCBD}/source/autodoc_abbr_options_plugins.rst
                ${CCBD}/source/autodoc_available_plugins.rst
                ${CCBD}/source/autodir_plugins/glossary__aointegrals.rst
                ${CCBD}/source/autodir_plugins/module__aointegrals.rst
                ${CCBD}/source/autodir_plugins/aointegrals__print.rst
        DEPENDS ${CCSD}/document_plugins.py
                ${psi4_SOURCE_DIR}/plugins/*/*.py
                ${psi4_SOURCE_DIR}/plugins/*/doc.rst
        COMMENT "Autodoc plugins")


    # Autodoc efp fragments
    if(${ENABLE_libefp})
        add_custom_command(
            COMMAND ${PYTHON_EXECUTABLE} ${CCSD}/document_efpfrag.py ${CCSD}
            OUTPUT  ${CCBD}/source/autodoc_available_efpfrag.rst
            DEPENDS ${CCSD}/document_efpfrag.py
                    ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/psi4/efpfrag/*efp
            COMMENT "Autodoc efp fragments")
    else()
        file(WRITE ${CCBD}/source/autodoc_available_efpfrag.rst "")
    endif()


    # Format git log for feed
    add_custom_command(
        COMMAND ${CMAKE_COMMAND} -E make_directory feed
        COMMAND ${PYTHON_EXECUTABLE} ${CCSD}/extract_changeset.py ${CCBD}
        OUTPUT  ${CCBD}/feed/history_trac_changeset.txt
                ${CCBD}/feed/latest_trac_changeset.txt
        DEPENDS ${CCSD}/extract_changeset.py
        COMMENT "Format git --log for psicode.org feed")

    # Establishing dependencies for docs targets

    # * mini github feed for psicode.org
    add_custom_target(ghfeed
        DEPENDS ${CCBD}/feed/history_trac_changeset.txt
                ${CCBD}/feed/latest_trac_changeset.txt
        COMMENT "Preparing GitHub feed")

    # * announcement
    if (LATEX_COMPILER AND PDFLATEX_COMPILER)
        message(STATUS "Documentation targets ${_jupy}available: sphinxman (html), sphinxmini (quick html), sphinxpdf (LaTeX --> PDF)")
    else()
        message(STATUS "Documentation targets ${_jupy}available: sphinxman (html), sphinxmini (quick html)")
    endif()

    # * quick mini-build for docs development
    add_custom_target(sphinxmini
        DEPENDS ${CCBD}/source/conf.py
        COMMENT "Building Sphinx abridged HTML documentation build")

    # * full, proper build
    add_custom_target(sphinxman
        DEPENDS ${CCBD}/source/autodoc_psi4bases.rst
                ${CCBD}/source/autodoc_psifiles.rst
                ${CCBD}/source/autodoc_physconst.rst
                ${CCBD}/source/autodoc_abbr_options_c.rst
                ${CCBD}/source/autodoc_glossary_options_c.rst
                ${CCBD}/source/autodoc_options_c_bymodule.rst
                ${CCBD}/source/autodir_options_c/module__scf.rst  # representative
                ${CCBD}/source/autodir_options_c/scf__scf_type.rst  # representative
                ${CCBD}/source/autodoc_psivariables_bymodule.rst
                ${CCBD}/source/autodir_psivariables/module__scf.rst  # representative
                ${CCBD}/source/autodir_psivariables/module__cfour.rst
                ${CCBD}/source/autodoc_available_databases.rst
                ${CCBD}/source/autodoc_testsuite_corepsi4.rst
                ${CCBD}/source/autodoc_testsuite_dftd3.rst
                ${CCBD}/source/autodoc_dft_energy.rst
                ${CCBD}/source/autodoc_dft_opt.rst
                ${CCBD}/source/autodoc_dft_all.rst
                ${CCBD}/source/autodoc_dft_gga.rst
                ${CCBD}/source/autodoc_dft_meta.rst
                ${CCBD}/source/autodoc_dft_lrc.rst
                ${CCBD}/source/autodoc_dft_hybrid.rst
                ${CCBD}/source/autodoc_dft_dhybrid.rst
                #${CCBD}/source/autodoc_dft_scsdhybrid.rst  # commented until populated again to appease Sphinx
                ${CCBD}/source/autodoc_dft_disp.rst
                ${CCBD}/source/autodoc_basissets_byfamily.rst
                ${CCBD}/source/autodoc_available_efpfrag.rst
        COMMENT "Preparing Sphinx HTML documentation build")
    add_dependencies(sphinxman psi4-core)

    install(
        DIRECTORY ${CCBD}/html
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DOCDIR}
        OPTIONAL)

    # * full, proper build with LaTeX
    if (LATEX_COMPILER AND PDFLATEX_COMPILER)
        add_custom_target(sphinxpdf
            DEPENDS sphinxman
            COMMENT "Preparing Sphinx HTML and PDF documentation build")

        # * full, proper build with latex builder
        add_custom_command(
            TARGET sphinxpdf POST_BUILD
            COMMAND ${SPHINX_EXECUTABLE}
                -b latex              # builder
                -d ${CCBD}/_doctrees  # cache w/ pickled reST
                -c ${CCBD}/source     # whereabouts of conf.py
                ${CCBD}/source        # build from
                ${CCBD}/latex         # build to
            COMMENT "Building latex documentation ...")
        add_custom_command(
            TARGET sphinxpdf POST_BUILD
            COMMAND ${PDFLATEX_COMPILER}
                -interaction=nonstopmode
                "Psi4.tex" > /dev/null 2>&1
            WORKING_DIRECTORY ${CCBD}/latex
            COMMENT "Building remarkably ugly PDF documentation from LaTeX ... (ignore the make exit error)")
    endif()


    # Establishing commands to build docs targets

    # * full, proper build with html builder
    add_custom_command(
        TARGET sphinxman POST_BUILD
        COMMAND ${SPHINX_EXECUTABLE}
            # -v  # Repeat for greater verbosity
            # -q  # Enable for quiet mode
            # -Q  # Enable for even more quiet mode
            # -E  # don't use a saved environment, always read all file
            # -a  # write all files; default is to only write new and changed files
            -b html               # builder
            -d ${CCBD}/_doctrees  # cache w/ pickled reST
            -c ${CCBD}/source     # whereabouts of conf.py
            ${CCBD}/source        # build from
            ${CCBD}/html          # build to
        COMMENT "Building html documentation ...")

    # * Special abridged-docs target for fast-compile when developing static docs
    #   * removing autodocs alleviates dependencies
    #   * thar will be warnings galore for unresolved links
    add_custom_command(
        TARGET sphinxmini
        COMMAND ${CMAKE_COMMAND} -E remove ${CCBD}/source/autodoc*
        COMMAND ${CMAKE_COMMAND} -E remove_directory ${CCBD}/source/autodir_options_c
        COMMAND ${CMAKE_COMMAND} -E remove_directory ${CCBD}/source/autodir_psivariables
        COMMAND ${CMAKE_COMMAND} -E remove_directory ${CCBD}/source/autodir_psivariables
        COMMAND ${CMAKE_COMMAND} -E remove_directory ${CCBD}/source/autodir_plugins
        COMMAND ${SPHINX_EXECUTABLE}
            -b html
            -d ${CCBD}/_doctrees
            -c ${CCBD}/source
            ${CCBD}/source
            ${CCBD}/html
        COMMENT "Building abridged html documentation ...")

endif(PERL_FOUND AND SPHINX_FOUND)

