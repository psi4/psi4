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

if(SPHINX_THEME STREQUAL sphinx_psi_theme)
    execute_process(COMMAND ${Python_EXECUTABLE} -c
        "import sphinx_psi_theme; import cloud_sptheme; import sphinx_automodapi; import graphviz"
        OUTPUT_QUIET ERROR_QUIET RESULT_VARIABLE _exec_python_err)
elseif(SPHINX_THEME STREQUAL basic)
    # these themes are built-in to Sphinx
    execute_process(COMMAND ${Python_EXECUTABLE} -c
        "import sphinx_automodapi; import graphviz"
        OUTPUT_QUIET ERROR_QUIET RESULT_VARIABLE _exec_python_err)
else()
    execute_process(COMMAND ${Python_EXECUTABLE} -c
        "import ${SPHINX_THEME}; import sphinx_automodapi; import graphviz"
        OUTPUT_QUIET ERROR_QUIET RESULT_VARIABLE _exec_python_err)
endif()
if(_exec_python_err)
    message(STATUS "No theme, no docs (try: conda install sphinx-psi-theme cloud_sptheme sphinx-automodapi sphinx-autodoc-typehints graphviz python-graphviz autodoc-pydantic sphinx-tabs -c conda-forge -c psi4)")
    set(SPHINX_STUFF_FOUND "OFF")

else()
    set(SPHINX_STUFF_FOUND "ON")
endif()

if(SPHINXMAN_STRICT)
    set(SPHINX_EXTRAS "-n;-W;--keep-going")
    #set(SPHINX_EXTRAS "-W;--keep-going")
    # return to above after qcel 0.25
    #set(SPHINX_EXTRAS "")
else()
    set(SPHINX_EXTRAS "")
endif()

# Setup and build Sphinx documentation
if(PERL_FOUND AND SPHINX_FOUND AND SPHINX_STUFF_FOUND)
    set(CCSD "${CMAKE_CURRENT_SOURCE_DIR}")
    set(CCBD "${CMAKE_CURRENT_BINARY_DIR}")

    # Static-doc reST files
    set(STATICDOC abbr_accents.rst adc.rst adcc.rst appendices.rst basissets.rst
    basissets_byelement.rst bibliography.rst cbs.rst cbs_eqn.rst ddx.rst
    nbody.rst customizing.rst db.rst detci.rst dfmp2.rst dlpnomp2.rst
    diatomic.rst documentation.rst energy.rst external.rst fchk.rst
    freq.rst gdma.rst glossary_psivariables.rst index.rst intercalls.rst
    interfacing.rst introduction.rst methods.rst notes_c.rst notes_py.rst
    oeprop.rst opt.rst optking.rst plugins.rst prog_corelibs.rst prop.rst
    psithonfunc.rst psithoninput.rst quickaddalias.rst mp2f12.rst
    basissets_tables.rst quickadddatabase.rst sapt.rst sowreap.rst
    testsuite.rst tutorial.rst dft_byfunctional.rst scf.rst mrcc.rst
    psimrcc.rst cc.rst dft.rst fnocc.rst optionshandling.rst
    basissets_byfamily.rst psipep.rst pep0001.rst dftd3.rst occ.rst
    molden.rst cubeprop.rst dct.rst relativistic.rst dkh.rst proc_py.rst
    psithonmol.rst cfour.rst mrcc_table_energy.rst cfour_table_energy.rst
    cfour_table_grad.rst libefp.rst cppe.rst conda.rst pcmsolver.rst thermo.rst
    chemps2.rst fisapt.rst plugin_v2rdm_casscf.rst psi4api.rst
    manage_addon.rst numpy.rst build_planning.rst build_faq.rst
    build_obtaining.rst libint.rst simint.rst gcp.rst
    index_tutorials.rst prog_faq.rst manage_index.rst manage_git.rst
    prog_style_c.rst prog_blas.rst add_tests.rst plugin_snsmp2.rst dftnl.rst
    libxc.rst gau2grid.rst prog_conduct.rst prog_intro.rst prog_newcode.rst
    prog_setup.rst prog_welcome.rst prog_ways_to_add.rst
    prog_full_integration.rst prog_style_py.rst prog_help.rst
    prog_debug_profile.rst manage_release.rst plugin_psi4fockci.rst
    brianqc.rst tdscf.rst prog_integrals.rst prog_tour.rst mcscf.rst
    code_style.rst nitpick-exceptions external_apis.rst ecpint.rst
    capabilities.rst variables.rst
    )

    # * compute relative path btwn top_srcdir and objdir/doc/sphinxman
    execute_process(COMMAND ${Python_EXECUTABLE} -c
        "import sys; import os; print('/' + os.path.relpath('${psi4_SOURCE_DIR}', '${CCBD}') + '/../')"
        OUTPUT_VARIABLE SFNX_INCLUDE OUTPUT_STRIP_TRAILING_WHITESPACE)
    # * add ipynb doc if build env can handle it
    execute_process(COMMAND ${Python_EXECUTABLE} -c
        "import nbsphinx; import jupyter_client; import ipykernel"
        OUTPUT_QUIET ERROR_QUIET RESULT_VARIABLE _exec_python_err)
    if(_exec_python_err)
        message(STATUS "Some ipynb docs require extra modules (try: conda install nbsphinx jupyter_client ipykernel -c conda-forge)")
        list(APPEND STATICDOC "psiapi.rst")
        set(_jupy "")
        set(_jupyconf "")
    else()
        list(APPEND STATICDOC "psiapi.ipynb")
        set(_jupy "w/Jupyter ")
        set(_jupyconf "'nbsphinx',")
    endif()
    # configure Sphinx build
    configure_file(
        "${CCSD}/source/conf.py.in"
        "${CCBD}/source/conf.py"
        @ONLY)
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
    set(STATICRES
      psi4square.png
      favicon-psi4.ico
      AlphaMOs.png
      gdma-2.3.pdf
      Total.png
      VA.png
      prflow.001.jpeg
      prflow.002.jpeg
      prflow.003.jpeg
      moxy.png
      )
    # * copy image and pdf files all over to objdir (configure_file corrupts file)
    foreach(fl ${STATICRES})
        execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CCSD}/source/resources/${fl} ${CCBD}/source)
    endforeach()

    # Extension files
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CCSD}/source/extensions ${CCBD}/source/extensions)

    # Autodoc basis sets
    add_custom_command(
        COMMAND ${PERL_EXECUTABLE} ${CCSD}/document_bases.pl ${CCSD}
        OUTPUT  ${CCBD}/source/autodoc_psi4bases.rst
        DEPENDS ${CCSD}/document_bases.pl
                $<IF:$<STREQUAL:${CMAKE_GENERATOR},"Unix Makefiles">,${psi4_SOURCE_DIR}/psi4/share/psi4/basis/*.gbs,${psi4_SOURCE_DIR}/psi4/share/psi4/basis/cc-pvdz.gbs>
        COMMENT "Autodoc basis sets")


    # Autodoc binary scratch files
    add_custom_command(
        COMMAND ${Python_EXECUTABLE} ${CCSD}/document_psifiles.py ${CCSD}
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
                ${CCBD}/source/autodir_options_c/scf__guess.rst  # representative
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
                $<IF:$<STREQUAL:${CMAKE_GENERATOR},"Unix Makefiles">,${psi4_SOURCE_DIR}/psi4/src/*/*/*.cc,${psi4_SOURCE_DIR}/psi4/src/psi4/fnocc/ccsd.cc>
        COMMENT "Autodoc PSI variables c-side")


    # Autodoc Cfour PSI variables
    add_custom_command(
        COMMAND ${CMAKE_COMMAND} -E make_directory source/autodir_psivariables
        COMMAND ${Python_EXECUTABLE} ${CCSD}/document_cfour.py ${CCSD}
        OUTPUT  ${CCBD}/source/autodir_psivariables/module__cfour.rst
        DEPENDS ${CCSD}/document_cfour.py
                ${psi4_SOURCE_DIR}/psi4/driver/qcdb/cfour.py
                ${CCBD}/source/autodoc_psivariables_bymodule.rst
        COMMENT "Autodoc Cfour PSI variables")


    # Autodoc databases
    add_custom_command(
        COMMAND ${Python_EXECUTABLE} ${CCSD}/document_databases.py ${CCSD}
        OUTPUT  ${CCBD}/source/autodoc_available_databases.rst
        DEPENDS ${CCSD}/document_databases.py
                $<IF:$<STREQUAL:${CMAKE_GENERATOR},"Unix Makefiles">,${psi4_SOURCE_DIR}/psi4/share/psi4/databases/*.py,${psi4_SOURCE_DIR}/psi4/share/psi4/databases/S22.py>
        COMMENT "Autodoc databases")


    # Autodoc test cases
    add_custom_command(
        COMMAND ${PERL_EXECUTABLE} ${CCSD}/document_tests.pl ${CCSD}
        OUTPUT  ${CCBD}/source/autodoc_testsuite_corepsi4.rst
                ${CCBD}/source/autodoc_testsuite_dftd3.rst
        DEPENDS ${CCSD}/document_tests.pl
                $<IF:$<STREQUAL:${CMAKE_GENERATOR},"Unix Makefiles">,${psi4_SOURCE_DIR}/tests/*/input.dat,${psi4_SOURCE_DIR}/tests/ghosts/input.dat>
                $<IF:$<STREQUAL:${CMAKE_GENERATOR},"Unix Makefiles">,${psi4_SOURCE_DIR}/tests/*/*/input.dat,${psi4_SOURCE_DIR}/tests/dftd3/energy/input.dat>
        COMMENT "Autodoc test cases")


    # Autodoc capabilities
    # * tables are generated in pairs: one to build directory for Sphinx to use and one for the repo for devs to track
    # * PYTHONPATH bit is because this is the only document_* script that uses a fn from Psi4
    add_custom_command(
        COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CCBD}/../../stage/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}:$PYTHONPATH ${Python_EXECUTABLE} ${CCSD}/document_capabilities.py --stdsuite ${CCSD}/../../samples/stdsuite_psi4.txt --sphinx --quiet --writefile         "source/autodoc_capabilities_details.rst"   --mode details
        COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CCBD}/../../stage/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}:$PYTHONPATH ${Python_EXECUTABLE} ${CCSD}/document_capabilities.py --stdsuite ${CCSD}/../../samples/stdsuite_psi4.txt          --quiet --writefile "${CCSD}/source/preview_capabilities_details.rst"   --mode details
        COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CCBD}/../../stage/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}:$PYTHONPATH ${Python_EXECUTABLE} ${CCSD}/document_capabilities.py --stdsuite ${CCSD}/../../samples/stdsuite_psi4.txt --sphinx --quiet --writefile         "source/autodoc_capabilities_summary.rst"   --mode summary
        COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CCBD}/../../stage/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}:$PYTHONPATH ${Python_EXECUTABLE} ${CCSD}/document_capabilities.py --stdsuite ${CCSD}/../../samples/stdsuite_psi4.txt          --quiet --writefile "${CCSD}/source/preview_capabilities_summary.rst"   --mode summary
        COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CCBD}/../../stage/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}:$PYTHONPATH ${Python_EXECUTABLE} ${CCSD}/document_capabilities.py --stdsuite ${CCSD}/../../samples/stdsuite_psi4.txt --sphinx --quiet --writefile         "source/autodoc_capabilities_ccenergy.rst"  --mode ccenergy
        COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CCBD}/../../stage/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}:$PYTHONPATH ${Python_EXECUTABLE} ${CCSD}/document_capabilities.py --stdsuite ${CCSD}/../../samples/stdsuite_psi4.txt          --quiet --writefile "${CCSD}/source/preview_capabilities_ccenergy.rst"  --mode ccenergy
        COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CCBD}/../../stage/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}:$PYTHONPATH ${Python_EXECUTABLE} ${CCSD}/document_capabilities.py --stdsuite ${CCSD}/../../samples/stdsuite_psi4.txt --sphinx --quiet --writefile         "source/autodoc_capabilities_fnocc.rst"     --mode fnocc --driver e
        COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CCBD}/../../stage/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}:$PYTHONPATH ${Python_EXECUTABLE} ${CCSD}/document_capabilities.py --stdsuite ${CCSD}/../../samples/stdsuite_psi4.txt          --quiet --writefile "${CCSD}/source/preview_capabilities_fnocc.rst"     --mode fnocc --driver e
        COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CCBD}/../../stage/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}:$PYTHONPATH ${Python_EXECUTABLE} ${CCSD}/document_capabilities.py --stdsuite ${CCSD}/../../samples/stdsuite_psi4.txt --sphinx --quiet --writefile         "source/autodoc_capabilities_dfmp2.rst"     --mode dfmp2
        COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CCBD}/../../stage/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}:$PYTHONPATH ${Python_EXECUTABLE} ${CCSD}/document_capabilities.py --stdsuite ${CCSD}/../../samples/stdsuite_psi4.txt          --quiet --writefile "${CCSD}/source/preview_capabilities_dfmp2.rst"     --mode dfmp2
        COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CCBD}/../../stage/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}:$PYTHONPATH ${Python_EXECUTABLE} ${CCSD}/document_capabilities.py --stdsuite ${CCSD}/../../samples/stdsuite_psi4.txt --sphinx --quiet --writefile         "source/autodoc_capabilities_occ_oo.rst"    --mode occ_oo
        COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CCBD}/../../stage/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}:$PYTHONPATH ${Python_EXECUTABLE} ${CCSD}/document_capabilities.py --stdsuite ${CCSD}/../../samples/stdsuite_psi4.txt          --quiet --writefile "${CCSD}/source/preview_capabilities_occ_oo.rst"    --mode occ_oo
        COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CCBD}/../../stage/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}:$PYTHONPATH ${Python_EXECUTABLE} ${CCSD}/document_capabilities.py --stdsuite ${CCSD}/../../samples/stdsuite_psi4.txt --sphinx --quiet --writefile         "source/autodoc_capabilities_occ_nonoo.rst" --mode occ_nonoo --second-in-rst
        COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CCBD}/../../stage/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}:$PYTHONPATH ${Python_EXECUTABLE} ${CCSD}/document_capabilities.py --stdsuite ${CCSD}/../../samples/stdsuite_psi4.txt          --quiet --writefile "${CCSD}/source/preview_capabilities_occ_nonoo.rst" --mode occ_nonoo --second-in-rst
        COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CCBD}/../../stage/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}:$PYTHONPATH ${Python_EXECUTABLE} ${CCSD}/document_capabilities.py --stdsuite ${CCSD}/../../samples/stdsuite_psi4.txt --sphinx --quiet --writefile         "source/autodoc_capabilities_scf.rst"       --mode scf --driver egh
        COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CCBD}/../../stage/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}:$PYTHONPATH ${Python_EXECUTABLE} ${CCSD}/document_capabilities.py --stdsuite ${CCSD}/../../samples/stdsuite_psi4.txt          --quiet --writefile "${CCSD}/source/preview_capabilities_scf.rst"       --mode scf --driver egh
        OUTPUT  ${CCBD}/source/autodoc_capabilities_details.rst
                ${CCBD}/source/autodoc_capabilities_summary.rst
                ${CCBD}/source/autodoc_capabilities_ccenergy.rst
                ${CCBD}/source/autodoc_capabilities_fnocc.rst
                ${CCBD}/source/autodoc_capabilities_dfmp2.rst
                ${CCBD}/source/autodoc_capabilities_occ_oo.rst
                ${CCBD}/source/autodoc_capabilities_occ_nonoo.rst
                ${CCBD}/source/autodoc_capabilities_scf.rst
        DEPENDS ${CCSD}/document_capabilities.py
                ${CCSD}/../../samples/stdsuite_psi4.txt
        COMMENT "Autodoc capabilities")

    # Autodoc stubs
    add_custom_command(
        COMMAND ${Python_EXECUTABLE} ${CCSD}/document_stubs.py ${CCSD}
        OUTPUT  ${CCBD}/source/api/psi4.core.set_variable.rst
        DEPENDS ${CCSD}/document_stubs.py
        COMMENT "Autodoc stubs")


#    # Autodoc driver
#    add_custom_command(
#        COMMAND ${Python_EXECUTABLE} ${CCSD}/document_driver.py ${CCSD}
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
        DEPENDS ${psi4_SOURCE_DIR}/psi4/driver/procrouting/dft/superfunctionals.py
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
    #add_custom_command(
    #    COMMAND ${CMAKE_COMMAND} -E make_directory source/autodir_plugins
    #    COMMAND ${Python_EXECUTABLE} ${CCSD}/document_plugins.py ${CCSD}
    #    OUTPUT  ${CCBD}/source/autodoc_abbr_options_plugins.rst
    #            ${CCBD}/source/autodoc_available_plugins.rst
    #            ${CCBD}/source/autodir_plugins/glossary__aointegrals.rst
    #            ${CCBD}/source/autodir_plugins/module__aointegrals.rst
    #            ${CCBD}/source/autodir_plugins/aointegrals__print.rst
    #    DEPENDS ${CCSD}/document_plugins.py
    #            ${psi4_SOURCE_DIR}/plugins/*/*.py
    #            ${psi4_SOURCE_DIR}/plugins/*/doc.rst
    #    COMMENT "Autodoc plugins")

    # Autodoc efp fragments
    if(${ENABLE_libefp})
        add_custom_command(
            COMMAND ${Python_EXECUTABLE} ${CCSD}/document_efpfrag.py ${CCSD}
            OUTPUT  ${CCBD}/source/autodoc_available_efpfrag.rst
            DEPENDS ${CCSD}/document_efpfrag.py
                    $<IF:$<STREQUAL:${CMAKE_GENERATOR},"Unix Makefiles">,${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/psi4/efpfrag/*efp,${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/psi4/efpfrag/h2o.efp>
            COMMENT "Autodoc efp fragments")
    else()
        file(WRITE ${CCBD}/source/autodoc_available_efpfrag.rst "")
    endif()


    # Format git log for feed
    add_custom_command(
        COMMAND ${CMAKE_COMMAND} -E make_directory feed
        COMMAND ${Python_EXECUTABLE} ${CCSD}/extract_changeset.py ${CCBD}
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
                ${CCBD}/source/autodir_options_c/scf__guess.rst  # representative
                ${CCBD}/source/autodoc_psivariables_bymodule.rst
                ${CCBD}/source/autodir_psivariables/module__scf.rst  # representative
                ${CCBD}/source/autodir_psivariables/module__cfour.rst
                ${CCBD}/source/autodoc_available_databases.rst
                ${CCBD}/source/autodoc_testsuite_corepsi4.rst
                ${CCBD}/source/autodoc_testsuite_dftd3.rst
                ${CCBD}/source/autodoc_capabilities_details.rst
                ${CCBD}/source/autodoc_capabilities_summary.rst
                ${CCBD}/source/autodoc_capabilities_ccenergy.rst
                ${CCBD}/source/autodoc_capabilities_fnocc.rst
                ${CCBD}/source/autodoc_capabilities_dfmp2.rst
                ${CCBD}/source/autodoc_capabilities_occ_oo.rst
                ${CCBD}/source/autodoc_capabilities_occ_nonoo.rst
                ${CCBD}/source/autodoc_capabilities_scf.rst
                ${CCBD}/source/api/psi4.core.set_variable.rst
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

    add_custom_target(linkcheck
        DEPENDS sphinxman
        COMMENT "Preparing Sphinx linkcheck documentation build")

    add_custom_command(
        TARGET linkcheck POST_BUILD
        COMMAND ${SPHINX_EXECUTABLE}
            -b linkcheck          # builder
            -d ${CCBD}/_doctrees  # cache w/ pickled reST
            -c ${CCBD}/source     # whereabouts of conf.py
            ${CCBD}/source        # build from
            ${CCBD}/latex         # build to
        COMMENT "Building linkcheck documentation ...")

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
            # -n                  # nit-picky mode, warn about all missing references
            # -W                  # turn warnings into errors
            # --keep-going        # with -W, keep going when getting warnings
            ${SPHINX_EXTRAS}
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
        TARGET sphinxmini POST_BUILD
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

endif(PERL_FOUND AND SPHINX_FOUND AND SPHINX_STUFF_FOUND)
