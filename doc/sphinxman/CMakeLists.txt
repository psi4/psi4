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
    set(STATICDOC abbr_accents.rst adc.rst appendices.rst basissets.rst basissets_byelement.rst bestpractices_py.rst bibliography.rst cbs.rst cbs_eqn.rst contributing.rst cp.rst customizing.rst db.rst detci.rst dfmp2.rst diatomic.rst documentation.rst energy.rst external.rst freq.rst glossary_psivariables.rst index.rst intercalls.rst interfacing.rst introduction.rst methods.rst notes_c.rst notes_py.rst oeprop.rst opt.rst optking.rst plugins.rst programming.rst prop.rst psithonfunc.rst psithoninput.rst quickaddalias.rst basissets_tables.rst quickadddatabase.rst sapt.rst sowreap.rst testsuite.rst tutorial.rst dft_byfunctional.rst scf.rst mrcc.rst psimrcc.rst cc.rst dft.rst fnocc.rst installfile.rst optionshandling.rst basissets_byfamily.rst psipep.rst pep0001.rst dftd3.rst occ.rst dfocc.rst molden.rst dcft.rst proc_py.rst psithonmol.rst cfour.rst mrcc_table_energy.rst cfour_table_energy.rst cfour_table_grad.rst prog_basissets.rst libefp.rst)
    # * compute relative path btwn top_srcdir and objdir/doc/sphinxman
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c 
        "import sys; import os; print '/' + os.path.relpath('${CMAKE_SOURCE_DIR}', '${CCBD}') + '/../'" 
        OUTPUT_VARIABLE SFNX_INCLUDE OUTPUT_STRIP_TRAILING_WHITESPACE)
    # * copy all over to objdir, substituting literalinclude along the way
    foreach(rst ${STATICDOC})
       configure_file(
        ${CCSD}/source/${rst}
        ${CCBD}/source/${rst}
        @ONLY)
    endforeach(rst)
    # * duplicate index.rst as contents.rst for toc link
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CCBD}/source/index.rst ${CCBD}/source/contents.rst)
    # * copy fake psi4.py over to objdir
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CCSD}/source/psi4.py ${CCBD}/source)
    # * copy logo image file (configure_file corrupts file)
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CCSD}/source/psi4banner.png ${CCBD}/source)
    

    # Autodoc basis sets
    add_custom_command(
        COMMAND ${PERL_EXECUTABLE} ${CCSD}/document_bases.pl ${CCSD}
        OUTPUT  ${CCBD}/source/autodoc_psi4bases.rst
        DEPENDS ${CCSD}/document_bases.pl
                ${CMAKE_SOURCE_DIR}/lib/basis/*.gbs
        COMMENT "Autodoc basis sets")
    

    # Autodoc binary scratch files
    add_custom_command(
        COMMAND ${PYTHON_EXECUTABLE} ${CCSD}/document_psifiles.py ${CCSD}
        OUTPUT  ${CCBD}/source/autodoc_psifiles.rst
        DEPENDS ${CCSD}/document_psifiles.py
                ${CMAKE_SOURCE_DIR}/include/psifiles.h
        COMMENT "Autodoc binary scratch files")


    # Autodoc physical constants
    add_custom_command(
        COMMAND ${PERL_EXECUTABLE} ${CCSD}/document_physconst.pl ${CCSD}
        OUTPUT  ${CCBD}/source/autodoc_physconst.rst
        DEPENDS ${CCSD}/document_physconst.pl
                ${CMAKE_SOURCE_DIR}/include/physconst.h
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
                ${CMAKE_SOURCE_DIR}/src/bin/psi4/read_options.cc
        COMMENT "Autodoc options c-side")


    # Autodoc PSI variables c-side
    add_custom_command(
        COMMAND ${CMAKE_COMMAND} -E make_directory source/autodir_psivariables
        COMMAND ${PERL_EXECUTABLE} ${CCSD}/document_psivariables.pl ${CCSD}
        OUTPUT  ${CCBD}/source/autodoc_psivariables_bymodule.rst
                ${CCBD}/source/autodir_psivariables/module__scf.rst  # representative
        DEPENDS ${CCSD}/document_psivariables.pl
                ${CMAKE_SOURCE_DIR}/src/*/*/*.cc
        COMMENT "Autodoc PSI variables c-side")


    # Autodoc Cfour PSI variables
    add_custom_command(
        COMMAND ${CMAKE_COMMAND} -E make_directory source/autodir_psivariables
        COMMAND ${PYTHON_EXECUTABLE} ${CCSD}/document_cfour.py ${CCSD}
        OUTPUT  ${CCBD}/source/autodir_psivariables/module__cfour.rst
        DEPENDS ${CCSD}/document_cfour.py
                ${CMAKE_SOURCE_DIR}/lib/python/qcdb/cfour.py
                ${CCBD}/source/autodoc_psivariables_bymodule.rst
        COMMENT "Autodoc Cfour PSI variables")


    # Autodoc databases
    add_custom_command(
        COMMAND ${PYTHON_EXECUTABLE} ${CCSD}/document_databases.py ${CCSD}
        OUTPUT  ${CCBD}/source/autodoc_available_databases.rst
        DEPENDS ${CCSD}/document_databases.py
                ${CMAKE_SOURCE_DIR}/lib/databases/*.py
        COMMENT "Autodoc databases")


    # Autodoc test cases
    add_custom_command(
        COMMAND ${PERL_EXECUTABLE} ${CCSD}/document_tests.pl ${CCSD}
        OUTPUT  ${CCBD}/source/autodoc_testsuite_corepsi4.rst
                ${CCBD}/source/autodoc_testsuite_dftd3.rst
        DEPENDS ${CCSD}/document_tests.pl
                ${CMAKE_SOURCE_DIR}/tests/*/input.dat
                ${CMAKE_SOURCE_DIR}/tests/*/*/input.dat
        COMMENT "Autodoc test cases")


    # Autodoc driver
    add_custom_command(
        COMMAND ${PYTHON_EXECUTABLE} ${CCSD}/document_driver.py ${CCSD}
        OUTPUT  ${CCBD}/source/autodoc_driver.rst
        DEPENDS ${CMAKE_SOURCE_DIR}/lib/python/*.py
                ${CCSD}/document_driver.py
        COMMENT "Autodoc driver")


    # Autodoc psi4 module
    add_custom_command(
        COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_BINARY_DIR} ctest -R docs-psimod
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/tests/docs-psimod/autodoc_psimod.rst ${CCBD}/source
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/tests/docs-psimod/autodoc_man_psi4.txt ${CCBD}/source
        OUTPUT  ${CCBD}/source/autodoc_psimod.rst
                ${CCBD}/source/autodoc_man_psi4.txt
        DEPENDS ${CMAKE_SOURCE_DIR}/src/*/*/*.cc
                ${CMAKE_SOURCE_DIR}/tests/docs-psimod/input.dat
                ${CMAKE_BINARY_DIR}/bin/psi4${CMAKE_EXECUTABLE_SUFFIX}
        COMMENT "Autodoc psi4 module")


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
        DEPENDS ${CMAKE_SOURCE_DIR}/lib/python/functional.py
                ${CMAKE_SOURCE_DIR}/tests/docs-dft/input.dat
                ${CMAKE_BINARY_DIR}/bin/psi4${CMAKE_EXECUTABLE_SUFFIX}
        COMMENT "Autodoc functionals")


    # Autodoc auxiliary basis sets
    add_custom_command(
        COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_BINARY_DIR} ctest -R docs-bases
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/tests/docs-bases/autodoc_basissets_byfamily.rst ${CCBD}/source
        OUTPUT ${CCBD}/source/autodoc_basissets_byfamily.rst
        DEPENDS ${CMAKE_SOURCE_DIR}/lib/python/p4util/basislist.py
                ${CMAKE_SOURCE_DIR}/lib/python/p4util/basislistdunning.py
                ${CMAKE_SOURCE_DIR}/lib/python/p4util/basislistother.py
                ${CMAKE_SOURCE_DIR}/tests/docs-bases/input.dat
                ${CMAKE_BINARY_DIR}/bin/psi4${CMAKE_EXECUTABLE_SUFFIX}
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
                ${CMAKE_SOURCE_DIR}/plugins/*/*.py
                ${CMAKE_SOURCE_DIR}/plugins/*/doc.rst
        COMMENT "Autodoc plugins")
    
    
    # Autodoc efp fragments
    add_custom_command(
        COMMAND ${PYTHON_EXECUTABLE} ${CCSD}/document_efpfrag.py ${CCSD}
        OUTPUT  ${CCBD}/source/autodoc_available_efpfrag.rst
        DEPENDS ${CCSD}/document_efpfrag.py
                ${CMAKE_SOURCE_DIR}/lib/efpfrag/*.efp
        COMMENT "Autodoc efp fragments")


    # Establishing dependencies for docs targets

    # * announcement
    if (LATEX_COMPILER AND PDFLATEX_COMPILER)
        message("-- Documentation targets available: sphinxman (html), sphinxmini (quick html), sphinxpdf (LaTeX --> PDF)")
    else()
        message("-- Documentation targets available: sphinxman (html), sphinxmini (quick html)")
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
                ${CCBD}/source/autodoc_driver.rst
                ${CCBD}/source/autodoc_psimod.rst
                ${CCBD}/source/autodoc_man_psi4.txt
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

    # * full, proper build with LaTeX
    if (LATEX_COMPILER AND PDFLATEX_COMPILER)
        add_custom_target(sphinxpdf
            DEPENDS sphinxman
            COMMENT "Preparing Sphinx HTML and PDF documentation build")
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


endif(PERL_FOUND AND SPHINX_FOUND)

