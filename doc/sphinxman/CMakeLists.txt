
# Settle docs-specific dependencies
find_package(Perl)
if(NOT PERL_FOUND)
    message(WARNING "No Perl, no docs. Pre-built documentation at http://sirius.chem.vt.edu/psi4manual/latest/index.html")
endif()
find_package(Sphinx)
if(NOT SPHINX_FOUND)
    message(WARNING "No Sphinx, no docs. Pre-built documentation at http://sirius.chem.vt.edu/psi4manual/latest/index.html")
endif()

# Proceed to build
if(PERL_FOUND AND SPHINX_FOUND)
#file(MAKE_DIRECTORY ${directory})
    set(CCSD "${CMAKE_CURRENT_SOURCE_DIR}")
    set(CCBD "${CMAKE_CURRENT_BINARY_DIR}")
    
    configure_file(
        "${CCSD}/source/conf.py.in"
        "${CCBD}/source/conf.py"
        @ONLY)
    
    
    set(STATICDOC abbr_accents.rst adc.rst appendices.rst basissets.rst basissets_byelement.rst bestpractices_py.rst bibliography.rst cbs.rst cbs_eqn.rst contributing.rst cp.rst customizing.rst db.rst detci.rst dfmp2.rst diatomic.rst documentation.rst energy.rst external.rst freq.rst glossary_psivariables.rst index.rst intercalls.rst interfacing.rst introduction.rst methods.rst notes_c.rst notes_py.rst oeprop.rst opt.rst optking.rst plugins.rst programming.rst prop.rst psithonfunc.rst psithoninput.rst quickaddalias.rst basissets_tables.rst quickadddatabase.rst sapt.rst sowreap.rst testsuite.rst tutorial.rst dft_byfunctional.rst scf.rst mrcc.rst psimrcc.rst cc.rst dft.rst fnocc.rst installfile.rst optionshandling.rst basissets_byfamily.rst psipep.rst pep0001.rst dftd3.rst occ.rst dfocc.rst molden.rst dcft.rst proc_py.rst psithonmol.rst cfour.rst mrcc_table_energy.rst cfour_table_energy.rst cfour_table_grad.rst)
    
    # TODO ok to add index, appendices to staticdoc?
    
    # Static docs reST files
    # * compute relative path btwn top_srcdir and objdir/doc/sphinxman
    # * copy all over to objdir, substituting literalinclude along the way
    # * duplicate index.rst as contents.rst for toc link
    # * copy fake psi4.py over to objdir
    # * copy logo image file (configure_file corrupts file)
    execute_process(COMMAND ${PYTHON_EXECUTABLE} ${CCSD}/extract_relpath.py ${CMAKE_SOURCE_DIR} ${CCBD} OUTPUT_VARIABLE SFNX_INCLUDE OUTPUT_STRIP_TRAILING_WHITESPACE)
    foreach(rst ${STATICDOC})
       configure_file(
        "${CCSD}/source/${rst}"
        "${CCBD}/source/${rst}"
        @ONLY)
    endforeach(rst)
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy "${CCBD}/source/index.rst" "${CCBD}/source/contents.rst")
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy "${CCSD}/source/psi4.py" "${CCBD}/source")
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy "${CCSD}/source/psi4banner.png" "${CCBD}/source")
    
    # Autodoc basis sets
    add_custom_command(
        COMMAND ${PERL_EXECUTABLE} ${CCSD}/document_bases.pl ${CCSD}
        OUTPUT  ${CCBD}/source/autodoc_psi4bases.rst
        DEPENDS ${CCSD}/document_bases.pl
        COMMENT "Autodoc basis sets"
        # TODO leftout dependency @top_srcdir@/lib/basis/*.gbs \
    )
    
    
    # Autodoc binary scratch files
    add_custom_command(
        COMMAND ${PYTHON_EXECUTABLE} ${CCSD}/document_psifiles.py ${CCSD}
        OUTPUT  ${CCBD}/source/autodoc_psifiles.rst
                ${top_srcdir}/lib/python/p4const/psifiles.py
        DEPENDS ${CCSD}/document_psifiles.py
                ${top_srcdir}/include/psifiles.h
        COMMENT "Autodoc binary scratch files"
    )

    # Autodoc physical constants
    add_custom_command(
        COMMAND ${PERL_EXECUTABLE} ${CCSD}/document_physconst.pl ${CCSD}
        OUTPUT  ${top_srcdir}/lib/python/p4const/physconst.py
        DEPENDS ${CCSD}/document_physconst.pl
                ${top_srcdir}/include/physconst.h
        COMMENT "Autodoc physical constants"
    )

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
                ${top_srcdir}/src/bin/psi4/read_options.cc
        COMMENT "Autodoc options c-side"
    )

    # Autodoc PSI variables c-side
    add_custom_command(
        COMMAND ${CMAKE_COMMAND} -E make_directory source/autodir_psivariables
        COMMAND ${PERL_EXECUTABLE} ${CCSD}/document_psivariables.pl ${CCSD}
        OUTPUT  ${CCBD}/source/autodoc_psivariables_bymodule.rst
                ${CCBD}/source/autodir_psivariables/module__scf.rst  # representative
        DEPENDS ${CCSD}/document_psivariables.pl
                ${top_srcdir}/src/*/*/*.cc
        COMMENT "Autodoc PSI variables c-side"
    )

    # Autodoc Cfour PSI variables
    add_custom_command(
        COMMAND ${CMAKE_COMMAND} -E make_directory source/autodir_psivariables
        COMMAND ${PYTHON_EXECUTABLE} ${CCSD}/document_cfour.py ${CCSD}
        OUTPUT  ${CCBD}/source/autodir_psivariables/module__cfour.rst
        DEPENDS ${CCSD}/document_cfour.py
                ${top_srcdir}/lib/python/qcdb/cfour.py
        COMMENT "Autodoc Cfour PSI variables"
    )



#    add_custom_command(
#        COMMAND ${CMAKE_COMMAND} -E make_directory source/autodir_plugins
#    )







    
    
    set(SPHINX_TARGETS html)


    # The doc target will build all documentation targets.
    add_custom_target(sphinxman
        DEPENDS ${CCBD}/source/autodoc_psi4bases.rst
                ${CCBD}/source/autodoc_psifiles.rst
                ${top_srcdir}/lib/python/p4const/physconst.py
                ${CCBD}/source/autodoc_abbr_options_c.rst
                ${CCBD}/source/autodoc_glossary_options_c.rst
                ${CCBD}/source/autodoc_options_c_bymodule.rst
                ${CCBD}/source/autodir_options_c/module__scf.rst  # representative
                ${CCBD}/source/autodir_options_c/scf__scf_type.rst  # representative
                ${CCBD}/source/autodoc_psivariables_bymodule.rst
                ${CCBD}/source/autodir_psivariables/module__scf.rst  # representative
                ${CCBD}/source/autodir_psivariables/module__cfour.rst
        COMMENT "Building documentation with Sphinx")

    foreach (SPHINX_TARGET ${SPHINX_TARGETS})
        add_custom_command(
            TARGET sphinxman POST_BUILD
            COMMAND ${SPHINX_EXECUTABLE} 
                # -v  # Repeat for greater verbosity
                # -q  # Enable for quiet mode 
                # -Q  # Enable for even more quiet mode
                # -E  # don't use a saved environment, always read all file
                # -a  # write all files; default is to only write new and changed files
                -b ${SPHINX_TARGET}
                #-d "${CCBD}/_doctrees"      # cache w/ pickled reST
                -d "${CCBD}/${SPHINX_TARGET}/doctrees"  # cache w/ pickled reST
                -c "${CCBD}/source"         # whereabouts of conf.py
                "${CCBD}/source"            # build from
                "${CCBD}/${SPHINX_TARGET}"  # build to
            COMMENT "Generating ${SPHINX_TARGET} documentation ...")
    endforeach()

endif(PERL_FOUND AND SPHINX_FOUND)



### Provide options for controlling different types of output
##option(SPHINX_OUTPUT_HTML "Output standalone HTML files" ON)
##option(SPHINX_OUTPUT_MAN "Output man pages" ON)
