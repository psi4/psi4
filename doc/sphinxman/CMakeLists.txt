

#file(MAKE_DIRECTORY ${directory})
#set(rstsrcdir "${CMAKE_CURRENT_SOURCE_DIR}/source")

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/source/conf.py.in"
    "${CMAKE_CURRENT_BINARY_DIR}/source/conf.py"
    @ONLY)


set(STATICDOC abbr_accents.rst adc.rst appendices.rst basissets.rst basissets_byelement.rst bestpractices_py.rst bibliography.rst cbs.rst cbs_eqn.rst contributing.rst cp.rst customizing.rst db.rst detci.rst dfmp2.rst diatomic.rst documentation.rst energy.rst external.rst freq.rst glossary_psivariables.rst index.rst intercalls.rst interfacing.rst introduction.rst methods.rst notes_c.rst notes_py.rst oeprop.rst opt.rst optking.rst plugins.rst programming.rst prop.rst psithonfunc.rst psithoninput.rst quickaddalias.rst basissets_tables.rst quickadddatabase.rst sapt.rst sowreap.rst testsuite.rst tutorial.rst dft_byfunctional.rst scf.rst mrcc.rst psimrcc.rst cc.rst dft.rst fnocc.rst installfile.rst optionshandling.rst basissets_byfamily.rst psipep.rst pep0001.rst dftd3.rst occ.rst dfocc.rst molden.rst dcft.rst proc_py.rst psithonmol.rst cfour.rst mrcc_table_energy.rst cfour_table_energy.rst cfour_table_grad.rst)

# TODO ok to add index, appendices to staticdoc?

# Process manual doc reST files
# * compute relative path btwn top_srcdir and objdir/doc/sphinxman
# * copy all over to objdir, substituting literalinclude along the way
# * duplicate index.rst as contents.rst for toc link
# * copy fake psi4.py over to objdir
# * copy logo image file (configure_file corrupts file)
execute_process(COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/relpath.py ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} OUTPUT_VARIABLE SFNX_INCLUDE OUTPUT_STRIP_TRAILING_WHITESPACE)
foreach(rst ${STATICDOC})
   configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/source/${rst}"
    "${CMAKE_CURRENT_BINARY_DIR}/source/${rst}"
    @ONLY)
endforeach(rst)
execute_process(COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_BINARY_DIR}/source/index.rst" "${CMAKE_CURRENT_BINARY_DIR}/source/contents.rst")
execute_process(COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/source/psi4.py" "${CMAKE_CURRENT_BINARY_DIR}/source")
execute_process(COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/source/psi4banner.png" "${CMAKE_CURRENT_BINARY_DIR}/source")





set(SPHINX_TARGETS html)

find_package(Sphinx)
if(SPHINX_FOUND)

    # The doc target will build all documentation targets.
    add_custom_target(sphinxman
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
                -d "${CMAKE_CURRENT_BINARY_DIR}/_doctrees"  # cache w/ pickled reST
                -c "${CMAKE_CURRENT_BINARY_DIR}/source" 
                "${CMAKE_CURRENT_BINARY_DIR}/source"
                "${CMAKE_CURRENT_BINARY_DIR}/${SPHINX_TARGET}"
            COMMENT "Generating ${SPHINX_TARGET} documentation ...")
    endforeach()
else()
    message(WARNING "Sphinx not found. Cannot generate documentation!")
endif(SPHINX_FOUND)




### Provide options for controlling different types of output
##option(SPHINX_OUTPUT_HTML "Output standalone HTML files" ON)
##option(SPHINX_OUTPUT_MAN "Output man pages" ON)
