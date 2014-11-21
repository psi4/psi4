set(headers_list "")
# List of headers
list(APPEND headers_list chkpt.h config.h )

# If you want to remove some headers specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM headers_list "")
else()
   list(REMOVE_ITEM headers_list "")
endif()
# Sort alphabetically
list(SORT headers_list)

set(sources_list "")
# List of sources
list(APPEND sources_list grad.cc max_am.cc exist.cc natom_per_fragment.cc num_unique_atom.cc nshell.cc puream.cc eccsd.cc zmat.cc nallatom.cc openpi.cc usotao.cc frzcpi.cc sprim.cc num_unique_shell.cc rottype.cc ncalcs.cc fgeom.cc ecorr.cc e_labeled.cc natom.cc atom_dummy.cc label.cc nso.cc escf.cc fock.cc clsdpi.cc disp_irrep.cc override_occ.cc frzvpi.cc cdsalc2cd.cc symoper.cc nfragment.cc emp2.cc sym_label.cc nmo.cc keyword.cc atom_pos.cc rotconst.cc scf.cc efzc.cc iopen.cc lagr.cc cdsalcpi.cc sopi.cc prefix.cc disp.cc rot_symm_num.cc us2s.cc nprim.cc geom.cc close.cc nsymhf.cc init.cc exps.cc phase_check.cc nfzv.cc eref.cc contr.cc snuc.cc enuc.cc shells_per_am.cc fragment_coeff.cc shell_transm.cc nallatom_per_fragment.cc nirreps.cc usotbf.cc stype.cc ua2a.cc sloc_new.cc contr_full.cc nfzc.cc zvals.cc orbspi.cc ccvecs.cc rref.cc sloc.cc cartrep.cc e_t.cc statespi.cc exist_add_prefix.cc irr_labs.cc vib_freqs.cc ict.cc evals.cc felement.cc nref_per_fragment.cc ref.cc snumg.cc etot.cc nao.cc am2canon_shell_order.cc )

# If you want to remove some sources specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM sources_list "vib_freqs.cc")
else()
   list(REMOVE_ITEM sources_list "vib_freqs.cc")
endif()

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}" "${headers_list}")
write_to_cloc_list("${to_count}")

# Build static library
add_library(chkpt STATIC ${sources_list})
# Specify dependencies for the library (if any)
#add_dependencies(chkpt )
set(libs_to_merge chkpt ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(chkpt custom_boost)
endif()
install(TARGETS chkpt ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libchkpt)
