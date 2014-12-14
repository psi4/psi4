set(headers_list "")
# List of headers
list(APPEND headers_list MCSCF_indpairs.h tpool.h globals.h civect.h ci_tol.h ciwave.h odometer.h structs.h slaterd.h )

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
list(APPEND sources_list MCSCF_cleanup.cc MCSCF_step.cc graphset.cc params.cc shift.cc MCSCF_detcas.cc MCSCF_thetas.cc h0block.cc printing.cc sigma.cc MCSCF_diis.cc b2brepl.cc import_vector.cc s1.cc slater.cc MCSCF_f_act.cc calc_d.cc ints.cc s1v.cc slaterd.cc MCSCF_get_mo_info.cc calc_hd_block.cc misc.cc s2.cc ssq.cc MCSCF_gradient.cc calc_pt_block.cc mitrush_iter.cc s2v.cc stringlist.cc MCSCF_hessian.cc check_energy.cc mpn.cc s3.cc time.cc MCSCF_indpairs.cc civect.cc odometer.cc s3_block_bz.cc tpdm.cc MCSCF_lagcalc.cc compute_cc.cc og_addr.cc s3v.cc tpool.cc MCSCF_read_dens.cc detci.cc olsengraph.cc sem.cc vector.cc MCSCF_ref_orbs.cc form_ov.cc olsenupdt.cc sem_test.cc MCSCF_rotate_orbs.cc get_mo_info.cc opdm.cc set_ciblks.cc )

# If you want to remove some sources specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM sources_list "calc_pt_block.cc" "graphset.cc")
else()
   list(REMOVE_ITEM sources_list "calc_pt_block.cc" "graphset.cc")
endif()

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}" "${headers_list}")
write_to_cloc_list("${to_count}")

# Build static library
add_library(detci STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(detci mints)
set(libs_to_merge detci ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(detci custom_boost)
endif()
install(TARGETS detci ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libdetci)
