set(sources_list params.cc h0block.cc printing.cc sigma.cc b2brepl.cc import_vector.cc slater.cc calc_d.cc ints.cc s1v.cc slaterd.cc misc.cc mitrush_iter.cc s2v.cc stringlist.cc diag_h.cc compute_mpn.cc civect.cc odometer.cc s3_block_bz.cc tpdm.cc compute_cc.cc og_addr.cc s3v.cc ciwave.cc detci.cc olsengraph.cc sem.cc vector.cc form_ov.cc olsenupdt.cc sem_test.cc get_mo_info.cc opdm.cc set_ciblks.cc slater_matel.cc)

SET_SOURCE_FILES_PROPERTIES(s1v.cc PROPERTIES COMPILE_FLAGS -O3)
SET_SOURCE_FILES_PROPERTIES(s2v.cc PROPERTIES COMPILE_FLAGS -O3)
SET_SOURCE_FILES_PROPERTIES(s3v.cc PROPERTIES COMPILE_FLAGS -O3)
SET_SOURCE_FILES_PROPERTIES(vector.cc PROPERTIES COMPILE_FLAGS -O3)
SET_SOURCE_FILES_PROPERTIES(tpdm.cc PROPERTIES COMPILE_FLAGS -O3)
psi4_add_module(bin detci sources_list mints thce trans)
