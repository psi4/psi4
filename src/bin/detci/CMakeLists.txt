set(headers_list "")
# List of headers
list(APPEND headers_list tpool.h globals.h civect.h ci_tol.h ciwave.h odometer.h structs.h slaterd.h )

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
list(APPEND sources_list graphset.cc sigma.cc tpool.cc mpn.cc check_energy.cc calc_pt_block.cc get_mo_info.cc calc_hd_block.cc s3_block_bz.cc ssq.cc sem.cc vector.cc printing.cc shift.cc opdm.cc olsenupdt.cc b2brepl.cc s2v.cc h0block.cc s1.cc sem_test.cc ints.cc s3v.cc detci.cc params.cc og_addr.cc compute_cc.cc tpdm.cc time.cc s2.cc slaterd.cc slater.cc olsengraph.cc form_ov.cc set_ciblks.cc stringlist.cc import_vector.cc odometer.cc misc.cc s3.cc calc_d.cc civect.cc s1v.cc mitrush_iter.cc )

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
