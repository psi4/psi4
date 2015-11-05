set(headers_list "")
# List of headers
list(APPEND headers_list tpool.h civect.h ci_tol.h ciwave.h odometer.h structs.h slaterd.h )

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
list(APPEND sources_list params.cc h0block.cc printing.cc sigma.cc b2brepl.cc import_vector.cc slater.cc calc_d.cc ints.cc s1v.cc slaterd.cc misc.cc mitrush_iter.cc s2v.cc stringlist.cc diag_h.cc compute_mpn.cc civect.cc odometer.cc s3_block_bz.cc tpdm.cc compute_cc.cc og_addr.cc s3v.cc tpool.cc ciwave.cc detci.cc olsengraph.cc sem.cc vector.cc form_ov.cc olsenupdt.cc sem_test.cc get_mo_info.cc opdm.cc set_ciblks.cc mcscf.cc slater_matel.cc)

# If you want to remove some sources specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM sources_list "graphset.cc")
else()
   list(REMOVE_ITEM sources_list "graphset.cc")
endif()

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}" "${headers_list}")
write_to_cloc_list("${to_count}")

# Build static library
add_library(detci STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(detci mints thce trans)
set_property(GLOBAL APPEND PROPERTY PSILIB detci)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(detci custom_boost)
endif()

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libdetci)
