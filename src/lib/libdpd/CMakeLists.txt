set(headers_list "")
# List of headers
list(APPEND headers_list dpd.h )

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
list(APPEND sources_list T3_AAA.cc buf4_dump.cc trace42_13.cc file2_dirprd.cc buf4_mat_irrep_rd_block.cc cc3_sigma_RHF_ic.cc file4_mat_irrep_init.cc file2_axpy.cc dot14.cc buf4_mat_irrep_close.cc file4_mat_irrep_row_init.cc buf4_sort.cc file2_close.cc T3_RHF.cc block_matrix.cc file4_mat_irrep_row_close.cc file4_print.cc buf4_mat_irrep_init.cc buf4_mat_irrep_wrt_block.cc file4_mat_irrep_rd_block.cc buf4_mat_irrep_shift31.cc file2_init.cc buf4_mat_irrep_row_wrt.cc buf4_scmcopy.cc trans4_mat_irrep_wrt.cc contract422.cc contract444.cc contract244.cc buf4_mat_irrep_row_init.cc set_default.cc buf4_close.cc buf4_mat_irrep_init_block.cc buf4_print.cc buf4_mat_irrep_wrt.cc contract424.cc T3_RHF_ic.cc file2_axpbycz.cc buf4_symm2.cc file4_mat_irrep_rd.cc trans4_mat_irrep_rd.cc dot24.cc file2_trace.cc file2_dot.cc buf4_axpy.cc cc3_sigma_UHF.cc file2_mat_init.cc buf4_dot.cc buf4_scm.cc buf4_sort_ooc.cc close.cc contract442.cc trans4_mat_irrep_init.cc file2_mat_rd.cc init.cc file2_dot_self.cc memfree.cc buf4_mat_irrep_shift13.cc dot13.cc file4_mat_irrep_row_rd.cc T3_AAB.cc file2_scm.cc dot23.cc buf4_axpbycz.cc buf4_mat_irrep_row_zero.cc file2_mat_close.cc trans4_mat_irrep_close.cc buf4_mat_irrep_rd.cc 4mat_irrep_print.cc 3d_sort.cc file4_mat_irrep_row_wrt.cc file4_mat_irrep_wrt.cc buf4_symm.cc file2_print.cc buf4_mat_irrep_close_block.cc error.cc file4_init.cc file4_close.cc file2_mat_print.cc file4_init_nocache.cc trans4_close.cc file4_mat_irrep_close.cc file4_mat_irrep_wrt_block.cc contract222.cc buf4_dirprd.cc buf4_dot_self.cc buf4_init.cc buf4_mat_irrep_row_rd.cc file2_copy.cc file4_cache.cc trans4_mat_irrep_shift31.cc file4_mat_irrep_row_zero.cc contract444_df.cc file2_mat_wrt.cc buf4_copy.cc trans4_init.cc buf4_mat_irrep_row_close.cc file2_cache.cc buf4_sort_axpy.cc trans4_mat_irrep_shift13.cc cc3_sigma_RHF.cc )

# If you want to remove some sources specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM sources_list "")
else()
   list(REMOVE_ITEM sources_list "")
endif()

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}" "${headers_list}")
write_to_cloc_list("${to_count}")

# Build static library
add_library(dpd STATIC ${sources_list})
# Specify dependencies for the library (if any)
#add_dependencies(dpd )
set(libs_to_merge dpd ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(dpd custom_boost)
endif()
install(TARGETS dpd ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libdpd)
