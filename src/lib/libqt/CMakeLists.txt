set(headers_list "")
# List of headers
list(APPEND headers_list lapack_intfc_mangle.h blas_intfc_mangle.h lapack_intfc.h slaterdset.h blas_intfc23_mangle.h qt.h )

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
list(APPEND sources_list zmat_point.cc lapack_intfc.cc slaterdset.cc strncpy.cc sort.cc mat_in.cc dx_write.cc orient_fragment.cc rotate_vecs.cc dirprd_block.cc cc_wfn.cc pople.cc schmidt_add.cc mat_print.cc filter.cc ras_set.cc fill_sym_matrix.cc rootfind.cc cc_excited.cc probabil.cc dot_block.cc timer.cc dx_read.cc blas_intfc.cc normalize.cc newmm_rking.cc 3d_array.cc blas_intfc23.cc reorder_qt.cc ci_wfn.cc schmidt.cc invert.cc solve_pep.cc v_3.cc david.cc )

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
add_library(qt STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(qt mints)
set(libs_to_merge qt ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(qt custom_boost)
endif()
install(TARGETS qt ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libqt)
