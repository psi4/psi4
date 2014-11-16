set(headers_list "")
# List of headers
list(APPEND headers_list MOInfo.h globals.h Params.h Local.h )

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
list(APPEND sources_list local.cc local_polar.cc transpert.cc cc_memcheck.cc classify.cc f_spinad.cc cache.cc build_F_RHF.cc local_magnetic.cc sort_tei.cc cphf_F.cc fock.cc get_params.cc scf_check.cc f_sort.cc cphf_B.cc idx_permute_multipass.cc d_sort.cc b_sort.cc build_B_RHF.cc build_abcd_packed.cc c_sort.cc distribute.cc sort_oei.cc idx_permute.cc e_spinad.cc file_build_multipass.cc get_moinfo.cc file_build.cc ccsort.cc d_spinad.cc file_build_presort.cc e_sort.cc idx_error.cc denom.cc sort_pert.cc domain_print.cc idx_permute_presort.cc )

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
add_library(ccsort STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(ccsort mints)
set(libs_to_merge ccsort ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(ccsort custom_boost)
endif()
install(TARGETS ccsort ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libccsort)
