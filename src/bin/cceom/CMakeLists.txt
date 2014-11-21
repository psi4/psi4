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
list(APPEND sources_list local.cc form_diagonal.cc WmnieSD.cc WabejDS.cc diagSS.cc check_sum.cc sigma_full.cc hbar_extra.cc sigmaSS.cc restart_with_root.cc cache.cc sigmaCC3.cc rzero.cc overlap.cc sigmaCC3_RHF.cc local_guess.cc get_params.cc read_guess.cc dgeev_eom.cc schmidt_add.cc norm_HC1.cc restart.cc sort_amps.cc cc3_HC1ET1.cc FSD.cc write_Rs.cc follow_root.cc WmnefDD.cc WamefSD.cc get_eom_params.cc sort_C.cc FDD.cc WmaijDS.cc precondition.cc WabefDD.cc get_moinfo.cc sigmaDD.cc WmbejDD.cc cc2_hbar_extra.cc cc3_HC1.cc c_clean.cc cceom.cc amp_write.cc cc2_sigma.cc WmnijDD.cc WbmfeDS.cc sigmaDS.cc sigmaSD.cc hbar_norms.cc norm.cc WnmjeDS.cc diag.cc )

# If you want to remove some sources specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM sources_list "read_guess.cc")
else()
   list(REMOVE_ITEM sources_list "read_guess.cc")
endif()

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}" "${headers_list}")
write_to_cloc_list("${to_count}")

# Build static library
add_library(cceom STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(cceom mints)
set(libs_to_merge cceom ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(cceom custom_boost)
endif()
install(TARGETS cceom ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libcceom)
