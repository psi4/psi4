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
list(APPEND sources_list local.cc lambda_residuals.cc LHX1Y1.cc roa.cc hbar_extra.cc analyze.cc cache.cc cc2_LHX1Y1.cc cc2_X1.cc optrot.cc get_params.cc cc2_sort_X.cc converged.cc X2.cc print_X.cc ccresponse.cc diis.cc X1.cc save_X.cc LHX1Y2.cc compute_X.cc polar.cc linresp.cc cc2_LHX1Y2.cc get_moinfo.cc cc2_hbar_extra.cc HXY.cc amp_write.cc init_X.cc update_X.cc pseudopolar.cc sort_lamps.cc pertbar.cc LHX2Y2.cc denom.cc cc2_X2.cc sort_X.cc sort_pert.cc LCX.cc preppert.cc )

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
add_library(ccresponse STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(ccresponse mints)
set(libs_to_merge ccresponse ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(ccresponse custom_boost)
endif()
install(TARGETS ccresponse ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libccresponse)
