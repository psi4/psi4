set(headers_list "")
# List of headers
list(APPEND headers_list defines.h dcft.h )

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
list(APPEND sources_list dcft_lambda.cc main.cc AO_contribute.cc dcft.cc dcft_gradient.cc dcft_density.cc dcft_integrals.cc dcft_tau.cc dcft_oo.cc dcft_qc.cc half_transform.cc dcft_compute.cc dcft_scf.cc dcft_memory.cc dcft_n_representability.cc dcft_energy.cc dcft_mp2.cc dcft_triples.cc dcft_intermediates.cc )

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
add_library(dcft STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(dcft mints)
set(libs_to_merge dcft ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(dcft custom_boost)
endif()
install(TARGETS dcft ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libdcft)
