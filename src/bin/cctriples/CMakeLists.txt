set(headers_list "")
# List of headers
list(APPEND headers_list MOInfo.h globals.h Params.h )

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
list(APPEND sources_list T3_UHF_ABC.cc triples.cc ET_ABB.cc cache.cc ET_RHF.cc count_ijk.cc T3_grad_UHF_AAA.cc ET_AAB.cc transpose_integrals.cc ET_AAA.cc test_abc_loops.cc T3_grad_UHF_AAB.cc ET_UHF_AAB.cc get_moinfo.cc ET_UHF_AAA.cc ET_UHF_ABB.cc T3_UHF_AAB.cc EaT_RHF.cc ET_BBB.cc ET_UHF_BBB.cc T3_UHF_AAA.cc T3_grad_UHF_BBA.cc T3_grad_UHF_BBB.cc T3_grad_RHF.cc )

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
add_library(cctriples STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(cctriples mints)
set(libs_to_merge cctriples ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(cctriples custom_boost)
endif()
install(TARGETS cctriples ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libcctriples)
