set(headers_list "")
# List of headers
list(APPEND headers_list MOInfo.h globals.h backsort.h Params.h yoshimine.h )

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
list(APPEND sources_list transform_two_uhf.cc transform_two_mp2r12a_t.cc transform_one.cc frozen_core.cc transform_two_mp2.cc transform_two_backtr_uhf.cc yoshimine.cc transqt.cc semicanonical.cc cleanup.cc transform_two.cc backsort.cc transform_two_mp2r12a_gr.cc )

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
add_library(transqt STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(transqt mints)
set(libs_to_merge transqt ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(transqt custom_boost)
endif()
install(TARGETS transqt ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libtransqt)
