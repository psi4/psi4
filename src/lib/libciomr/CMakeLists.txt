set(headers_list "")
# List of headers
list(APPEND headers_list libciomr.h )

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
list(APPEND sources_list mxmb.cc eigout.cc init_array.cc init_matrix.cc block_matrix.cc tqli.cc flin.cc sq_rsp.cc add_arr.cc int_array.cc ffile.cc zero.cc print_mat.cc tri_to_sq.cc rsp.cc add_mat.cc mmult.cc dot.cc tred2.cc tstart.cc eivout.cc eigsort.cc ludcmp.cc print_array.cc long_int_array.cc sq_to_tri.cc lubksb.cc )

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
add_library(ciomr STATIC ${sources_list})
# Specify dependencies for the library (if any)
#add_dependencies(ciomr )
set(libs_to_merge ciomr ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(ciomr custom_boost)
endif()
install(TARGETS ciomr ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libciomr)
