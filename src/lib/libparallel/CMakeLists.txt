include(GlobUtils)

set(headers_list "")
# Get headers with .h extension
get_h_headers(headers_list)
# If you want to remove some headers specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM headers_list "")
else()
   list(REMOVE_ITEM headers_list "")
endif()
if(NOT MPI_FOUND)
   list(REMOVE_ITEM headers_list "mpi_wrapper.h")
endif()
# Sort alphabetically
list(SORT headers_list)

set(sources_list "")
# Get sources with .cc extension
get_cc_sources(sources_list)
# If you want to remove some sources specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM sources_list "")
else()
   list(REMOVE_ITEM sources_list "")
endif()
if(NOT MPI_FOUND)
   list(REMOVE_ITEM sources_list "mpi_wrapper.cc")
endif()

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}" "${headers_list}")
write_to_cloc_list("${to_count}")

# Build static library
add_library(parallel STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(parallel mints)
set(libs_to_merge parallel ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(parallel custom_boost)
endif()
install(TARGETS parallel ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libparallel)
