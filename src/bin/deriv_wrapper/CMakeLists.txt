set(headers_list "")

set(sources_list "")
# List of sources
list(APPEND sources_list main.cc )

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
add_library(deriv_wrapper STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(deriv_wrapper mints)
set(libs_to_merge deriv_wrapper ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(deriv_wrapper custom_boost)
endif()
install(TARGETS deriv_wrapper ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libderiv_wrapper)
