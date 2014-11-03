include(GlobUtils)

set(sources_list "")
# Get sources with Fortran extension
get_fortran_sources(sources_list)
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
add_library(erd STATIC ${sources_list})
# Specify dependencies for the library (if any)
#add_dependencies(erd )
set(libs_to_merge erd ${libs_to_merge} PARENT_SCOPE)
install(TARGETS erd ARCHIVE DESTINATION lib)
