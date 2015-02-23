set(headers_list "")
# List of headers
list(APPEND headers_list quad.h memory_manager.h class_macros.h ref.h
libpsi4util.h )

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
list(APPEND sources_list quad.cc stl_string.cc memory_manager.cc exception.cc combinations.cc )

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
add_library(psi4util STATIC ${sources_list})
# Specify dependencies for the library (if any)
#add_dependencies(psi4util )
set_property(GLOBAL APPEND PROPERTY LIBLIST psi4util)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(psi4util custom_boost)
endif()
install(TARGETS psi4util ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libpsi4util)
