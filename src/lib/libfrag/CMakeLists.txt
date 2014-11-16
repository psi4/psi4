set(headers_list "")
# List of headers
list(APPEND headers_list AtomSet.h MBE.h libfrag.h Set.h Connections.h Capper.h Fragmenter.h Embedder.h MBEFrag.h MBEFragSet.h LibFragTypes.h FragOptions.h Properties.h GMBE.h BSSEer.h )

# If you want to remove some headers specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM headers_list "Properties.h")
else()
   list(REMOVE_ITEM headers_list "Properties.h")
endif()
# Sort alphabetically
list(SORT headers_list)

set(sources_list "")
# List of sources
list(APPEND sources_list BSSEer.cc MBEFragSet.cc MBE.cc Fragmenter.cc MBEFrag.cc libfrag.cc Capper.cc Embedder.cc GMBE.cc AtomSet.cc Properties.cc FragOptions.cc Connections.cc )

# If you want to remove some sources specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM sources_list "Properties.cc")
else()
   list(REMOVE_ITEM sources_list "Properties.cc")
endif()

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}" "${headers_list}")
write_to_cloc_list("${to_count}")

# Build static library
add_library(frag STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(frag mints parallel)
set(libs_to_merge frag ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(frag custom_boost)
endif()
install(TARGETS frag ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libfrag)
