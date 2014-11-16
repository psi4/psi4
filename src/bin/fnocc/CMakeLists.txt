set(headers_list "")
# List of headers
list(APPEND headers_list ccsd.h blas_mangle.h blas.h frozen_natural_orbitals.h )

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
list(APPEND sources_list frozen_natural_orbitals.cc triples.cc ccsd.cc lowmemory_triples.cc sortintegrals.cc coupled_pair.cc mp2.cc blas.cc df_ccsd.cc opdm.cc quadratic.cc diis.cc local_triples.cc fnocc.cc linear.cc )

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
add_library(fnocc STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(fnocc mints)
set(libs_to_merge fnocc ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(fnocc custom_boost)
endif()
install(TARGETS fnocc ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libfnocc)
