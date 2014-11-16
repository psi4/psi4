set(headers_list "")
# List of headers
list(APPEND headers_list cholesky.h qr.h 3index.h schwarz.h pstensor.h fitter.h dftensor.h denominator.h )

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
list(APPEND sources_list dealias.cc dftensor.cc pstensor.cc denominator.cc pseudotrial.cc fitter.cc schwarz.cc fittingmetric.cc qr.cc cholesky.cc )

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
add_library(3index STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(3index mints)
set(libs_to_merge 3index ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(3index custom_boost)
endif()
install(TARGETS 3index ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/lib3index)
