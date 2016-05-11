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
list(APPEND sources_list transtwo_rhf.cc cache.cc get_params.cc transtwo_uhf.cc transqt.cc transone.cc get_moinfo.cc file_build_presort.cc idx_error.cc idx_permute_presort.cc )

# If you want to remove some sources specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM sources_list "")
else()
   list(REMOVE_ITEM sources_list "")
endif()


# Build static library
add_library(transqt2 STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(transqt2 mints)
set_property(GLOBAL APPEND PROPERTY PSILIB transqt2)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(transqt2 custom_boost)
endif()

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/psi4/transqt2)
