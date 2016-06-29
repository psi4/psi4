set(headers_list "")
# List of headers
list(APPEND headers_list libmoinfo.h model_space.h moinfo_scf.h moinfo_base.h slater_determinant.h moinfo.h )

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
list(APPEND sources_list model_space_build.cc slater_determinant.cc moinfo_model_space.cc moinfo_mappings.cc moinfo_slaterdeterminant.cc moinfo.cc model_space.cc moinfo_scf.cc moinfo_base.cc )

# If you want to remove some sources specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM sources_list "")
else()
   list(REMOVE_ITEM sources_list "")
endif()


# Build static library
add_library(moinfo STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(moinfo mints)
set_property(GLOBAL APPEND PROPERTY LIBLIST moinfo)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(moinfo custom_boost)
endif()

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/psi4/lib/libmoinfo)
