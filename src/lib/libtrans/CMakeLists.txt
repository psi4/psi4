set(headers_list "")
# List of headers
list(APPEND headers_list mospace.h integraltransform_functors.h integraltransform.h )

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
list(APPEND sources_list integraltransform_moinfo.cc integraltransform_tpdm_unrestricted.cc integraltransform_sort_so_tei.cc integraltransform_tpdm.cc integraltransform_oei.cc integraltransform_sort_mo_tpdm.cc integraltransform_tei_2nd_half.cc integraltransform_tei_1st_half.cc integraltransform_sort_so_tpdm.cc integraltransform.cc integraltransform_tpdm_restricted.cc integraltransform_dpd_id.cc mospace.cc integraltransform_tei.cc )

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
add_library(trans STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(trans mints)
set(libs_to_merge trans ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(trans custom_boost)
endif()
install(TARGETS trans ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libtrans)
