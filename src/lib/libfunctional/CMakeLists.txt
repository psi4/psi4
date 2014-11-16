set(headers_list "")
# List of headers
list(APPEND headers_list PW91_Cfunctional.h PZ81_Cfunctional.h PBE_Cfunctional.h LYP_Cfunctional.h superfunctional.h P86_Cfunctional.h wpbec_functional.h xfunctional.h wpbex_functional.h VWN5_Cfunctional.h PW92_Cfunctional.h functional.h FT97_Cfunctional.h utility.h VWN3_Cfunctional.h FT97B_Xfunctional.h cfunctional.h )

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
list(APPEND sources_list PBE_Cfunctional.cc superfactory.cc PZ81_Cfunctional.cc PW91_Cfunctional.cc wpbec_functional.cc VWN3_Cfunctional.cc superfunctional.cc LYP_Cfunctional.cc FT97B_Xfunctional.cc wpbex_functional.cc FT97_Cfunctional.cc VWN5_Cfunctional.cc PW92_Cfunctional.cc P86_Cfunctional.cc cfunctional.cc factory.cc utility.cc functional.cc xfunctional.cc )

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
add_library(functional STATIC ${sources_list})
# Specify dependencies for the library (if any)
#add_dependencies(functional )
set(libs_to_merge functional ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(functional custom_boost)
endif()
install(TARGETS functional ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libfunctional)
