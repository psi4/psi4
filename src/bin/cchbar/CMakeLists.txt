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
list(APPEND sources_list Wabei_RHF.cc status.cc Wabij.cc cache.cc cc3_HET1.cc Fai.cc get_params.cc reference.cc cchbar.cc sort_amps.cc cc2_Zmbej.cc Wabei_BBBB_UHF.cc Wmbej.cc Wmnie.cc cc2_Wmbej.cc Wamef.cc get_moinfo.cc Wabei_RHF_FT2_a.cc Wmbij.cc F.cc Wabei_ABAB_UHF.cc HET1_Wabef.cc purge.cc tau.cc Wabei.cc norm_HET1.cc cc2_Wmbij.cc Wabei_ROHF.cc taut.cc Wabei_AAAA_UHF.cc cc2_Wabei.cc Wabei_BABA_UHF.cc )

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
add_library(cchbar STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(cchbar mints)
set(libs_to_merge cchbar ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(cchbar custom_boost)
endif()
install(TARGETS cchbar ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libcchbar)
