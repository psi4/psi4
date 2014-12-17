set(headers_list "")
# List of headers
list(APPEND headers_list iwl.h config.h )

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
list(APPEND sources_list buf_rd_arr2.cc buf_wrt_mp2.cc buf_wrt_mp2r12a.cc buf_rd_all_mp2r12a.cc rdtwo.cc buf_fetch.cc rdone.cc buf_wrt_val.cc buf_wrt_arr2.cc buf_rd_arr.cc buf_wrt_val_SI.cc buf_wrt_arr_SI_nocut.cc buf_rd_all_act.cc buf_close.cc buf_toend.cc buf_wrt_mat.cc buf_wrt_all.cc buf_wrt_arr.cc buf_wrt.cc buf_flush.cc buf_wrt_arr_SI.cc buf_rd.cc sortbuf.cc wrtone.cc buf_init.cc wrttwo.cc buf_put.cc buf_rd_all.cc )

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
add_library(iwl STATIC ${sources_list})
# Specify dependencies for the library (if any)
#add_dependencies(iwl )
set(libs_to_merge iwl ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(iwl custom_boost)
endif()
install(TARGETS iwl ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libiwl)
