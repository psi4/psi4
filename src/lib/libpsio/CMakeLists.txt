set(headers_list "")
# List of headers
list(APPEND headers_list aiohandler.h MOFile.h BinaryFile.h psio.h config.h )

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
list(APPEND sources_list rw.cc getpid.cc filemanager.cc tocwrite.cc write_entry.cc tocclean.cc read_entry.cc rename_file.cc tocscan.cc get_numvols.cc BinaryFile.cc change_namespace.cc tocdel.cc done.cc MOFile.cc get_volpath.cc toclen.cc get_address.cc close.cc init.cc read.cc get_filename.cc volseek.cc write.cc get_global_address.cc open_check.cc zero_disk.cc error.cc aio_handler.cc open.cc toclast.cc tocprint.cc get_length.cc tocread.cc filescfg.cc )

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
add_library(psio STATIC ${sources_list})
# Specify dependencies for the library (if any)
#add_dependencies(psio )
set(libs_to_merge psio ${libs_to_merge} PARENT_SCOPE)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(psio custom_boost)
endif()
install(TARGETS psio ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libpsio)
