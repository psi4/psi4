set(headers_list "")
# List of headers
list(APPEND headers_list StringArray.h aiohandler.h MOFile.h BinaryFile.h psio.h config.h psio.hpp )

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
   list(REMOVE_ITEM sources_list "tocdel.cc")
else()
   list(REMOVE_ITEM sources_list "tocdel.cc")
endif()


# Build static library
add_library(psio STATIC ${sources_list})
# Specify dependencies for the library (if any)
#add_dependencies(psio )
set_property(GLOBAL APPEND PROPERTY LIBLIST psio)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(psio custom_boost)
endif()

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/psi4/libpsio)
