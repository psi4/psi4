set(headers_list "")
# List of headers
list(APPEND headers_list LibBabel.h Psi4OBUnitCell.h )

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
list(APPEND sources_list Psi4OBUnitCell.cc LibBabel.cc )

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
add_library(babel STATIC ${sources_list})
if(NOT BABEL_FOUND)
  add_dependencies(babel OpenBabel)
endif()
include_directories(${BABEL_INCLUDE_DIR})
if(APPLE)
  target_link_libraries(babel ${BABEL_LIBRARY_DIR}/libopenbabel.dylib)
else()
  target_link_libraries(babel ${BABEL_LIBRARY_DIR}/libopenbabel.so)
endif()
# Specify dependencies for the library (if any)
#add_dependencies(babel )
set_property(GLOBAL APPEND PROPERTY LIBLIST babel)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(babel custom_boost)
endif()

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libbabel)
