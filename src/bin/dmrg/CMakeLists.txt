set(headers_list "")
# List of headers
list(APPEND headers_list )

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
list(APPEND sources_list dmrgscf.cc )

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

find_package(GSL REQUIRED)
find_package(HDF5 REQUIRED)

# Build static library
add_library(dmrg STATIC ${sources_list})
include_directories(${CMAKE_BINARY_DIR}/interfaces/chemps2/chemps2/install/include)
target_link_libraries(dmrg ${GSL_LIBRARIES} ${HDF5_LIBRARIES})

# Specify dependencies for the library (if any)
add_dependencies(dmrg mints)
set_property(GLOBAL APPEND PROPERTY PSILIB dmrg)
add_dependencies(dmrg chemps2)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(dmrg custom_boost)
endif()

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/libdmrg)
