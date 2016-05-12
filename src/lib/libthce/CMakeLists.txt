set(headers_list "")
# List of headers
list(APPEND headers_list laplace.h thce.h lreri.h thcew.h )

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
list(APPEND sources_list lreri.cc thcew.cc laplace.cc thce.cc )

# If you want to remove some sources specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM sources_list "")
else()
   list(REMOVE_ITEM sources_list "")
endif()


# Build static library
add_library(thce STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(thce mints)
set_property(GLOBAL APPEND PROPERTY LIBLIST thce)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(thce custom_boost)
endif()

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/psi4/lib/libthce)
