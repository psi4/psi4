# List of headers (needed to avoid an "empty list error" from CMake)
set(headers_list "")

set(sources_list "")
# List of sources
list(APPEND sources_list wrapper.cc atom_grids.f90 dma.f90 gdma.f90 input.F90  timing.F90)

# If you want to remove some sources specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM sources_list "")
else()
   list(REMOVE_ITEM sources_list "")
endif()


# Build static library
add_library(gdma STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(gdma mints)
set_property(GLOBAL APPEND PROPERTY LIBLIST gdma)

