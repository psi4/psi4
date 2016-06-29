set(headers_list "")
# List of headers
list(APPEND headers_list ccsd.h blas_mangle.h blas.h frozen_natural_orbitals.h )

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
list(APPEND sources_list frozen_natural_orbitals.cc triples.cc ccsd.cc lowmemory_triples.cc sortintegrals.cc coupled_pair.cc mp2.cc blas.cc df_cc_residual.cc df_t1_transformation.cc df_ccsd.cc opdm.cc quadratic.cc diis.cc df_scs.cc fnocc.cc linear.cc )

# If you want to remove some sources specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM sources_list "")
else()
   list(REMOVE_ITEM sources_list "")
endif()


# Build static library
add_library(fnocc STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(fnocc mints)
set_property(GLOBAL APPEND PROPERTY PSILIB fnocc)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(fnocc custom_boost)
endif()

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/psi4/bin/fnocc)
