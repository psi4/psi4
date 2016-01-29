set(headers_list "")
# List of headers
list(APPEND headers_list points.h direct_screening.h
            apps.h hamiltonian.h link.h cubature.h v.h solver.h
            jk_independent.h jk.h jk_independent_impl.h gridblocker.h soscf.h)

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
list(APPEND sources_list apps.cc v.cc hamiltonian.cc points.cc
            cubature.cc solver.cc link.cc direct_screening.cc
            wrapper.cc jk.cc DiskJK.cc PKJK.cc DirectJK.cc DFJK.cc
            CDJK.cc FastDFJK.cc PSJK.cc GTFockJK.cc soscf.cc)

# If you want to remove some sources specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM sources_list "")
else()
   list(REMOVE_ITEM sources_list "")
endif()


# Build static library
add_library(fock STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(fock mints)
set_property(GLOBAL APPEND PROPERTY LIBLIST fock)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(fock custom_boost)
endif()

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/psi4/libfock)
