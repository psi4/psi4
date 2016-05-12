set(headers_list "")

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
list(APPEND sources_list c_sort.cc d_sort.cc e_sort.cc fock.cc scf_check.cc a_spinad.cc cache.cc d_spinad.cc e_spinad.cc memcheck.cc sort_tei_rhf.cc b_spinad.cc cctransort.cc denom.cc f_sort.cc pitzer2qt.cc sort_tei_uhf.cc)

# If you want to remove some sources specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM sources_list "")
else()
   list(REMOVE_ITEM sources_list "")
endif()


# Build static library
add_library(cctransort STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(cctransort dpd iwl trans options qt psio mints)
set_property(GLOBAL APPEND PROPERTY PSILIB cctransort)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(cctransort custom_boost)
endif()

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/psi4/bin/cctransort)
