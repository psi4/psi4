set(headers_list "")
# List of headers
list(APPEND headers_list MOInfo.h globals.h Params.h )

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
list(APPEND sources_list T3_UHF_ABC.cc triples.cc ET_ABB.cc cache.cc ET_RHF.cc count_ijk.cc T3_grad_UHF_AAA.cc ET_AAB.cc transpose_integrals.cc ET_AAA.cc test_abc_loops.cc T3_grad_UHF_AAB.cc ET_UHF_AAB.cc get_moinfo.cc ET_UHF_AAA.cc ET_UHF_ABB.cc T3_UHF_AAB.cc EaT_RHF.cc ET_BBB.cc ET_UHF_BBB.cc T3_UHF_AAA.cc T3_grad_UHF_BBA.cc T3_grad_UHF_BBB.cc T3_grad_RHF.cc )

# If you want to remove some sources specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM sources_list "")
else()
   list(REMOVE_ITEM sources_list "")
endif()


# Build static library
add_library(cctriples STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(cctriples mints)
set_property(GLOBAL APPEND PROPERTY PSILIB cctriples)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(cctriples custom_boost)
endif()

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/psi4/bin/cctriples)
