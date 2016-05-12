set(headers_list "")
# List of headers
list(APPEND headers_list defines.h dcft.h )

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
list(APPEND sources_list dcft_lambda_UHF.cc main.cc AO_contribute.cc dcft.cc dcft_gradient.cc dcft_density_UHF.cc dcft_integrals_UHF.cc dcft_tau_UHF.cc dcft_oo_UHF.cc dcft_qc.cc half_transform.cc dcft_compute.cc dcft_compute_UHF.cc dcft_compute_RHF.cc dcft_scf_UHF.cc dcft_memory.cc dcft_n_representability.cc dcft_mp2_UHF.cc dcft_triples.cc dcft_intermediates_UHF.cc dcft_energy_RHF.cc dcft_energy_UHF.cc dcft_gradient_RHF.cc dcft_gradient_UHF.cc dcft_integrals_RHF.cc dcft_intermediates_RHF.cc dcft_lambda_RHF.cc dcft_mp2_RHF.cc dcft_oo_RHF.cc dcft_scf_RHF.cc dcft_tau_RHF.cc dcft_density_RHF.cc dcft_sort_mo_tpdm.cc dcft_df_tensor.cc)

# If you want to remove some sources specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM sources_list "")
else()
   list(REMOVE_ITEM sources_list "")
endif()


# Build static library
add_library(dcft STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(dcft mints)
set_property(GLOBAL APPEND PROPERTY PSILIB dcft)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(dcft custom_boost)
endif()

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/psi4/bin/dcft)
