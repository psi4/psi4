list(APPEND sources
  AO_contribute.cc
  dcft.cc
  dcft_compute.cc
  dcft_compute_RHF.cc
  dcft_compute_UHF.cc
  dcft_density_RHF.cc
  dcft_density_UHF.cc
  dcft_df_tensor.cc
  dcft_energy_RHF.cc
  dcft_energy_UHF.cc
  dcft_gradient.cc
  dcft_gradient_RHF.cc
  dcft_gradient_UHF.cc
  dcft_integrals_RHF.cc
  dcft_integrals_UHF.cc
  dcft_intermediates_RHF.cc
  dcft_intermediates_UHF.cc
  dcft_lambda_RHF.cc
  dcft_lambda_UHF.cc
  dcft_memory.cc
  dcft_mp2_RHF.cc
  dcft_mp2_UHF.cc
  dcft_n_representability.cc
  dcft_oo_RHF.cc
  dcft_oo_UHF.cc
  dcft_qc.cc
  dcft_scf_RHF.cc
  dcft_scf_UHF.cc
  dcft_sort_mo_tpdm.cc
  dcft_tau_RHF.cc
  dcft_tau_UHF.cc
  dcft_triples.cc
  half_transform.cc
  main.cc
  )

# See https://github.com/psi4/psi4/issues/923
if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18.0.3 AND
       (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 18.0 OR
        CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 18.0))
           message(WARNING "Intel C++ (ICC) 18.0 fails to compile dcft_gradient_UHF.cc with aggressive optimization flags")
           set_source_files_properties(dcft_gradient_UHF.cc PROPERTIES
           COMPILE_FLAGS -O1)
   endif ()
endif ()

psi4_add_module(bin dcft sources)
