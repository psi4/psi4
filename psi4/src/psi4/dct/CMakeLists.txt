list(APPEND sources
  AO_contribute.cc
  dct.cc
  dct_compute.cc
  dct_compute_RHF.cc
  dct_compute_UHF.cc
  dct_density_RHF.cc
  dct_density_UHF.cc
  dct_df_operations.cc
  dct_df_tensor.cc
  dct_energy_RHF.cc
  dct_energy_UHF.cc
  dct_gradient.cc
  dct_gradient_RHF.cc
  dct_gradient_UHF.cc
  dct_integrals_RHF.cc
  dct_integrals_UHF.cc
  dct_intermediates_RHF.cc
  dct_intermediates_UHF.cc
  dct_lambda_RHF.cc
  dct_lambda_UHF.cc
  dct_memory.cc
  dct_mp2_RHF.cc
  dct_mp2_UHF.cc
  dct_oo_RHF.cc
  dct_oo_UHF.cc
  dct_qc.cc
  dct_relaxed_density_UHF.cc
  dct_scf_RHF.cc
  dct_scf_UHF.cc
  dct_sort_mo_tpdm.cc
  dct_tau_RHF.cc
  dct_tau_UHF.cc
  dct_triples.cc
  dct_validate.cc
  half_transform.cc
  main.cc
  )

# See https://github.com/psi4/psi4/issues/923
if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18.0.3 AND
       (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 18.0 OR
        CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 18.0))
           message(WARNING "Intel C++ (ICC) 18.0 fails to compile dct_gradient_UHF.cc with aggressive optimization flags")
           set_source_files_properties(dct_gradient_UHF.cc PROPERTIES
           COMPILE_FLAGS -O1)
   endif ()
endif ()

psi4_add_module(bin dct sources)

if(TARGET gauxc::gauxc)
  target_include_directories(dct
    PUBLIC
      $<TARGET_PROPERTY:gauxc::gauxc,INTERFACE_INCLUDE_DIRECTORIES>
    )
endif()

target_link_libraries(dct PUBLIC diis)
