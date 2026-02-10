list(APPEND sources
  algebra_interface.cc
  block_matrix.cc
  block_vector.cc
  matrix_base.cc
  mcscf.cc
  sblock_matrix.cc
  sblock_vector.cc
  scf.cc
  scf_F.cc
  scf_Feff.cc
  scf_G.cc
  scf_S_inverse_sqrt.cc
  scf_canonicalize_MO.cc
  scf_check_orthonormality.cc
  scf_compute_energy.cc
  scf_density_matrix.cc
  scf_diis.cc
  scf_energy.cc
  scf_guess_occupation.cc
  scf_initial_guess.cc
  scf_iterate_scf_equations.cc
  scf_pairs.cc
  scf_print_eigenvectors_and_MO.cc
  scf_read_so_oei.cc
  scf_read_so_tei.cc
  scf_save_info.cc
  vector_base.cc
  )
psi4_add_module(bin mcscf sources)

if(TARGET gauxc::gauxc)
  target_include_directories(mcscf
    PUBLIC
      $<TARGET_PROPERTY:gauxc::gauxc,INTERFACE_INCLUDE_DIRECTORIES>
    )
endif()
