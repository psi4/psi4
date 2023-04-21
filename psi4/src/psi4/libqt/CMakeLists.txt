list(APPEND sources
  3d_array.cc
  blas_intfc.cc
  blas_intfc23.cc
  cc_excited.cc
  dirprd_block.cc
  dot_block.cc
  dx_read.cc
  dx_write.cc
  fill_sym_matrix.cc
  invert.cc
  lapack_intfc.cc
  mat_print.cc
  newmm_rking.cc
  normalize.cc
  pople.cc
  probabil.cc
  ras_set.cc
  reorder_qt.cc
  schmidt.cc
  schmidt_add.cc
  solve_pep.cc
  timer.cc
  )

psi4_add_module(lib qt sources)
target_compile_definitions(qt
  PRIVATE
    FC_SYMBOL=${FC_SYMBOL}
  )
