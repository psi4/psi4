list(APPEND sources
  frozen_natural_orbitals.cc
  triples.cc
  ccsd.cc
  lowmemory_triples.cc
  sortintegrals.cc
  coupled_pair.cc
  mp2.cc
  blas.cc
  df_cc_residual.cc
  df_t1_transformation.cc
  df_ccsd.cc
  opdm.cc
  quadratic.cc
  diis.cc
  df_scs.cc
  fnocc.cc
  linear.cc
  )
psi4_add_module(bin fnocc sources)

