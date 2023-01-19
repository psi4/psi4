list(APPEND sources
  block_matrix.cc
  dsyev_ascending.cc
  dsyev_descending.cc
  eigsort.cc
  eivout.cc
  flin.cc
  init_array.cc
  init_matrix.cc
  int_array.cc
  long_int_array.cc
  lubksb.cc
  ludcmp.cc
  print_array.cc
  print_mat.cc
  sq_rsp.cc
  sq_to_tri.cc
  tri_to_sq.cc
  tstart.cc
  zero.cc
  )
psi4_add_module(lib ciomr sources)
