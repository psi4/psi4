list(APPEND sources
  buf_close.cc
  buf_fetch.cc
  buf_flush.cc
  buf_init.cc
  buf_put.cc
  buf_wrt.cc
  buf_wrt_mat.cc
  buf_wrt_val.cc
  rdone.cc
  wrtone.cc
  )
psi4_add_module(lib iwl sources)
