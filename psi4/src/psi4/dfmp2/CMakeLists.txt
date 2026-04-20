list(APPEND sources
  mp2.cc
  corr_grad.cc
  wrapper.cc
  )
psi4_add_module(bin dfmp2 sources)
target_precompile_headers(dfmp2 REUSE_FROM mints)
