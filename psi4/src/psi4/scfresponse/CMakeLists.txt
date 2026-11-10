list(APPEND sources
  diis.cc
  ks_response.cc
  rks_response.cc
  scf_response.cc
  uks_response.cc
  )
psi4_add_module(bin scfresponse sources)
