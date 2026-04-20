list(APPEND sources
  dlpno.cc
  mp2.cc
  ccsd.cc
  triples.cc
  wrapper.cc
  sparse.cc
  )
psi4_add_module(bin dlpno sources)
target_precompile_headers(dlpno REUSE_FROM mints)

target_link_libraries(dlpno PUBLIC diis)
