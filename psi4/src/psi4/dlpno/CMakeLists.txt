list(APPEND sources
  dlpno.cc
  mp2.cc
  ccsd.cc
  triples.cc
  quads.cc
  wrapper.cc
  sparse.cc
  )
psi4_add_module(bin dlpno sources)

target_link_libraries(dlpno PUBLIC diis)
target_link_libraries(dlpno PRIVATE Einsums::Einsums)
