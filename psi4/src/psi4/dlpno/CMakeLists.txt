list(APPEND sources
  dlpno.cc
  mp2.cc
  wrapper.cc
  sparse.cc
  )
psi4_add_module(bin dlpno sources)

target_link_libraries(dlpno PUBLIC diis)
