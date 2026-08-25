list(APPEND sources
  dlpno.cc
  mp2.cc
  ccsd.cc
  triples.cc
  wrapper.cc
  sparse.cc
  )

if(TARGET Einsums::Einsums)
    list(APPEND sources quads.cc)
endif()

psi4_add_module(bin dlpno sources)
if (psi4_ENABLE_PRECOMPILE_HEADERS)
    target_precompile_headers(dlpno REUSE_FROM psio)
endif()

target_link_libraries(dlpno PUBLIC diis)

if(TARGET Einsums::Einsums)
    target_link_libraries(dlpno PRIVATE Einsums::Einsums)
endif()
