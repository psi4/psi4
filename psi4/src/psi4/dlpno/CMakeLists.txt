list(APPEND sources
  dlpno.cc
  mp2.cc
  ccsd.cc
  triples.cc
  wrapper.cc
  sparse.cc
  )
psi4_add_module(bin dlpno sources)
if (psi4_ENABLE_PRECOMPILE_HEADERS)
    target_precompile_headers(dlpno REUSE_FROM psio)
endif()

target_link_libraries(dlpno PUBLIC diis)

# test/ contains a standalone unit test for invert_map. It links no psi4
# library, so unlike the libpsio and libiwl tests it is built on every platform.
add_subdirectory(test)
