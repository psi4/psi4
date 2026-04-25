list(APPEND sources
  wrapper.cc
  )
psi4_add_module(bin sapt sources)
if (${psi4_ENABLE_PRECOMPILE_HEADERS})
    target_precompile_headers(sapt REUSE_FROM psio)
endif()
