list(APPEND sources
  wrapper.cc
  )
psi4_add_module(bin sapt sources)
target_precompile_headers(sapt REUSE_FROM mints)
