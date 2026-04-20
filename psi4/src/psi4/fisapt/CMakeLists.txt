list(APPEND sources
  fisapt.cc
  local2.cc
  )
psi4_add_module(bin fisapt sources)
target_precompile_headers(fisapt REUSE_FROM mints)

target_link_libraries(fisapt PUBLIC diis)
