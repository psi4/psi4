list(APPEND sources
  fisapt.cc
  local2.cc
  )
psi4_add_module(bin fisapt sources)
if (${psi4_ENABLE_PRECOMPILE_HEADERS})
    target_precompile_headers(fisapt REUSE_FROM psio)
endif()

target_link_libraries(fisapt PUBLIC diis)
