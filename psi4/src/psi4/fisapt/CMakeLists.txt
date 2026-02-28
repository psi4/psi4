list(APPEND sources
  fisapt.cc
  local2.cc
  )
psi4_add_module(bin fisapt sources)

target_link_libraries(fisapt PUBLIC diis)
