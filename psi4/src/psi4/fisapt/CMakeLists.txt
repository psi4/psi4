list(APPEND sources
  fisapt.cc
  )
psi4_add_module(bin fisapt sources)

target_link_libraries(fisapt PUBLIC diis)
