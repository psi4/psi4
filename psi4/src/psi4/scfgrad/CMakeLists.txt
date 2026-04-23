list(APPEND sources
  jk_grad.cc
  response.cc
  scf_grad.cc
  wrapper.cc
  )
psi4_add_module(bin scfgrad sources)
if (${psi4_ENABLE_PRECOMPILE_HEADERS})
    target_precompile_headers(scfgrad REUSE_FROM psio)
endif()

target_link_libraries(scfgrad PUBLIC scf_solver)

# Link to common optional dependencies interface library
target_link_libraries(scfgrad PUBLIC psi4_optional_deps)
