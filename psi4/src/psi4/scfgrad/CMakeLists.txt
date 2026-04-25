list(APPEND sources
  jk_grad.cc
  response.cc
  scf_grad.cc
  wrapper.cc
  )
psi4_add_module(bin scfgrad sources)
if (${psi4_ENABLE_PRECOMPILE_HEADERS})
    # Skip PCH when OpenOrbitalOptimizer is enabled to avoid conflicts with armadillo
    # scfgrad links to scf_solver which has USING_OpenOrbitalOptimizer as PUBLIC,
    # so scfgrad inherits the definition (see libmints/CMakeLists.txt for full explanation)
    if(NOT TARGET OpenOrbitalOptimizer::OpenOrbitalOptimizer)
        target_precompile_headers(scfgrad REUSE_FROM psio)
    endif()
endif()

target_link_libraries(scfgrad PUBLIC scf_solver)
