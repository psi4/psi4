list(APPEND sources
  cuhf.cc
  frac.cc
  hf.cc
  mom.cc
  rhf.cc
  rohf.cc
  sad.cc
  uhf.cc
  )
psi4_add_module(lib scf_solver sources)
if (psi4_ENABLE_PRECOMPILE_HEADERS)
    # Skip PCH when OpenOrbitalOptimizer is enabled to avoid conflicts with armadillo
    if(NOT TARGET OpenOrbitalOptimizer::OpenOrbitalOptimizer)
        target_precompile_headers(scf_solver REUSE_FROM psio)
    endif()
endif()
if(ENABLE_GTFOCK)
    add_dependencies(scf_solver INTERFACE GTFockInt GTFockCInt)
endif()

target_link_libraries(scf_solver PUBLIC diis)

# OpenOrbitalOptimizer support for SCF solver
# PUBLIC because scf_solver headers may expose OpenOrbitalOptimizer types
if(TARGET OpenOrbitalOptimizer::OpenOrbitalOptimizer)
  target_compile_definitions(scf_solver
    PUBLIC
      USING_OpenOrbitalOptimizer
    )
  target_link_libraries(scf_solver
    PUBLIC
      OpenOrbitalOptimizer::OpenOrbitalOptimizer
    )
endif()
