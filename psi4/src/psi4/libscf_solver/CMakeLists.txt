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
target_precompile_headers(scf_solver REUSE_FROM mints)
if(ENABLE_GTFOCK)
    add_dependencies(scf_solver INTERFACE GTFockInt GTFockCInt)
endif()

target_link_libraries(scf_solver PUBLIC diis)

# Link to common optional dependencies interface library
target_link_libraries(scf_solver PUBLIC psi4_optional_deps)
