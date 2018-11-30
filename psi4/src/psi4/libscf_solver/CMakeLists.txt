list(APPEND sources
  cuhf.cc
  frac.cc
  hf.cc
  mom.cc
  rhf.cc
  rohf.cc
  sad.cc
  stability.cc
  uhf.cc
  )
psi4_add_module(lib scf_solver sources)
if(ENABLE_GTFOCK)
    add_dependencies(scf_solver INTERFACE GTFockInt GTFockCInt)
endif()
