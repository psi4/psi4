list(APPEND sources
  LibXCfunctional.cc
  factory.cc
  functional.cc
  superfunctional.cc
  )
psi4_add_module(lib functional sources)
if (${psi4_ENABLE_PRECOMPILE_HEADERS})
    target_precompile_headers(functional REUSE_FROM psio)
endif()
