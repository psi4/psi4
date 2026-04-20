list(APPEND sources
  LibXCfunctional.cc
  factory.cc
  functional.cc
  superfunctional.cc
  )
psi4_add_module(lib functional sources)
target_precompile_headers(functional REUSE_FROM mints)
target_link_libraries(functional
  PUBLIC
    Libxc::xc
  )
