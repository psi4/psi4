list(APPEND sources
  LibXCfunctional.cc
  GauXCfunctional.cc
  factory.cc
  functional.cc
  superfunctional.cc
  )
psi4_add_module(lib functional sources)
target_link_libraries(functional
  PUBLIC
    Libxc::xc
    gauxc::gauxc
  )
