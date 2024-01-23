list(APPEND sources
  LibXCfunctional.cc
  GauXCfunctional.cc
  factory.cc
  functional.cc
  superfunctional.cc
  )

if(TARGET gauxc::gauxc)
  list(APPEND sources
    GauXCfunctional.cc
  )
endif()


psi4_add_module(lib functional sources)

target_link_libraries(functional
  PUBLIC
    Libxc::xc
  )

if(TARGET gauxc::gauxc)
  target_link_libraries(fock
    PUBLIC 
      gauxc::gauxc 
    )
endif()


