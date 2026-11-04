list(APPEND sources
  convert.cc
  integrals.cc
  )

if(TARGET ambit::ambit)
  psi4_add_module(lib ambit_interface sources)
endif()
