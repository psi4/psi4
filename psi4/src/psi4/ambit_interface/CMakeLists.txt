list(APPEND sources
  convert.cc
  integrals.cc
  )

if(TARGET ambit::ambit)
  psi4_add_module(lib ambit_interface sources)
  target_link_libraries(ambit_interface
    PUBLIC
      ambit::ambit
    )
endif()
