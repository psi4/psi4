if(TARGET cppe::cppe)
  list(APPEND sources
    psipe.cc
    )
  psi4_add_module(lib psipe sources)
  target_link_libraries(psipe
    PUBLIC
      cppe::cppe
    )
endif()
