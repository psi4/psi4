list(APPEND sources
  dfep2.cc
  )
psi4_add_module(bin dfep2 sources)

if(TARGET cuEST::cuEST)
  target_link_libraries(dfep2
    PUBLIC
      cuEST::cuEST
    )
endif()

