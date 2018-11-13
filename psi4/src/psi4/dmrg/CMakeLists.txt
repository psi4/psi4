list(APPEND sources
  dmrgscf.cc
  )
psi4_add_module(bin dmrg sources)
if(TARGET CheMPS2::chemps2)
  target_link_libraries(dmrg
    PUBLIC
      CheMPS2::chemps2
    )
endif()
