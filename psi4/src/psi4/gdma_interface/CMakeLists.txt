list(APPEND sources
  wrapper.cc
  )
psi4_add_module(bin gdma_interface sources)
if(TARGET gdma::gdma)
  target_link_libraries(gdma_interface
    PUBLIC
      gdma::gdma
    )
endif()
