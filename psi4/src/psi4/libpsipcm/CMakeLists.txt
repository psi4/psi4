if(TARGET PCMSolver::pcm)
  list(APPEND sources
    psipcm.cc
    )
  psi4_add_module(lib psipcm sources)
endif()

if(TARGET cuEST::cuEST)
  set(cuest_sources cuestpcm.cc)
  psi4_add_module(lib cuestpcm cuest_sources)
  target_link_libraries(cuestpcm
    PUBLIC
      mints
      cuEST::cuEST
    )
endif()
