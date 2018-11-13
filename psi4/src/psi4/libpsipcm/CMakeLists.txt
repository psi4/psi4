if(TARGET PCMSolver::pcm)
  list(APPEND sources
    psipcm.cc
    )
  psi4_add_module(lib psipcm sources)
  target_link_libraries(psipcm
    PUBLIC
      PCMSolver::pcm
    )
endif()
