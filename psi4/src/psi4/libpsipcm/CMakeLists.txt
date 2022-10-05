if(TARGET PCMSolver::pcm)
  list(APPEND sources
    psipcm.cc
    )
  psi4_add_module(lib psipcm sources)
  target_link_libraries(psipcm
    PUBLIC
      PCMSolver::pcm
      mints
      Libint2::cxx
    )
    # Because psipcm.cc includes potentialint.h, which includes libint2's engine.h
    # we have to make sure we reduce the optimization level to match libmints
    set_property(TARGET psipcm PROPERTY CXX_STANDARD 14)
endif()
