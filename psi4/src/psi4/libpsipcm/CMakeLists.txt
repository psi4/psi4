if (ENABLE_PCMSOLVER)
    set(sources_list psipcm.cc host_writer.cpp)
    psi4_add_library(psipcm sources_list mints PCMSolver)
endif ()
