set(sources_list psipcm.cc )
psi4_add_library(psipcm sources_list)
target_link_libraries(psipcm mints PCMSolver::PCMSolver)
